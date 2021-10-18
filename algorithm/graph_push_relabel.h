#pragma once

#include <vector>
#include <queue>

#include <tbb/scalable_allocator.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/buffered_vector.h"
#include "../datastructure/queue.h"
#include "../datastructure/stack.h"

#include <boost/circular_buffer.hpp>

namespace whfc {

template<typename T> using vec = std::vector<T, tbb::scalable_allocator<T> >;

class GraphPushRelabel {
public:

	static constexpr bool log = false;
	GraphPushRelabel(FlowHypergraph& hg, bool restrict_capacity_to_in_nodes = false) {
		first_out.assign(hg.numNodes() + 2 * hg.numHyperedges() + 1, 0);
		arcs.assign(4 * hg.numPins() + 2 * hg.numHyperedges(), Arc());

		level.assign(numNodes(), 0);
		excess.assign(numNodes(), 0);

		size_t in_node_offset = hg.numNodes(), out_node_offset = hg.numNodes() + hg.numHyperedges();

		for (Node u : hg.nodeIDs()) {
			first_out[u + 1] = first_out[u] + 2 * hg.degree(u);
		}
		assert(first_out[in_node_offset] == 2 * hg.numPins());
		first_out[out_node_offset] = 3 * hg.numPins() + hg.numHyperedges();
		for (Hyperedge e : hg.hyperedgeIDs()) {
			first_out[e + in_node_offset + 1] = first_out[e + in_node_offset] + hg.pinCount(e) + 1;
			first_out[e + out_node_offset + 1] = first_out[e + out_node_offset] + hg.pinCount(e) + 1;
		}
		assert(first_out[out_node_offset] == 3 * hg.numPins() + hg.numHyperedges());
		assert(first_out[numNodes()] == numArcs());

		vec<size_t> offsets(numNodes(), 0);

		// bridge arcs
		for (Hyperedge e : hg.hyperedgeIDs()) {
			Node e_in(e + in_node_offset);
			Node e_out(e + out_node_offset);
			size_t i = first_out[e_in], j = first_out[e_out];
			arcs[i].head = e_out;
			arcs[i].capacity = hg.capacity(e);
			arcs[i].back_arc = j;
			arcs[j].head = e_in;
			arcs[j].capacity = 0;
			arcs[j].back_arc = i;
			offsets[e_in] = 1;
			offsets[e_out] = 1;
		}

		for (Hyperedge e : hg.hyperedgeIDs()) {
			Node e_in(e + in_node_offset);
			Node e_out(e + out_node_offset);
			for (auto& p : hg.pinsOf(e)) {
				Node v = p.pin;
				// (v, e_in))
				size_t i = first_out[v] + offsets[v];
				size_t j = first_out[e_in] + offsets[e_in];
				arcs[i].head = e_in;
				arcs[i].capacity = restrict_capacity_to_in_nodes ? hg.capacity(e) : std::numeric_limits<Flow>::max();
				arcs[i].back_arc = j;
				arcs[j].head = v;
				arcs[j].capacity = 0;
				arcs[j].back_arc = i;

				// (e_out, v)
				i = i + hg.degree(v);
				j = first_out[e_out] + offsets[e_out];
				arcs[i].head = e_out;
				arcs[i].capacity = 0;
				arcs[i].back_arc = j;
				arcs[j].head = v;
				arcs[j].capacity = std::numeric_limits<Flow>::max();	// could put hg.capacity(e) as well, but want so see what happens here.
				arcs[j].back_arc = i;

				offsets[v]++;
				offsets[e_in]++;
				offsets[e_out]++;
			}
		}

		for (size_t i = 0; i < arcs.size(); ++i) {
			assert(i == arcs[arcs[i].back_arc].back_arc);
		}
	}

	Flow computeFlow(Node s, Node t) {
		source = s; target = t;

		level[source] = numNodes();
		queue.reserve(numNodes());
		while (!active_nodes.empty()) active_nodes.pop();
		// active_nodes.clear();
		// active_nodes.set_capacity(numNodes());
		for (Arc& a : arcsOf(source)) {
			Flow d = a.rcap();
			if (d > 0) {
				if (excess[a.head] == 0 && a.head != target) {
					// active_nodes.push_back(a.head);
					active_nodes.push(a.head);
				}
				push(source, a, d);
			}
		}
		LOGGER << V(active_nodes.size());
		size_t num_discharges = 0;
		while (!active_nodes.empty()) {
			if (work > 1.2 * double(numNodes() + numArcs())) {
				work = 0;
				globalRelabel();
			}
			const Node x = active_nodes.front();
			// active_nodes.pop_front();
			active_nodes.pop();
			if (num_discharges % 100000 == 0)
				LOGGER << V(num_discharges) << V(x) << V(excess[target]);
			num_discharges++;
			discharge(x);
		}
		return excess[target];
	}

	Flow dinitz() {
		queue.reserve(numNodes());
		FixedCapacityStack<Node> stack(numNodes());
		vec<size_t> current_arc(numNodes(), 0);

		auto bfs = [&] {
			queue.clear();
			int inf = numNodes() + 1;
			level.assign(numNodes(), inf);
			level[source] = 0;
			current_arc[source] = first_out[source];
			queue.push(source);
			bool target_found = false;
			while (!queue.empty()) {
				Node u = queue.pop();
				assert(queue.size() <= numNodes());
				for (Arc& a : arcsOf(u)) {
					if (level[a.head] == inf && a.rcap() > 0) {
						if (a.head != target) {
							current_arc[a.head] = first_out[a.head];
							level[a.head] = level[u] + 1;
							queue.push(a.head);
						} else {
							target_found = true;
						}
					}
				}
			}
			return target_found;
		};

		auto dfs = [&] {
			Flow f_dfs = 0;
			stack.clear();
			stack.push(source);
			while (!stack.empty()) {
				Node u = stack.top();
				Node v = invalidNode;
				for ( ; current_arc[u] < first_out[u+1]; current_arc[u]++) {
					const Arc& a = arcs[current_arc[u]];
					if ((a.head == target || level[a.head] == level[u] + 1) && a.rcap() > 0) {
						v = a.head;
						break;
					}
				}
				if (v == invalidNode) {
					stack.pop();
					level[u] = std::numeric_limits<int>::max();
				} else if (v == target) {
					int64_t lowest_bottleneck = std::numeric_limits<int64_t>::max();
					Flow d = std::numeric_limits<Flow>::max();
					for (int64_t stack_pointer = stack.size() - 1; stack_pointer >= 0; --stack_pointer) {
						Flow res = arcs[current_arc[stack.at(stack_pointer)]].rcap();
						if (res <= d) {
							d = res;
							lowest_bottleneck = stack_pointer;
						}
					}
					for (int64_t stack_pointer = stack.size() - 1; stack_pointer >= 0; --stack_pointer) {
						Arc& a = arcs[current_arc[stack.at(stack_pointer)]];
						a.flow += d;
						arcs[a.back_arc].flow -= d;
					}
					stack.popDownTo(lowest_bottleneck);
					assert(stack.size() == lowest_bottleneck + 1);
					f_dfs += d;
				} else {
					stack.push(v);
				}
			}
			return f_dfs;
		};

		Flow f = 0;
		while (bfs()) {
			Flow d = dfs();
			f += d;
		}
		return f;
	}

	void discharge(Node u) {
		while (excess[u] > 0 && level[u] < int(numNodes())) {
			int min_level = numNodes();
			for (Arc& a : arcsOf(u)) {
				Node v = a.head;
				if (a.rcap() > 0) {
					assert(level[u] <= level[v] + 1);
					if (level[u] == level[v] + 1) {
						if (excess[v] == 0 && v != target) {
							// active_nodes.push_back(v);
							active_nodes.push(v);
						}
						push(u, a, std::min(excess[u], a.rcap()));
					} else {
						min_level = std::min(min_level, level[v] + 1);
					}
				}
			}
			if (excess[u] > 0) {
				level[u] = min_level;
			}
			work += arcsOf(u).size() + 5;
		}
	}

	void globalRelabel() {
		LOGGER << "global relabel";
		queue.clear();
		int inf = numNodes();
		level.assign(numNodes(), inf);
		level[target] = 0;
		queue.push(target);
		while (!queue.empty()) {
			assert(queue.size() <= numNodes());
			Node v = queue.pop();
			for (Arc& a : arcsOf(v)) {
				Node u = a.head;
				if (level[u] == inf && arcs[a.back_arc].rcap() > 0) {
					level[u] = level[v] + 1;
					queue.push(u);
				}
			}
		}
	}

	size_t work = std::numeric_limits<size_t>::max();
	vec<int> level;
	vec<Flow> excess;
	Node source, target;
	// boost::circular_buffer<Node> active_nodes;
	std::queue<Node> active_nodes;
	LayeredQueue<Node> queue;

	struct Arc {
		Node head;
		Flow flow = 0;
		Flow capacity = 0;
		size_t back_arc = 0;
		Flow rcap() const { return capacity - flow; }
	};

	void push(Node u, Arc& a, Flow d) {
		excess[u] -= d;
		excess[a.head] += d;
		a.flow += d;
		arcs[a.back_arc].flow -= d;
	}

	vec<size_t> first_out;
	vec<Arc> arcs;
	size_t numNodes() const { return first_out.size() - 1; }
	size_t numArcs() const { return arcs.size(); }
	mutable_range<vec<Arc>> arcsOf(Node u) {
		assert(u < numNodes());
		return mutable_range<vec<Arc>>(arcs.begin() + first_out[u], arcs.begin() + first_out[u+1]);
	}
};
}
