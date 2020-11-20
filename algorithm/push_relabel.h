#include "cutter_state.h"

#include "../datastructure/queue.h"
#include "../datastructure/stack.h"

#include <boost/circular_buffer.hpp>

#include "../datastructure/bitset_reachable_sets.h"

namespace whfc {

	class PushRelabel{
	public:
		static constexpr bool same_traversal_as_grow_assimilated = false;
		static constexpr bool grow_reachable_marks_flow_sending_pins_when_marking_all_pins = true;

		static constexpr bool relabel_to_front = false;

		using Type = PushRelabel;
		using ScanList = LayeredQueue<Node>;
		using Pin = FlowHypergraph::Pin;
		using InHe = FlowHypergraph::InHe;

		using ReachableNodes = BitsetReachableNodes;			// TODO temporary so CutterState won't complain.
		using ReachableHyperedges = BitsetReachableHyperedges;

		FlowHypergraph& hg;

		LayeredQueue<Node> bfs_queue;		// for global relabeling
		FixedCapacityStack<Pin> stack;		// for flow decomposition
		int direction = 0, previous_cutter_state_direction = 0;

		std::vector<PinIndex> current_pin;
		std::vector<InHeIndex> current_hyperedge;

		std::vector<Flow> excess;
		boost::circular_buffer<Node> active_vertices_and_edges;
		std::vector<int> level;

		int max_level;
		Node source, target;

		Flow upperFlowBound = std::numeric_limits<Flow>::max();

		PushRelabel(FlowHypergraph& hg) :
				hg(hg), bfs_queue(hg.numNodes()), stack(hg.numNodes()),
				current_pin(hg.numHyperedges(), PinIndex::Invalid()),
				current_hyperedge(hg.numNodes(), InHeIndex::Invalid()),
				excess(hg.numNodes() + hg.numHyperedges(), 0),
				level(hg.numNodes() + hg.numHyperedges(), 0),
				timer("Push Relabel")
		{

		}

		void flipViewDirection() {
			direction = 1 - direction;
		}


		void reset() {

		}

		ScanList& getScanList() {
			return bfs_queue;
		}

		static constexpr bool log = false;

		bool exhaustFlow(CutterState<Type>& cs) {
			prepare(cs);
			max_level = hg.numNodes() + hg.numHyperedges();
			if (active_vertices_and_edges.capacity() < size_t(max_level)) {
				active_vertices_and_edges.set_capacity(max_level);
			}
			source = cs.sourcePiercingNodes.front().node;
			target = cs.targetPiercingNodes.front().node;
			Flow old_excess = excess[target];

			timer.start("Levels and Initial Excesses");
			// do level assignment before excesses are set --> global relabeling doesn't add anything to the active_vertices queue
			level.assign(hg.numNodes() + hg.numHyperedges(), 0);		// still TODO global relabeling instead of this
			level[source] = max_level;

			for (Hyperedge e : hg.hyperedgeIDs()) {
				current_pin[e] = hg.beginIndexPins(e);
			}
			for (Node u : hg.nodeIDs()) {
				current_hyperedge[u] = hg.beginIndexHyperedges(u);
			}

			for (InHe& in_he : hg.hyperedgesOf(source)) {
				LOGGER << "initial push" << V(in_he.e) << V(hg.residualCapacity(in_he.e) + hg.absoluteFlowReceived(in_he));
				pushToHyperedge(source, Node(in_he.e + hg.numNodes()), in_he, hg.residualCapacity(in_he.e) + hg.absoluteFlowReceived(in_he));
			}
			timer.stop("Levels and Initial Excesses");

			LOGGER << V(max_level) << V(source) << V(target) << V(excess[source]) << V(excess[target]) << V(active_vertices_and_edges.size());

			if constexpr (relabel_to_front) {
				// must insert all non-terminal vertices into the queue, and cannot insert new excess vertices during pushes

				std::vector<Node> front;
				while (!active_vertices_and_edges.empty()) {
					const Node x = active_vertices_and_edges.front();
					active_vertices_and_edges.pop_front();
					int old_level = level[x];
					discharge(x);

					if (old_level < level[x]) {
						// was relabeled --> move to front
						for (auto it = front.crbegin(); it != front.crend(); ++it) {
							active_vertices_and_edges.push_front(*it);
						}
						front.clear();
						active_vertices_and_edges.push_front(x);	// will be removed straight away but then kept in front
					} else if (level[x] < max_level) {
						front.push_back(x);
					}
				}
			} else {
				timer.start("Discharge");
				LOGGER << V(active_vertices_and_edges.size());
				print();
				while (/*excess[target] <= upperFlowBound && */ !active_vertices_and_edges.empty()) {
					const Node x = active_vertices_and_edges.front();
					active_vertices_and_edges.pop_front();
					if (x < hg.numNodes()) {
						LOGGER << "discharge node" << x;
					} else {
						LOGGER << "discharge hyperedge" << (x - hg.numNodes());
					}
					discharge(x);
					print();
				}
				timer.stop("Discharge");
			}


			active_vertices_and_edges.clear();

			finish(cs);
			timer.report(std::cout);

			Flow flow_delta = excess[target] - old_excess;
			LOGGER << V(flow_delta) << V(excess[target]);
			assert(flow_delta >= 0);
			cs.flowValue += flow_delta;
			return flow_delta > 0;
		}

		Flow growFlowOrSourceReachable(CutterState<Type>& cs) {
			prepare(cs);
			Flow f = 0;
			finish(cs);
			return f;
		}

		void growReachable(CutterState<Type>& cs) {
			prepare(cs);
			finish(cs);
		}

		Flow recycleDatastructuresFromGrowReachablePhase(CutterState<Type> &cs) {
			return 0;		// can't do anything here
		}

	private:

		void print() {
			std::cout << "----\nexcesses:\t";
			size_t i = 0;
			for (Flow& exc : excess) {
				if (i == hg.numNodes()) { std::cout << "|| "; }
				std::cout << exc << " ";
				++i;
			}
			std::cout << "\nlevels:\t\t";
			i = 0;
			for (int l : level) {
				if (i == hg.numNodes()) { std::cout << "|| "; }
				std::cout << l << " ";
				++i;
			}
			std::cout << "\nflow:\t\t";
			for (Hyperedge e : hg.hyperedgeIDs()) {
				std::cout << V(e) << " " << V(hg.flow(e)) << " ";
				for (Pin& pin : hg.pinsOf(e)) {
					std::cout << ", " << pin.pin << " : " << hg.getInHe(pin).flow;
				}
				std::cout << " | ";
			}
			std::cout << std::endl;
		}

		void discharge(Node u) {
			if (excess[u] > 0 && level[u] < max_level) {
				if (u < hg.numNodes()) {
					dischargeNode(u);
				} else {
					dischargeHyperedge(Hyperedge(u - hg.numNodes()), u);
				}
			}
		}

		void pushToHyperedge(Node u, Node e_node, InHe& inc_he, Flow residual) {
			LOGGER << "push to hyperedge" << V(inc_he.e) << V(inc_he.flow) << V(residual) << V(u) ;
			if constexpr (!relabel_to_front) {
				if (excess[e_node] == 0) {
					active_vertices_and_edges.push_back(e_node);
				}
			}
			excess[e_node] += residual;
			excess[u] -= residual;
			hg.flow(inc_he.e) -= std::min(residual, hg.absoluteFlowReceived(inc_he));
			inc_he.flow += residual;

			assert(hg.flow(inc_he.e) == out_flow(e_node));
			assert(excess[e_node] == in_flow(e_node) - out_flow(e_node));	// breaks running time in debug mode
		}

		void dischargeNode(Node u) {
			assert(u != source && u != target);
			assert(hg.degree(u) > 0);
			assert(level[u] < max_level);

			while (excess[u] > 0) {
				InHe& inc_he = hg.getInHe(current_hyperedge[u]);
				const Node e_node(inc_he.e + hg.numNodes());
				//Flow residual = hg.residualCapacity(inc_he.e) + hg.absoluteFlowReceived(inc_he);
				//if (residual > 0 && level[e_node] + 1 == level[u]) {
				//	pushToHyperedge(u, e_node, inc_he, std::min(excess[u], residual));
				if (level[e_node] + 1 == level[u]) {
					// dump it all to edge_in. never makes sense to push more than cap(e)
					// can push cap(e)-flow(e) to edge_out, and push at most flow(e) back to other pins sending flow in.
					pushToHyperedge(u, e_node, inc_he, excess[u] /* std::min(excess[u], hg.capacity(e) */ );
				} else {
					// don't advance iterator if pushed
					if (++current_hyperedge[u] == hg.endIndexHyperedges(u)) {
						// relabel
						int min_level = std::numeric_limits<int>::max();
						for (const InHe& inc_he2 : hg.hyperedgesOf(u)){
							//if (hg.residualCapacity(inc_he2.e) + hg.absoluteFlowReceived(inc_he2) > 0) {
								assert(level[u] <= level[inc_he2.e + hg.numNodes()] + 1);	// how is this assert ever gonna hold ?
								min_level = std::min(min_level, level[inc_he2.e + hg.numNodes()]);
							//}
						}
						assert(min_level >= level[u]);
						LOGGER << "Relabel node" << V(u) << V(level[u]) << V(min_level + 1);
						level[u] = min_level + 1;
						current_hyperedge[u] = hg.beginIndexHyperedges(u);
					}
				}
			}
		}

		void dischargeHyperedge(Hyperedge e, Node e_node) {
			assert(!hg.pinsOf(e).empty());
			assert(level[e_node] < max_level);

			while (excess[e_node] > 0) {
				Pin& pin = hg.getPin(current_pin[e]);
				Flow residual = std::min(hg.residualCapacity(e) + hg.absoluteFlowSent(pin), excess[e_node]);
				if (level[pin.pin] + 1 == level[e_node] && residual > 0) {
					LOGGER << "push to node" << V(pin.pin) << V(hg.getInHe(pin).flow) << V(e) << V(excess[e_node]) << V(residual);
					if constexpr (!relabel_to_front) {
						if (pin.pin != target && excess[pin.pin] == 0) {
							active_vertices_and_edges.push_back(pin.pin);
						}
					}

					hg.flow(e) += std::max(0, residual - hg.absoluteFlowSent(pin));		// drain flow sent first, then push the rest
					hg.getInHe(pin).flow -= residual;
					excess[e_node] -= residual;
					excess[pin.pin] += residual;

					assert(hg.flow(e) == out_flow(e_node));
					assert(excess[e_node] == in_flow(e_node) - out_flow(e_node));
				} else {
					if (++current_pin[e] == hg.endIndexPins(e)) {
						// relabel
						int min_level = std::numeric_limits<int>::max();
						for (const Pin& pin2 : hg.pinsOf(e)) {
							if (hg.residualCapacity(e) + hg.absoluteFlowSent(pin2) > 0) {
								assert(level[e_node] <= level[pin2.pin] + 1);
								min_level = std::min(min_level, level[pin2.pin]);
							}
						}
						if (min_level == std::numeric_limits<int>::max()) {
							LOGGER << V(e) << V(hg.flow(e)) << V(in_flow(e_node)) << V(out_flow(e_node)) << V(excess[e_node]) << V(hg.capacity(e));
							for (Pin& p2 : hg.pinsOf(e)) {
								LOGGER << V(p2.pin) << V(hg.getInHe(p2).flow) << V(level[p2.pin]);
							}
							std::abort();
						}
						assert(min_level >= level[e_node]);
						LOGGER << "Relabel hyperedge" << V(e) << V(level[e_node]) << V(min_level + 1);
						level[e_node] = min_level + 1;
						current_pin[e] = hg.beginIndexPins(e);
					}
				}
			}
		}

		void prepare(CutterState<Type>& cs) {
			previous_cutter_state_direction = cs.currentViewDirection();
			if (previous_cutter_state_direction != 0) {
				cs.flipViewDirection();
			}
			assert(cs.currentViewDirection() == 0);
		}

		void finish(CutterState<Type>& cs) {

		}

		std::pair<Flow,Flow> local_flow(const Node u) const {
			Flow in_flow = 0, out_flow = 0;
			if (u < hg.numNodes()) {
				for (InHe& inc_he : hg.hyperedgesOf(u)) {
					if (inc_he.flow > 0) {
						out_flow += inc_he.flow;
					} else {
						in_flow -= inc_he.flow;
					}
				}
			} else {
				for (Pin& pin : hg.pinsOf(Hyperedge(u - hg.numNodes()))) {
					const Flow f = hg.getInHe(pin).flow;
					if (f > 0) {
						in_flow += f;
					} else {
						out_flow -= f;
					}
				}
			}
			return std::make_pair(in_flow, out_flow);
		}

		Flow in_flow(const Node u) const {
			return local_flow(u).first;
		}

		Flow out_flow(const Node u) const {
			return local_flow(u).second;
		}

		TimeReporter timer;
	};
}