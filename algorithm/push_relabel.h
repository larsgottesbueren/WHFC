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
			level.assign(hg.numNodes(), 0);		// still TODO global relabeling instead of this
			level[source] = max_level;

			for (InHe& in_he : hg.hyperedgesOf(source)) {
				pushToHyperedge(source, Node(in_he.e + hg.numNodes()), in_he, hg.capacity(in_he.e));
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
				while (excess[target] <= upperFlowBound && !active_vertices_and_edges.empty()) {
					const Node x = active_vertices_and_edges.front();
					active_vertices_and_edges.pop_front();
					discharge(x);
				}
				timer.stop("Discharge");
			}


			active_vertices_and_edges.clear();

			finish(cs);
			timer.report(std::cout);

			Flow flow_delta = excess[target] - old_excess;
			LOGGER << V(flow_delta);
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
			//LOGGER << "push to hyperedge" << V(u) << V(e_node) << V(residual);
			if constexpr (!relabel_to_front) {
				if (excess[e_node] == 0) {
					active_vertices_and_edges.push_back(e_node);
				}
			}
			// push p4 and p2
			residual = std::min(residual, excess[u]);
			excess[e_node] += residual;
			excess[u] -= residual;
			inc_he.flow += residual;
			// update flow value of hyperedge due to p4. First consume p2 completely, then p4
			hg.flow(inc_he.e) += std::max(0, residual - hg.absoluteFlowReceived(inc_he));	// TODO actually hg.flowSent(max(0, res - flow_rec(inc_he)) but multiplier is 1 for now?

			assert(hg.flow(inc_he.e) == in_flow(e_node));
			assert(excess[e_node] == in_flow(e_node) - out_flow(e_node));	// these break running time guarantees in debug mode
		}

		void dischargeNode(Node u) {
			assert(hg.degree(u) > 0);
			assert(level[u] < max_level);

			const Hyperedge num_nodes(hg.numNodes());	// TaggedInteger conversion...
			while (excess[u] > 0) {
				InHe& inc_he = hg.getInHe(current_hyperedge[u]);
				const Node e_node(inc_he.e + num_nodes);
				Flow residual = hg.residualCapacity(inc_he.e) + hg.absoluteFlowReceived(inc_he);
				if (residual > 0 && level[e_node] + 1 == level[u]) {
					pushToHyperedge(u, e_node, inc_he, residual);
				} else {
					// don't advance iterator if pushed
					if (++current_hyperedge[u] == hg.endIndexHyperedges(u)) {
						// relabel
						int min_level = std::numeric_limits<int>::max();
						for (const InHe& inc_he2 : hg.hyperedgesOf(u)){
							if (hg.residualCapacity(inc_he2.e) + hg.absoluteFlowReceived(inc_he2) > 0) {
								min_level = std::min(min_level, level[inc_he2.e + num_nodes]);
							}
						}
						level[u] = min_level + 1;
						current_hyperedge[u] = hg.beginIndexHyperedges(u);
					}
				}
			}
		}

		void dischargeHyperedge(Hyperedge e, Node e_node) {
			assert(!hg.pinsOf(e).empty());
			assert(level[e_node] < max_level);
			// if hyperedge has excess we always hit the first pin with appropriate level. --> put it right into push for node?
			// is it so smart to always push everything out? just leads to lot of relabels to get it back. find out if anyone ever pushes less than possible residual capacity.

			while (excess[e_node] > 0) {
				Pin& pin = hg.getPin(current_pin[e]);
				if (level[pin.pin] + 1 == level[e_node]) {
					LOGGER << "push to node" << V(pin.pin) << V(e) << V(excess[e_node]);
					if constexpr (!relabel_to_front) {
						if (excess[pin.pin] == 0) {
							active_vertices_and_edges.push_back(e_node);
						}
					}

					hg.flow(e) -= hg.absoluteFlowSent(pin);		// first push back on p1
					hg.getInHe(pin).flow -= excess[e_node];		// then dump it all via p2/p4
					excess[pin.pin] += excess[e_node];
					excess[e_node] = 0;
				} else {
					if (++current_pin[e] == hg.endIndexPins(e)) {
						// relabel
						int min_level = std::numeric_limits<int>::max();
						for (const Pin& pin2 : hg.pinsOf(e)) {
							min_level = std::min(min_level, level[pin2.pin]);
						}
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
						out_flow += f;
					} else {
						in_flow -= f;
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