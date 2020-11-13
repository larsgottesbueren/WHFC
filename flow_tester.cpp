#include <queue>
#include "datastructure/flow_hypergraph.h"
#include "io/hmetis_io.h"
#include "io/whfc_io.h"
#include "datastructure/flow_hypergraph_builder.h"
#include "algorithm/dinic_bidirectional.h"
#include "algorithm/dinic.h"
#include "algorithm/dinic_scaling.h"


namespace whfc {
	
	template<typename FlowAlgorithm>
	int runSnapshotTester(const std::string& filename) {
		static constexpr bool log = true;
		TimeReporter timer;
		WHFC_IO::WHFCInformation info = WHFC_IO::readAdditionalInformation(filename);
		Node s = info.s;
		Node t = info.t;
		//LOGGER << s << t << info.maxBlockWeight[0] << info.maxBlockWeight[1] << info.upperFlowBound;
		FlowHypergraph hg = HMetisIO::readFlowHypergraph(filename);

		if (s >= hg.numNodes() || t >= hg.numNodes())
			throw std::runtime_error("s or t not within node id range");
		
		LOGGER << V(hg.numNodes()) << V(hg.numHyperedges()) << V(hg.numPins());

		FlowAlgorithm flow_algo(hg);
		CutterState<FlowAlgorithm> cs(hg, timer);
		for (int i = 0; i < 2; ++i)
			cs.setMaxBlockWeight(i, info.maxBlockWeight[i]);
		cs.initialize(s, t);

		flow_algo.exhaustFlow(cs);
		//LOGGER << V(cs.flowValue);

		#ifndef NDEBUG
		{
			// test if flow is maximal
			BitVector vis(hg.numNodes()), he_vis(hg.numHyperedges()), he_flow_sending_vis(hg.numHyperedges());
			std::queue<Node> q;
			
			size_t n_vis = 0;
			auto visit = [&](Node v) {
				assert(!cs.n.isTarget(v));
				if (!vis[v]) {
					vis.set(v);
					q.push(v);
					n_vis++;
				}
			};
			
			for (auto& sp : cs.sourcePiercingNodes)
				visit(sp.node);
			
			while (!q.empty()) {
				const Node u = q.front(); q.pop();
				for (InHeIndex inc_u_iter : hg.incidentHyperedgeIndices(u)) {
					const auto& inc_u = hg.getInHe(inc_u_iter);
					const Hyperedge e = inc_u.e;
					if (!he_vis[e]) {
						if (!hg.isSaturated(e) || hg.flowReceived(inc_u) > 0) {
							he_vis.set(e);
							for (auto& pv : hg.pinsOf(e)) {
								visit(pv.pin);
							}
						} else if(!he_flow_sending_vis[e]) {
							he_flow_sending_vis.set(e);
							for (auto& pv : hg.pinsOf(e)) {
								if (hg.flowSent(pv) > 0) {
									visit(pv.pin);
								}
							}
						}
					}
				}
			}
			assert(n_vis < hg.numNodes());
		}
		#endif
		

		return cs.flowValue;
	}
}

int main(int argc, const char* argv[]) {
	if (argc < 2 || argc > 3)
		throw std::runtime_error("Usage: ./FlowTester hypergraphfile");
	std::string hgfile = argv[1];
	std::cout << "Bidir Dinic" << std::endl;
	int f1 = whfc::runSnapshotTester<whfc::BidirectionalDinic>(hgfile); (void)(f1);
	std::cout << "Plain Dinic" << std::endl;
	int f2 = whfc::runSnapshotTester<whfc::Dinic>(hgfile); (void)(f2);
	//assert(f1 == f2);
	//int f3 = whfc::runSnapshotTester<whfc::ScalingDinic>(hgfile);
	//assert(f2 == f3);
	return 0;
}