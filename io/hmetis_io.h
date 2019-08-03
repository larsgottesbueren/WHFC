#pragma once

#include <fstream>
#include "../datastructure/flow_hypergraph.h"

namespace whfc {
	class HMetisIO {
		inline void mgetline(std::ifstream& f, std::string& line) {
			std::getline(f, line);
			while (line[0] == '%') {
				std::getline(f,line);
			}
		}

		using FlowHypergraph::HyperedgeData;
		using FlowHypergraph::InHe;
		using FlowHypergraph::NodeData;
		using FlowHypergraph::Pin;

		enum class HGType : uint8_t {
			Unweighted = 0,
			EdgeWeights = 1,
			NodeWeights = 10,
			EdgeAndNodeWeights = 11,
		};

		FlowHypergraph readFlowHypergraph(std::string& filename) {

			std::vector<NodeWeight> nodeWeights;
			std::vector<HyperedgeWeight> hyperedgeWeights;
			std::vector<Node> pins;
			std::vector<PinIndex> hyperedgeSizes;
			size_t num_hes, numNodes;
			HGType hg_type = HGType::Unweighted;

			std::ifstream f(filename);
			if (!f)
				throw std::runtime_error("File: " + filename + " not found.");

			std::string line;
			{
				//read header
				mgetline(f,line);
				std::istringstream iss(line);
				iss >> num_hes >> numNodes;
				uint8_t type;
				if (iss >> type)
					hg_type = static_cast<HGType>(type);
			}

			bool hasHyperedgeWeights = hg_type == HGType::EdgeAndNodeWeights || hg_type == HGType ::EdgeWeights;
			if (!hasHyperedgeWeights)
				hyperedgeWeights.resize(num_hes, HyperedgeWeight(1));
			bool hasNodeWeights = hg_type == HGType::EdgeAndNodeWeights || hg_type == HGType::NodeWeights;
			if (!hasNodeWeights)
				nodeWeights.resize(numNodes, NodeWeight(1));

			for (size_t e = 0; e < num_hes; ++e) {
				mgetline(f, line);
				std::istringstream iss(line);
				uint32_t pin;
				uint32_t he_size = 0;
				if (hasHyperedgeWeights) {
					uint32_t we;
					iss >> we;
					hyperedgeWeights.emplace_back(we);
				}
				while (iss >> pin) {
					if (pin < 1)
						throw std::runtime_error("File: " + filename + " has pin id < 1 (in one-based ids).");
					if (pin > numNodes)
						throw std::runtime_error("File: " + filename + " has pin id > number of nodes.");
					he_size++;
					pins.emplace_back(pin-1);
				}
				if (he_size > numNodes)
					throw std::runtime_error("File: " + filename + " has hyperedge with more pins than nodes in the hypergraph.");
				if (he_size == 0)
					throw std::runtime_error("File: " + filename + " has hyperedge with zero pins.");

				if (he_size == 1) //ignore single pin hyperedges
					pins.pop_back();
				else
					hyperedgeSizes.emplace_back(he_size);
			}

			if (hasNodeWeights) {
				for (size_t u = 0; u < numNodes; ++u) {
					uint32_t nw;
					mgetline(f, line);
					std::istringstream iss(line);
					iss >> nw;
					nodeWeights.emplace_back(nw);
				}
			}

			f.close();

			return FlowHypergraph(nodeWeights, hyperedgeWeights, hyperedgeSizes, pins);
		}

		void writeFlowHypergraph(FlowHypergraph& hg, std::string& filename) {
			if (filename.empty())
				throw std::runtime_error("No filename for Flow Hypergraph specified");
			std::ofstream f(filename);
			if (!f)
				throw std::runtime_error("Failed at creating Flow Hypergraph file " + filename);

			bool hasNodeWeights = hg.hasNodeWeights();
			bool hasHyperedgeWeights = hg.hasHyperedgeWeights();

			{
				//write header
				f << hg.numHyperedges() << " " << hg.numNodes() << " ";
				if (hasNodeWeights)
					if (hasHyperedgeWeights)
						f << static_cast<uint8_t>(HGType::EdgeAndNodeWeights);
					else
						f << static_cast<uint8_t>(HGType::NodeWeights);
				else if (hasHyperedgeWeights)
					f << static_cast<uint8_t>(HGType::EdgeWeights);
				f << "\n";
			}

			for (Hyperedge e : hg.hyperedgeIDs()) {
				if (hasHyperedgeWeights)
					f << hg.capacity(e) << " ";
				for (const FlowHypergraph::Pin& p : hg.pinsOf(e))
					f << (p.pin + 1) << " ";		//yes... hMetis insists on 1-based IDs -.-
				f << "\n";
			}

			if (hasNodeWeights)
				for (Node u : hg.nodeIDs())
					f << hg.nodeWeight(u) << "\n";


			f << std::flush;
			f.close();
		}

	};
}
