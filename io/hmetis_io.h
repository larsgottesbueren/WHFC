#pragma once

#include <fstream>
#include "../datastructure/flow_hypergraph.h"
#include "../datastructure/flow_hypergraph_builder.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

namespace whfc {
	class HMetisIO {
	private:
		inline static void mgetline(std::ifstream& f, std::string& line) {
			std::getline(f, line);
			while (line[0] == '%') {
				std::getline(f,line);
			}
		}

		static int open_file(const std::string& filename) {
			int fd = open(filename.c_str(), O_RDONLY);
			if ( fd == -1 ) {
				throw std::runtime_error("Could not open:" + filename);
			}
			return fd;
		}

		static size_t file_size(int fd) {
			struct stat file_info;
			if ( fstat(fd, &file_info) == -1 ) {
				throw std::runtime_error("Error while getting file stats");
			}
			return static_cast<size_t>(file_info.st_size);
		}

		static char* mmap_file(int fd, const size_t length) {
			char* mapped_file = (char*) mmap(0, length, PROT_READ, MAP_SHARED, fd, 0);
			if ( mapped_file == MAP_FAILED ) {
				close(fd);
				throw std::runtime_error("Error while mapping file to memory");
			}
			return mapped_file;
		}

		static void munmap_file(char* mapped_file, int fd, const size_t length) {
			if ( munmap(mapped_file, length) == -1 ) {
				close(fd);
				throw std::runtime_error("Error while unmapping file from memory");
			}
		}

		static inline void goto_next_line(char* mapped_file, size_t& pos, const size_t length) {
			for ( ; ; ++pos ) {
				if ( pos == length || mapped_file[pos] == '\n' ) {
					++pos;
					break;
				}
			}
		}

		static void skip_comments(char* mapped_file, size_t& pos, size_t length) {
			while ( mapped_file[pos] == '%' ) {
				goto_next_line(mapped_file, pos, length);
				assert(pos < length);
			}
		}

		static inline int64_t read_number(char* mapped_file, size_t& pos, const size_t length) {
			int64_t number = 0;
			for ( ; pos < length; ++pos ) {
				if ( mapped_file[pos] == ' ' || mapped_file[pos] == '\n' ) {
					while ( mapped_file[pos] == ' ' || mapped_file[pos] == '\n' ) {
						++pos;
					}
					break;
				}
				assert(mapped_file[pos] >= '0' && mapped_file[pos] <= '9');
				number = number * 10 + (mapped_file[pos] - '0');
			}
			return number;
		}


		static auto readHeader(char* mapped_file, size_t& pos, const size_t length) {
			skip_comments(mapped_file, pos, length);
			int64_t num_edges = read_number(mapped_file, pos, length);
			int64_t num_nodes = read_number(mapped_file, pos, length);
			HGType hg_type = HGType::Unweighted;
			if ( mapped_file[pos - 1] != '\n' ) {
				hg_type = static_cast<HGType>(read_number(mapped_file, pos, length));
			}
			//std::cout << std::endl;
			assert(mapped_file[pos - 1] == '\n');
			return std::make_tuple(num_nodes, num_edges, hg_type);
		}
		
	public:
		
		enum class HGType : uint8_t {
			Unweighted = 0,
			EdgeWeights = 1,
			NodeWeights = 10,
			EdgeAndNodeWeights = 11,
		};
		
		
		static auto readHeader(std::ifstream& f) {
			std::string line;
			size_t numHEs, numNodes;
			HGType hg_type = HGType::Unweighted;
			{
				//read header
				mgetline(f,line);
				std::istringstream iss(line);
				iss >> numHEs >> numNodes;
				uint32_t type = 0;
				if (iss >> type) {
					hg_type = static_cast<HGType>(type);
				}
			}
			return std::make_tuple(numNodes, numHEs, hg_type);
		}
		
		static FlowHypergraphBuilder readFlowHypergraphWithBuilder(const std::string& filename) {
			FlowHypergraphBuilder hgb;
			return readFlowHypergraphWithBuilder(hgb, filename);
		}
		
		static FlowHypergraphBuilder& readFlowHypergraphWithBuilder(FlowHypergraphBuilder& hgb, const std::string& filename) {
			std::ifstream f(filename);
			if (!f)
				throw std::runtime_error("File: " + filename + " not found.");
			
			auto [numNodes, numHEs, hg_type] = readHeader(f);
			hgb.clear();
			hgb.reinitialize(numNodes);
			
			std::string line;
			
			bool hasHyperedgeWeights = hg_type == HGType::EdgeAndNodeWeights || hg_type == HGType ::EdgeWeights;
			bool hasNodeWeights = hg_type == HGType::EdgeAndNodeWeights || hg_type == HGType::NodeWeights;


			for (size_t e = 0; e < numHEs; ++e) {
				mgetline(f, line);
				std::istringstream iss(line);
				uint32_t pin;
				uint32_t he_weight = 1;

				if (hasHyperedgeWeights)
					iss >> he_weight;
				
				hgb.startHyperedge(he_weight);
				size_t he_size = 0;
				while (iss >> pin) {
					if (pin < 1)
						throw std::runtime_error("File: " + filename + " has pin id < 1 (in one-based ids).");
					if (pin > numNodes)
						throw std::runtime_error("File: " + filename + " has pin id > number of nodes.");
					hgb.addPin(Node(pin-1));
					he_size++;
				}
				if (he_size <= 1)
					throw std::runtime_error("File: " + filename + " has pin with zero or one pins.");
			}

			for (Node u(0); u < numNodes; ++u) {
				NodeWeight nw(1);
				if (hasNodeWeights) {
					mgetline(f, line);
					std::istringstream iss(line);
					iss >> nw;
				}
				hgb.nodeWeight(u) = nw;
			}
			
			hgb.finalize();
			
			f.close();
			return hgb;
		}

		static FlowHypergraph readFlowHypergraph(const std::string& filename) {
			int fd = open_file(filename);
			const size_t length = file_size(fd);
			char* mapped_file = mmap_file(fd, length);
			size_t pos = 0;

			auto [num_nodes, num_edges, hg_type] = readHeader(mapped_file, pos, length);
			bool has_edge_weights = hg_type == HGType::EdgeAndNodeWeights || hg_type == HGType ::EdgeWeights;
			bool has_node_weights = hg_type == HGType::EdgeAndNodeWeights || hg_type == HGType::NodeWeights;

			std::vector<NodeWeight> node_weights;
			std::vector<HyperedgeWeight> edge_weights;
			std::vector<PinIndex> edge_sizes;
			std::vector<Node> pins;

			edge_weights.reserve(num_edges);
			edge_sizes.reserve(num_edges);
			for (int64_t e = 0; e < num_edges; ++e) {
				skip_comments(mapped_file, pos, length);

				uint32_t he_weight = 1;
				if (has_edge_weights) {
					he_weight = read_number(mapped_file, pos, length);
				}
				edge_weights.emplace_back(he_weight);

				uint32_t he_size = 0, pin;
				do {
					pin = read_number(mapped_file, pos, length);
					he_size++;
					pins.emplace_back(pin - 1);
					if (pin == 0) { throw std::runtime_error("read number didn't give a number"); }
				} while (mapped_file[pos-1] != '\n');

				edge_sizes.emplace_back(he_size);
				if (he_size <= 1) { throw std::runtime_error("File: " + filename + " has pin with zero or one pins."); }
			}

			if (has_node_weights) {
				node_weights.reserve(num_nodes);
				for (int64_t u = 0; u < num_nodes; ++u) {
					node_weights.emplace_back( read_number(mapped_file, pos, length) );
				}
			} else {
				node_weights.resize(num_nodes, 1);
			}

			skip_comments(mapped_file, pos, length);
			assert(pos == length);

			return FlowHypergraph(node_weights, edge_weights, edge_sizes, pins);
		}


		static FlowHypergraph readFlowHypergraphOld(const std::string& filename) {

			std::vector<NodeWeight> nodeWeights;
			std::vector<HyperedgeWeight> hyperedgeWeights;
			std::vector<Node> pins;
			std::vector<PinIndex> hyperedgeSizes;

			std::ifstream f(filename);
			if (!f)
				throw std::runtime_error("File: " + filename + " not found.");

			auto [numNodes, numHEs, hg_type] = readHeader(f);
			
			std::string line;
			
			bool hasHyperedgeWeights = hg_type == HGType::EdgeAndNodeWeights || hg_type == HGType ::EdgeWeights;
			
			bool hasNodeWeights = hg_type == HGType::EdgeAndNodeWeights || hg_type == HGType::NodeWeights;
			if (!hasNodeWeights)
				nodeWeights.resize(numNodes, NodeWeight(1));

			for (size_t e = 0; e < numHEs; ++e) {
				mgetline(f, line);
				std::istringstream iss(line);
				uint32_t pin;
				uint32_t he_size = 0;
				uint32_t he_weight = 1;
				
				if (hasHyperedgeWeights)
					iss >> he_weight;
				hyperedgeWeights.emplace_back(he_weight);
				
				while (iss >> pin) {
					if (pin < 1)
						throw std::runtime_error("File: " + filename + " has pin id < 1 (in one-based ids).");
					if (pin > numNodes)
						throw std::runtime_error("File: " + filename + " has pin id > number of nodes.");
					he_size++;
					pins.emplace_back(pin-1);
				}
				hyperedgeSizes.emplace_back(he_size);
				
				if (he_size > numNodes)
					throw std::runtime_error("File: " + filename + " has hyperedge with more pins than nodes in the hypergraph.");
				if (he_size == 0)
					throw std::runtime_error("File: " + filename + " has hyperedge with zero pins.");
				if (he_size == 1) {
					//ignore single pin hyperedges
					pins.pop_back();
					hyperedgeWeights.pop_back();
					hyperedgeSizes.pop_back();
				}
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

		static void writeFlowHypergraph(FlowHypergraph& hg, std::string& filename) {
			if (filename.empty())
				throw std::runtime_error("No filename for Flow Hypergraph specified");
			std::ofstream f(filename);
			if (!f)
				throw std::runtime_error("Failed at creating Flow Hypergraph file " + filename);

			bool hasNodeWeights = hg.hasNodeWeights();
			bool hasHyperedgeWeights = hg.hasHyperedgeWeights();

			{
				//write header
				f << hg.numHyperedges() << " " << hg.numNodes();
				if (hasNodeWeights)
					if (hasHyperedgeWeights)
						f << " " << static_cast<uint32_t>(HGType::EdgeAndNodeWeights);
					else
						f << " " << static_cast<uint32_t>(HGType::NodeWeights);
				else if (hasHyperedgeWeights)
					f << " " << static_cast<uint32_t>(HGType::EdgeWeights);
				f << "\n";
			}

			for (Hyperedge e : hg.hyperedgeIDs()) {
				auto pinsOfE = hg.pinsOf(e);
				if (pinsOfE.empty())
					throw std::runtime_error("Hypergraph has hyperedge with zero pins");
				if (hasHyperedgeWeights)
					f << hg.capacity(e) << " ";
				
				f << (pinsOfE.begin()->pin + 1); pinsOfE.advance_begin();	//special case first pin since we have |e|-1 white spaces
				for (const FlowHypergraph::Pin& p : pinsOfE)
					f << " " << (p.pin + 1);		//yes... hMetis insists on 1-based IDs -.-
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
