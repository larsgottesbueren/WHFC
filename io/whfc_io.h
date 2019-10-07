#pragma once

#include <fstream>
#include "../definitions.h"

namespace whfc {
	class WHFC_IO {
	public:
		struct WHFCInformation {
			NodeWeight maxBlockWeight;
			Flow upperFlowBound;
			Node s, t;
		};
		static WHFCInformation readAdditionalInformation(const std::string& hgpath) {
			std::string fileSuffix = ".whfc";
			std::ifstream f(hgpath + fileSuffix);
			WHFCInformation i;
			f >> i.maxBlockWeight
			  >> i.upperFlowBound
			  >> i.s
			  >> i.t;
			f.close();
			return i;
		}

		static void writeAdditionalInformation(std::string& hgpath, WHFCInformation& i) {
			std::string fileSuffix = ".whfc";
			std::ofstream f(hgpath + fileSuffix);
			f << i.maxBlockWeight << " " << i.upperFlowBound << " " << i.s << " " << i.t << std::endl;
			f.close();
		}
	};
}