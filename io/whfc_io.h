#pragma once

#include <fstream>
#include "../definitions.h"

namespace whfc {
	class WHFC_IO {
		static constexpr std::string fileSuffix = ".whfc";
		struct WHFCInformation {
			NodeWeight maxBlockWeight;
			Flow currentCutSize;
			Node s, t;
		};
	public:
		static WHFCInformation readAdditionalInformation(std::string& hgpath) {
			std::ifstream f(hgpath + fileSuffix);
			WHFCInformation i;
			f >> i.maxBlockWeight
			  >> i.currentCutSize
			  >> i.s
			  >> i.t;
			f.close();
			return i;
		}

		static void writeAdditionalInformation(std::string& hgpath, WHFCInformation& i) {
			std::ofstream f(hgpath + fileSuffix);
			f << i.maxBlockWeight << " " << i.currentCutSize << " " << i.s << " " << i.t << std::endl;
			f.close();
		}
	};
}