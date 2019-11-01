#pragma once

#include <fstream>
#include "../util/random.h"
#include "../definitions.h"

namespace whfc {
	class WHFC_IO {
	public:
		struct WHFCInformation {
			NodeWeight maxBlockWeight;
			Flow upperFlowBound;
			Node s, t;
		};
		
		static void readRandomGeneratorState(const std::string& hgpath) {
			std::ifstream df(hgpath + ".distribution");
			if (df) {
				df >> Random::instance().get64BitUintDistribution();
			}
			df.close();
			
			std::ifstream genf(hgpath + ".generator");
			if (genf) {
				genf >> Random::instance().getGenerator();
			}
			genf.close();
		}
		
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
			
			std::ofstream df(hgpath + ".distribution");
			df << Random::instance().get64BitUintDistribution();
			df.close();
			
			std::ofstream genf(hgpath + ".generator");
			genf << Random::instance().getGenerator();
			genf.close();
		}
	};
}