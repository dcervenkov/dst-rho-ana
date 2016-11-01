/*
 * Continuum.cpp
 *
 *  Created on: May 25, 2015
 *      Author: cervenkov
 */

#include <string>

// For some reason ksfwmoments.h can't be included without these two
#include "panther/panther.h"
#include "eid/eid.h"

#include "ksfwmoments.h"
#include "Continuum.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

Continuum::Continuum() {
	reader = new TMVA::Reader("!Color:!Silent");

	reader->AddVariable("mbc", &vars[0]);
	reader->AddVariable("shcosthb", &vars[1]);
	reader->AddVariable("shcostht", &vars[2]);
	reader->AddVariable("k1mm2", &vars[3]);
	reader->AddVariable("k1et", &vars[4]);
	reader->AddVariable("k1hoo0", &vars[16]);
	reader->AddVariable("k1hoo1", &vars[17]);
	reader->AddVariable("k1hoo2", &vars[18]);
	reader->AddVariable("k1hoo3", &vars[19]);
	reader->AddVariable("k1hoo4", &vars[20]);
	reader->AddVariable("k1hso00", &vars[5]);
	reader->AddVariable("k1hso01", &vars[6]);
	reader->AddVariable("k1hso02", &vars[7]);
	reader->AddVariable("k1hso03", &vars[8]);
	reader->AddVariable("k1hso04", &vars[9]);
	reader->AddVariable("k1hso10", &vars[10]);
	reader->AddVariable("k1hso12", &vars[11]);
	reader->AddVariable("k1hso14", &vars[12]);
	reader->AddVariable("k1hso20", &vars[13]);
	reader->AddVariable("k1hso22", &vars[14]);
	reader->AddVariable("k1hso24", &vars[15]);
}

Continuum::~Continuum() {
}

// Too lazy to do this properly :)
void Continuum::readInWeights(const char* channel, const char* svd) {
	const char* baseFilePath = "weights/TMVAClassification_";
	std::string weightFilePath;

	weightFilePath = baseFilePath;
	weightFilePath += channel;
	weightFilePath += "_";
	weightFilePath += "MLPBNN";
	weightFilePath += ".weights.xml";
	reader->BookMVA("MLPBNN method", weightFilePath.c_str());

	weightFilePath = baseFilePath;
	weightFilePath += channel;
	weightFilePath += "_";
	weightFilePath += "BDTG";
	weightFilePath += ".weights.xml";
	reader->BookMVA("BDTG method", weightFilePath.c_str());
}

void Continuum::readInVars(const double& mbc, const double& cosThetaThrust,
		const double& cosThetaB, ksfwmoments& km) {

	vars[0] = (float) mbc;
	// CS was trained with abs(cosThetaThrust)
	// See journal entry <deid:10205 CS and shcosthb>
	vars[1] = fabs((float) cosThetaThrust);
	vars[2] = (float) cosThetaB;

	km.usefinal(1);
	vars[3] = (float) km.mm2();
    vars[4] = (float) km.et();
    vars[5] = (float) km.Hso(0, 0);
    vars[6] = (float) km.Hso(0, 1);
    vars[7] = (float) km.Hso(0, 2);
    vars[8] = (float) km.Hso(0, 3);
    vars[9] = (float) km.Hso(0, 4);
    vars[10] = (float) km.Hso(1, 0);
    vars[11] = (float) km.Hso(1, 2);
    vars[12] = (float) km.Hso(1, 4);
    vars[13] = (float) km.Hso(2, 0);
    vars[14] = (float) km.Hso(2, 2);
    vars[15] = (float) km.Hso(2, 4);
    vars[16] = (float) km.Hoo(0);
    vars[17] = (float) km.Hoo(1);
    vars[18] = (float) km.Hoo(2);
    vars[19] = (float) km.Hoo(3);
    vars[20] = (float) km.Hoo(4);
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

