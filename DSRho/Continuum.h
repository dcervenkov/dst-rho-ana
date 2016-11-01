/*
 * Continuum.h
 *
 *  Created on: May 25, 2015
 *      Author: cervenkov
 */

#ifndef CONTINUUM_H_
#define CONTINUUM_H_

#include "TMVA/Reader.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class ksfwmoments;

class Continuum {
public:
	Continuum();
	virtual ~Continuum();

	void readInWeights(const char* channel, const char* svd);
	void readInVars(const double& mbc, const double& cosThetaThrust,
			const double& cosThetaB, ksfwmoments& km);
	double evaluateBDTG() {return reader->EvaluateMVA("BDTG method");}
	double evaluateMLP() {return reader->EvaluateMVA("MLPBNN method");}

private:
	TMVA::Reader* reader;
	float vars[21];
};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif /* CONTINUUM_H_ */
