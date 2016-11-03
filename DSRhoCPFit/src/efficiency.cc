/**
 *  @file    efficiency.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-01-22
 *
 *  @brief Class that calculates efficiency (detector acceptance) from a model for a set of parameters
 *
 */

#include "efficiency.h"

Efficiency::Efficiency() {
	// TODO Auto-generated constructor stub
}

Efficiency::~Efficiency() {
	// TODO Auto-generated destructor stub
}

/**
 * Return the actual efficiency value for the supplied parameters
 */
double Efficiency::GetEfficiency(double thetat, double thetab, double phit) const {
	thetat_->setVal(thetat);
	thetab_->setVal(thetab);
	phit_->setVal(phit);

//	return GetModel1Efficiency();
//	return GetModel2Efficiency();
	return GetModel3Efficiency();
}
