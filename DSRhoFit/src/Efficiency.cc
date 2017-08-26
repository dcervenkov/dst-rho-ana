/*
 * Efficiency.cc
 *
 *  Created on: Jan 22, 2016
 *      Author: cervenkov
 */

#include "Efficiency.h"

Efficiency::Efficiency() {
	// TODO Auto-generated constructor stub
}

Efficiency::~Efficiency() {
	// TODO Auto-generated destructor stub
}

double Efficiency::GetEfficiency(double thetat, double thetab, double phit) const {
	thetat_->setVal(thetat);
	thetab_->setVal(thetab);
	phit_->setVal(phit);

	return GetModel1Efficiency();
//	return GetModel2Efficiency();
//	return GetModel3Efficiency();
//	return GetModel4Efficiency();
}
