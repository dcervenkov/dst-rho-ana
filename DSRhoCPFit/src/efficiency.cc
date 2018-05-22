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
	binned_efficiency = new BinnedDensity("binned_efficiency", &phasespace, "efficiency");
}

Efficiency::~Efficiency() {
	// TODO Auto-generated destructor stub
}

/**
 * Return the actual efficiency value for the supplied parameters
 */
double Efficiency::GetEfficiency(double thetat, double thetab, double phit, int efficiency_model) const {
	thetat_->setVal(thetat);
	thetab_->setVal(thetab);
	phit_->setVal(phit);

	switch (efficiency_model) {
		case 1:
			return GetModel1Efficiency();
		case 2:
			return GetModel2Efficiency();
		case 3:
			return GetModel3Efficiency();
		case 4:
			return GetModel4Efficiency();
		case 5:
			// RescaleVars(thetat, thetab, phit, 0.1);
			std::vector<Double_t> coords(3);
			coords[0] = thetat;
			coords[1] = thetab;
			coords[2] = phit;
			// printf("EFFDBG: model4 = %f, model5 = %f\n", GetModel4Efficiency(), binned_efficiency->density(coords));
			double eff = binned_efficiency->density(coords);
			if (eff == 0) {
				eff = GetModel4Efficiency();
			}
			return eff;
	}
	return 0;
}

// double Efficiency::EfficiencyInterface(double* vars, double* pars) const {
double Efficiency::EfficiencyInterface(double* x, double* p) const {
	GetEfficiency(x[0], x[1], x[2], 5);
}

/**
 * Rescale variables to account for margin-mirrored phase space from DSRhoEfficiency.
 */
void Efficiency::RescaleVars(double& thetat, double& thetab, double& phit, const double margin) const {
	const double vars[3] = {thetat, thetab, phit};
	const double min[3] = {0, 0.6, -constants::pi};
	const double max[3] = {constants::pi, 2.95, constants::pi};

	double center[3];
	double new_vars[3];

	for (int i = 0; i < 3; i++) {
		center[i] = (min[i] + max[i]) / 2;
		new_vars[i] = center[i] + (vars[i] - center[i])/(1 + 2 * margin);
	}

	thetat = new_vars[0];
	thetab = new_vars[1];
	phit = new_vars[2];
}
