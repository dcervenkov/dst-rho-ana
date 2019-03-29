/**
 *  @file    efficiency.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-01-22
 *
 *  @brief Class that calculates efficiency (detector acceptance) from a model for a set of
 * parameters
 *
 */

#include "efficiency.h"

// Boost includes
#include <boost/algorithm/string/predicate.hpp>

// ROOT includes
#include "TFile.h"

// Local includes
#include "log.h"

Efficiency::Efficiency() {
}

Efficiency::Efficiency(const char* filename) {
    ReadInFile(filename);
}

void Efficiency::ReadInFile(const char* filename) {
    if (boost::algorithm::ends_with(filename, ".root")) {
        TFile efficiency_file(filename, "read");
        histo_efficiency = (TH3F*)efficiency_file.Get("eff_histo");
        histo_efficiency->SetDirectory(0);
        efficiency_file.Close();
    } else {
        binned_efficiency = new BinnedDensity(filename, &phasespace, filename);
    }
    Log::print(Log::info, "Reading in efficiency file '%s'\n", filename);
    if (histo_efficiency && binned_efficiency) {
        histo_normalization = GetNormalization();
    }
}

Efficiency::~Efficiency() {
}

/**
 * Return the actual efficiency value for the supplied parameters
 */
double Efficiency::GetEfficiency(double thetat, double thetab, double phit,
                                 int efficiency_model) const {
    thetat_->setVal(thetat);
    thetab_->setVal(thetab);
    phit_->setVal(phit);

    switch (efficiency_model) {
        case 0:
            // Uniform efficiency for generator level fits
            return 1;
        case 1:
            return GetModel1Efficiency();
        case 2:
            return GetModel2Efficiency();
        case 3:
            return GetModel3Efficiency();
        case 4:
            return GetModel4Efficiency();
        case 5: {
            assert(binned_efficiency != nullptr);
            return GetKDEEfficiency(thetat, thetab, phit);
        }
        case 6: {
            assert(histo_efficiency != nullptr);
            return GetHistogramEfficiency(thetat, thetab, phit);
        }
        case 7: {
            assert(binned_efficiency != nullptr);
            assert(histo_efficiency != nullptr);
            std::vector<Double_t> coords(3);
            coords[0] = thetat;
            coords[1] = thetab;
            coords[2] = phit;
            if (CloseToEdge(coords, 0.05)) {
                return GetHistogramEfficiency(thetat, thetab, phit) * histo_normalization;
            } else {
                return GetKDEEfficiency(thetat, thetab, phit);
            }
        }

    }
    return 0;
}

/**
 * Check whether any of the vars is too close to the histogram edge to use
 * interpolation.
 * 
 * From the ROOT documentation: The given values (x,y,z) must be between first
 * bin center and last bin center for each coordinate.
 */
bool Efficiency::CanUseInterpolation(const double& phit, const double& transtht,
                               const double& transthb) const {
    double vars[3] = {phit, transtht, transthb};
    TAxis* axes[3] = {histo_efficiency->GetXaxis(), histo_efficiency->GetYaxis(),
                      histo_efficiency->GetZaxis()};

    for (int i = 0; i < 3; i++) {
        int last_bin = axes[i]->GetNbins();
        double low_center = axes[i]->GetBinCenter(1);
        double high_center = axes[i]->GetBinCenter(last_bin);
        if (vars[i] < low_center || vars[i] > high_center) {
            return false;
        }
    }
    return true;
}

double Efficiency::EfficiencyInterface(double* x, double* p) const {
	return GetEfficiency(x[0], x[1], x[2], 5);
}

/**
 * Rescale variables to account for margin-mirrored phase space from DSRhoEfficiency.
 */
// void Efficiency::RescaleVars(double& thetat, double& thetab, double& phit, const double margin) const {
// 	const double vars[3] = {thetat, thetab, phit};
// const double min[3] = {constants::cuts::thetat_low, constants::cuts::phit_low, constants::cuts::phit_low};
// const double max[3] = {constants::cuts::thetat_high, constants::cuts::phit_high, constants::cuts::phit_high};

// 	double center[3];
// 	double new_vars[3];

// 	for (int i = 0; i < 3; i++) {
// 		center[i] = (min[i] + max[i]) / 2;
// 		new_vars[i] = center[i] + (vars[i] - center[i])/(1 + 2 * margin);
// 	}

// 	thetat = new_vars[0];
// 	thetab = new_vars[1];
// 	phit = new_vars[2];
// }

int Efficiency::CloseToEdge(const std::vector<Double_t> vals, const double margin) const {
    RooRealVar* vars[3] = {thetat_, thetab_, phit_};
    for (int var_num = 0; var_num < 3; var_num++) {
        const double min = vars[var_num]->getMin();
        const double max = vars[var_num]->getMax();
        const double range = max - min;
        if (vals[var_num] < min + range * margin) {
            return 1;
        } else if (vals[var_num] > max - range * margin) {
            return 2;
        }
    }
    return 0;
}

double Efficiency::GetHistogramEfficiency(double thetat, double thetab, double phit) const {
    double eff = 0;
    // double transtht = thetat / TMath::Pi() * 2 - 1;
    // double transthb = thetab / TMath::Pi() * 2 - 1;
    double transtht = thetat;
    double transthb = thetab;
    if (CanUseInterpolation(transtht, transthb, phit)) {
        eff = histo_efficiency->Interpolate(transtht, transthb, phit);
    } else {
        int binx = histo_efficiency->GetXaxis()->FindBin(transtht);
        int biny = histo_efficiency->GetYaxis()->FindBin(transthb);
        int binz = histo_efficiency->GetZaxis()->FindBin(phit);
        int bin = histo_efficiency->GetBin(binx, biny, binz);
        eff = histo_efficiency->GetBinContent(bin);
    }
    // Add a tiny number to avoid PDF being 0, which causes trouble
    if (eff == 0) eff += 0.000001;
    return eff;
}

double Efficiency::GetKDEEfficiency(double thetat, double thetab, double phit) const {
    std::vector<Double_t> coords(3);
    coords[0] = thetat;
    coords[1] = thetab;
    coords[2] = phit;
    double eff = binned_efficiency->density(coords);
    return eff;
}

double Efficiency::GetNormalization() {
    TH3F* kde_efficiency = dynamic_cast<TH3F*>(histo_efficiency->Clone());
    kde_efficiency->Reset();
    assert(kde_efficiency->GetSumOfWeights() == 0);

    double thetat, thetab, phit;
    for (int x = 1; x <= histo_efficiency->GetXaxis()->GetNbins(); x++) {
        for (int y = 1; y <= kde_efficiency->GetYaxis()->GetNbins(); y++) {
            for (int z = 1; z <= kde_efficiency->GetZaxis()->GetNbins(); z++) {
                thetat = kde_efficiency->GetXaxis()->GetBinCenter(x);
                thetab = kde_efficiency->GetYaxis()->GetBinCenter(y);
                phit = kde_efficiency->GetZaxis()->GetBinCenter(z);
                kde_efficiency->SetBinContent(kde_efficiency->GetBin(x, y, z),
                                     GetKDEEfficiency(thetat, thetab, phit));
            }
        }
    }

    return kde_efficiency->GetSumOfWeights() / histo_efficiency->GetSumOfWeights();
}