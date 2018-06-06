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

// ROOT includes
#include "TFile.h"

Efficiency::Efficiency() {
    binned_efficiency = new BinnedDensity("binned_efficiency", &phasespace, "efficiency");

    TFile efficiency_file("efficiency.root", "read");
    histo_efficiency = (TH3F*)efficiency_file.Get("eff_histo");
    histo_efficiency->SetDirectory(0);
    efficiency_file.Close();
}

Efficiency::~Efficiency() {
    // TODO Auto-generated destructor stub
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
        case 1:
            return GetModel1Efficiency();
        case 2:
            return GetModel2Efficiency();
        case 3:
            return GetModel3Efficiency();
        case 4:
            return GetModel4Efficiency();
        case 5: {
            std::vector<Double_t> coords(3);
            coords[0] = thetat;
            coords[1] = thetab;
            coords[2] = phit;
            // printf("EFFDBG: model4 = %f, model5 = %f\n", GetModel4Efficiency(),
            // binned_efficiency->density(coords));
            double eff = binned_efficiency->density(coords);
            if (eff == 0) {
                eff = GetModel4Efficiency();
            }
            return eff;
        }
        case 6:
            double eff = 0;
            double transtht = thetat / TMath::Pi() * 2 - 1;
            double transthb = thetab / TMath::Pi() * 2 - 1;
            if (CanUseInterpolation(phit, transtht, transthb)) {
                eff = histo_efficiency->Interpolate(phit, transtht, transthb);
            } else {
                int binx = histo_efficiency->GetXaxis()->FindBin(phit);
                int biny = histo_efficiency->GetYaxis()->FindBin(transtht);
                int binz = histo_efficiency->GetZaxis()->FindBin(transthb);
                int bin = histo_efficiency->GetBin(binx, biny, binz);
                eff = histo_efficiency->GetBinContent(bin);
            }
            return eff;
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
