/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "angularpdf.h"

// ROOT includes
#include "Math/Functor.h"
#include "Math/IntegratorMultiDim.h"

// Local includes
#include "log.h"

//ClassImp(AngularPDF)

AngularPDF::AngularPDF(const char *name, const char *title, bool _B_bar, int _efficiency_model, std::vector<std::string> _efficiency_files,
                   RooAbsReal& _tht,
                   RooAbsReal& _thb,
                   RooAbsReal& _phit,
                   RooAbsReal& _ap,
                   RooAbsReal& _apa,
                   RooAbsReal& _a0,
                   RooAbsReal& _ata) :
    RooAbsPdf(name,title),
    tht("tht","tht",this,_tht),
    thb("thb","thb",this,_thb),
    phit("phit","phit",this,_phit),
    ap("ap","ap",this,_ap),
    apa("apa","apa",this,_apa),
    a0("a0","a0",this,_a0),
    ata("ata","ata",this,_ata),
    B_bar(_B_bar),
    efficiency_model(_efficiency_model)
{
    for (auto file : _efficiency_files) {
        eff.ReadInFile(file.c_str());
    }

    // The rest of this constructor computes angular integration
    // of certain terms of the PDF. This is used to speed up
    // computation of normalization; see AngularPDF::analyticalIntegral
    CalculateBinnedIntegralFunctors();

    ROOT::Math::Functor wf1(this, &AngularPDF::f1e, 3);
    ROOT::Math::Functor wf2(this, &AngularPDF::f2e, 3);
    ROOT::Math::Functor wf3(this, &AngularPDF::f3e, 3);
    ROOT::Math::Functor wf4(this, &AngularPDF::f4e, 3);
    ROOT::Math::Functor wf5(this, &AngularPDF::f5e, 3);
    ROOT::Math::Functor wf6(this, &AngularPDF::f6e, 3);

    double a[] = {tht.min(), thb.min(), phit.min()};
    double b[] = {tht.max(), thb.max(), phit.max()};

    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(wf1);
    int_tht_thb_phit[0] = ig.Integral(a,b);
    ig.SetFunction(wf2);
    int_tht_thb_phit[1] = ig.Integral(a,b);
    ig.SetFunction(wf3);
    int_tht_thb_phit[2] = ig.Integral(a,b);
    ig.SetFunction(wf4);
    int_tht_thb_phit[3] = ig.Integral(a,b);
    ig.SetFunction(wf5);
    int_tht_thb_phit[4] = ig.Integral(a,b);
    ig.SetFunction(wf6);
    int_tht_thb_phit[5] = ig.Integral(a,b);

}


AngularPDF::AngularPDF(const AngularPDF& other, const char* name) :
    RooAbsPdf(other,name),
    tht("tht",this,other.tht),
    thb("thb",this,other.thb),
    phit("phit",this,other.phit),
    ap("ap",this,other.ap),
    apa("apa",this,other.apa),
    a0("a0",this,other.a0),
    ata("ata",this,other.ata),
    B_bar(other.B_bar),
    efficiency_model(other.efficiency_model),
    eff(other.eff)
{
    for (int i = 0; i < 6; i++) {
        int_tht_thb_phit[i] = other.int_tht_thb_phit[i];
        fxeff_tht[i] = other.fxeff_tht[i];
        fxeff_thb[i] = other.fxeff_thb[i];
        fxeff_phit[i] = other.fxeff_phit[i];
    }
}


Double_t AngularPDF::evaluate() const {

    Double_t a0a = 0;
    Double_t at = sqrt(1-ap*ap-a0*a0);

    Double_t ap0r = ap*a0*cos(-apa+a0a);
    Double_t a0ti = a0*at*sin(-a0a+ata);
    Double_t apti = ap*at*sin(-apa+ata);

    // Add pi to at phase for antiparticles
    if (B_bar) {
        a0ti = -a0ti;
        apti = -apti;
    }

    Double_t value = ap*ap*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
                     at*at*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
                     a0*a0*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
                     sqrt(2)*ap0r*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
                     sqrt(2)*a0ti*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
                     2*apti*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit);

    // Log::print(Log::debug, "eval: %f\n", value * eff.GetEfficiency(tht, thb, phit, efficiency_model));

    return value * eff.GetEfficiency(tht, thb, phit, efficiency_model);
}



Int_t AngularPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {
    // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED,
    // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS
    // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
    // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs
    // EXPRESSION MULTIPLE TIMES

    // if (matchArgs(allVars,analVars,x)) return 1 ;

    if (matchArgs(allVars, analVars, tht, thb, phit)) return 1;
    if (matchArgs(allVars, analVars, tht, thb)) return 2;
    if (matchArgs(allVars, analVars, tht, phit)) return 3;
    if (matchArgs(allVars, analVars, thb, phit)) return 4;
    //    if(matchArgs(allVars,analVars,tht)) return 5;
    //    if(matchArgs(allVars,analVars,thb)) return 6;
    //    if(matchArgs(allVars,analVars,phit)) return 7;

    return 0 ;
}



Double_t AngularPDF::analyticalIntegral(Int_t code, const char* rangeName) const {
    // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
    // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
    // BOUNDARIES FOR EACH OBSERVABLE x

    // assert(code==1) ;
    // return (x.max(rangeName)-x.min(rangeName)) ;

    Double_t a0a = 0;
    Double_t at = sqrt(1-ap*ap-a0*a0);

    Double_t ap0r = ap*a0*cos(-apa+a0a);
//    Double_t ap0i = ap*a0*sin(-apa+a0a);
//    Double_t a0tr = a0*at*cos(-a0a+ata);
    Double_t a0ti = a0*at*sin(-a0a+ata);
//    Double_t aptr = ap*at*cos(-apa+ata);
    Double_t apti = ap*at*sin(-apa+ata);

    switch(code) {
    case 1: // Int[g,{tht,thb,phit}]
        // return 32.*TMath::Pi()/9;
        return ap*ap * int_tht_thb_phit[0] +
               at*at * int_tht_thb_phit[1] +
               a0*a0 * int_tht_thb_phit[2] +
               ap0r * int_tht_thb_phit[3] -
               a0ti * int_tht_thb_phit[4] -
               apti * int_tht_thb_phit[5];

    case 2: // Int[g,{tht,thb}]
        // return 16./9.*(at*at + 2*a0*a0*cos(phit)*cos(phit) + 2*ap*ap*sin(phit)*sin(phit));
        return ap*ap * fxeff_phit[0]->Interpolate(phit) +
               at*at * fxeff_phit[1]->Interpolate(phit) +
               a0*a0 * fxeff_phit[2]->Interpolate(phit) +
               ap0r  * fxeff_phit[3]->Interpolate(phit) -
               a0ti  * fxeff_phit[4]->Interpolate(phit) -
               apti  * fxeff_phit[5]->Interpolate(phit);

    case 3: // Int[g,{tht,phit}]
        // return 8.*TMath::Pi()/3*((ap*ap+at*at)*sin(thb)*sin(thb) + 2*a0*a0*cos(thb)*cos(thb))*sin(thb);
        return ap*ap * fxeff_thb[0]->Interpolate(thb) +
               at*at * fxeff_thb[1]->Interpolate(thb) +
               a0*a0 * fxeff_thb[2]->Interpolate(thb) +
               ap0r  * fxeff_thb[3]->Interpolate(thb) -
               a0ti  * fxeff_thb[4]->Interpolate(thb) -
               apti  * fxeff_thb[5]->Interpolate(thb);

    case 4: // Int[g,{thb,phit}]
        // return 8.*TMath::Pi()/3*((ap*ap+a0*a0)*sin(tht)*sin(tht) + 2*at*at*cos(tht)*cos(tht))*sin(tht);
        return ap*ap * fxeff_tht[0]->Interpolate(tht) +
               at*at * fxeff_tht[1]->Interpolate(tht) +
               a0*a0 * fxeff_tht[2]->Interpolate(tht) +
               ap0r  * fxeff_tht[3]->Interpolate(tht) -
               a0ti  * fxeff_tht[4]->Interpolate(tht) -
               apti  * fxeff_tht[5]->Interpolate(tht);

    case 5: // Int[g,{tht}]
        return 4./3.*(4*a0*a0*cos(phit)*cos(phit)*cos(thb)*cos(thb) + \
       (at*at + 2*ap*ap*sin(phit)*sin(phit))*sin(thb)*sin(thb) + sqrt(2)*ap0r*sin(2*phit)*sin(2*thb))*sin(thb);

    case 6: // Int[g,{thb}]
        return 8./3.*(at*at*cos(tht)*cos(tht) + a0*a0*cos(phit)*cos(phit)*sin(tht)*sin(tht) + \
       sin(phit)*(ap*ap*sin(phit)*sin(tht)*sin(tht) - apti*sin(2*tht)))*sin(tht);

    case 7: // Int[g,{phit}]
        return 2*TMath::Pi()*(2*at*at*cos(tht)*cos(tht)*sin(thb)*sin(thb) +
       (2*a0*a0*cos(thb)*cos(thb) + ap*ap*sin(thb)*sin(thb))*sin(tht)*sin(tht))*(sin(tht)*sin(thb));

    default:
        return 0;
    }

}

Double_t AngularPDF::f1(const double * vars) {
    return 2 * sin(vars[0]) * sin(vars[0]) * sin(vars[0]) * sin(vars[1]) * sin(vars[1]) * sin(vars[1]) * sin(vars[2]) * sin(vars[2]);
}

Double_t AngularPDF::f2(const double * vars) {
    return 2 * cos(vars[0]) * cos(vars[0]) * sin(vars[0]) * sin(vars[1]) * sin(vars[1]) * sin(vars[1]);
}

Double_t AngularPDF::f3(const double * vars) {
    return 4 * sin(vars[0]) * sin(vars[0]) * sin(vars[0]) * cos(vars[1]) * cos(vars[1]) * sin(vars[1]) * cos(vars[2]) * cos(vars[2]);
}

Double_t AngularPDF::f4(const double * vars) {
    return sqrt(2) * sin(vars[0]) * sin(vars[0]) * sin(vars[0]) * sin(2 * vars[1]) * sin(vars[1]) * sin(2 * vars[2]);
}

Double_t AngularPDF::f5(const double * vars) {
    return sqrt(2) * sin(2 * vars[0]) * sin(vars[0]) * sin(2 * vars[1]) * sin(vars[1]) * cos(vars[2]);
}

Double_t AngularPDF::f6(const double * vars) {
    return 2 * sin(2 * vars[0]) * sin(vars[0]) * sin(vars[1]) * sin(vars[1]) * sin(vars[1]) * sin(vars[2]);
}


Double_t AngularPDF::f1e(const double * vars) {
    return f1(vars) * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t AngularPDF::f2e(const double * vars) {
    return f2(vars) * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t AngularPDF::f3e(const double * vars) {
    return f3(vars) * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t AngularPDF::f4e(const double * vars) {
    return f4(vars) * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t AngularPDF::f5e(const double * vars) {
    return f5(vars) * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t AngularPDF::f6e(const double * vars) {
    return f6(vars) * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}
// Double_t AngularPDF::

TH3F* AngularPDF::GetBinnedTrigXEfficiency(Double_t (*trig)(const double * vars), int model, int nbins, double limits[6]) {
    TH3F* binned = new TH3F("binned_eff", "binned_eff", nbins, limits[0], limits[1], nbins,
                            limits[2], limits[3], nbins, limits[4], limits[5]);

    // double thetat, thetab, phit;
    double vars[3];
    for (int x = 1; x <= nbins; x++) {
        for (int y = 1; y <= nbins; y++) {
            for (int z = 1; z <= nbins; z++) {
                vars[0] = binned->GetXaxis()->GetBinCenter(x);
                vars[1] = binned->GetYaxis()->GetBinCenter(y);
                vars[2] = binned->GetZaxis()->GetBinCenter(z);
                binned->SetBinContent(binned->GetBin(x, y, z),
                                      trig(vars) * eff.GetEfficiency(vars[0], vars[1], vars[2], model));
            }
        }
    }
    return binned;
}

void AngularPDF::CalculateBinnedIntegralFunctors() {
    double limits[6] = {tht.min(), tht.max(), thb.min(), thb.max(), phit.min(), phit.max()};

    Double_t (* func [6])(const double*);
    func[0] = f1;
    func[1] = f2;
    func[2] = f3;
    func[3] = f4;
    func[4] = f5;
    func[5] = f6;

    const int nbins = 100;
    double tht_bin_width = (tht.max() - tht.min())/nbins;
    double thb_bin_width = (thb.max() - thb.min())/nbins;
    double phit_bin_width = (phit.max() - phit.min())/nbins;

    Log::LogLine(Log::info) << "Computing binned numerical integral functors...";
    for (int i = 0; i < 6; i++) {
        TH3F* fxeff = GetBinnedTrigXEfficiency(func[i], efficiency_model, nbins, limits);
        TString name_base = "f" + i;
        fxeff_tht[i] = fxeff->ProjectionX(name_base + "_tht");
        fxeff_thb[i] = fxeff->ProjectionY(name_base + "_thb");
        fxeff_phit[i] = fxeff->ProjectionZ(name_base + "_phit");

        fxeff_tht[i]->Scale(thb_bin_width * phit_bin_width);
        fxeff_thb[i]->Scale(tht_bin_width * phit_bin_width);
        fxeff_phit[i]->Scale(tht_bin_width * thb_bin_width);
        delete fxeff;
    }
}
