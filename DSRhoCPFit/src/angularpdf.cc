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

//ClassImp(AngularPDF)

AngularPDF::AngularPDF(const char *name, const char *title, bool _B_bar, int _efficiency_model,
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
    // The rest of this constructor computes angular integration
    // of certain terms of the PDF. This is used to speed up
    // computation of normalization; see AngularPDF::analyticalIntegral
    ROOT::Math::Functor wf1(this, &AngularPDF::f1, 3);
    ROOT::Math::Functor wf2(this, &AngularPDF::f2, 3);
    ROOT::Math::Functor wf3(this, &AngularPDF::f3, 3);
    ROOT::Math::Functor wf4(this, &AngularPDF::f4, 3);
    ROOT::Math::Functor wf5(this, &AngularPDF::f5, 3);
    ROOT::Math::Functor wf6(this, &AngularPDF::f6, 3);

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
    efficiency_model(other.efficiency_model)
{
    for (int i = 0; i < 6; i++) {
        int_tht_thb_phit[i] = other.int_tht_thb_phit[i];
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

    // printf("eval: %f\n", value * eff.GetEfficiency(tht, thb, phit, efficiency_model));

    return value * eff.GetEfficiency(tht, thb, phit, efficiency_model);
}



Int_t AngularPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {
    // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED,
    // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS
    // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
    // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs
    // EXPRESSION MULTIPLE TIMES

    // if (matchArgs(allVars,analVars,x)) return 1 ;

    if(matchArgs(allVars,analVars,tht,thb,phit)) return 1;
//    if(matchArgs(allVars,analVars,tht,thb)) return 2;
//    if(matchArgs(allVars,analVars,tht,phit)) return 3;
//    if(matchArgs(allVars,analVars,thb,phit)) return 4;
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
//        return 32.*TMath::Pi()/9;
        // printf("int: %f\n", ap*ap * int_tht_thb_phit[0] +
        //        at*at * int_tht_thb_phit[1] +
        //        a0*a0 * int_tht_thb_phit[2] +
        //        ap0r * int_tht_thb_phit[3] -
        //        a0ti * int_tht_thb_phit[4] -
        //        apti * int_tht_thb_phit[5]);

        return ap*ap * int_tht_thb_phit[0] +
               at*at * int_tht_thb_phit[1] +
               a0*a0 * int_tht_thb_phit[2] +
               ap0r * int_tht_thb_phit[3] -
               a0ti * int_tht_thb_phit[4] -
               apti * int_tht_thb_phit[5];

    case 2: // Int[g,{tht,thb}]
        return 16./9.*(at*at + 2*a0*a0*cos(phit)*cos(phit) + 2*ap*ap*sin(phit)*sin(phit));

    case 3: // Int[g,{tht,phit}]
        return 8.*TMath::Pi()/3*((ap*ap+at*at)*sin(thb)*sin(thb) + 2*a0*a0*cos(thb)*cos(thb))*sin(thb);

    case 4: // Int[g,{thb,phit}]
        return 8.*TMath::Pi()/3*((ap*ap+a0*a0)*sin(tht)*sin(tht) + 2*at*at*cos(tht)*cos(tht))*sin(tht);

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
    Double_t val = 2 * sin(vars[0]) * sin(vars[0]) * sin(vars[0]) * sin(vars[1]) * sin(vars[1]) * sin(vars[1]) * sin(vars[2]) * sin(vars[2]);
    return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t AngularPDF::f2(const double * vars) {
    Double_t val = 2 * cos(vars[0]) * cos(vars[0]) * sin(vars[0]) * sin(vars[1]) * sin(vars[1]) * sin(vars[1]);
    return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t AngularPDF::f3(const double * vars) {
    Double_t val = 4 * sin(vars[0]) * sin(vars[0]) * sin(vars[0]) * cos(vars[1]) * cos(vars[1]) * sin(vars[1]) * cos(vars[2]) * cos(vars[2]);
    return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t AngularPDF::f4(const double * vars) {
    Double_t val = sqrt(2) * sin(vars[0]) * sin(vars[0]) * sin(vars[0]) * sin(2 * vars[1]) * sin(vars[1]) * sin(2 * vars[2]);
    return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t AngularPDF::f5(const double * vars) {
    Double_t val = sqrt(2) * sin(2 * vars[0]) * sin(vars[0]) * sin(2 * vars[1]) * sin(vars[1]) * cos(vars[2]);
    return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t AngularPDF::f6(const double * vars) {
    Double_t val = 2 * sin(2 * vars[0]) * sin(vars[0]) * sin(vars[1]) * sin(vars[1]) * sin(vars[1]) * sin(vars[2]);
    return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}