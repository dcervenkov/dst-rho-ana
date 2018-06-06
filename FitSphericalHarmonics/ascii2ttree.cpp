#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

int ascii2ttree(const char* file) {

    RooAbsData::setDefaultStorageType(RooAbsData::Tree);

    RooRealVar thetat("thetat", "thetat", 0, TMath::Pi());
    RooRealVar thetab("thetab", "thetab", 0, TMath::Pi());
    RooRealVar phit("phit", "phit", -TMath::Pi(), TMath::Pi());
    RooRealVar dt("dt", "dt", -1000, 1000);

    RooCategory decaytype("decaytype", "decaytype");
    decaytype.defineType("a", 1);
    decaytype.defineType("ab", 2);
    decaytype.defineType("b", 3);
    decaytype.defineType("bb", 4);

    RooDataSet* dataset = RooDataSet::read(file, RooArgList(thetat, thetab, phit, dt, decaytype));
    
    TFile new_file("from_ascii.root", "RECREATE");
    TTree* tree = (TTree*)dataset->tree();
    tree->Write();
    new_file.Close();

    return 0;
}
