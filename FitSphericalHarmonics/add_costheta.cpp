#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

int add_costheta_gsim(const char* filename) {

    TFile* f = new TFile(filename, "update");
    TTree* T = (TTree*)f->Get("h2000");
    float thetat, thetab;
    float costhetat, costhetab;
    TBranch* bcosthetat = T->Branch("costht", &costhetat, "costht/F");
    TBranch* bcosthetab = T->Branch("costhb", &costhetab, "costhb/F");
    T->SetBranchAddress("thetat", &thetat);
    T->SetBranchAddress("thetab", &thetab);
    Long64_t nentries = T->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        T->GetEntry(i);
        costhetat = TMath::Cos(thetat);
        costhetab = TMath::Cos(thetab);
        bcosthetat->Fill();
        bcosthetab->Fill();
    }
    T->Print();
    T->Write();
    delete f;

    return 0;
}

int add_costheta(const char* filename) {

    TFile* f = new TFile(filename, "update");
    TTree* T = (TTree*)f->Get("dataset");
    double thetat, thetab;
    double costhetat, costhetab;
    TBranch* bcosthetat = T->Branch("costht", &costhetat, "costht/D");
    TBranch* bcosthetab = T->Branch("costhb", &costhetab, "costhb/D");
    T->SetBranchAddress("thetat", &thetat);
    T->SetBranchAddress("thetab", &thetab);
    Long64_t nentries = T->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        T->GetEntry(i);
        costhetat = TMath::Cos(thetat);
        costhetab = TMath::Cos(thetab);
        bcosthetat->Fill();
        bcosthetab->Fill();
    }
    T->Print();
    T->Write();
    delete f;

    return 0;
}

int transform_thetab_gsim(const char* filename) {

    TFile* f = new TFile(filename, "update");
    TTree* T = (TTree*)f->Get("h2000");
    float thetat, thetab;
    float transthetat, transthetab;
    TBranch* btransthetat = T->Branch("transtht", &transthetat, "transtht/F");
    TBranch* btransthetab = T->Branch("transthb", &transthetab, "transthb/F");
    T->SetBranchAddress("thetat", &thetat);
    T->SetBranchAddress("thetab", &thetab);
    Long64_t nentries = T->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        T->GetEntry(i);
        transthetat = thetat/TMath::Pi() * 2 - 1;
        transthetab = thetab/TMath::Pi() * 2 - 1;
        btransthetat->Fill();
        btransthetab->Fill();
    }
    T->Print();
    T->Write();
    delete f;

    return 0;
}

int transform_thetab(const char* filename) {

    TFile* f = new TFile(filename, "update");
    TTree* T = (TTree*)f->Get("dataset");
    double thetat, thetab;
    double transthetat, transthetab;
    TBranch* btransthetat = T->Branch("transtht", &transthetat, "transtht/D");
    TBranch* btransthetab = T->Branch("transthb", &transthetab, "transthb/D");
    T->SetBranchAddress("thetat", &thetat);
    T->SetBranchAddress("thetab", &thetab);
    Long64_t nentries = T->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        T->GetEntry(i);
        transthetat = thetat/TMath::Pi() * 2 - 1;
        transthetab = thetab/TMath::Pi() * 2 - 1;
        btransthetat->Fill();
        btransthetab->Fill();
    }
    T->Print();
    T->Write();
    delete f;

    return 0;
}