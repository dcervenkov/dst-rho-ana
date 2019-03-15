/*
 * DSRhoPeek.cc
 *
 *  Created on: Aug 20, 2014
 *      Author: D. Cervenkov
 */

#include "DSRhoPeek.h"

// Standard includes
#include <string>
#include <cstdio>

// ROOT includes
#include "TROOT.h"
#include "TEnv.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"
#include "THStack.h"
#include "TH1.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TColor.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooPlot.h"

// Local includes
#include "constants.h"
#include "colors.h"
#include "log.h"
#include "tools.h"

//#define GRAPHIC

int main(int argc, char* argv[]) {

#ifdef GRAPHIC
	TApplication* rootapp = new TApplication("DSRhoPeek",&argc,argv); //For graphic-output apps only
	/**
	 * When using TApplication it removes arguments it "handles" from
	 * the argument array. E.g. -b, -x, -q, --help, <dir>, <file>, etc.
	 * For more info read TApplication's GetOptions function help.
	 * The solution is to use rootapp->Argc() and rootapp->Argv(i).
	 * The next few lines are for compatibility of GRAPHIC vs. non-GRAPHIC -
	 * they recreate the original argc and argv when using GRAPHIC.
	 **/
	argc = rootapp->Argc();
	for(int i = 0; i < argc; i++) {
		argv[i] = rootapp->Argv(i);
	}
#endif // GRAPHIC

	if (argc != 2 && argc != 3){
		printf("ERROR: Wrong number of arguments.\n");
		printf("Usage: DSRhoPeek DATA-FILE [WITH_TRUE_CANDIDATE]\n");
		return 2;
	}
	const char* inputFileName = argv[1];
	std::string strOutputFileName(inputFileName);
	std::string nameAppendix = "_results";

	int withTrueCandidate = -1;
	if (argc == 3) {
		withTrueCandidate = atoi(argv[2]);
		nameAppendix += "_tc";
		nameAppendix += argv[2];
	}

	const char* outputFileName =
			(strOutputFileName.insert(strOutputFileName.rfind(".root"), nameAppendix.c_str())).c_str();

	tools::SetPlotDir("plots");
	colors::setColors();

	// Open file for reading
	Log::print(Log::debug, "Opening file %s...\n",inputFileName);
	TFile* inputFile = new TFile(inputFileName);
	if(inputFile->IsOpen() == kTRUE){
		Log::print(Log::info, "File %s opened.\n",inputFileName);
	} else {
		Log::print(Log::error, "Couldn't open file %s for reading. Exiting...\n", inputFileName);
		return 2;
	}

	TTree* Bcand = getBASFTree(inputFile,BASF_HIST_ID_BEST);

	// Open file for writing
	Log::print(Log::debug, "Opening file %s...\n",outputFileName);
	TFile* outputFile = new TFile(outputFileName,"RECREATE");
	if(outputFile->IsOpen() == kTRUE) {
		Log::print(Log::info, "File %s opened.\n",outputFileName);
	} else {
		Log::print(Log::error, "Couldn't open file %s for writing. Exiting...\n", inputFileName);
		return 2;
	}

	TCanvas* c_stacked = new TCanvas("c_stacked","B candidates stacked",1280,800);

	TObjArray* listOfBranches = Bcand->GetListOfBranches();
	int numBranches = listOfBranches->GetEntries();
	THStack* h_stacks[numBranches];
	TH1F* h_stacked[numBranches][NUM_MC_FLAGS];
	for (int branchNo = 0; branchNo < numBranches; ++branchNo) {
		const char* branchName = listOfBranches->At(branchNo)->GetName();
		h_stacks[branchNo] = new THStack(branchName,branchName);
		char name[200];
		double min = 0;
		double max = 0;
		getOptimalHistogramScope(Bcand, branchName, min, max);
		for(int flagNo = 0; flagNo < NUM_MC_FLAGS; ++flagNo) {
			sprintf(name,"h_stacked_%i_%i",branchNo, flagNo);
			h_stacked[branchNo][flagNo] = new TH1F(name,name,100,min,max);
			getHistogramFromTree(h_stacked[branchNo][flagNo],Bcand,branchName,MC_FLAGS_ASSIGNMENT[flagNo],withTrueCandidate);
			// Skip histograms with less than MIN_ENTRIES events to avoid clutter
			if(h_stacked[branchNo][flagNo]->GetEntries() < MIN_ENTRIES) continue;
			h_stacked[branchNo][flagNo]->SetLineColor(COLOR_MAP[flagNo]);
			h_stacked[branchNo][flagNo]->SetLineWidth(2);
			h_stacks[branchNo]->Add(h_stacked[branchNo][flagNo]);
		}
		h_stacks[branchNo]->Draw("nostack");
		c_stacked->SetName(branchName);
		//c_stacked->SaveAs(OUTPUT_FORMAT);
		h_stacks[branchNo]->Write();
		c_stacked->Write();
	}

	Log::print(Log::info, "# MC truth candidates:\n");
	for(int i = 0; i < NUM_MC_FLAGS; ++i) {
		Log::print(Log::info, "%3i: %.0f\n", MC_FLAGS_ASSIGNMENT[i],h_stacked[0][i]->GetEntries());
	}

	double numSignalEvts = 0;
	double numBackgroundEvts = 0;
	// Count number of signal and background events for figure-of-merit calculation
	// TODO: This has to be changed
	for(int i = 0; i < NUM_MC_FLAGS; ++i) {
		if (MC_FLAGS_ASSIGNMENT[i] == MC_FLAG_SIGNAL || MC_FLAGS_ASSIGNMENT[i] == MC_FLAG_SIGNAL_WITHOUT_FSR) {
			numSignalEvts += h_stacked[0][i]->GetEntries();
		} else {
			numBackgroundEvts += h_stacked[0][i]->GetEntries();
		}
	}
	Log::print(Log::info, "FOM = %f\n", numSignalEvts/sqrt(numSignalEvts + numBackgroundEvts));

	outputFile->Close();

#ifdef GRAPHIC
	printf("\nProgram execution has finished.\n");
	rootapp->Run();
#endif
	return 0;
}

TTree* getBASFTree(TFile* file, const int treeNumber){
	char numStr[21]; // all numbers up to 64-bits
	sprintf(numStr, "%i", treeNumber);
	std::string histName("h");
	histName += numStr;
	if(!file->GetListOfKeys()->Contains(histName.c_str())){
		Log::print(Log::error, "Object '%s' doesn't exist in '%s'! Exiting...\n",histName.c_str(),file->GetName());
		exit(3);
	}
	return (TTree*)file->Get(histName.c_str());
}

unsigned int getHistogramFromTree(TH1F* histogram, TTree* tree, const char* branchName, const int flag, const int withTrueCandidate){
	float datapoint;
	float eventCandidate;
	float mcflag;

	TBranch* MCFlagBranch = tree->GetBranch("evmcflag");
	TBranch* branch = tree->GetBranch(branchName);
	TBranch* candidateBranch = tree->GetBranch("candsel");
	branch->SetAddress(&datapoint);
	candidateBranch->SetAddress(&eventCandidate);
	MCFlagBranch->SetAddress(&mcflag);

	int nEvents = tree->GetEntries();
	for (int i = 0; i < nEvents; ++i) {
		candidateBranch->GetEntry(i);
		MCFlagBranch->GetEntry(i);

		// Separate events into histos based on mcflag
		if (int(mcflag) != flag) continue;
		// If true candidate requirement given as an argument, skip events with/without true candidate.
		// mcflag 8 is without the true candidate, but it has a signal-like peak, so we move it to
		// the same plot as signal for now.
		if (withTrueCandidate != -1) {
			if (withTrueCandidate == 0 && !(int(eventCandidate) == 0 &&	(int(mcflag) != 1 && int(mcflag) != 8)))
				continue;
			if (withTrueCandidate == 1 &&  (int(eventCandidate) == 0 &&	(int(mcflag) != 1 && int(mcflag) != 8)))
				continue;
		}
		// If branch == xBranch the first SetAddress is overridden by the second
		if (branch == candidateBranch) datapoint = eventCandidate;
		else if (branch == MCFlagBranch) datapoint = mcflag;
		else branch->GetEntry(i);

		histogram->Fill(datapoint);
	}
	// Required to avoid segfault when the datapoint, eventFlag and eventCandidate go out of scope
	branch->ResetAddress();
	candidateBranch->ResetAddress();
	MCFlagBranch->ResetAddress();

	return nEvents;
}

void getOptimalHistogramScope(TTree* tree, const char* branchName, double& min, double& max){
	double dataMin = tree->GetMinimum(branchName);
	double dataMax = tree->GetMaximum(branchName);
	double range = dataMax - dataMin;
	// The histogram range should be larger than the data range,
	// otherwise last bin (and possibly first bin) may be clipped
	min = dataMin - range * 0.05;
	max = dataMax + range * 0.05;

	if (min == max){
		min -= 1;
		max += 1;
	}
}



