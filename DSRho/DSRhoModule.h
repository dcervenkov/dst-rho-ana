/*
 Module recostructs B- -> D*- rho+ decays and c.c.,
 where D0 decays are reconstructed in the following final
 states:
 D0 -> K- pi+
 D0 -> K- pi+ pi0
 D0 -> Kshort pi+ pi-

 The module saves the relevant info to ntuple.
 */

#ifndef DSRHOMODULE_H_
#define DSRHOMODULE_H_
// BASF headers
#include "belle.h"
#include "basf/module.h"
#include "particle/Particle.h"
#include "tuple/BelleTupleManager.h"
#include "kfitter/kmassfitter.h"

// Local includes
#include "constants.h"
#include "printTree.h"
#include "Continuum.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class DSRhoModule: public Module {
public:
	DSRhoModule(void);
	~DSRhoModule(void);
	void init(int*);
	void term(void);
	void disp_stat(const char *cmd);
	void hist_def(void);
	void event(BelleEvent*, int*);
	void begin_run(BelleEvent *ptr, int*);
	void end_run(BelleEvent *ptr, int*);
	void other(int *rec, BelleEvent *ptr, int*);

	// Parameter set from outside; used to select which channel to reconstruct
	char channel[128];
	char svd[10];
	char sidebands[10];

private:
	void saveToTuple(Particle Bcand, BelleTuple* tupleB);
	void printEvent(Particle Bcand);
	void fitD0(std::vector<Particle> &particles);
	void fitDS(std::vector<Particle> &particles);

	// ntuples for saving the B candidates
	BelleTuple *tupleB;

	// ntuple for saving metadata for the whole run
	// (not run, but the data being analyzed, actually)
	BelleTuple* tupleRun;
	// Metadata variables
	unsigned eventsWithRightDecay; // num of events where the right decay was generated
	unsigned eventsWithOutRightDecay; // num of events where the right decay was NOT generated
	unsigned eventsWithRightCandidate; // num of events where the right candidate was present
	unsigned eventsWithOutRightCandidate; // num of events where the right candidate was NOT present
	unsigned numCorrectlySelectedCandidates; // number of times BCS selected the right candidate when present

	unsigned eventsDiscardedByNumCandidatesCut;
	unsigned eventsDiscardedByDEandMbcCuts;
	unsigned eventsDiscardedByCSCut;

	int expMC, detVer, expNo, runNo, evtNo;

	// beam energy
	double myHER;
	double myLER;

	// boost to CMS vector
	Hep3Vector myCMSBoost;

	// CMS energy
	double myCMSE;

	int dataType; // 0 = DATA, 1 = MC
	int ipUsable;

	int channelNum;
	bool useSidebands;

	printTree* pt;
	Continuum continuumSupression;

	int numPrintedEvents;
};

// Register module with BASF
extern "C" Module_descr *mdcl_dsrhomodule() {
	DSRhoModule *module = new DSRhoModule;
	Module_descr *dscr = new Module_descr("dsrhomodule", module);
	dscr->define_param("channel" ,"Decay channel to be reconstructed", sizeof(module->channel)/sizeof(char), module->channel);
	dscr->define_param("svd" ,"Weight file for which SVD to be used", sizeof(module->svd)/sizeof(char), module->svd);
	dscr->define_param("sidebands" ,"Whether to process sidebands or not", sizeof(module->sidebands)/sizeof(char), module->sidebands);
	basf_save_mdsc(dscr);
	return dscr;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif /* DSRHOMODULE_H_ */

