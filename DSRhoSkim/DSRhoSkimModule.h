/*
 Module recostructs B- -> D*- rho+ decays and c.c.,
 where D0 decays are reconstructed in the following final
 states:
 D0 -> K- pi+
 D0 -> K- pi+ pi0
 D0 -> Kshort pi+ pi-

 The module saves the relevant info to ntuple.
 */

#ifndef DSRHOSKIMMODULE_H_
#define DSRHOSKIMMODULE_H_
// BASF headers
#include "belle.h"
#include "basf/module.h"
#include "particle/Particle.h"
#include "tuple/BelleTupleManager.h"

// Local includes
#include "constants.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class DSRhoSkimModule: public Module {
public:
	DSRhoSkimModule(void);
	~DSRhoSkimModule(void);
	void init(int*);
	void term(void);
	void disp_stat(const char *cmd);
	void hist_def(void);
	void event(BelleEvent*, int*);
	void begin_run(BelleEvent *ptr, int*);
	void end_run(BelleEvent *ptr, int*);
	void other(int *rec, BelleEvent *ptr, int*);

private:
	int expNo, runNo, evtNo;

	// beam energy
	double myHER;
	double myLER;

	// boost to CMS vector
	Hep3Vector myCMSBoost;

	// CMS energy
	double myCMSE;

	int dataType; // 0 = DATA, 1 = MC
	int ipUsable;
};

// Register module with BASF
extern "C" Module_descr *mdcl_dsrhoskimmodule() {
	DSRhoSkimModule *module = new DSRhoSkimModule;
	Module_descr *dscr = new Module_descr("dsrhoskimmodule", module);
	basf_save_mdsc(dscr);
	return dscr;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif /* DSRHOSKIMMODULE_H_ */

