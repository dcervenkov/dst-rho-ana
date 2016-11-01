//  C++ Standard Library
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <unistd.h>

// BASF headers
#include "basf/basfout.h"
#include "basf/basfshm.h"
#include "belle.h"
#include "benergy/BeamEnergy.h"
#include "eid/eid.h"
#include "event/BelleEvent.h"
#include "ip/IpProfile.h"
#include "mdst/Muid_mdst.h"
#include "particle/combination.h"
#include "particle/utility.h"
#include "tables/belletdf.h"
#include "tables/evtcls.h"
#include "tables/hepevt.h"
#include "tables/level4.h"
#include "tables/mctype.h"
#include "tables/mdst.h"
#include "tuple/BelleTupleManager.h"

// Local includes
#include "constants.h"
#include "DSRhoSkimModule.h"
#include "geninfo.h"
#include "myutils.h"
#include "UserInfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

int expNo, runNo, evtNo;

DSRhoSkimModule::DSRhoSkimModule(void) {
	// Module constructor
	// initialize module members
	expNo = 0;
	runNo = 0;
	evtNo = 0;

	myHER = 0.0;
	myLER = 0.0;

	myCMSBoost = Hep3Vector();
	myCMSE = 0.0;
	dataType = -1;
	ipUsable = -1;

}

DSRhoSkimModule::~DSRhoSkimModule(void) {
}

void DSRhoSkimModule::disp_stat(const char*) {
}

void DSRhoSkimModule::end_run(BelleEvent*, int*) {
}

void DSRhoSkimModule::other(int*, BelleEvent*, int*) {
}

void DSRhoSkimModule::init(int*) {
	Ptype ptype_dummy("VPHO");
}

void DSRhoSkimModule::begin_run(BelleEvent*, int *status) {
	// Used for skimming; set to 1 if you want to save the index of an event
	*status = 0;

	// Get IP profile data from $BELLE_POSTGRES_SERVER
	IpProfile::begin_run();
	bool usableIP = IpProfile::usable();
	if (usableIP)
		ipUsable = 1;
	else
		ipUsable = 0;

	// get the Beam Energy and boost to CMS Vector
	BeamEnergy::begin_run();
	myHER = BeamEnergy::E_HER();
	myLER = BeamEnergy::E_LER();
	myCMSBoost = -BeamEnergy::CMBoost(); // note the "-" sign
	myCMSE = BeamEnergy::E_beam_corr();

	//Running MC = 2 or DATA = 1 // TODO: What is this?
	Belle_runhead_Manager &runhead_m = Belle_runhead_Manager::get_manager();
	Belle_runhead &runhead = runhead_m((Panther_ID) 1);

	if (runhead) {
		if (runhead.ExpMC() == 1) {
			dataType = 0; // running on real data
		} else {
			dataType = 1; // running on MC
		}
	}

	// init eid
	eid::init_data();

}

void DSRhoSkimModule::term(void) {
}

void DSRhoSkimModule::hist_def(void) {
}

void DSRhoSkimModule::event(BelleEvent*, int *status) {
	// Used for skimming; default to 0, change to 1 if we find a candidate
	*status = 0;

	// Obtain the ExpnNo, runNo and evtNo
	expNo = runNo = evtNo = 0;
	Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
	if (evtMgr.count()) {
		expNo = evtMgr[0].ExpNo();
		runNo = evtMgr[0].RunNo();
		evtNo = evtMgr[0].EvtNo();
	}

	// obtain final state particles
	// Charged kaons and pions
	std::vector<Particle> k_m;
	std::vector<Particle> k_p;
	std::vector<Particle> pi_m;
	std::vector<Particle> pi_p;

	// neutral pions
	std::vector<Particle> pi0;

	// use all tracks to build kaons and pions
	// there is no PID selection applied
	// function is implemented in myutils.(h,cc)
	fillKaonsAndPions(k_p, k_m, pi_p, pi_m);

	// PID cuts
	withKaonId(k_p, KAON_ID, 3, 1, 5, 3, 2);
	withKaonId(k_m, KAON_ID, 3, 1, 5, 3, 2);

	withPionId(pi_p, PION_ID, 3, 1, 5, 3, 2);
	withPionId(pi_m, PION_ID, 3, 1, 5, 3, 2);

	// Impact parameter cuts
	withImpactPars(k_p, IMPACT_PAR_DR, IMPACT_PAR_DZ);
	withImpactPars(k_m, IMPACT_PAR_DR, IMPACT_PAR_DZ);
	withImpactPars(pi_p, IMPACT_PAR_DR, IMPACT_PAR_DZ);
	withImpactPars(pi_m, IMPACT_PAR_DR, IMPACT_PAR_DZ);

	// Get Pi0 candidates; certain selection criteria used
	fillAllPi0(pi0);

	// Reconstruction of neutral D modes
	std::vector<Particle> D0_ch1;
	std::vector<Particle> antiD0_ch1;
//	combination(D0_ch1, Ptype(421), k_m, pi_p, D0_MASS_WINDOW);
//	combination(antiD0_ch1, Ptype(-421), k_p, pi_m, D0_MASS_WINDOW);
//	setDecayMode(D0_ch1, 1);
//	setDecayMode(antiD0_ch1, 1);

	std::vector<Particle> D0_ch2;
	std::vector<Particle> antiD0_ch2;
	combination(D0_ch2, Ptype(421), k_m, pi_p, pi0, D0_MASS_WINDOW);
	combination(antiD0_ch2, Ptype(-421), k_p, pi_m, pi0, D0_MASS_WINDOW);
	setDecayMode(D0_ch2, 2);
	setDecayMode(antiD0_ch2, 2);

	std::vector<Particle> D0_ch3;
	std::vector<Particle> antiD0_ch3;
//	combination(D0_ch3, Ptype(421), k_m, pi_p, pi_p, pi_m, D0_MASS_WINDOW);
//	combination(antiD0_ch3, Ptype(-421), k_p, pi_m, pi_p, pi_m, D0_MASS_WINDOW);
//	setDecayMode(D0_ch3, 3);
//	setDecayMode(antiD0_ch3, 3);

	std::vector<Particle> D0;
	std::vector<Particle> antiD0;
	PlistCopy(D0_ch1, D0);
	PlistCopy(D0_ch2, D0);
	PlistCopy(D0_ch3, D0);
	PlistCopy(antiD0_ch1, antiD0);
	PlistCopy(antiD0_ch2, antiD0);
	PlistCopy(antiD0_ch3, antiD0);


	std::vector<Particle> DS_p;
	std::vector<Particle> DS_m;

	// combination, but without any mass cut
	// (we will make |M_DS - M_D| cuts later)
	combination(DS_p, Ptype(413), D0, pi_p);
	combination(DS_m, Ptype(-413), antiD0, pi_m);

	// Cut on the mass difference between D* and its first daughter - D0
	withMassDifCut(DS_p, DSD0_MASS_DIFF_LOW, DSD0_MASS_DIFF_HIGH, 0);
	withMassDifCut(DS_m, DSD0_MASS_DIFF_LOW, DSD0_MASS_DIFF_HIGH, 0);

	std::vector<Particle> rho_p;
	std::vector<Particle> rho_m;

	combination(rho_p, Ptype(213), pi0, pi_p, RHO_MASS_WINDOW);
	combination(rho_m, Ptype(-213), pi0, pi_m, RHO_MASS_WINDOW);

	std::vector<Particle> B0;
	std::vector<Particle> antiB0;

	// Combination, but without any mass cut
	// (we will make deltaE, Mbc cuts later)
	combination(B0, Ptype(511), DS_m, rho_p);
	combination(antiB0, Ptype(-511), DS_p, rho_m);

	// Copy the B0 and antiB0 lists into a common B list
	std::vector<Particle> B;
	PlistCopy(B0, B);
	PlistCopy(antiB0, B);

	// Loop over all B candidates, calculate DeltaE and Mbc,
	// keep only candidates passing Mbc and DeltaE cuts
	for (unsigned i = 0; i < B.size(); i++) {
		Particle Bcand = B[i];

		HepLorentzVector B0mom(Bcand.p());
		B0mom.boost(myCMSBoost);
		double mbc = sqrt((myCMSE * myCMSE) - B0mom.vect().mag2());

		double Energy = B0mom.e();
		double deltaE = Energy - myCMSE;

		// save Mbc & DeltaE to UserInfo
		if (!&B[i].userInfo())
			setUserInfo(B[i]);
		dynamic_cast<UserInfo&>(B[i].userInfo()).setMbc(mbc);
		dynamic_cast<UserInfo&>(B[i].userInfo()).setDeltaE(deltaE);

		if (mbc < MBC_CUT || fabs(deltaE) > DELTA_E_CUT) {
			// bad candidate, remove it
			B.erase(B.begin() + i);
			--i;
			continue;
		}
	}

	if (B.size()){
		*status = 1;
	}
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
