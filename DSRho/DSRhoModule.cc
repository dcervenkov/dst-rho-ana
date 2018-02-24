//  C++ Standard Library
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>
#include <unistd.h>

// BASF headers
#include "basf/basfout.h"
#include "basf/basfshm.h"
#include "belle.h"
#include "benergy/BeamEnergy.h"
#include "eid/eid.h"
#include "event/BelleEvent.h"
#include "hamlet/Hamlet.h"
#include "ip/IpProfile.h"
#include "kfitter/kmassfitter.h"
#include "mdst/mdst.h"
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
#include "DSRhoModule.h"
#include "geninfo.h"
#include "myutils.h"
#include "UserInfo.h"
#include "ksfwmoments.h"
#include "printTree.h"
//#include "Continuum.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

int expNo, runNo, evtNo;

DSRhoModule::DSRhoModule(void) {
	// Module constructor
	// initialize module members
	expMC = 0;
	detVer = 0;
	expNo = 0;
	runNo = 0;
	evtNo = 0;

	myHER = 0.0;
	myLER = 0.0;

	myCMSBoost = Hep3Vector();
	myCMSE = 0.0;
	dataType = -1;
	ipUsable = -1;
	tupleB = 0;

	tupleRun = 0;
	eventsWithRightCandidate = 0;
	eventsWithOutRightCandidate = 0;
	eventsWithRightDecay = 0;
	eventsWithOutRightDecay = 0;
	numCorrectlySelectedCandidates = 0;

	eventsDiscardedByNumCandidatesCut = 0;
	eventsDiscardedByDEandMbcCuts = 0;
	eventsDiscardedByCSCut = 0;

	sprintf(channel, "none");
	channelNum = 0;

	pt = 0;
}

DSRhoModule::~DSRhoModule(void) {
	delete pt;
}

void DSRhoModule::disp_stat(const char*) {
}

void DSRhoModule::end_run(BelleEvent*, int*) {
}

void DSRhoModule::other(int*, BelleEvent*, int*) {
}

void DSRhoModule::init(int*) {
	Hamlet::init();
	Ptype ptype_dummy("VPHO");
	printf ("*** DC: Init - channel = %s, svd = %s\n", channel, svd);
	pt = new printTree(channel);
	continuumSupression.readInWeights(channel, svd);
}

void DSRhoModule::begin_run(BelleEvent*, int *status) {
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

	Hamlet::begin_run(Hamlet::MULT_DIM_LH);

	if (strcmp(channel, "D0Kpi") == 0) channelNum = 1;
	else if (strcmp(channel, "D0Kpipi0") == 0) channelNum = 2;
	else if (strcmp(channel, "D0K3pi") == 0) channelNum = 3;
}

void DSRhoModule::term(void) {
	if (dataType) {
		printf("\n*** DC: Candidate Selection Metadata ***\n");

		printf("Evts discarded by DE and Mbc cuts: %i\n", eventsDiscardedByDEandMbcCuts);
		printf("Evts discarded by num candidates cut: %i\n", eventsDiscardedByNumCandidatesCut);
		printf("Evts discarded by CS cut: %i\n", eventsDiscardedByCSCut);

		printf("Evts without right decay: %i\nEvts with right decay: %i\n",\
				eventsWithOutRightDecay, eventsWithRightDecay);

		printf("--- Evts having at least one candidate ---\n");

		printf("Evts without right candidate: %i\nEvts with right candidate: %i\n",\
				eventsWithOutRightCandidate, eventsWithRightCandidate);

		printf("Correct candidate selected: %i\nWrong candidate selected: %i\n\n",\
				numCorrectlySelectedCandidates, eventsWithRightCandidate - numCorrectlySelectedCandidates);

		printf("Correct/PotentiallyCorrect = %.1f%%\n",
				100*float(numCorrectlySelectedCandidates)/eventsWithRightCandidate);

		printf("Selection efficiency = %.1f%%\n\n", 100*float(numCorrectlySelectedCandidates)/eventsWithRightDecay);
	}
}

void DSRhoModule::hist_def(void) {
	extern BelleTupleManager* BASF_Histogram;

	// Tuples to be created in the HBOOK file
	std::string tuples("benergy brecflavor btagmclink candsel csbdtg csmlp d0chi2 d0mass d0massbf d0mcflag d0mclink d0pcms "
			"d0pi0chi de detver dpi0mcflag dsd0diff dspiinv dsmass dsmcflag dsmclink dspcms "
			"e9oe25g1 e9oe25g2 energyg1 energyg2 evmcflag evtno expmc expno mbc mcflag mclink mdspi mdspi0 "
			"nocand phit phitg pi0chi2 pi0decangle pi0mass pi0mcflag pi0pcms rhomass rhomcflag rhomclink "
			"rhopcms rhopidr rhopidz rhopimcflag rhopinrf rhopinz rhopipcms runno "
			"shcosthb shcostht shcosttb spidr spidz spimcflag spimclink spinrf spinz spipcms "
			"thetab thetabg thetag1 thetag2 thetat thetatg "
			"vrusable vrvtxx vrvtxy vrvtxz vrerr1 vrerr2 vrerr3 vrerr4 vrerr5 vrerr6 vrchi2 vreffxi vrndf vreffndf vrntrk "
			"vtusable vtvtxx vtvtxy vtvtxz vterr1 vterr2 vterr3 vterr4 vterr5 vterr6 vtchi2 vtndf vtntrk "
			"vrgexist vrgvtxx vrgvtxy vrgvtxz vtgexist vtgvtxx vtgvtxy vtgvtxz vtistagl "
			"tagflavor tagqr tagwtag ");

	tuples +=
			"k0mm2 k0et "
			"k0hoo0 k0hoo1 k0hoo2 k0hoo3 k0hoo4 "
			"k0hso00 k0hso01 k0hso02 k0hso03 k0hso04 "
			"k0hso10 k0hso12 k0hso14 k0hso20 k0hso22 k0hso24 "
			"k1mm2 k1et "
			"k1hoo0 k1hoo1 k1hoo2 k1hoo3 k1hoo4 "
			"k1hso00 k1hso01 k1hso02 k1hso03 k1hso04 "
			"k1hso10 k1hso12 k1hso14 k1hso20 k1hso22 k1hso24 ";

	// All candidates will be saved to a single ntuple
	tupleB = (*BASF_Histogram).ntuple("Bcand", tuples.c_str(), BASF_HIST_ID_BEST);
}

void DSRhoModule::event(BelleEvent*, int *status) {
	// Used for skimming; default to 0, change to 1 if we find a candidate
	*status = 0;

	// Obtain the ExpnNo, runNo and evtNo
	expMC = detVer = expNo = runNo = evtNo = 0;
	Belle_event_Manager &evtMgr = Belle_event_Manager::get_manager();
	if (evtMgr.count()) {
		expMC = evtMgr[0].ExpMC();
		detVer = evtMgr[0].DetVer();
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

	// store pis without SVD cuts for slow pi candidates
	std::vector<Particle> pi_m_wo_SVD_cuts = pi_m;
	std::vector<Particle> pi_p_wo_SVD_cuts = pi_p;

	// Cuts on the minimal number of hits in the r-phi a z planes
	withSVD2(k_p, N_HITS_SVD_R, N_HITS_SVD_Z);
	withSVD2(k_m, N_HITS_SVD_R, N_HITS_SVD_Z);
	withSVD2(pi_p, N_HITS_SVD_R, N_HITS_SVD_Z);
	withSVD2(pi_m, N_HITS_SVD_R, N_HITS_SVD_Z);

	// Get Pi0 candidates; certain selection criteria used, then set the error matrix
	fillAllPi0(pi0);
	setPi0Error(pi0);

	// Reconstruction of neutral D modes
	std::vector<Particle> D0_ch1;
	std::vector<Particle> antiD0_ch1;
	if (channelNum == 1) {
		combination(D0_ch1, Ptype(421), k_m, pi_p, D0_MASS_WINDOW);
		combination(antiD0_ch1, Ptype(-421), k_p, pi_m, D0_MASS_WINDOW);
		setDecayMode(D0_ch1, 1);
		setDecayMode(antiD0_ch1, 1);
	}

	std::vector<Particle> D0_ch2;
	std::vector<Particle> antiD0_ch2;
	if (channelNum == 2) {
		combination(D0_ch2, Ptype(421), k_m, pi_p, pi0, D0_PI0_MASS_WINDOW_LOW, D0_PI0_MASS_WINDOW_HIGH);
		combination(antiD0_ch2, Ptype(-421), k_p, pi_m, pi0, D0_PI0_MASS_WINDOW_LOW, D0_PI0_MASS_WINDOW_HIGH);
		setDecayMode(D0_ch2, 2);
		setDecayMode(antiD0_ch2, 2);
	}

	std::vector<Particle> D0_ch3;
	std::vector<Particle> antiD0_ch3;
	if (channelNum == 3) {
		combination(D0_ch3, Ptype(421), k_m, pi_p, pi_p, pi_m, D0_3PI_MASS_WINDOW);
		combination(antiD0_ch3, Ptype(-421), k_p, pi_m, pi_p, pi_m, D0_3PI_MASS_WINDOW);
		setDecayMode(D0_ch3, 3);
		setDecayMode(antiD0_ch3, 3);
	}

	std::vector<Particle> D0;
	std::vector<Particle> antiD0;
	PlistCopy(D0_ch1, D0);
	PlistCopy(D0_ch2, D0);
	PlistCopy(D0_ch3, D0);
	PlistCopy(antiD0_ch1, antiD0);
	PlistCopy(antiD0_ch2, antiD0);
	PlistCopy(antiD0_ch3, antiD0);

	fitD0(D0);
	fitD0(antiD0);

	std::vector<Particle> DS_p;
	std::vector<Particle> DS_m;

	// combination, but without any mass cut
	// (we will make |M_DS - M_D| cuts later)
	combination(DS_p, Ptype(413), D0, pi_p_wo_SVD_cuts);
	combination(DS_m, Ptype(-413), antiD0, pi_m_wo_SVD_cuts);

	fitDS(DS_p);
	fitDS(DS_m);

	// Cut on the mass difference between D* and its first daughter - D0
	if (channelNum == 1) {
		withMassDifCut(DS_p, DSD0_MASS_DIFF_LOW, DSD0_MASS_DIFF_HIGH, 0);
		withMassDifCut(DS_m, DSD0_MASS_DIFF_LOW, DSD0_MASS_DIFF_HIGH, 0);
	} else {
		withMassDifCut(DS_p, DSD0_TIGHT_MASS_DIFF_LOW, DSD0_TIGHT_MASS_DIFF_HIGH, 0);
		withMassDifCut(DS_m, DSD0_TIGHT_MASS_DIFF_LOW, DSD0_TIGHT_MASS_DIFF_HIGH, 0);
	}

	std::vector<Particle> rho_p;
	std::vector<Particle> rho_m;

	combination(rho_p, Ptype(213), pi0, pi_p, RHO_MASS_WINDOW_LOW, RHO_MASS_WINDOW_HIGH);
	combination(rho_m, Ptype(-213), pi0, pi_m, RHO_MASS_WINDOW_LOW, RHO_MASS_WINDOW_HIGH);

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

	int trueDecayWasGenerated;

	// Set MC truth if running on MC
	if (dataType) {
		setMCtruth(B);
		trueDecayWasGenerated = trueDecayGenerated(channelNum);
		if (trueDecayWasGenerated == 1) {
			eventsWithRightDecay++;
		} else {
			//printf("*** DBG: Wrong decay in Exp %i, Run %i, Evt %i\n", expNo, runNo, evtNo);
			eventsWithOutRightDecay++;
		}
	}
	// Loop over all B candidates, keep only candidates passing Mbc, DeltaE, rho pi chi^2
	// and continuum suppression cuts
	applyCandidateCuts(B, continuumSupression, eventsDiscardedByDEandMbcCuts, eventsDiscardedByCSCut);

	// Cut on the number of candidates, because the chance to pick the
	// right one decreases with their number
	if (B.size() > MAX_NUM_B_CANDIDATES) {
		eventsDiscardedByNumCandidatesCut++;
		return;
	}

	// Select the best B candidate, store its index number in the candidate list
	// to a variable and save it to the best candidate ntuple
	if (B.size() > 0) {
		unsigned indexOfBestCandidate;
		Particle& Bcand = selectBestCandidate2(indexOfBestCandidate, B);
		if (dataType) {
			setBCSMetadata(trueDecayWasGenerated, indexOfBestCandidate, B,
					eventsWithRightCandidate, eventsWithOutRightCandidate, numCorrectlySelectedCandidates);
		}

		dynamic_cast<UserInfo&>(Bcand.userInfo()).setNumCandidates(B.size());

//		setBRecVertex(Bcand);
		setBRecVertexKFitter(Bcand);
		setBAscVertex(Bcand);

		tag(Bcand);

		saveToTuple(Bcand, tupleB);

//		int DSFlag = getMCtruthFlag(Bcand.child(0));
//		int D0Flag = getMCtruthFlag(Bcand.child(0).child(0));
//		int rhoFlag = getMCtruthFlag(Bcand.child(1));
//		int candSel = dynamic_cast<const UserInfo&>(Bcand.userInfo()).getCandidateSelection();
//
//		if ((int(DSFlag) == 1 || int(DSFlag) == 10) &&
//			(int(D0Flag == 1) || int(D0Flag == 10)) &&
//			(int(rhoFlag) == 3) && (candSel == 0) &&
//			dynamic_cast<const UserInfo&>(Bcand.userInfo()).getDeltaE() < -0.13) {
//				printEvent(Bcand);
//		}

	}
}

// Particle Bcand cannot be const, because Belle functions like IDHEP can't handle it
void DSRhoModule::saveToTuple(Particle Bcand, BelleTuple* tuple) {
	tuple->column("expmc", expMC);
	tuple->column("detver", detVer);
	tuple->column("expno", expNo);
	tuple->column("evtno", evtNo);
	tuple->column("runno", runNo);

	tuple->column("benergy", myCMSE);

	tuple->column("brecflavor", (Bcand.lund() > 0 ? 1 : -1));

	tuple->column("nocand", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getNumCandidates());

	tuple->column("mbc", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getMbc());
	tuple->column("de", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getDeltaE());

	tuple->column("dsmass", Bcand.child(0).mass());
	tuple->column("dspcms", pStar(Bcand.child(0)).vect().mag());

	tuple->column("rhomass", Bcand.child(1).mass());
	tuple->column("rhopcms", pStar(Bcand.child(1)).vect().mag());
	tuple->column("rhoplab", Bcand.child(1).ptot());

	tuple->column("d0chi2", dynamic_cast<const UserInfo&>(Bcand.child(0).child(0).userInfo()).getMassConstrainedChi2());
	tuple->column("d0mass", Bcand.child(0).child(0).mass());
	tuple->column("d0massbf", dynamic_cast<const UserInfo&>(Bcand.child(0).child(0).userInfo()).getMassBeforeVertexFit());
	tuple->column("d0pcms", pStar(Bcand.child(0).child(0)).vect().mag());

	if (dynamic_cast<const UserInfo&>(Bcand.child(0).child(0).userInfo()).decayMode() == 2) {
		tuple->column("d0pi0chi", dynamic_cast<const UserInfo&>(Bcand.child(0).child(0).child(2).userInfo()).getMassConstrainedChi2());
	} else {
		tuple->column("d0pi0chi", -1);
	}

	tuple->column("spipcms", pStar(Bcand.child(0).child(1)).vect().mag());
	tuple->column("spiplab", Bcand.child(0).child(1).ptot());

	tuple->column("dsd0diff", Bcand.child(0).mass() - Bcand.child(0).child(0).mass());

	tuple->column("dspiinv", invariantMass(Bcand.child(0), Bcand.child(1).child(1)));

	tuple->column("pi0mass", dynamic_cast<const UserInfo&>(Bcand.child(1).child(0).userInfo()).getMassBeforeVertexFit());
	tuple->column("pi0chi2", dynamic_cast<const UserInfo&>(Bcand.child(1).child(0).userInfo()).getMassConstrainedChi2());
	tuple->column("pi0decangle", dynamic_cast<const UserInfo&>(Bcand.child(1).child(0).userInfo()).getDecayAngle());
	tuple->column("e9oe25g1", dynamic_cast<const UserInfo&>(Bcand.child(1).child(0).userInfo()).getMinE9E25());
	tuple->column("e9oe25g2", dynamic_cast<const UserInfo&>(Bcand.child(1).child(0).userInfo()).getMaxE9E25());
	tuple->column("energyg1", dynamic_cast<const UserInfo&>(Bcand.child(1).child(0).userInfo()).getMinEnergy());
	tuple->column("energyg2", dynamic_cast<const UserInfo&>(Bcand.child(1).child(0).userInfo()).getMaxEnergy());
	tuple->column("thetag1", dynamic_cast<const UserInfo&>(Bcand.child(1).child(0).userInfo()).getMinTheta());
	tuple->column("thetag2", dynamic_cast<const UserInfo&>(Bcand.child(1).child(0).userInfo()).getMaxTheta());

	tuple->column("rhopidr", dynamic_cast<const UserInfo&>(Bcand.child(1).child(1).userInfo()).getDr());
	tuple->column("rhopidz", dynamic_cast<const UserInfo&>(Bcand.child(1).child(1).userInfo()).getDz());

	tuple->column("spidr", dynamic_cast<const UserInfo&>(Bcand.child(0).child(1).userInfo()).getDr());
	tuple->column("spidz", dynamic_cast<const UserInfo&>(Bcand.child(0).child(1).userInfo()).getDz());

	tuple->column("rhopipcms", pStar(Bcand.child(1).child(1)).vect().mag());
	tuple->column("pi0pcms", pStar(Bcand.child(1).child(0)).vect().mag());

	tuple->column("rhopinrf", Bcand.child(1).child(1).mdstCharged().trk().mhyp(2).nhits(3));
	tuple->column("rhopinz", Bcand.child(1).child(1).mdstCharged().trk().mhyp(2).nhits(4));

	tuple->column("spinrf", Bcand.child(0).child(1).mdstCharged().trk().mhyp(2).nhits(3));
	tuple->column("spinz", Bcand.child(0).child(1).mdstCharged().trk().mhyp(2).nhits(4));

	tuple->column("thetat", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getThetaT());
	tuple->column("phit", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getPhiT());
	tuple->column("thetab", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getThetaB());

	tuple->column("thetatg", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getThetaTG());
	tuple->column("phitg", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getPhiTG());
	tuple->column("thetabg", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getThetaBG());

	tuple->column("shcosthb", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getCosThetaB());
	tuple->column("shcosttb", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getCosThetaTB());
	tuple->column("shcostht", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getCosThetaT());

	tuple->column("mdspi", twoParticleInvariantMass(Bcand.child(0), Bcand.child(1).child(1)));
	tuple->column("mdspi0", twoParticleInvariantMass(Bcand.child(0), Bcand.child(1).child(0)));

	const vfit_info vtxRec = dynamic_cast<UserInfo&>(Bcand.userInfo()).getVertexRec();
	tuple->column("vrusable", vtxRec.m_usable);
	tuple->column("vrvtxx", vtxRec.m_pos.x());
	tuple->column("vrvtxy", vtxRec.m_pos.y());
	tuple->column("vrvtxz", vtxRec.m_pos.z());
	tuple->column("vrerr1", vtxRec.m_err[0][0]);
	tuple->column("vrerr2", vtxRec.m_err[1][0]);
	tuple->column("vrerr3", vtxRec.m_err[1][1]);
	tuple->column("vrerr4", vtxRec.m_err[2][0]);
	tuple->column("vrerr5", vtxRec.m_err[2][1]);
	tuple->column("vrerr6", vtxRec.m_err[2][2]);
	tuple->column("vrchi2", vtxRec.m_chi2);
	tuple->column("vreffxi", vtxRec.m_effxi);
	tuple->column("vrndf", vtxRec.m_ndf);
	tuple->column("vreffndf", vtxRec.m_effndf);
	tuple->column("vrntrk", vtxRec.m_ntrk);

	// Save generated vertices, if genHepevt link exists
	if (Bcand.genHepevt()) {
		if (abs(Bcand.genHepevt().idhep()) != 511) {
			tuple->column("vrgexist", 0);
			tuple->column("vrgvtxx", 100);
			tuple->column("vrgvtxy", 100);
			tuple->column("vrgvtxz", 100);
		} else {
			tuple->column("vrgexist", 1);
			tuple->column("vrgvtxx", Bcand.child(0).genHepevt().VX());
			tuple->column("vrgvtxy", Bcand.child(0).genHepevt().VY());
			tuple->column("vrgvtxz", Bcand.child(0).genHepevt().VZ());
		}
	} else {
		tuple->column("vrgexist", 0);
		tuple->column("vrgvtxx", 100);
		tuple->column("vrgvtxy", 100);
		tuple->column("vrgvtxz", 100);
	}

	// Save generated vertices for tag-side, if genHepevt link exists
	// Also save the btag mclink
	if (Bcand.genHepevt()) {
		Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
		Gen_hepevt* genBtag;
		const int BcandID = Bcand.genHepevt().get_ID();
		if (Bcand.genHepevt().mother()) {
			const int UpsilonFirstDaughterID = Bcand.genHepevt().mother().daFirst();
			const int UpsilonLastDaughterID = Bcand.genHepevt().mother().daLast();
			if (UpsilonLastDaughterID - UpsilonFirstDaughterID != 1) {
				printf("WARNING: Bcand's mother has more than two children!\n");
			}
			if (BcandID == UpsilonFirstDaughterID) {
				genBtag = &GenMgr(Panther_ID(UpsilonLastDaughterID));
			} else {
				genBtag = &GenMgr(Panther_ID(UpsilonFirstDaughterID));
			}

			Gen_hepevt* genBTagDaughter = &GenMgr(Panther_ID(genBtag->daFirst()));

			tuple->column("vtgexist", 1);
			tuple->column("vtgvtxx", genBTagDaughter->VX());
			tuple->column("vtgvtxy", genBTagDaughter->VY());
			tuple->column("vtgvtxz", genBTagDaughter->VZ());
			tuple->column("btagmclink", genBtag->idhep());
		} else {
			tuple->column("vtgexist", 0);
			tuple->column("vtgvtxx", 100);
			tuple->column("vtgvtxy", 100);
			tuple->column("vtgvtxz", 100);
			tuple->column("btagmclink", 0);
		}
	} else {
		tuple->column("vtgexist", 0);
		tuple->column("vtgvtxx", 100);
		tuple->column("vtgvtxy", 100);
		tuple->column("vtgvtxz", 100);
		tuple->column("btagmclink", 0);

	}


	vfit_info vtxTag = dynamic_cast<UserInfo&>(Bcand.userInfo()).getVertexTag();
	tuple->column("vtusable", vtxTag.m_usable);
	tuple->column("vtvtxx", vtxTag.m_pos.x());
	tuple->column("vtvtxy", vtxTag.m_pos.y());
	tuple->column("vtvtxz", vtxTag.m_pos.z());
	tuple->column("vterr1", vtxTag.m_err[0][0]);
	tuple->column("vterr2", vtxTag.m_err[1][0]);
	tuple->column("vterr3", vtxTag.m_err[1][1]);
	tuple->column("vterr4", vtxTag.m_err[2][0]);
	tuple->column("vterr5", vtxTag.m_err[2][1]);
	tuple->column("vterr6", vtxTag.m_err[2][2]);
	tuple->column("vtchi2", vtxTag.m_chi2);
	tuple->column("vtndf", vtxTag.m_ndf);
	tuple->column("vtntrk", vtxTag.m_ntrk);
	tuple->column("vtistagl", vtxTag.m_tmp[1]);

	tuple->column("tagflavor", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getFlavor());
	tuple->column("tagqr", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getQr());
	tuple->column("tagwtag", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getWtag());

	(dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW()).usefinal(0);
	tuple->column("k0mm2",   dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().mm2());
	tuple->column("k0et",    dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().et());
	tuple->column("k0hso00", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(0, 0));
	tuple->column("k0hso01", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(0, 1));
	tuple->column("k0hso02", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(0, 2));
	tuple->column("k0hso03", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(0, 3));
	tuple->column("k0hso04", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(0, 4));
	tuple->column("k0hso10", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(1, 0));
	tuple->column("k0hso12", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(1, 2));
	tuple->column("k0hso14", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(1, 4));
	tuple->column("k0hso20", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(2, 0));
	tuple->column("k0hso22", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(2, 2));
	tuple->column("k0hso24", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(2, 4));
	tuple->column("k0hoo0",  dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hoo(0));
	tuple->column("k0hoo1",  dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hoo(1));
	tuple->column("k0hoo2",  dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hoo(2));
	tuple->column("k0hoo3",  dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hoo(3));
	tuple->column("k0hoo4",  dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hoo(4));

	(dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW()).usefinal(1);
	tuple->column("k1mm2",   dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().mm2());
	tuple->column("k1et",    dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().et());
	tuple->column("k1hso00", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(0, 0));
	tuple->column("k1hso01", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(0, 1));
	tuple->column("k1hso02", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(0, 2));
	tuple->column("k1hso03", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(0, 3));
	tuple->column("k1hso04", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(0, 4));
	tuple->column("k1hso10", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(1, 0));
	tuple->column("k1hso12", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(1, 2));
	tuple->column("k1hso14", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(1, 4));
	tuple->column("k1hso20", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(2, 0));
	tuple->column("k1hso22", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(2, 2));
	tuple->column("k1hso24", dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hso(2, 4));
	tuple->column("k1hoo0",  dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hoo(0));
	tuple->column("k1hoo1",  dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hoo(1));
	tuple->column("k1hoo2",  dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hoo(2));
	tuple->column("k1hoo3",  dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hoo(3));
	tuple->column("k1hoo4",  dynamic_cast<UserInfo&>(Bcand.userInfo()).getKSFW().Hoo(4));

	tuple->column("csbdtg", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getContSuppressionBDTG());
	tuple->column("csmlp", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getContSuppressionMLP());

	//		// D0 decay mode
	//		tuple->column("d0dm", dynamic_cast<UserInfo&>(Bcand.child(0).child(0).userInfo()).decayMode());
	//
	//		// save pion related info (pi from B decay)
	//		// total momentum in lab frame
	//		tuple->column("rhoptot", Bcand.child(1).ptot());
	//		tuple->column("rhopid", dynamic_cast<UserInfo&>(Bcand.child(1).userInfo()).getPid());
	//		tuple->column("rhoeid", dynamic_cast<UserInfo&>(Bcand.child(1).userInfo()).getEid());
	//		tuple->column("rhomuid", dynamic_cast<UserInfo&>(Bcand.child(1).userInfo()).getMuid());
	//		tuple->column("rhodr", dynamic_cast<UserInfo&>(Bcand.child(1).userInfo()).getDr());
	//		tuple->column("rhodz", dynamic_cast<UserInfo&>(Bcand.child(1).userInfo()).getDz());
	//
	//		tuple->column("pi0mass", dynamic_cast<UserInfo&>(Bcand.child(1).child(0).userInfo()).getMassBeforeVertexFit());

	// save MC truth info
	if (dataType) {
		tuple->column("candsel", dynamic_cast<const UserInfo&>(Bcand.userInfo()).getCandidateSelection());
		tuple->column("mclink", IDhep(Bcand));
		tuple->column("mcflag", getMCtruthFlag(Bcand));
		tuple->column("dsmclink", IDhep(Bcand.child(0)));
		tuple->column("dsmcflag", getMCtruthFlag(Bcand.child(0)));
		tuple->column("rhomclink", IDhep(Bcand.child(1)));
		tuple->column("rhomcflag", getMCtruthFlag(Bcand.child(1)));
		tuple->column("d0mclink", IDhep(Bcand.child(0).child(0)));
		tuple->column("d0mcflag", getMCtruthFlag(Bcand.child(0).child(0)));
		tuple->column("spimclink", IDhep(Bcand.child(0).child(1)));
		tuple->column("spimcflag", getMCtruthFlag(Bcand.child(0).child(1)));
		tuple->column("rhopimcflag", getMCtruthFlag(Bcand.child(1).child(1)));
		tuple->column("pi0mcflag", getMCtruthFlag(Bcand.child(1).child(0)));
		if (Bcand.child(0).child(0).nChildren() == 3) {
			if (Bcand.child(0).child(0).child(2).lund() == 111) {
				tuple->column("dpi0mcflag", getMCtruthFlag(Bcand.child(0).child(0).child(2)));
			} else {
				tuple->column("dpi0mcflag", -999);
			}
		} else {
			tuple->column("dpi0mcflag", -999);
		}
		tuple->column("evmcflag", getEventMCFlag(dynamic_cast<const UserInfo&>(Bcand.userInfo()).getCandidateSelection(),
				getMCtruthFlag(Bcand), getMCtruthFlag(Bcand.child(0)), getMCtruthFlag(Bcand.child(1)),
				getMCtruthFlag(Bcand.child(0).child(0))));

	} else {
		tuple->column("candsel", 0.0);
		tuple->column("mclink", 0.0);
		tuple->column("mcflag", 0.0);
		tuple->column("dsmclink", 0.0);
		tuple->column("dsmcflag", 0.0);
		tuple->column("rhomclink", 0.0);
		tuple->column("rhomcflag", 0.0);
		tuple->column("d0mclink", 0.0);
		tuple->column("d0mcflag", 0.0);
		tuple->column("spimclink", 0.0);
		tuple->column("spimcflag", 0.0);
		tuple->column("rhopimcflag", 0.0);
		tuple->column("pi0mcflag", 0.0);
		tuple->column("dpi0mcflag", 0.0);
		tuple->column("evmcflag", 0.0);
	}
	tuple->dumpData();
}

void DSRhoModule::printEvent(Particle Bcand) {
	Gen_hepevt_Manager& EvtMgr = Gen_hepevt_Manager::get_manager();
	HepAList<Gen_hepevt> PartList;
	for (std::vector<Gen_hepevt>::iterator it = EvtMgr.begin();
			it != EvtMgr.end(); it++) {
		Gen_hepevt& event = *it;
		PartList.append(event);
	}
	if (PartList.last()) {
		pt->print(PartList, expNo, runNo, evtNo);
	}
	pt->printParticleChildren(Bcand);

//	for (std::vector<Particle>::iterator it = pi_m.begin(); it != pi_m.end(); it++) {
//		pt.printParticle(*it);
//		pt.printPiParticle(*it);
//	}
}

void DSRhoModule::fitD0(std::vector<Particle> &particles) {
	for (unsigned i = 0; i < particles.size(); ++i) {
		kmassvertexfitter kf;
		for (unsigned j = 0; j < particles[i].nChildren(); j++){
			addTrack2fit(kf, particles[i].child(j));
		}
		kf.invariantMass(Ptype("D0").mass());
		unsigned err = kf.fit();
		if (err == 0) {
			dynamic_cast<UserInfo&>(particles[i].userInfo()).setMassBeforeVertexFit(particles[i].mass());
			makeMother(kf, particles[i]);
			dynamic_cast<UserInfo&>(particles[i].userInfo()).setMassConstrainedChi2(kf.chisq());
		} else {
			HepPoint3D vtx(999., 999., 999.);
			HepSymMatrix errVtx(3, 0);
			particles[i].momentum().decayVertex(vtx, errVtx);
		}
	}

}

void DSRhoModule::fitDS(std::vector<Particle> &particles) {
	for (unsigned i = 0; i < particles.size(); ++i) {
		kvertexfitter kf;
		for (unsigned j = 0; j < particles[i].nChildren(); j++){
			addTrack2fit(kf, particles[i].child(j));
		}
		unsigned err = kf.fit();
		if (err == 0) {
			makeMother(kf, particles[i]);
		} else {
			HepPoint3D vtx(999., 999., 999.);
			HepSymMatrix errVtx(3, 0);
			particles[i].momentum().decayVertex(vtx, errVtx);
		}
	}

}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
