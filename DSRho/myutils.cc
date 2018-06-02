// ******************************************************************
// Utilty functions
// A. Zupanc jun-06
//
// ******************************************************************
#include <stdio.h>

#include "belle.h"
#include "belleCLHEP/Vector/Rotation.h"
#include "benergy/BeamEnergy.h"
#include "eid/eid.h"
#include "exkfitter/ExKFitter.h"
#include "hamlet/Hamlet.h"
#include "helix/Helix.h"
#include HEPEVT_H
#include "ip/IpProfile.h"
#include "kid/atc_pid.h"
#include "mdst/mdst.h"
#include "mdst/Muid_mdst.h"
#include "mdst/findKs.h"
#include MDST_H
#include "panther/panther.h"
#include "particle/utility.h"
#include "tagv/TagV.h"
#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"

#include "cstdlib"
#include "cmath"

#include "myutils.h"
#include "constants.h"
#include "UserInfo.h"
#include "geninfo.h"
#include "Continuum.h"


#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


//#include <algorithm>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//======
// Set proper error matrix for gamma. Note errCart length is 4.
// Copied from FindGamma in icpv_skim package.
//======
void setGammaError(Particle & gamma) {
	HepSymMatrix errEcl(3, 0); // 3x3 initialize to zero
	errEcl[0][0] = gamma.mdstGamma().ecl().error(0); // Energy
	errEcl[1][0] = gamma.mdstGamma().ecl().error(1);
	errEcl[1][1] = gamma.mdstGamma().ecl().error(2); // Phi
	errEcl[2][0] = gamma.mdstGamma().ecl().error(3);
	errEcl[2][1] = gamma.mdstGamma().ecl().error(4);
	errEcl[2][2] = gamma.mdstGamma().ecl().error(5); // Theta

	double cp = cos(gamma.mdstGamma().ecl().phi());
	double sp = sin(gamma.mdstGamma().ecl().phi());
	double ct = cos(gamma.mdstGamma().ecl().theta());
	double st = sin(gamma.mdstGamma().ecl().theta());
	double E = gamma.mdstGamma().ecl().energy();

	HepMatrix jacobian(4, 3, 0);
	jacobian[0][0] = cp * st;
	jacobian[0][1] = -E * sp * st;
	jacobian[0][2] = E * cp * ct;
	jacobian[1][0] = sp * st;
	jacobian[1][1] = E * cp * st;
	jacobian[1][2] = E * sp * ct;
	jacobian[2][0] = ct;
	jacobian[2][1] = 0.0;
	jacobian[2][2] = -E * st;
	jacobian[3][0] = 1.0;
	jacobian[3][1] = 0.0;
	jacobian[3][2] = 0.0;

	HepSymMatrix errCart = errEcl.similarity(jacobian);

	const HepPoint3D origin; // Default constructor sets (0,0,0).
	static double largeError = 1.0; // 1.0*1.0cm^2.
	HepSymMatrix dx(3, 0);
	dx[0][0] = largeError;
	dx[1][1] = largeError;
	dx[2][2] = largeError;
	// Convert Mdst_gamma into Particle object.

	gamma.momentum().momentum(gamma.p(), errCart);
	gamma.momentum().position(origin, dx);
	return;
}

void setGammaError(std::vector<Particle> &p) {
	for (std::vector<Particle>::iterator i = p.begin(); i != p.end(); ++i)
		setGammaError(*i);
}

/************ Set error matrix (dpx) for pi0 ***************************/
void setPi0Error(Particle &p) {
	if (p.nChildren() != 2)
		return;
	if (!p.child(0).mdstGamma() || !p.child(1).mdstGamma())
		return;
	for (unsigned i = 0; i < 2; i++)
		setGammaError(p.child(i));

	kmassvertexfitter kmv;
	kmv.invariantMass(p.pType().mass());
	for (unsigned i = 0; i < 2; ++i)
		addTrack2fit(kmv, p.child(i));
	int err = kmv.fit();
	if (!err) {
		makeMother(kmv, p);
		p.momentum().vertex(kmv.vertex(), kmv.errVertex());
	}
}

void setPi0Error(std::vector<Particle> &p) {
	for (std::vector<Particle>::iterator i = p.begin(); i != p.end(); ++i)
		setPi0Error(*i);
}



// fills K+/K- and Pi+/Pi- vectors with +/- tracks
// no PID cuts
// saves the impact parameters dr & dz in userInfo every particle (see UserInfo.h)
// saves eid pid and muid in userInfo of every particle
// (the PID info is saved to userInfo so that it doesn't need to be recalculated later)
void fillKaonsAndPions(std::vector<Particle> &kp,
		std::vector<Particle> &km,
		std::vector<Particle> &pi_p,
		std::vector<Particle> &pi_m) {

	// Loop over all charged tracks
	Mdst_charged_Manager& chgMgr = Mdst_charged_Manager::get_manager();
	for(std::vector<Mdst_charged>::iterator it=chgMgr.begin(); it!=chgMgr.end(); it++) {

		eid sel_e(*it);
		Muid_mdst Muid(*it);
		double eid = sel_e.prob(3,-1,5);
		double muid= Muid.Muon_likelihood();

		double im_dr, im_dz, pid, protonIDkaon, protonIDpion;

		// prob(K:pi)
		pid = atc_pid(3,1,5,3,2).prob(*it);
		// prob(K:p)
		protonIDkaon = atc_pid(3,1,5,3,4).prob(*it);
		// prob(K:p)
		protonIDpion = atc_pid(3,1,5,2,4).prob(*it);

		Ptype ptype_use;
		if(it->charge() == 1) {
			// PI+
			ptype_use = Ptype("PI+");

			Particle pi(*it,ptype_use);

			impact_param(pi, im_dr, im_dz);
			saveTrackInfo(pi, im_dr, im_dz, eid, pid, muid, protonIDpion);
			pi_p.push_back(pi);

			// K+
			ptype_use = Ptype("K+");

			Particle k(*it,ptype_use);

			impact_param(k, im_dr, im_dz);
			saveTrackInfo(k, im_dr, im_dz, eid, pid, muid, protonIDkaon);
			kp.push_back(k);

		} else if(it->charge() == -1) {
			// PI-
			ptype_use = Ptype("PI-");

			Particle pi(*it,ptype_use);

			impact_param(pi, im_dr, im_dz);
			saveTrackInfo(pi, im_dr, im_dz, eid, pid, muid, protonIDpion);
			pi_m.push_back(pi);

			// K-
			ptype_use = Ptype("K-");

			Particle k(*it,ptype_use);

			impact_param(k, im_dr, im_dz);
			saveTrackInfo(k, im_dr, im_dz, eid, pid, muid, protonIDkaon);
			km.push_back(k);
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// Ks candidates are filled from MdstVee2 bank
// "good Kshort" selection is applied
void fillGoodKShort(std::vector<Particle> &K0short) {
	Mdst_vee2_Manager& VeeMgr = Mdst_vee2_Manager::get_manager();

	for(std::vector<Mdst_vee2>::const_iterator iv=VeeMgr.begin(); iv!=VeeMgr.end(); iv++) {
		// Is it a K_short, or another V-particle ?
		if (iv->kind() != 1)
			continue;

		// Check for daughters
		if(!iv->chgd(0) || !iv->chgd(1))
			continue;

		Particle Kshort(*iv);

		// use Fang Fang's Kshort finder to check K-short quality
		FindKs KSfinder;
		// 1 as an argument in IpProfile::position = EVENT_DEPENDENT_IP
		KSfinder.candidates(*iv, IpProfile::position(1));
		int goodKsFlag = KSfinder.goodKs();

		if(!goodKsFlag)
			continue; // not a good Ks candidate so we skip it

		K0short.push_back(Kshort);
	}
}


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Loops over Mdst_pi0 panther table and selects Pi0 candidates, with
// mass between 0.115 and 0.150 MeV and energy of both photons being larger than 50 MeV
// it also saves few addition information to the userInfo
void fillAllPi0(std::vector<Particle> &pi0) {

	Mdst_pi0_Manager &mgrPi0 = Mdst_pi0_Manager::get_manager();

	for (Mdst_pi0_Manager::iterator it = mgrPi0.begin(); it != mgrPi0.end(); it++) {
		Particle pizero(*it);

		// make the selection
		if(it->mass() < PI0_MASS_CUT_LOW || it->mass() > PI0_MASS_CUT_HIGH ||
				pizero.child(0).e() < GAMMA_E_CUT_BARREL || pizero.child(1).e() < GAMMA_E_CUT_BARREL)
			continue;


		// Min and Max E9/E25
		float e9oe25_g1 = -999;
		float e9oe25_g2 = -999;

		// Min and Max theta
		float theta_g1 = -999;
		float theta_g2 = -999;

		Mdst_ecl_aux_Manager& aux = Mdst_ecl_aux_Manager::get_manager();
		if(pizero.child(0).mdstGamma()) {
			Mdst_ecl_aux & shower_aux(aux(Panther_ID(pizero.child(0).mdstGamma().ecl().get_ID())));
			e9oe25_g1 = (float)shower_aux.e9oe25();

			theta_g1 = pizero.child(0).mdstGamma().ecl().theta();
		}

		if(pizero.child(1).mdstGamma()) {
			Mdst_ecl_aux & shower_aux(aux(Panther_ID(pizero.child(1).mdstGamma().ecl().get_ID())));
			e9oe25_g2 = (float)shower_aux.e9oe25();

			theta_g2 = pizero.child(1).mdstGamma().ecl().theta();
		}

		// Min and Max energy
		float energy_g1 = pizero.child(0).e();
		float energy_g2 = pizero.child(1).e();

		// Further selection; it's here, because you need theta from mdstGamma
		if((theta_g1 < FWD_ENDCAP_ANGLE || theta_g1 > BKW_ENDCAP_ANGLE) && energy_g1 < GAMMA_E_CUT_ENDCAP)
			continue;

		if((theta_g2 < FWD_ENDCAP_ANGLE || theta_g2 > BKW_ENDCAP_ANGLE) && energy_g2 < GAMMA_E_CUT_ENDCAP)
			continue;

		// Angle between mother momentum in lab frame and child momentum in mother's rest frame
		float decAng = -999;

		float beta = pizero.momentum().p().boostVector().mag();
		decAng = fabs(1./beta * (energy_g1 - energy_g2)/(energy_g1 + energy_g2));

		// creates the userinfor for this pi0 candidate
		setUserInfo(pizero);

		// saves the pi0 mass before mass constrained fit
		dynamic_cast<UserInfo&>(pizero.userInfo()).setMassBeforeVertexFit(it->mass());

		// saves the pi0 mass-constrained chi2
		dynamic_cast<UserInfo&>(pizero.userInfo()).setMassConstrainedChi2(it->chisq());

		savePi0Info(pizero, e9oe25_g1, e9oe25_g2, energy_g1, energy_g2, theta_g1, theta_g2, decAng);

		pi0.push_back(pizero);
	}
}


/////////////////////////////////////////////////////////////////////////////////
// save basic track info in particle userInfo
void saveTrackInfo(Particle &p, double dr, double dz, double eid, double pid, double muid, double protonID) {
	// check if the particleUserinfo already exists
	if(! &p.userInfo() ) setUserInfo(p);

	dynamic_cast<UserInfo&>(p.userInfo()).setType(p.lund());
	dynamic_cast<UserInfo&>(p.userInfo()).setDr(dr);
	dynamic_cast<UserInfo&>(p.userInfo()).setDz(dz);
	dynamic_cast<UserInfo&>(p.userInfo()).setEid(eid);
	dynamic_cast<UserInfo&>(p.userInfo()).setPid(pid);
	dynamic_cast<UserInfo&>(p.userInfo()).setMuid(muid);
	dynamic_cast<UserInfo&>(p.userInfo()).setProtonID(protonID);
}

// save basic pi0 info in particle userInfo
void savePi0Info(Particle &p, double min_e9oe25, double max_e9oe25, double min_e, double max_e,
		double min_theta, double max_theta, double decAng) {
	if(! &p.userInfo() ) setUserInfo(p);

	dynamic_cast<UserInfo&>(p.userInfo()).setType(p.lund());
	dynamic_cast<UserInfo&>(p.userInfo()).setMinE9E25(min_e9oe25);
	dynamic_cast<UserInfo&>(p.userInfo()).setMaxE9E25(max_e9oe25);
	dynamic_cast<UserInfo&>(p.userInfo()).setMinEnergy(min_e);
	dynamic_cast<UserInfo&>(p.userInfo()).setMaxEnergy(max_e);
	dynamic_cast<UserInfo&>(p.userInfo()).setMinTheta(min_theta);
	dynamic_cast<UserInfo&>(p.userInfo()).setMaxTheta(max_theta);
	dynamic_cast<UserInfo&>(p.userInfo()).setDecayAngle(decAng);
}



// saves decay mode identifier (integer) to the userInfo
void setDecayMode(std::vector<Particle> &list, int dMode) {
	for (int i = 0; i < (int)list.size(); ++i) {
		if(! &list[i].userInfo() ) setUserInfo(list[i]);
		dynamic_cast<UserInfo&>(list[i].userInfo()).decayMode(dMode);
	}
}

// copies Particles from plist to copy_plist
void PlistCopy(std::vector<Particle> &plist, std::vector<Particle> &copy_plist)
{

	for(int i=0;i<(int)plist.size();++i) {
		copy_plist.push_back(plist[i]);
	}

}

void impact_param(Particle &p, double &dr, double &dz){

	const Mdst_charged charged(p.mdstCharged());

	double thisMass = p.mass();

	int hyp = 4;
	if(thisMass < 0.005){ // e = 0.000511
		hyp = 0;
	}else if(thisMass < 0.110){ // mu = 0.1056
		hyp = 1;
	}else if(thisMass < 0.200){ // pi = 0.13956
		hyp = 2;
	}else if(thisMass < 0.5){ // K = 0.4936
		hyp = 3;
	}

	const HepPoint3D pivot(charged.trk().mhyp(hyp).pivot_x(),
			charged.trk().mhyp(hyp).pivot_y(),
			charged.trk().mhyp(hyp).pivot_z());

	HepVector  a(5);
	a[0] = charged.trk().mhyp(hyp).helix(0);
	a[1] = charged.trk().mhyp(hyp).helix(1);
	a[2] = charged.trk().mhyp(hyp).helix(2);
	a[3] = charged.trk().mhyp(hyp).helix(3);
	a[4] = charged.trk().mhyp(hyp).helix(4);
	HepSymMatrix Ea(5,0);
	Ea[0][0] = charged.trk().mhyp(hyp).error(0);
	Ea[1][0] = charged.trk().mhyp(hyp).error(1);
	Ea[1][1] = charged.trk().mhyp(hyp).error(2);
	Ea[2][0] = charged.trk().mhyp(hyp).error(3);
	Ea[2][1] = charged.trk().mhyp(hyp).error(4);
	Ea[2][2] = charged.trk().mhyp(hyp).error(5);
	Ea[3][0] = charged.trk().mhyp(hyp).error(6);
	Ea[3][1] = charged.trk().mhyp(hyp).error(7);
	Ea[3][2] = charged.trk().mhyp(hyp).error(8);
	Ea[3][3] = charged.trk().mhyp(hyp).error(9);
	Ea[4][0] = charged.trk().mhyp(hyp).error(10);
	Ea[4][1] = charged.trk().mhyp(hyp).error(11);
	Ea[4][2] = charged.trk().mhyp(hyp).error(12);
	Ea[4][3] = charged.trk().mhyp(hyp).error(13);
	Ea[4][4] = charged.trk().mhyp(hyp).error(14);
	Helix helix(pivot, a, Ea);

	const Hep3Vector&   IP     = IpProfile::position(1);
	if (IP.mag())
		helix.pivot(IP);

	dr = helix.dr();
	dz = helix.dz();
}

void withImpactPars(std::vector<Particle> &list, const double max_dr, const double max_dz){
	for(int i=0; i<(int)list.size(); ++i){
		if(fabs(dynamic_cast<UserInfo&>(list[i].userInfo()).getDr()) > max_dr ||
				fabs(dynamic_cast<UserInfo&>(list[i].userInfo()).getDz()) > max_dz){
			list.erase(list.begin()+i);
			--i;
		}
	}
}

double invariantMass(const Particle& p1, const Particle& p2) {
	HepLorentzVector sumOfMomenta = p1.p() + p2.p();
	return sumOfMomenta.mag();
}

Particle& selectBestCandidate(unsigned& index, std::vector<Particle>& B){
	double smallestMbcDelta = 0;
	unsigned indexOfSmallest = 0;
	double currentMbcDelta = 0;

	for(unsigned i=0; i<B.size(); i++) {
		currentMbcDelta = fabs(MBC - dynamic_cast<const UserInfo&>(B[i].userInfo()).getMbc());
		if(i == 0){
			smallestMbcDelta = currentMbcDelta;
		} else if(smallestMbcDelta > currentMbcDelta){
			smallestMbcDelta = currentMbcDelta;
			indexOfSmallest = i;
		}
	}
	index = indexOfSmallest;
	return B[indexOfSmallest];
}

Particle& selectBestCandidate2(unsigned& index, std::vector<Particle>& B){
	double MbcChi2 = 0;
	double Mpi0Chi2 = 0;
	double smallestVal = 0;
	unsigned indexOfSmallest = 0;
	double currentVal = 0;

	for(unsigned i=0; i<B.size(); i++) {
		MbcChi2 = pow((MBC - dynamic_cast<const UserInfo&>(B[i].userInfo()).getMbc())/MBCSIGMA, 2);
		Mpi0Chi2 = pow((PI0MASS - dynamic_cast<UserInfo&>(B[i].child(1).child(0).userInfo()).getMassBeforeVertexFit())/PI0WIDTH, 2);
		currentVal = MbcChi2 + Mpi0Chi2;
		if(i == 0){
			smallestVal = currentVal;
		} else if(smallestVal > currentVal){
			smallestVal = currentVal;
			indexOfSmallest = i;
		}
	}
	index = indexOfSmallest;
	return B[indexOfSmallest];
}

Particle& selectRandomBestCandidate(unsigned& index, std::vector<Particle>& B){
	index = (rand()/RAND_MAX) * (B.size() - 1);
	return B[index];
}

bool isSignal(const Particle& B) {
	std::vector<Particle> particles;
	particles.push_back(B);
	particles.push_back(B.child(0)); //DS
	particles.push_back(B.child(0).child(0)); //D0
	particles.push_back(B.child(1)); //rho
	particles.push_back(B.child(1).child(0)); //pi0

	if (B.child(0).child(0).nChildren() == 3) {
		if (B.child(0).child(0).child(2).lund() == 111) {
			particles.push_back(B.child(0).child(0).child(2)); //pi0
		}
	}

	for (std::vector<Particle>::iterator it = particles.begin(); it != particles.end(); it++) {
		if (getMCtruthFlag(*it) != 1 && getMCtruthFlag(*it) != 10) return 0;
	}

	return 1;
}

bool findRightCandidate(unsigned& index, std::vector<Particle>& B){
	bool foundRightCandidate = false;
	for(unsigned i=0; i<B.size(); i++) {
		if (isSignal(B[i])) {
			index = i;
			foundRightCandidate = true;
			break;
		}
	}
	return foundRightCandidate;
}

int trueDecayGenerated(int channelNum) {
	Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
	Gen_hepevt* U4S = 0;
	Gen_hepevt* B01 = 0;
	Gen_hepevt* B02 = 0;
	Gen_hepevt* DS = 0;
	Gen_hepevt* rho = 0;
	Gen_hepevt* DSD0 = 0;
	Gen_hepevt* DSpi = 0;
	Gen_hepevt* D0K = 0;
	Gen_hepevt* D0pi1 = 0;
	Gen_hepevt* D0pi2 = 0;
	Gen_hepevt* D0pi3 = 0;
	Gen_hepevt* D0pi0 = 0;
	Gen_hepevt* rhopi = 0;
	Gen_hepevt* rhopi0 = 0;

	std::vector<Gen_hepevt**> U4SChildren;
	std::vector<int> U4SChildrenIdhep;
	U4SChildren.push_back(&B01);	U4SChildrenIdhep.push_back(IDHEP_B0);
	U4SChildren.push_back(&B02);	U4SChildrenIdhep.push_back(IDHEP_B0);

	std::vector<Gen_hepevt**> B0Children;
	std::vector<int> B0ChildrenIdhep;
	B0Children.push_back(&DS);	B0ChildrenIdhep.push_back(IDHEP_DS);
	B0Children.push_back(&rho);	B0ChildrenIdhep.push_back(IDHEP_RHO);

	std::vector<Gen_hepevt**> DSChildren;
	std::vector<int> DSChildrenIdhep;
	DSChildren.push_back(&DSD0);	DSChildrenIdhep.push_back(IDHEP_D0);
	DSChildren.push_back(&DSpi);	DSChildrenIdhep.push_back(IDHEP_PI);

	std::vector<Gen_hepevt**> D0Children;
	std::vector<int> D0ChildrenIdhep;
	if (channelNum == 1) {
		D0Children.push_back(&D0K);		D0ChildrenIdhep.push_back(IDHEP_K);
		D0Children.push_back(&D0pi1);	D0ChildrenIdhep.push_back(IDHEP_PI);
	} else if (channelNum == 2) {
		D0Children.push_back(&D0K);		D0ChildrenIdhep.push_back(IDHEP_K);
		D0Children.push_back(&D0pi1);	D0ChildrenIdhep.push_back(IDHEP_PI);
		D0Children.push_back(&D0pi0);	D0ChildrenIdhep.push_back(IDHEP_PI0);
	} else if (channelNum == 3) {
		D0Children.push_back(&D0K);		D0ChildrenIdhep.push_back(IDHEP_K);
		D0Children.push_back(&D0pi1);	D0ChildrenIdhep.push_back(IDHEP_PI);
		D0Children.push_back(&D0pi2);	D0ChildrenIdhep.push_back(IDHEP_PI);
		D0Children.push_back(&D0pi3);	D0ChildrenIdhep.push_back(IDHEP_PI);
	}

	std::vector<Gen_hepevt**> rhoChildren;
	std::vector<int> rhoChildrenIdhep;
	rhoChildren.push_back(&rhopi);	rhoChildrenIdhep.push_back(IDHEP_PI);
	rhoChildren.push_back(&rhopi0);	rhoChildrenIdhep.push_back(IDHEP_PI0);

	// TODO: Explanation
	std::vector<Gen_hepevt**> B0ChildrenNonResonant;
	std::vector<int> B0ChildrenIdhepNonResonant;
	B0ChildrenNonResonant.push_back(&DS);		B0ChildrenIdhepNonResonant.push_back(IDHEP_DS);
	B0ChildrenNonResonant.push_back(&rhopi);	B0ChildrenIdhepNonResonant.push_back(IDHEP_PI);
	B0ChildrenNonResonant.push_back(&rhopi0);	B0ChildrenIdhepNonResonant.push_back(IDHEP_PI0);


	if (GenMgr(Panther_ID(1)).idhep() == IDHEP_Y4S) {
		U4S = & GenMgr(Panther_ID(1));
	}

	// TODO: Rewrite this ugly thing; remove redundant and excessive calls to findChildren
	if (findChildren(U4S, U4SChildren, U4SChildrenIdhep) &&
			findChildren(B01, B0Children, B0ChildrenIdhep) &&
			findChildren(DS, DSChildren, DSChildrenIdhep) &&
			findChildren(DSD0, D0Children, D0ChildrenIdhep, true) &&
			findChildren(rho, rhoChildren, rhoChildrenIdhep)) {
		//printf("*** DBG: Exiting trueDecayGenerated with true (B01)\n");
		return 1;
	} else if (findChildren(U4S, U4SChildren, U4SChildrenIdhep) &&
			findChildren(B02, B0Children, B0ChildrenIdhep) &&
			findChildren(DS, DSChildren, DSChildrenIdhep) &&
			findChildren(DSD0, D0Children, D0ChildrenIdhep, true) &&
			findChildren(rho, rhoChildren, rhoChildrenIdhep)) {
		//	printf("*** DBG: Exiting trueDecayGenerated with true (B02)\n");
		return 1;
	} else if (findChildren(U4S, U4SChildren, U4SChildrenIdhep) &&
			findChildren(B01, B0ChildrenNonResonant, B0ChildrenIdhepNonResonant) &&
			findChildren(DS, DSChildren, DSChildrenIdhep) &&
			findChildren(DSD0, D0Children, D0ChildrenIdhep, true)) {
		//	printf("*** DBG: Exiting trueDecayGenerated with true (B02)\n");
		return 2;
	} else if (findChildren(U4S, U4SChildren, U4SChildrenIdhep) &&
			findChildren(B02, B0ChildrenNonResonant, B0ChildrenIdhepNonResonant) &&
			findChildren(DS, DSChildren, DSChildrenIdhep) &&
			findChildren(DSD0, D0Children, D0ChildrenIdhep, true)) {
		//	printf("*** DBG: Exiting trueDecayGenerated with true (B02)\n");
		return 2;
	} else {
		//printf("*** DBG: Exiting trueDecayGenerated with false because\n");
		//printf("*** DBG: findChildren(U4S, U4SChildren) = %i\n", findChildren(U4S, U4SChildren, U4SChildrenIdhep));
		//printf("*** DBG: findChildren(B01, B0Children) = %i\n", findChildren(B01, B0Children, B0ChildrenIdhep));
		//printf("*** DBG: findChildren(DS, DSChildren) = %i\n", findChildren(DS, DSChildren, DSChildrenIdhep));
		//printf("*** DBG: findChildren(DSD0, D0Children) = %i\n", findChildren(DSD0, D0Children, D0ChildrenIdhep));
		//printf("*** DBG: findChildren(rho, rhoChildren) = %i\n", findChildren(rho, rhoChildren, rhoChildrenIdhep));
		return 0;
	}
}

// Given a Gen_hepevt particle check that it has exactly the daughters given by idheps vector.
// It assigns the found children to the daughters vector elements.
// In principle, the next level - calling findChildren on the just found children should be done
// only when findChildren returns true, but you can call it even if it doesn't as it can deal
// with null pointers.
bool findChildren(Gen_hepevt* particle, std::vector<Gen_hepevt**>& daughters, std::vector<int>& idheps, bool lookAtFSPChildren) {
	if (particle == NULL) return false;
	Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();

	std::vector<Gen_hepevt*> actualDaughters;
	if (lookAtFSPChildren) {
		appendGenFSP(*particle, actualDaughters);
	} else {
		for (int i = particle->daFirst(); i <= particle->daLast(); i++) {
			actualDaughters.push_back(& GenMgr(Panther_ID(i)));
		}
	}

	bool* usedDaughters = new bool [daughters.size()];
	for (unsigned int i = 0; i < daughters.size(); i++) {
		usedDaughters[i] = false;
	}

	int foundGammas = 0;

	for (std::vector<Gen_hepevt*>::iterator child = actualDaughters.begin(); child != actualDaughters.end(); child++) {
		bool foundChild = false;

		if ((*child)->idhep() == IDHEP_GAMMA) {
			foundGammas++;
			if (lookAtFSPChildren) {
				// If this is the first gamma we found among the FSP, ignore it
				// If not, then it's not FSR and this is not signal
				if (foundGammas == 1) continue;
			// If we are not looking only at FSP, the gammas will be FSR
			} else {
				continue;
			}
		}

		for (unsigned int j = 0; j < daughters.size(); j++) {
			if (usedDaughters[j]) continue;

			if (idheps[j] == abs((*child)->idhep())) {
				*daughters[j] =  *child;
				usedDaughters[j] = true;
				foundChild = true;
				break;
			}
		}

		if (!foundChild) {
			// Found an extra particle
			return false;
		}
	}

	for (unsigned int i = 0; i < daughters.size(); i++) {
		if (usedDaughters[i] == false) {
			// Didn't find all children
			return false;
		}
	}

	// Successfully found all children
	return true;
}

// Save best candidate selection metadata to userInfo and to vars used for printing
// summary information in Module::term
void setBCSMetadata(int trueCandidateWasGenerated, const unsigned indexOfBestCandidate, std::vector<Particle>& B,
		unsigned& eventsWithRightCandidate, unsigned& eventsWithOutRightCandidate,
		unsigned& numCorrectlySelectedCandidates) {

	if (trueCandidateWasGenerated == 1) {
		unsigned int indexOfRightCandidate;
		// If the right (= true) candidate is present in the list of candidates
		if (findRightCandidate(indexOfRightCandidate, B)) {
			eventsWithRightCandidate++;
			if (indexOfBestCandidate == indexOfRightCandidate) {
				numCorrectlySelectedCandidates++;
				// True candidate present and selected
				dynamic_cast<UserInfo&>(B[indexOfBestCandidate].userInfo()).setCandidateSelection(4);
			} else {
				// True candidate present but not selected
				dynamic_cast<UserInfo&>(B[indexOfBestCandidate].userInfo()).setCandidateSelection(3);
			}
		} else {
			eventsWithOutRightCandidate++;
			// True candidate not present but true decay generated
			dynamic_cast<UserInfo&>(B[indexOfBestCandidate].userInfo()).setCandidateSelection(2);
		}
	} else if (trueCandidateWasGenerated == 2) {
		// Almost true decay generated, but with non-resonant pi pi0 (different angular PDF)
		dynamic_cast<UserInfo&>(B[indexOfBestCandidate].userInfo()).setCandidateSelection(1);
	} else {
		// True decay not generated
		dynamic_cast<UserInfo&>(B[indexOfBestCandidate].userInfo()).setCandidateSelection(0);
	}
}

void getTransversityAngles(const Particle& mother, const bool generatedValues, double& thetaT, double& phiT, double& thetaB){
	bool genPresent;

	// If all the relevant particles have genHepevt links, we can save the generated values
	if (mother.genHepevt() &&
			mother.child(0).genHepevt() && mother.child(0).child(0).genHepevt() && mother.child(0).child(1).genHepevt() &&
			mother.child(1).genHepevt() && mother.child(1).child(0).genHepevt() && mother.child(1).child(1).genHepevt()) {
		genPresent = true;
	} else {
		genPresent = false;
	}

	HepLorentzVector branchA[3];
	HepLorentzVector branchB[3];
	Hep3Vector CMSboost;

	if (generatedValues && genPresent) {
		branchA[0] = HepLorentzVector(
				mother.child(0).genHepevt().PX(),
				mother.child(0).genHepevt().PY(),
				mother.child(0).genHepevt().PZ(),
				mother.child(0).genHepevt().E());			// DS
		branchA[1] = HepLorentzVector(
				mother.child(0).child(0).genHepevt().PX(),
				mother.child(0).child(0).genHepevt().PY(),
				mother.child(0).child(0).genHepevt().PZ(),
				mother.child(0).child(0).genHepevt().E());	// DSD0
		branchA[2] = HepLorentzVector(
				mother.child(0).child(1).genHepevt().PX(),
				mother.child(0).child(1).genHepevt().PY(),
				mother.child(0).child(1).genHepevt().PZ(),
				mother.child(0).child(1).genHepevt().E());	// DSPi
		branchB[0] = HepLorentzVector(
				mother.child(1).genHepevt().PX(),
				mother.child(1).genHepevt().PY(),
				mother.child(1).genHepevt().PZ(),
				mother.child(1).genHepevt().E());			// Rho
		branchB[1] = HepLorentzVector(
				mother.child(1).child(0).genHepevt().PX(),
				mother.child(1).child(0).genHepevt().PY(),
				mother.child(1).child(0).genHepevt().PZ(),
				mother.child(1).child(0).genHepevt().E());	// RhoPi0
		branchB[2] = HepLorentzVector(
				mother.child(1).child(1).genHepevt().PX(),
				mother.child(1).child(1).genHepevt().PY(),
				mother.child(1).child(1).genHepevt().PZ(),
				mother.child(1).child(1).genHepevt().E());	// RhoPi

		HepLorentzVector motherP(
				mother.genHepevt().PX(),
				mother.genHepevt().PY(),
				mother.genHepevt().PZ(),
				mother.genHepevt().E());
		CMSboost = motherP.boostVector();
	} else {
		branchA[0] = mother.child(0).p();			// DS
		branchA[1] = mother.child(0).child(0).p();	// DSD0
		branchA[2] = mother.child(0).child(1).p();	// DSPi
		branchB[0] = mother.child(1).p();			// Rho
		branchB[1] = mother.child(1).child(0).p();	// RhoPi0
		branchB[2] = mother.child(1).child(1).p();	// RhoPi
		CMSboost = mother.p().boostVector();
	}

	if ((generatedValues && genPresent) || !generatedValues) {
		for (int i = 0; i < 3; i++) {
			branchA[i].boost(-CMSboost);
			branchB[i].boost(-CMSboost);
		}

		HepRotation rotationToZ = getRotationToZ(branchA[0]);
		for (int i = 0; i < 3; i++) {
			branchA[i] *= rotationToZ;
			branchB[i] *= rotationToZ;
		}

		HepRotation rotationZToX = getZRotationToX(branchB[2]);
		for (int i = 0; i < 3; i++) {
			branchA[i] *= rotationZToX;
			branchB[i] *= rotationZToX;
		}

		HepRotation rotationTransversity;
		rotationTransversity.rotateX(PI/2);
		rotationTransversity.rotateZ(PI/2);
		for (int i = 0; i < 3; i++) {
			branchA[i] *= rotationTransversity;
			branchB[i] *= rotationTransversity;
		}

		Hep3Vector DSCMSboost = branchA[0].boostVector();
		Hep3Vector RhoCMSboost = branchB[0].boostVector();
		for (int i = 0; i < 3; i++) {
			branchA[i].boost(-DSCMSboost);
			branchB[i].boost(-RhoCMSboost);
		}

		thetaT = branchA[1].theta();
		phiT = branchA[1].phi();
		thetaB = PI - branchB[2].phi();
	} else {
		thetaT = -1;
		phiT = -1;
		thetaB = -1;
	}

}

/// This function returns a rotation needed to transform the given momentum
/// so that it is in the direction of the Z axis
HepRotation getRotationToZ(Hep3Vector momentum){
	HepRotation rotation;

	double phi = acos(momentum.z()/(sqrt((momentum.y())*(momentum.y())+(momentum.z())*(momentum.z()))));
	if (momentum.y() < 0) {
		rotation.rotateX(-phi);
	} else {
		rotation.rotateX(phi);
	}

	momentum = rotation*momentum;

	double theta = acos(momentum.z()/(sqrt((momentum.z())*(momentum.z())+(momentum.x())*(momentum.x()))));
	if (momentum.x() < 0) {
		rotation.rotateY(theta);
	} else {
		rotation.rotateY(-theta);
	}

	return rotation;
}

/// This function returns a rotation needed to transform the given momentum
/// so that it is in the of the XZ plane
HepRotation getZRotationToX(Hep3Vector momentum){
	HepRotation rotation;

	double phi = acos(momentum.x()/(sqrt((momentum.x())*(momentum.x())+(momentum.y())*(momentum.y()))));
	if(momentum.y() < 0) {
		rotation.rotateZ(phi);
	} else {
		rotation.rotateZ(-phi);
	}

	return rotation;
}

double getCosThetaTB(std::vector<Hep3Vector>& momenta) {
	Hep3Vector thrustVector = thrust(momenta.begin(), momenta.end(), SelfFunc(Hep3Vector()));
	return fabs(thrustVector.cosTheta());
}

double getCosThetaT(std::vector<Hep3Vector>& signalMomenta, std::vector<Hep3Vector>& otherMomenta) {
	Hep3Vector signalThrust = thrust(signalMomenta.begin(), signalMomenta.end(), SelfFunc(Hep3Vector()));
	Hep3Vector otherThrust = thrust(otherMomenta.begin(), otherMomenta.end(), SelfFunc(Hep3Vector()));
	return fabs(cos(otherThrust.angle(signalThrust)));
}

void applyCandidateCuts(std::vector<Particle> &B, Continuum& continuumSupression,
		unsigned& eventsDiscardedByDEandMbcCuts, unsigned& eventsDiscardedByCSCut, const bool useSidebands) {
	Hep3Vector myCMSBoost = -BeamEnergy::CMBoost(); // note the "-" sign
	double myCMSE = BeamEnergy::E_beam_corr();
	// Loop over all B candidates, calculate DeltaE and Mbc,
	// keep only candidates passing Mbc, DeltaE and continuum suppression cuts
	for (unsigned i = 0; i < B.size(); i++) {
		Particle& Bcand = B[i];
		HepLorentzVector B0mom(Bcand.p());
		B0mom.boost(myCMSBoost);
		double mbc = sqrt((myCMSE * myCMSE) - B0mom.vect().mag2());
		double Energy = B0mom.e();
		double deltaE = Energy - myCMSE;

		if (((!useSidebands && mbc < MBC_CUT) || 
			 (useSidebands && (mbc < MBC_SIDEBAND_CUT || mbc > MBC_CUT))) || 
			 deltaE > DELTA_E_CUT_HIGH || deltaE < DELTA_E_CUT_LOW) {
			// bad candidate, remove it
			B.erase(B.begin() + i);
			--i;
			if (B.size() == 0) {
				eventsDiscardedByDEandMbcCuts++;
			}
			continue;
		}

		// If no user info exists yet, create it
		if (!&B[i].userInfo())
			setUserInfo(B[i]);

		// save Mbc & DeltaE to UserInfo
		dynamic_cast<UserInfo&>(B[i].userInfo()).setMbc(mbc);
		dynamic_cast<UserInfo&>(B[i].userInfo()).setDeltaE(deltaE);

		if (dynamic_cast<const UserInfo&>(Bcand.child(1).child(0).userInfo()).getMassConstrainedChi2() > RHO_PI0_CHI2_CUT) {
			// bad candidate, remove it
			B.erase(B.begin() + i);
			--i;
			continue;
		}

		// Calculate and save transversity angles (both reco'ed and generated) to userInfo
		double thetaT = 0;
		double phiT = 0;
		double thetaB = 0;
		double thetaTG = 0;
		double phiTG = 0;
		double thetaBG = 0;
		getTransversityAngles(Bcand, false, thetaT, phiT, thetaB);
		getTransversityAngles(Bcand, true, thetaTG, phiTG, thetaBG);
		dynamic_cast<UserInfo&>(Bcand.userInfo()).setTransversityAngles(thetaT, phiT, thetaB);
		dynamic_cast<UserInfo&>(Bcand.userInfo()).setGeneratedTransversityAngles(thetaTG, phiTG, thetaBG);

		ksfwmoments km(Bcand, BeamEnergy::Ecm()/2, -BeamEnergy::CMBoost());
		dynamic_cast<UserInfo&>(Bcand.userInfo()).setKSFW(km);

		// p_cms_sigB has use_finalstate_for_sig == 1; meaning that thrust is computed
		// from signal final state particles rather than from B's immediate daughters;
		// see rooksfw README for more info
		std::vector<Hep3Vector> signalMomenta = km.get_p_cms_sigB();
		std::vector<Hep3Vector> otherMomenta = km.get_p_cms_other();

		dynamic_cast<UserInfo&>(Bcand.userInfo()).setCosThetaB(pStar(Bcand).vect().cosTheta());
		dynamic_cast<UserInfo&>(Bcand.userInfo()).setCosThetaTB(getCosThetaTB(signalMomenta));
		dynamic_cast<UserInfo&>(Bcand.userInfo()).setCosThetaT(getCosThetaT(signalMomenta, otherMomenta));

		continuumSupression.readInVars(
				dynamic_cast<const UserInfo&>(Bcand.userInfo()).getMbc(),
				dynamic_cast<const UserInfo&>(Bcand.userInfo()).getCosThetaB(),
				dynamic_cast<const UserInfo&>(Bcand.userInfo()).getCosThetaT(),
				km);

		dynamic_cast<UserInfo&>(Bcand.userInfo()).setContSuppressionBDTG(continuumSupression.evaluateBDTG());
		dynamic_cast<UserInfo&>(Bcand.userInfo()).setContSuppressionMLP(continuumSupression.evaluateMLP());

		// Continuum suppression cut
		if (dynamic_cast<const UserInfo&>(Bcand.userInfo()).getContSuppressionBDTG() < CS_CUT) {
			B.erase(B.begin() + i);
			--i;
			eventsDiscardedByCSCut++;
			continue;
		}
	}
}

int getExKFitterPID(const Particle& particle) {
	int exkfitter_pid;
	if(abs(particle.lund()) == 321) { // If K+-
		exkfitter_pid = 3;
	} else if(abs(particle.lund()) == 211) { // If pi+-
		exkfitter_pid = 2;
	} else {
		exkfitter_pid = -1;
	}
	return exkfitter_pid;
}

void setBRecVertex(Particle& Bcand) {
	// Create vectors to hold the children used for the vertex fit as the number
	// of D children changes between channels and more importantly only children
	// passing SVD hit cuts are added to these vectors and in turn to the vertex fit
	std::vector<ExKFitterParticle> D_children;
	std::vector<ExKFitterParticle> B_children;

	Particle& D_particle = Bcand.child(0).child(0);
	for (unsigned i = 0; i < D_particle.nChildren(); i++) {
		if (enoughSVDHits(D_particle.child(i))) {
			D_children.push_back(ExKFitterParticle(D_particle.child(i).mdstCharged(), getExKFitterPID(D_particle.child(i))));
		}
	}
	Particle& slow_Pi_particle = Bcand.child(0).child(1);
	if (enoughSVDHits(slow_Pi_particle)) {
		B_children.push_back(ExKFitterParticle(slow_Pi_particle.mdstCharged(), getExKFitterPID(slow_Pi_particle)));
	}
	Particle& Rho_Pi_particle = Bcand.child(1).child(1);
	if (enoughSVDHits(Rho_Pi_particle)) {
		B_children.push_back(ExKFitterParticle(Rho_Pi_particle.mdstCharged(), getExKFitterPID(Rho_Pi_particle)));
	}

	// Remember it's better to set this Initial_D_Vertex
	// for the case of multi-vertex simultaneous fit
//	ExKFitterVertex D_Vertex(Initial_D_Vertex);
	ExKFitterVertex D_Vertex;

	// Assign IP and IPerr(3x3 error matrix) to B_Vertex = IP constrained fit
	// Even for IPtube contrained fit, the vertex has to be initialized like this
	ExKFitterVertex B_Vertex(IpProfile::position(), IpProfile::position_err());

	ExKFitterParticle D; // Create a virtual particle D
	for(std::vector<ExKFitterParticle>::iterator it = D_children.begin(); it != D_children.end(); ++it) {
		D.LinkParticle(&*it);
	}
	D.LinkVertex(&D_Vertex);

	ExKFitterConstrain con1; // D vertex constraint
	con1.SetVertexConstrain();
	for(std::vector<ExKFitterParticle>::iterator it = D_children.begin(); it != D_children.end(); ++it) {
		con1.LinkParticle(&*it);
	}
	con1.LinkVertex(&D_Vertex);

	B_children.push_back(D);

	ExKFitterConstrain con2; // B vertex constraint
	ExKFitterParticle IPtube(IpProfile::ip_tube()); // This is the way to add IP *tube* constraint
	con2.SetVertexConstrain();
	con2.LinkParticle(&IPtube);
	for(std::vector<ExKFitterParticle>::iterator it = B_children.begin(); it != B_children.end(); ++it) {
		con2.LinkParticle(&*it);
	}
	con2.LinkVertex(&B_Vertex);

	ExKFitterMass D_Mass(1.8645);

	ExKFitterConstrain con3; // D mass constraint
	con3.SetMassConstrain();
	for(std::vector<ExKFitterParticle>::iterator it = D_children.begin(); it != D_children.end(); ++it) {
		con3.LinkParticle(&*it);
	}
	con3.LinkVertex(&D_Vertex);
	con3.LinkMass(&D_Mass);

	ExKFitter Core;
	Core.LinkConstrain(&con1);
	Core.LinkConstrain(&con2);
	Core.LinkConstrain(&con3);

	vfit_info vtx_rec;

	const int vtx_fit_error = Core.Minimize();
	if (!vtx_fit_error) {
		vtx_rec.m_usable   = true;
		vtx_rec.m_pos      = B_Vertex.Vertex();
		vtx_rec.m_err      = B_Vertex.ErrVertex();
		vtx_rec.m_chi2     = Core.Chisq();
		vtx_rec.m_effxi	   = Core.EffectiveXi(&con2);
		vtx_rec.m_ndf      = Core.N_DegreeOfFreedom();
		vtx_rec.m_effndf   = Core.N_EffectiveDOF(&con2);
		vtx_rec.m_ntrk     = Core.N_VertexingTracks();
	}

	dynamic_cast<UserInfo&>(Bcand.userInfo()).setVertexRec(vtx_rec);
}

void setBRecVertexKFitter(Particle& Bcand) {
	kvertexfitter kf;

//	int num_tracks = 0;
//	if (enoughSVDHits(Bcand.child(0).child(0).child(0))) {
//		addTrack2fit(kf, Bcand.child(0).child(0).child(0));
//		num_tracks++;
//	}
//	if (enoughSVDHits(Bcand.child(0).child(0).child(1))) {
//		addTrack2fit(kf, Bcand.child(0).child(0).child(1));
//		num_tracks++;
//	}

	int num_tracks = 2;
	addTrack2fit(kf, Bcand.child(0).child(0));

	// No SVD hits cut on the slow pi
	addTrack2fit(kf, Bcand.child(0).child(1));
	num_tracks++;

//	addTrack2fit(kf, Bcand.child(0)); // Chi2 too large
	if (enoughSVDHits(Bcand.child(1).child(1))) {
		addTrack2fit(kf, Bcand.child(1).child(1));
		num_tracks++;
	}
	addTube2fit(kf);

	vfit_info vtx_rec;
	unsigned vtx_fit_error = kf.fit();
	if (!vtx_fit_error) {
//		makeMother(kf, p); // change Momentum Info.
//		myMakeMother(kf, p); // change Momentum Info.
		vtx_rec.m_usable = true;
		vtx_rec.m_pos = kf.vertex();
		vtx_rec.m_err = kf.errVertex();
		vtx_rec.m_chi2 = kf.chisq_tracks();
		vtx_rec.m_ndf = kf.dgf_tracks();
		vtx_rec.m_ntrk = num_tracks; //TODO: Generalize
	} else {
		HepPoint3D vtx(999., 999., 999.);
		HepSymMatrix errVtx(3, 0);
		Bcand.momentum().decayVertex(vtx, errVtx);
	}
	dynamic_cast<UserInfo&>(Bcand.userInfo()).setVertexRec(vtx_rec);
}

int setDVertex(Particle& Dcand) {
	kvertexfitter kf;
	int num_tracks = 0;
	for (unsigned i = 0; i < Dcand.nChildren(); i++) {
		if (enoughSVDHits(Dcand.child(i))) {
			addTrack2fit(kf, Dcand.child(i));
			num_tracks++;
		}
	}
	unsigned vtx_fit_error = kf.fit();
	if (vtx_fit_error == 0) {
		makeMother(kf, Dcand);
	} else {
		HepPoint3D vtx(999., 999., 999.);
		HepSymMatrix errVtx(3, 0);
		Dcand.momentum().decayVertex(vtx, errVtx);
	}
	return num_tracks;
}

void setBAscVertex(Particle& Bcand) {
	vfit_info vtx_asc;
	const vfit_info vtx_rec = dynamic_cast<UserInfo&>(Bcand.userInfo()).getVertexRec();
	TagVK tagv_asc;

	std::vector<Particle> taglist;
	std::vector<Particle> k_p, k_m, pi_p, pi_m;
	makeKPi(k_p, k_m, pi_p, pi_m);
	taglist = pi_p;
	taglist.insert(taglist.end(), pi_m.begin(), pi_m.end());

	Particle bcand_copy(Bcand);
//	tagv_asc.setdefault(bcand_copy, vtx_rec.m_pos, TAG_TRACK_DR_CUT, TAG_TRACK_DZ_CUT, TAG_TRACK_DZ_ERROR_CUT, 1);
	tagv_asc.setdefault(bcand_copy, vtx_rec.m_pos);

	for (unsigned i = 0; i < Bcand.relation().nFinalStateParticles(); i++) {
		removeParticle(taglist, Bcand.relation().finalStateParticle(i));
	}

	for (std::vector<Particle>::iterator it = taglist.begin(); it != taglist.end(); ++it) {
		tagv_asc.push_back(&(*it));
	}

	const int fit_error = tagv_asc.fit();
	if (!fit_error) {
		vtx_asc.m_usable = true;
		vtx_asc.m_pos = tagv_asc.vtx();
		vtx_asc.m_err = tagv_asc.errVtx();
		vtx_asc.m_chi2 = tagv_asc.chisq_tracks();
		vtx_asc.m_ndf = (int) tagv_asc.ndf_tracks();
		vtx_asc.m_ntrk = tagv_asc.used_particles().size();
		vtx_asc.m_tmp[1] = tagv_asc.isTagLeptonVertex() ? (int) tagv_asc.VertexTagLepton().get_ID() : 0;
	}

//	vtx_asc.m_tmp[2] = 0;
//	for( unsigned i=0; i<tagv_asc.used_particles().size(); i++ ){
//		Particle p(*(tagv_asc.used_particles())[i]);
//		if( p.mdstVee2() ) vtx_asc.m_tmp[2]++;
//	}

	dynamic_cast<UserInfo&>(Bcand.userInfo()).setVertexTag(vtx_asc);
}

bool enoughSVDHits(const Particle& p) {
	if (!p.mdstCharged()) return false;

	Mdst_trk_fit& trk_fit = p.mdstCharged().trk().mhyp(2);
	if (trk_fit.nhits(3) < TRACK_NHITS_RPHI) return false;
	if (trk_fit.nhits(4) < TRACK_NHITS_Z) return false;

	return true;
}

void tag(Particle& Bcand) {
	Hamlet flavor_tagger;
	flavor_tagger.setBcp(Bcand, 0);
	flavor_tagger.setTagMethod(Hamlet::MULT_DIM_LH);

	dynamic_cast<UserInfo&>(Bcand.userInfo()).setFlavor(flavor_tagger.flavor());
	// Hamlet's q is what is usually refered to as q*r; q is the flavor, r
	// (a.k.a. D), the dilution factor r = D = 1 - 2*wtag
	dynamic_cast<UserInfo&>(Bcand.userInfo()).setQr(flavor_tagger.q());
	dynamic_cast<UserInfo&>(Bcand.userInfo()).setWtag((1.0 - flavor_tagger.q()*flavor_tagger.flavor())*0.5);
}

double twoParticleInvariantMass(const Particle& part1, const Particle& part2) {
	HepLorentzVector p1(part1.p());
	HepLorentzVector p2(part2.p());
	p1 += p2;
	return p1.m();
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
