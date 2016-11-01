// ******************************************************************
// Utilty functions
// A. Zupanc jun-06
//
// ******************************************************************
#include "belle.h"
#include "panther/panther.h"
#include "eid/eid.h"
#include "kid/atc_pid.h"
#include "mdst/mdst.h"
#include "mdst/Muid_mdst.h"
#include "mdst/findKs.h"
#include MDST_H
#include HEPEVT_H
#include "particle/utility.h"
#include "ip/IpProfile.h"
#include "helix/Helix.h"
#include "belleCLHEP/Vector/Rotation.h"

#include "cstdlib"

#include "myutils.h"
#include "constants.h"
#include "UserInfo.h"

#include <stdio.h>

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


//#include <algorithm>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

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
		if(it->mass()<0.115 || it->mass()>0.15 || pizero.child(0).e()<0.05 || pizero.child(1).e()<0.05)
			continue;

		// creates the userinfor for this pi0 candidate
		setUserInfo(pizero);

		// saves the pi0 mass before mass constrained fit
		dynamic_cast<UserInfo&>(pizero.userInfo()).setMassBeforeVertexFit(it->mass());

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

		// Angle between mother momentum in lab frame and child momentum in mother's rest frame
		float decAng = -999;

		HepLorentzVector mother_mom = pizero.momentum().p();
		Hep3Vector mother_boost = -(mother_mom.boostVector());

		HepLorentzVector child0_mom = pizero.child(0).momentum().p();
		child0_mom.boost(mother_boost);

		decAng = (float)cos(child0_mom.angle(mother_mom.vect()));

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

Particle selectBestCandidate(unsigned& index, std::vector<Particle>& B){
	double smallestMbcDelta = fabs(MBC - dynamic_cast<const UserInfo&>(B[0].userInfo()).getMbc());
	unsigned indexOfSmallest = 0;
	double currentMbcDelta = 0;

	for(unsigned i=1; i<B.size(); i++) {
		currentMbcDelta = fabs(MBC - dynamic_cast<const UserInfo&>(B[i].userInfo()).getMbc());
		if(smallestMbcDelta > currentMbcDelta){
			smallestMbcDelta = currentMbcDelta;
			indexOfSmallest = i;
		}
	}
	index = indexOfSmallest;
	return B[indexOfSmallest];
}

Particle selectRandomBestCandidate(unsigned& index, std::vector<Particle>& B){
	index = (rand()/RAND_MAX) * (B.size() - 1);
	return B[index];
}

void GetTransversityAngles(const Particle& mother, double& thetaT, double& phiT, double& thetaB){

	HepLorentzVector branchA[] = {
		mother.child(0).p(),			// DS
		mother.child(0).child(0).p(),	// DSD0
		mother.child(0).child(1).p()	// DSPi
	};

	HepLorentzVector branchB[] = {
		mother.child(1).p(),			// Rho
		mother.child(1).child(0).p(),	// RhoPi0
		mother.child(1).child(1).p()	// RhoPi
	};

	Hep3Vector CMSboost = mother.p().boostVector();
	for (int i = 0; i < 3; i++) {
		branchA[i].boost(-CMSboost);
		branchB[i].boost(-CMSboost);
	}

	HepRotation rotationToZ = GetRotationToZ(branchA[0]);
    for (int i = 0; i < 3; i++) {
		branchA[i] *= rotationToZ;
		branchB[i] *= rotationToZ;
	}

    HepRotation rotationZToX = GetZRotationToX(branchB[2]);
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

}

/// This function returns a rotation needed to transform the given momentum
/// so that it is in the direction of the Z axis
HepRotation GetRotationToZ(Hep3Vector momentum){
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
HepRotation GetZRotationToX(Hep3Vector momentum){
    HepRotation rotation;

    double phi = acos(momentum.x()/(sqrt((momentum.x())*(momentum.x())+(momentum.y())*(momentum.y()))));
    if(momentum.y() < 0)
        rotation.rotateZ(phi);
    else
        rotation.rotateZ(-phi);

    return rotation;
}



#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
