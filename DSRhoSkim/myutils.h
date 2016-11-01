#ifndef MYUTILS_H_
#define MYUTILS_H_

#include "belle.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

void impact_param(Particle &p, double &dr, double &dz);
void fillKaonsAndPions(std::vector<Particle> &kp,
		std::vector<Particle> &km,
		std::vector<Particle> &pi_p,
		std::vector<Particle> &pi_m);
void fillGoodKShort(std::vector<Particle> &K0short);
void fillAllPi0(std::vector<Particle> &pi0);
void saveTrackInfo(Particle &p, double dr, double dz, double eid, double pid, double muid, double protonID);
void savePi0Info(Particle &p, double min_e9oe25, double max_e9oe25, double min_e, double max_e, double min_theta, double max_theta, double decAng);
void setDecayMode(std::vector<Particle> &list, int dMode) ;
void impact_param(Particle &p, double &dr, double &dz);
void PlistCopy(std::vector<Particle> &plist, std::vector<Particle> &copy_plist);
void withImpactPars(std::vector<Particle> &list, const double max_dr, const double max_dz);
Particle selectBestCandidate(unsigned& index, std::vector<Particle>& list);
Particle selectRandomBestCandidate(unsigned& index, std::vector<Particle>& list);
void GetTransversityAngles(const Particle& mother, double& thetaT, double& phiT, double& thetaB);
HepRotation GetRotationToZ(Hep3Vector momentum);
HepRotation GetZRotationToX(Hep3Vector momentum);

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif /* MYUTILS_H_ */
