#ifndef MYUTILS_H_
#define MYUTILS_H_

#include "belle.h"
#include "Continuum.h"

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
void setGammaError(Particle & gamma);
void setGammaError(std::vector<Particle> &p);
void setPi0Error(Particle &p);
void setPi0Error(std::vector<Particle> &p);
void saveTrackInfo(Particle &p, double dr, double dz, double eid, double pid, double muid, double protonID);
void savePi0Info(Particle &p, double min_e9oe25, double max_e9oe25, double min_e, double max_e, double min_theta, double max_theta, double decAng);
void setDecayMode(std::vector<Particle> &list, int dMode) ;
void impact_param(Particle &p, double &dr, double &dz);
void PlistCopy(std::vector<Particle> &plist, std::vector<Particle> &copy_plist);
void withImpactPars(std::vector<Particle> &list, const double max_dr, const double max_dz);
Particle& selectBestCandidate(unsigned& index, std::vector<Particle>& list);
Particle& selectBestCandidate2(unsigned& index, std::vector<Particle>& list);
Particle& selectRandomBestCandidate(unsigned& index, std::vector<Particle>& list);
bool findRightCandidate(unsigned& index, std::vector<Particle>& list);
void getTransversityAngles(const Particle& mother, const bool generatedValues, double& thetaT, double& phiT, double& thetaB);
HepRotation getRotationToZ(Hep3Vector momentum);
HepRotation getZRotationToX(Hep3Vector momentum);
double getCosThetaTB(std::vector<Hep3Vector>& momenta);
double getCosThetaT(std::vector<Hep3Vector>& signalMomenta, std::vector<Hep3Vector>& otherMomenta);
bool isSignal(const Particle& B);
int trueDecayGenerated(int channelNum);
bool findChildren(Gen_hepevt* particle, std::vector<Gen_hepevt**>& daughters, std::vector<int>& idheps, bool lookAtFSPChildren = false);
void setBCSMetadata(int trueCandidateWasGenerated, const unsigned indexOfBestCandidate, std::vector<Particle>& B,
		unsigned& eventsWithRightCandidate, unsigned& eventsWithOutRightCandidate, unsigned& numCorrectlySelectedCandidates);
double invariantMass(const Particle& p1, const Particle& p2);
void applyCandidateCuts(std::vector<Particle> &B, Continuum& continuumSupression,
		unsigned& eventsDiscardedByDEandMbcCuts, unsigned& eventsDiscardedByCSCut);
void setBRecVertex(Particle& Bcand);
void setBRecVertexKFitter(Particle& Bcand);
int setDVertex(Particle& Dcand);
void setBAscVertex(Particle& Bcand);
bool enoughSVDHits(const Particle& p);
int getExKFitterPID(const Particle& particle);
void tag(Particle& Bcand);

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif /* MYUTILS_H_ */
