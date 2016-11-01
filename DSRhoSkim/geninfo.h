#ifndef GENINFO_H_
#define GENINFO_H_

#include "belle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

void setMCtruth(std::vector<Particle> &plist);
void setMCtruth(Particle &p);
void setMCtruthPi0(Particle &p);
int isFSP(Gen_hepevt P);
void appendRecFSP(Particle p, std::vector<Particle> &children);
int appendGenFSP(const Gen_hepevt &gen, std::vector<Gen_hepevt> &children);
int dumpGenFSP(const Gen_hepevt &gen);
int commonMother(std::vector<int> &mothers);
int findCommonMother(Gen_hepevt pThis, Gen_hepevt pOther, int level);
int compareFinalStates(std::vector<Gen_hepevt> reconstructed, std::vector<Gen_hepevt> generated);
void fillMothers(Particle &p, std::vector<int> &thisMothers, std::vector<int> &otherMothers);
void fillMothers(Particle &A, Particle &B, std::vector<int> &thisMothers, std::vector<int> &otherMothers);
void fillMothers(Particle &A, Particle &B, Particle &C, std::vector<int> &thisMothers, std::vector<int> &otherMothers);

int getCommonMother(Particle &A, Particle &B);
int getCommonMother(Particle &A, Particle &B, Particle &C);

int findCommonMother(int nChildren, std::vector<int> thisMothers, std::vector<int> otherMothers);

int getMCtruthPi0Flag(Particle &p);
int getMCtruthChargedTrackFlag(Particle &p);
int getMCtruthFlag(Particle &p);

void genDecayChain(Particle p, int* dChain);

int IDhep(Particle &part);
int NdecayProd(Particle &part);


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif /* GENINFO_H_ */
