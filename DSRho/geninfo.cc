// ******************************************************************
// MC Matching
// author: A. Zupanc (anze.zupanc@ijs.si)
//
// Description: TODO
//
// ******************************************************************
#include "belle.h"
#include "panther/panther.h"
#include "mdst/mdst.h"
#include MDST_H
#include HEPEVT_H

#include "particle/utility.h"

#include "geninfo.h"
#include "constants.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

int isFSP(Gen_hepevt P) {
	switch (abs(P.idhep())) {
	case 211:	// pi
		return 1;
	case 321:	// K
		return 1;
	case 11:	// e
		return 1;
	case 13:	// mu
		return 1;
	case 22:	// gamma
		return 1;
	case 2212:	// p
		return 1;
	case 111:	// pi0
		return 1;
//	case 310:	// K_S
//		return 1;
	case 130:	// K_L
		return 1;
	case 2112:	// n
		return 1;
	default:
		return 0;
	}
}
void appendRecFSP(Particle p, std::vector<Particle> &children) {
	for (int i = 0; i < (int)p.nChildren(); ++i) {
		if(p.child(i).nChildren() && p.child(i).lund()!=111 && p.child(i).lund()!=310) {
			appendRecFSP(p.child(i),children);
		} else {
			children.push_back(p.child(i));
		}
	}
}

int appendGenFSP(const Gen_hepevt &gen, std::vector<Gen_hepevt> &children) {
	Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();

	for (int i = gen.daFirst(); i <= gen.daLast(); ++i) {
		if(i==0) {
			std::cout << "[GenInfo] appendGenFSP: requesting Particle with ID = 0! Exiting." << std::endl;
			return 0;
		}

		const Gen_hepevt& child = GenMgr(Panther_ID(i));
		int ndaug2 = (child.daLast()-child.daFirst()) + 1;
		if(ndaug2 && !isFSP(child)) {
			appendGenFSP(child,children);
		} else {
			children.push_back(child);
		}
	}
	return 1;
}

int appendGenFSP(const Gen_hepevt &gen, std::vector<Gen_hepevt*> &children) {
	Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();

	for (int i = gen.daFirst(); i <= gen.daLast(); ++i) {
		if(i==0) {
			std::cout << "[GenInfo] appendGenFSP: requesting Particle with ID = 0! Exiting." << std::endl;
			return 0;
		}

		Gen_hepevt* child = & GenMgr(Panther_ID(i));
		int ndaug2 = (child->daLast()-child->daFirst()) + 1;
		if(ndaug2 && !isFSP(*child)) {
			appendGenFSP(*child,children);
		} else {
			children.push_back(child);
		}
	}
	return 1;
}

int dumpGenFSP(const Gen_hepevt &gen) {
	Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();

	int ndaug = (gen.daLast()-gen.daFirst()) + 1;
	std::cout << "Appending FS children for " << gen.idhep() << " with nChildren = " << ndaug << std::endl;
	for (int i = gen.daFirst(); i <= gen.daLast(); ++i) {
		std::cout << "searching for particle with ID = " << i << std::endl;
		if(i==0) {
			std::cout << "[GenInfo] appendGenFSP: requesting Particle with ID = 0! Exiting." << std::endl;
			return 0;
		}

		const Gen_hepevt& child = GenMgr(Panther_ID(i));
		if(!child)
			std::cout << "Particle with this PantherID does not exist." << std::endl;
		int ndaug2 = (child.daLast()-child.daFirst()) + 1;
		std::cout << " -> child " << i << ": ID = " << i << "with nChildren = " << ndaug2 << std::endl;
		if(ndaug2 && !isFSP(child)) {
			std::cout << "not FSP: " << child.idhep() << std::endl;
			dumpGenFSP(child);
		} else {
			std::cout << "Appending FSP: " << child.idhep() << std::endl;
		}
	}
	return 1;
}

int commonMother(std::vector<int> &mothers) {
	if(mothers.size()==0) {
		return 0;
	} else if(mothers.size()==1) {
		return mothers[0];
	}

	int motherID = mothers[0];
	for (int i = 1; i < (int)mothers.size(); ++i) {
		if(motherID!=mothers[i]) {
			return 0;
		}
	}
	return motherID;
}

void fillMothers(Particle &p, std::vector<int> &thisMothers, std::vector<int> &otherMothers) {
	Gen_hepevt motherThis(p.child(0).genHepevt());

	thisMothers.push_back(motherThis.get_ID());

	while(motherThis.mother()) {
		motherThis = motherThis.mother();
		thisMothers.push_back(motherThis.get_ID());
	}

	for(int i = 1; i<(int)p.nChildren(); ++i) {
		motherThis = p.child(i).genHepevt();
		// new
		otherMothers.push_back(motherThis.get_ID());
		while(motherThis.mother()) {
			motherThis = motherThis.mother();
			otherMothers.push_back(motherThis.get_ID());
		}
	}
}

int findCommonMother(int nChildren, std::vector<int> thisMothers, std::vector<int> otherMothers) {
	for (int i = 0; i < (int)thisMothers.size(); ++i) {
		int counter = 0;
		for (int j = 0; j < (int)otherMothers.size(); ++j) {
			if(thisMothers[i]==otherMothers[j])
				counter++;
		}
		if(counter==nChildren-1)
			return thisMothers[i];
	}
	return 0;
}
int findCommonMother(Gen_hepevt pThis, Gen_hepevt pOther, int level) {
	Gen_hepevt motherThis(pThis);

	int i = 1;
	while(motherThis.mother()) {
		motherThis = motherThis.mother();
		i++;

		if(i>level) {
			Gen_hepevt motherOther(pOther);
			while(motherOther.mother()) {
				motherOther = motherOther.mother();
				if(motherThis.get_ID()==motherOther.get_ID())
					return motherThis.get_ID();
			}
		}
	}
	return 0;
}

// TODO
int compareFinalStates(std::vector<Particle> reconstructed, std::vector<Gen_hepevt> generated) {
	if(reconstructed.size() == generated.size()) {
		int missID = 0;
		int missPi0 = 0;
		for (int i = 0; i < (int)reconstructed.size(); ++i) {
			if(reconstructed[i].genHepevt()) {
				int link = 0;
				for (int j = 0; j < (int)generated.size(); ++j) {
					if(reconstructed[i].genHepevt().get_ID()==generated[j].get_ID()) {
						link = 1;
						if(reconstructed[i].lund()!=generated[j].idhep())
							missID++;
						if(reconstructed[i].lund()==111 && getMCtruthPi0Flag(reconstructed[i])!=1) {
							missPi0++;
						}
						break;
					}
				}
				if(!link) {
					std::cout << "[GenInfo] compareFinalStates: Particle (" << reconstructed[i].lund() << " with hepevt = " << reconstructed[i].genHepevt().idhep() << ") not found in list of FSP (gen)!" << std::endl;
					return -11;
				}
			} else {
				std::cout << "[GenInfo] compareFinalStates: Particle without link to genHepevt! [-10]" << std::endl;
				return -10;
			}
		}
		if(missID && missPi0)
			return 6;
		else if(missID)
			return 2;
		else if(missPi0)
			return 5;

		return 1;
	} else if(reconstructed.size() < generated.size()) { // missing particle
		int missing = 0;
		std::vector<int> missP;
		for (int i = 0; i < (int)generated.size(); ++i) {
			int link = 0;
			for (int j = 0; j < (int)reconstructed.size(); ++j) {
				if(reconstructed[j].genHepevt()) {
					if(reconstructed[j].genHepevt().get_ID()==generated[i].get_ID()) {
						link = 1;
						break;
					}
				} else {
					std::cout << "[GenInfo] compareFinalStates: Particle without link to genHepevt! [-9]" << std::endl;
					return -9;
				}
			}
			if(!link) {
				missing++;
				missP.push_back(i);
			}
		}
		if(missing) {
			int masslessOnly = 1;
			int missNu = 0;
			int missGfsr = 0;
			int missGrad = 0;
			for (int i = 0; i < (int)missP.size(); ++i) {
				// new way : distinguish between missing neutrino and FSR gamma and radiative gamma (the latter separation is not 100% precise)
				if(abs(generated[missP[i]].idhep())==12 || abs(generated[missP[i]].idhep())==14 || abs(generated[missP[i]].idhep())==16) {
					// neutrino is missing
					missNu = 1;
				} else if(generated[missP[i]].idhep()==22) {
					// photon is missing (FSR or radiative?)
					if(generated[missP[i]].mother()) {
						int ndaug = (generated[missP[i]].mother().daLast()-generated[missP[i]].mother().daFirst()) + 1;
						if(ndaug==2)
							missGrad = 1;
						else
							missGfsr = 1;
					} else
						missGfsr = 1;
				} else
					masslessOnly = 0;
			}
			if(missNu) {
				if(missGrad && !masslessOnly)
					return 21;
				else if(missGrad && masslessOnly)
					return 24;
				else if(!missGrad && !masslessOnly)
					return 23;
				else
					return 20;
			} else {
				if(missGfsr && !missGrad && masslessOnly)
					return 10;
				else if(missGrad && !masslessOnly)
					return 11;
				else if(missGrad && masslessOnly)
					return 4;
				else if(!missGrad && !masslessOnly)
					return 3;
				else
					return -20;
			}
		} else {
			std::cout << "[GenInfo] compareFinalStates: At least one particle should be missing!" << std::endl;
			return -8;
		}
	} else {
		std::cout << "[GenInfo] compareFinalStates: More reconstructed than generated particle in final state!" << std::endl;
		return -5;
	}
}

//TODO
// Flags
// -1 - none of the below
//  0 - neither of the two gammas point to pi0
//  1 - correctly reconstructed pi0
//  2 - only one gamma has mother (CAUTION) - link is made to that pi0
//  3 - only one gamma has link to pi0 (CAUTION), the other gamma links to something else - link is made to that pi0
//  4 - mothers of gammas are different pi0s (CAUTION) - link is made to mother of gamma with largest energy
//  5 - only one gamma has link to gen_ehepvt - link is made to mother of gamma, if exists
int getMCtruthPi0Flag(Particle &p) {
	// if link is not to pi0
	if(IDhep(p)!=111)
		return 0;

	//Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
	if(p.child(0).genHepevt() && p.child(1).genHepevt()) {
		const Gen_hepevt & mother1(p.child(0).genHepevt());
		const Gen_hepevt & mother2(p.child(1).genHepevt());
		if(mother1.mother() && mother2.mother()) {
			if(mother1.mother().get_ID()==mother2.mother().get_ID()) {
				if(mother1.mother().get_ID()==0) {
					std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
					return 0;
				}
				return 1;
			} else if(mother1.mother().idhep()==111 && mother2.mother().idhep()==111) {
				if(p.child(0).e()>p.child(1).e()) {
					if(mother1.mother().get_ID()==0) {
						std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
						return 0;
					}
					return 4;
				} else {
					if(mother2.mother().get_ID()==0) {
						std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
						return 0;
					}
					return 4;
				}
			} else if(mother1.mother().idhep()==111 || mother2.mother().idhep()==111) {
				if(mother1.mother().idhep()==111) {
					if(mother1.mother().get_ID()==0) {
						std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
						return 0;
					}
					return 3;
				} else {
					if(mother2.mother().get_ID()==0) {
						std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
						return 0;
					}
					return 3;
				}
			}
		} else {
			if(mother1.mother()) {
				if(mother1.mother().get_ID()==0) {
					std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
					return 0;
				}
				return 2;
			} else if(mother2.mother()) {
				if(mother2.mother().get_ID()==0) {
					std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
					return 0;
				}
				return 2;
			}
		}
	} else if(p.child(0).genHepevt()) {
		const Gen_hepevt & mother(p.child(0).genHepevt());
		if(mother.mother()) {
			if(mother.mother().idhep()==111) {
				if(mother.mother().get_ID()==0) {
					std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
					return 0;
				}
				return 5;
			}
			return 0;
		}
	} else if(p.child(1).genHepevt()) {
		const Gen_hepevt & mother(p.child(1).genHepevt());
		if(mother.mother()) {
			if(mother.mother().idhep()==111) {
				if(mother.mother().get_ID()==0) {
					std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
					return 0;
				}
				return 5;
			}
			return 0;
		}
	} else {
		std::cout << "[GenInfo] setMCtruthPi0: Neither of two photons has link to gen_hepevt!" << std::endl;
	}

	return -1;
}


void setMCtruthPi0(Particle &p) {
	// MC truth already set
	if(IDhep(p)!=0) {
		return;
	}

	for(int i=0; i<(int)p.nChildren(); ++i) {
		if(p.child(i).mdstGamma()) {
			const Gen_hepevt & hep(gen_level(get_hepevt(p.child(i).mdstGamma())));
			if (hep){
				p.child(i).relation().genHepevt(hep);
			} else {
				std::cout << " [GenInfo] setMCtruthPi0: child " << i << " has no link! " << std::endl;
			}
		} else {
			std::cout << " [GenInfo] setMCtruthPi0: pi0 child is not mdstGamma! Exit. " << std::endl;
		}
	}

	Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
	if(p.child(0).genHepevt() && p.child(1).genHepevt()) {
		const Gen_hepevt & mother1(p.child(0).genHepevt());
		const Gen_hepevt & mother2(p.child(1).genHepevt());
		if(mother1.mother() && mother2.mother()) {
			if(mother1.mother().get_ID()==mother2.mother().get_ID()) {
				if(mother1.mother().get_ID()==0) {
					std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
					return;
				}
				const Gen_hepevt& thisGen = GenMgr(Panther_ID(mother1.mother().get_ID()));
				p.relation().genHepevt(thisGen);
				return;
			} else if(mother1.mother().idhep()==111 && mother2.mother().idhep()==111) {
				if(p.child(0).e()>p.child(1).e()) {
					if(mother1.mother().get_ID()==0) {
						std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
						return;
					}
					const Gen_hepevt& thisGen = GenMgr(Panther_ID(mother1.mother().get_ID()));
					p.relation().genHepevt(thisGen);
					return;
				} else {
					if(mother2.mother().get_ID()==0) {
						std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
						return;
					}
					const Gen_hepevt& thisGen = GenMgr(Panther_ID(mother2.mother().get_ID()));
					p.relation().genHepevt(thisGen);
					return;
				}
			} else if(mother1.mother().idhep()==111 || mother2.mother().idhep()==111) {
				if(mother1.mother().idhep()==111) {
					if(mother1.mother().get_ID()==0) {
						return;
					}
					const Gen_hepevt& thisGen = GenMgr(Panther_ID(mother1.mother().get_ID()));
					p.relation().genHepevt(thisGen);
					return;
				} else {
					if(mother2.mother().get_ID()==0) {
						std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
						return;
					}
					const Gen_hepevt& thisGen = GenMgr(Panther_ID(mother2.mother().get_ID()));
					p.relation().genHepevt(thisGen);
					return;
				}
			}
		} else {
			if(mother1.mother()) {
				if(mother1.mother().get_ID()==0) {
					std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
					return;
				}
				const Gen_hepevt& thisGen = GenMgr(Panther_ID(mother1.mother().get_ID()));
				p.relation().genHepevt(thisGen);
				return;
			} else if(mother2.mother()) {
				if(mother2.mother().get_ID()==0) {
					std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
					return;
				}
				const Gen_hepevt& thisGen = GenMgr(Panther_ID(mother2.mother().get_ID()));
				p.relation().genHepevt(thisGen);
				return;
			}
		}
	} else if(p.child(0).genHepevt()) {
		std::cout << "[GenInfo] setMCtruthPi0: Only first photon with link to gen_hepevt: " << std::endl;
		const Gen_hepevt & mother(p.child(0).genHepevt());
		std::cout << " -> idhep = " << mother.idhep() << std::endl;
		if(mother.mother()) {
			std::cout << " -> mother ID = " << mother.mother().idhep() << std::endl;
			if(mother.mother().idhep()==111) {
				if(mother.mother().get_ID()==0) {
					std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
					return;
				}
				const Gen_hepevt& thisGen = GenMgr(Panther_ID(mother.mother().get_ID()));
				p.relation().genHepevt(thisGen);
				return;
			}
			return;
		}
	} else if(p.child(1).genHepevt()) {
		std::cout << "[GenInfo] setMCtruthPi0: Only second photon with link to gen_hepevt: " << std::endl;
		const Gen_hepevt & mother(p.child(1).genHepevt());
		std::cout << " -> idhep = " << mother.idhep() << std::endl;
		if(mother.mother()) {
			std::cout << " -> mother ID = " << mother.mother().idhep() << std::endl;
			if(mother.mother().idhep()==111) {
				if(mother.mother().get_ID()==0) {
					std::cout << "[GenInfo] setMCtruthPi0: requesting particle with ID = 0! Exiting." << std::endl;
					return;
				}
				const Gen_hepevt& thisGen = GenMgr(Panther_ID(mother.mother().get_ID()));
				p.relation().genHepevt(thisGen);
				return;
			}
			return;
		}
	} else {
		std::cout << "[GenInfo] setMCtruthPi0: Neither of two photons has link to gen_hepevt!" << std::endl;
	}


	return;
}

// Flags
// 0 - link doesn't exist
// 1 - correctly identified charged track
// 2 - misidentified charged track
int getMCtruthChargedTrackFlag(Particle &p) {
	if(IDhep(p)==0)
		return 0;


	if(p.lund()==IDhep(p)) {
		return 1;
	} else {
		return 2;
	}
}

void setMCtruth(Particle &p){
	// MC truth already set
	if(IDhep(p)!=0) {
		return;
	}

	if(p.lund()==111) {
		setMCtruthPi0(p);
		return;
	}

	if(p.mdstCharged()) {
		const Gen_hepevt & hep(gen_level(get_hepevt(p.mdstCharged())));
		if (hep){
			p.relation().genHepevt(hep);
			return;
		}
	} else if(p.mdstGamma()) {
		const Gen_hepevt & hep(gen_level(get_hepevt(p.mdstGamma())));
		if (hep){
			p.relation().genHepevt(hep);
			return;
		}
	}

	// it is not mdstCharged or mdstGamma (combined particle)
	int nChildren = p.relation().nChildren();
	if(nChildren<2)
		return;

	// special treatment for pi0
	// in case that one of the gammas doesn't have link to gen_hepevt but the otherone has
	// set link for pi0 to the mother of that gamma

	// check that all child particles have MC truth
	for(int i=0; i<nChildren; ++i) {
		if(!p.relation().child(i).genHepevt())
			setMCtruth(p.relation().child(i));

		if(!p.relation().child(i).genHepevt()) {
			return;
		}
		if(p.relation().child(i).genHepevt().idhep() == 0) {
			return;
		}
	}

	// check that there are no clones
	for(int i=0; i<nChildren; ++i) {
		for(int j=i+1; j<nChildren; ++j) {
			if(p.relation().child(i).genHepevt().get_ID() == p.relation().child(j).genHepevt().get_ID()) {
				return;
			}
		}
	}

	// find common mother, if exists (p has at least 2 children)
	std::vector<int> thisMothers;
	std::vector<int> otherMothers;

	fillMothers(p, thisMothers, otherMothers);
	if(thisMothers.size()==0 || otherMothers.size()==0) {
		return;
	}
	int motherID =  findCommonMother(p.nChildren(), thisMothers, otherMothers);

	if(!motherID) {
		std::cout << p.lund() <<"[GenInfo] setMCtruth: requesting particle with ID = 0! Exiting." << std::endl;
		return;
	}

	Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
	const Gen_hepevt& thisGen = GenMgr(Panther_ID(motherID));
	// setting the relation
	p.relation().genHepevt(thisGen);
}

// number and type of final state particles should be the same
// Flags
// -11 : gen_hepevt link of one of children not found in list of daughters of matched gen_hepevt (something wrong)
// -10 : one of children has not link to gen_hepevt (something wrong)
//  -9 : as -10 but in addition size of reconstructed < size of generated
//  -5 : more reconstructed than generated particle in final state (something wrong)
//  -2 : gen_hepvt particle doesn't have any daughters (something wrong in maching)
//  -1 : random combination (common mother is virtual gamma (ccbar, uds mc), or Upsilon(4S), Upsilon(5S))
//   0 : no link to gen_hepevt particle
//   1 : particle is correctly reconstructed; including correct particle id of final state particles (SIGNAL)
//   2 : one or more FSP are misidentified, but have common mother
//   3 : FSP have common mother, but at least one massive particle is missing
//   4 : FSP have common mother, but at least one massless particle is missing (radiative photon)
//   5 : final state includes pi0 without perfect match to gen_hepevt
//   6 : ID = 2 and 5 are true
//  10 : particle is correctly reconstructed; including correct particle id of final state particles, but FSR photon is missing (SIGNAL)
//  11 : in addition to FSR photon missing one more radiative photon is missing
//  20 : missing neutrino
//  21 : missing neutrino and radiative photon
//  23 : missing neutrino and massive particle
//  24 : missing neutrino and another massles particle (FSR photon)

int getMCtruthFlag(Particle &p) {
	if(IDhep(p)==0)
		return 0;

	if(p.lund()==111) 
		return getMCtruthPi0Flag(p);

	if(p.mdstCharged())
		return getMCtruthChargedTrackFlag(p);

	Gen_hepevt thisGen = p.relation().genHepevt();

	int motherIDhep = thisGen.idhep();
	if(motherIDhep==10022 || motherIDhep==300553 || motherIDhep==9000553) {
		return -1;
	}

	// Ks doesn't have daughters in gen_hepevt table
	if(motherIDhep==310) {
		return 1;
	}

	std::vector<Particle>   reconstructed;
	std::vector<Gen_hepevt> generated;

	appendRecFSP(p, reconstructed);
	int a = appendGenFSP(thisGen, generated);
	if(a==0) {
		return -2;
	}

	int truth = compareFinalStates(reconstructed, generated);
	return truth;
}



void setMCtruth(std::vector<Particle> &plist){
	if(plist.size() == 0) return;

	for(std::vector<Particle>::iterator i = plist.begin();
			i != plist.end(); ++i){
		setMCtruth(*i);
	}
}
void fillMothers(Particle &A, Particle &B, std::vector<int> &thisMothers, std::vector<int> &otherMothers) {
	Gen_hepevt motherThis(A.genHepevt());

	while(motherThis.mother()) {
		motherThis = motherThis.mother();
		thisMothers.push_back(motherThis.get_ID());
	}

	motherThis = B.genHepevt();
	while(motherThis.mother()) {
		motherThis = motherThis.mother();
		otherMothers.push_back(motherThis.get_ID());
	}
}

void fillMothers(Particle &A, Particle &B, Particle &C, std::vector<int> &thisMothers, std::vector<int> &otherMothers) {
	Gen_hepevt motherThis(A.genHepevt());

	while(motherThis.mother()) {
		motherThis = motherThis.mother();
		thisMothers.push_back(motherThis.get_ID());
	}

	motherThis = B.genHepevt();
	while(motherThis.mother()) {
		motherThis = motherThis.mother();
		otherMothers.push_back(motherThis.get_ID());
	}

	motherThis = C.genHepevt();
	while(motherThis.mother()) {
		motherThis = motherThis.mother();
		otherMothers.push_back(motherThis.get_ID());
	}
}

int getCommonMother(Particle &A, Particle &B) {
	if(!A.genHepevt())
		return 0;
	if(!B.genHepevt())
		return 0;

	std::vector<int> thisMothers;
	std::vector<int> otherMothers;

	fillMothers(A, B, thisMothers, otherMothers);
	if(thisMothers.size()==0 || otherMothers.size()==0) {
		return 0;
	}

	int motherID = findCommonMother(2, thisMothers, otherMothers);

	if(!motherID) {
		std::cout << "[GenInfo] setMCtruth: requesting particle with ID = 0! Exiting." << std::endl;
		return 0;
	}

	Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
	const Gen_hepevt& thisGen = GenMgr(Panther_ID(motherID));

	return thisGen.idhep();
}

int getCommonMother(Particle &A, Particle &B, Particle &C) {
	if(!A.genHepevt())
		return 0;
	if(!B.genHepevt())
		return 0;
	if(!C.genHepevt())
		return 0;

	std::vector<int> thisMothers;
	std::vector<int> otherMothers;

	fillMothers(A, B, C, thisMothers, otherMothers);
	if(thisMothers.size()==0 || otherMothers.size()==0) {
		return 0;
	}

	int motherID =  findCommonMother(3, thisMothers, otherMothers);

	if(!motherID) {
		std::cout << "[GenInfo] setMCtruth: requesting particle with ID = 0! Exiting." << std::endl;
		return 0;
	}

	Gen_hepevt_Manager& GenMgr = Gen_hepevt_Manager::get_manager();
	const Gen_hepevt& thisGen = GenMgr(Panther_ID(motherID));

	return thisGen.idhep();
}

////////////////////////////////////////////////////////////////////////////////
// get decay chain for final stat particle

void genDecayChain(Particle p, int* dChain) {
	for(int i=0; i<=8; i++) dChain[i] = -1;

	if(p.relation().genHepevt()) {
		Gen_hepevt igen = p.relation().genHepevt();
		dChain[0] = igen.idhep();

		Gen_hepevt imot = igen.mother();
		if(imot) {
			dChain[1] = imot.idhep();
			dChain[2] = imot.daLast()-imot.daFirst()+1;

			Gen_hepevt immot = imot.mother();
			if(immot) {
				dChain[3] = immot.idhep();
				dChain[4] = immot.daLast()-immot.daFirst()+1;

				Gen_hepevt rg_mmmot = immot.mother();
				if(rg_mmmot) {
					dChain[5] = rg_mmmot.idhep();
					dChain[6] = rg_mmmot.daLast()-rg_mmmot.daFirst()+1;

					Gen_hepevt mrg_mmmot = rg_mmmot.mother();
					if(mrg_mmmot) {
						dChain[7] = mrg_mmmot.idhep();
						dChain[8] = mrg_mmmot.daLast()-mrg_mmmot.daFirst()+1;
					}
				}
			}
		}
	}
}

int IDhep(Particle &part) {

	if(! part.genHepevt()) return 0;
	else if (part.genHepevt().idhep() == IDHEP_Y4S) return IDHEP_Y4S_REDEFINED;
	return part.genHepevt().idhep();
}

int NdecayProd(Particle &part) {

	if(! part.genHepevt()) return 0;
	return part.genHepevt().daLast() - part.genHepevt().daFirst() +1;
}

int getEventMCFlag(const int& candSel, const int& B0Flag, const int& DSFlag, const int& rhoFlag, const int& D0Flag) {

	if ((int(B0Flag) == 1 || int(B0Flag) == 10) &&
		(int(DSFlag) == 1 || int(DSFlag) == 10) &&
		(int(rhoFlag) == 1 || int(rhoFlag) == 10) &&
		(int(D0Flag == 1) || int(D0Flag == 10)))
		return 1;

	if ((int(DSFlag) == 1 || int(DSFlag) == 10) &&
		!(int(rhoFlag) == 1 || int(rhoFlag) == 10) &&
		!(int(rhoFlag) == 3) &&
		(int(D0Flag == 1) || int(D0Flag == 10)))
		return 2;

	if ((int(rhoFlag) == 1 || int(rhoFlag) == 10) &&
		!(int(D0Flag == 1) || int(D0Flag == 10)))
		return 3;

	if (!(int(DSFlag) == 1 || int(DSFlag) == 10) &&
		(int(rhoFlag) == 1 || int(rhoFlag) == 10) &&
		(int(D0Flag == 1) || int(D0Flag == 10)))
		return 4;

	if (!(int(DSFlag) == 1 || int(DSFlag) == 10) &&
		!(int(rhoFlag) == 1 || int(rhoFlag) == 10) &&
		!(int(D0Flag == 1) || int(D0Flag == 10)))
		return 5;

	if (!(int(B0Flag) == 1 || int(B0Flag) == 10) &&
		(int(DSFlag) == 1 || int(DSFlag) == 10) &&
		(int(D0Flag == 1) || int(D0Flag == 10)) &&
		(int(rhoFlag) == 3))
		return 6;

	if (!(int(candSel) == 0) &&
		(int(B0Flag) == 1 || int(B0Flag) == 10) &&
		(int(DSFlag) == 1 || int(DSFlag) == 10) &&
		(int(D0Flag == 1) || int(D0Flag == 10)) &&
		(int(rhoFlag) == 3))
		return 7;

	if ((int(candSel) == 0) &&
		(int(B0Flag) == 1 || int(B0Flag) == 10) &&
		(int(DSFlag) == 1 || int(DSFlag) == 10) &&
		(int(D0Flag == 1) || int(D0Flag == 10)) &&
		(int(rhoFlag) == 3))
		return 8;

	return -1;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
