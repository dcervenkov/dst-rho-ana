#include "belle.h"
#include "event/BelleEvent.h"

#include "basf/module.h"
#include "basf/module_descr.h"

#include "mdst/mdst.h"

#include <cstring>
#include <iostream>
#include <iomanip>
#include <strings.h>
#include <bitset>

// Belle libraries
#include "tables/belletdf.h"
#include "tables/hepevt.h"
#include "tables/mdst.h"

// CLHEP libraries
#include "CLHEP/Alist/AList.h"

// Particle
#include "particle/utility.h"

#include "printTree.h"
#include "geninfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

printTree::printTree(char* channel) {
	sprintf(outFileName, "%s.txt", channel);
	sprintf(strVars, "MXYZEi");

	qqOnly = 1;
	whichB = 0;
	decDepth = 999;
	runNo = 0;
	evtNo = 338;
	nEvts = 1;

	level = 0;

	outFile.open(outFileName);
	if (!outFile) {
		std::cerr << __FILE__ << " " << __LINE__ << ": ";
		std::cerr << "printTree: Couldn't open output file " << outFileName << std::endl;
	}

	if (qqOnly)
		outFile << " Will print QQ particles only" << std::endl;
	else
		outFile << " Will print all particles    " << std::endl;

	outFile << " Variables to print are:     " << strVars << std::endl;

	if (whichB == 0) {
		outFile << " Print both B mezons         " << std::endl;
	} else {
		outFile << " Print B #" << whichB << std::endl;
	}

	outFile << " Will show decay depth of:   " << decDepth << std::endl;

	outFile << "\n********************************\n\n" << std::endl;
}

printTree::~printTree() {
	outFile.close();
}

void printTree::print(HepAList<Gen_hepevt> list, int iExp, int iRun, int iEvt) {

//	outFile << "DBG: Entering print\n";
	int i=0;
	int converge=0;
	int level=0;
	//int mom;
	int varlen=0;

	int  iB=0;
	bool dooutput=true;

	std::string leaders("");
	std::string blank("    ");
	std::string filled("|   ");

	outFile.precision(4);

	//
	// Header
	outFile << std::endl << std::setiosflags(std::ios::left) << std::setw(49+13*std::strlen(strVars)+1) << std::setfill('=') << "" << std::setfill(' ') << std::endl;
	outFile << " Decay Tree - Exp #" << std::setw(2) << iExp;
	outFile << " Run #"              << std::setw(4) << iRun;
	outFile << " Event #"            << std::setw(5) << iEvt;
	outFile << std::setw(5)          << "";

	if (strchr(strVars,ID)) { outFile << "|" << std::setw(12) << "    ID"   ; varlen++; }
	if (strchr(strVars,PX)) { outFile << "|" << std::setw(12) << "    PX"   ; varlen++; }
	if (strchr(strVars,PY)) { outFile << "|" << std::setw(12) << "    PY"   ; varlen++; }
	if (strchr(strVars,PZ)) { outFile << "|" << std::setw(12) << "    PZ"   ; varlen++; }
	if (strchr(strVars,VX)) { outFile << "|" << std::setw(12) << "    VX"   ; varlen++; }
	if (strchr(strVars,VY)) { outFile << "|" << std::setw(12) << "    VY"   ; varlen++; }
	if (strchr(strVars,VZ)) { outFile << "|" << std::setw(12) << "    VZ"   ; varlen++; }
	if (strchr(strVars,E))  { outFile << "|" << std::setw(12) << "     E"   ; varlen++; }
	if (strchr(strVars,M))  { outFile << "|" << std::setw(12) << "     M"   ; varlen++; }
	if (strchr(strVars,ST)) { outFile << "|" << std::setw(12) << "   ISTHEP"; varlen++; }
	if (strchr(strVars,PD)) { outFile << "|" << std::setw(12) << "     PID" ; varlen++; }
	outFile << "|" << std::endl;
	outFile << std::setw(49+13*varlen+1) << std::setfill('-') << "" << std::setfill(' ') << std::endl;

	//
	// Particles
//	outFile << "DBG: Before while\n";
	while (converge<list.index(list.last())) {
//		outFile << "DBG: Looping...\n";

		converge++;
		Gen_hepevt& event = *list[i];

		//if event = B+/- or B0/bar
		if ((whichB!=0) && ((fabs(event.idhep())==521) || (fabs(event.idhep())==511))) {

			iB++;
			if (iB!=whichB) dooutput=false;
			else            dooutput=true;
		}

		//Modified by jwicht (plot particles if level<depth)
//		outFile << "DBG: Looping point 1\n";
//		outFile.flush();
		if ((level<decDepth)&&(dooutput)) {


			outFile << std::setw(12)         << setiosflags(std::ios::left) << getRealName(event.idhep()); // Tree
			outFile << std::setw(40-4*level) << "";                                                        // Some space
//			outFile << "DBG: Looping point 1.1\n";
//			outFile.flush();

			if (strchr(strVars,ID)) outFile << std::setw(10) << i+1            <<  std::setw(3) << "";
			if (strchr(strVars,PX)) outFile << std::setw(10) << event.PX()     <<  std::setw(3) << "";
			if (strchr(strVars,PY)) outFile << std::setw(10) << event.PY()     <<  std::setw(3) << "";
			if (strchr(strVars,PZ)) outFile << std::setw(10) << event.PZ()     <<  std::setw(3) << "";
			if (strchr(strVars,VX)) outFile << std::setw(10) << event.VX()     <<  std::setw(3) << "";
			if (strchr(strVars,VY)) outFile << std::setw(10) << event.VY()     <<  std::setw(3) << "";
			if (strchr(strVars,VZ)) outFile << std::setw(10) << event.VZ()     <<  std::setw(3) << "";
			if (strchr(strVars,E))  outFile << std::setw(10) << event.E()      <<  std::setw(3) << "";
			if (strchr(strVars,M))  outFile << std::setw(10) << event.M()      <<  std::setw(3) << "";
//			outFile << "DBG: Looping point 1.2\n";
//			outFile.flush();
			if (strchr(strVars,ST)) outFile << std::setw(10) << event.isthep() <<  std::setw(3) << "";
//			outFile << "DBG: Looping point 1.3\n";
//			outFile.flush();
			if (strchr(strVars,PD)) outFile << std::setw(10) << event.idhep()  <<  std::setw(3) << "";
//			outFile << "DBG: Looping point 1.4\n";
//			outFile.flush();
			outFile << std::endl;
//			outFile << "DBG: Looping point 1.5\n";
//			outFile.flush();
		}
//		outFile << "DBG: Looping point 2\n";
//		outFile.flush();
		//
		// Find next particle to draw

		// Has daughter (qqOnly?? or created by EvtGen)
		if (event.daFirst()	&& (!qqOnly || list[event.daFirst() - 1]->isthep() > 0)) {
//			outFile << "DBG: Looping point 3A\n";
//			outFile.flush();

			// Set leaders according to mother
			if (event.moFirst()) {

				if (list[event.moFirst() - 1]->daLast() - 1 != i)
					leaders += filled;
				else
					leaders += blank;
			}

			i = event.daFirst() - 1;
			level++;
		}
		// Has mother
		else if (event.moFirst()) {
//			outFile << "DBG: Looping point 3B\n";
//			outFile.flush();

			// Not last daughter
			if (list[event.moFirst() - 1]->daLast() - 1 != i) {
				i++;
//				outFile << "DBG: Looping point 3B.A\n";
//				outFile.flush();
			}
			// Last daughter: find next mother
			else {
//				outFile << "DBG: Looping point 3B.B1\n";
//				outFile.flush();

				i = event.moFirst() - 1; // Rewind to mother
				level--;

//				outFile << "DBG: Looping point 3B.B.info\n";
//				outFile << "DBG: i = " << i;
//				outFile.flush();
//				outFile << ", list[i]->moFirst() = " << list[i]->moFirst();
//				outFile.flush();

				if (list[i]->moFirst() == 0) return;

				while (list[list[i]->moFirst() - 1]->daLast() - 1 == i) {
//					outFile << "DBG: Looping point 3B.B.Looping\n";
//					outFile.flush();

					i = list[i]->moFirst() - 1; // Find mother with daughters
					level--;

					if (!i)
						return; // It happens... (when there's an orphan in the event)
				}
//				outFile << "DBG: Looping point 3B.B2\n";
//				outFile.flush();

				i++;
				leaders = leaders.substr(0, (level - 1) * filled.length());

//				outFile << "DBG: Looping point 3B.B3\n";
//				outFile.flush();

			}
		} else {
//			outFile << "DBG: Looping point 3C\n";
//			outFile.flush();
			return; // No mother, no daughter! Poor girl...
		}
		// Prepare for next particle (modified by jwicht)
		if ((level < decDepth) && (dooutput))
			outFile << std::setw(0) << leaders << "+-> ";

//		outFile << "DBG: Looping point 4\n";
//		outFile.flush();

	} // While
//	outFile << "DBG: After while\n";
//	outFile.flush();

	//
	// Last particle has a separate treatment
	Gen_hepevt& event = *list[i];

	outFile << std::setw(10)         << setiosflags(std::ios::left) << getRealName(event.idhep()); // Tree
	outFile << std::setw(40-4*level) << "";                                                        // Some space
	if (strchr(strVars,ID)) outFile << std::setw(10) << i+1            <<  std::setw(3) << "";
	if (strchr(strVars,PX)) outFile << std::setw(10) << event.PX()     <<  std::setw(3) << "";
	if (strchr(strVars,PY)) outFile << std::setw(10) << event.PY()     <<  std::setw(3) << "";
	if (strchr(strVars,PZ)) outFile << std::setw(10) << event.PZ()     <<  std::setw(3) << "";
	if (strchr(strVars,VX)) outFile << std::setw(10) << event.VX()     <<  std::setw(3) << "";
	if (strchr(strVars,VY)) outFile << std::setw(10) << event.VY()     <<  std::setw(3) << "";
	if (strchr(strVars,VZ)) outFile << std::setw(10) << event.VZ()     <<  std::setw(3) << "";
	if (strchr(strVars,E))  outFile << std::setw(10) << event.E()      <<  std::setw(3) << "";
	if (strchr(strVars,M))  outFile << std::setw(10) << event.M()      <<  std::setw(3) << "";
	if (strchr(strVars,ST)) outFile << std::setw(10) << event.isthep() <<  std::setw(3) << "";
	if (strchr(strVars,PD)) outFile << std::setw(10) << event.idhep()  <<  std::setw(3) << "";
	outFile << std::endl;

//	outFile << "DBG: Exiting print\n";
//	outFile.flush();
}

void printTree::printParticleChildren(Particle& p) {

	printParticle(p);

	level++;

	//printParticleChildren(p);
	for (int i = 0; i < (int)p.nChildren(); ++i) {
		if(p.child(i).nChildren()){// && p.child(i).lund()!=111 && p.child(i).lund()!=310) {
			printParticleChildren(p.child(i));
		} else {
			printParticle(p.child(i));
		}
	}
	outFile << std::endl;
	level--;
}

void printTree::printParticle(Particle& p) {
	outFile << std::setw(2*level) << "" << "name = " << getRealName(p.lund());
	outFile << ", ID = ";
	if (p.genHepevt()) {
		outFile << p.genHepevt().get_ID();
	} else {
		outFile << "N/A";
	}
	//outFile << ", px = " << p.px() << ", pz = " << p.pz()
	outFile << ", MCFlag = " << getMCtruthFlag(p);

//	if (p.genHepevt()) {
//		outFile  << ", gpx = " << p.genHepevt().PX()  << ", gpz = " << p.genHepevt().PZ();
//	} else {
//		outFile  << ", gpx = " << "N/A"  << ", gpz = " << "N/A";
//	}

	outFile << std::endl;
}

void printTree::printPiParticle(Particle& p) {
	outFile << "good_charged = " << good_charged(p.mdstCharged());
	outFile << std::endl;
}

void printTree::printText(const char* text) {
	outFile << text;
}

// Return realName
std::string printTree::getRealName(int idhep) {

	std::string partName;

	// Got it from http://root.cern.ch/root/html225/src/TAttParticle.cxx.html
	switch (idhep) {

	case 1 : partName =  "d"; break;
	case -1 : partName =  "anti-d"; break;
	case 2 : partName =  "u"; break;
	case -2 : partName =  "anti-u"; break;
	case 3 : partName =  "s"; break;
	case -3 : partName =  "anti-s"; break;
	case 4 : partName =  "c"; break;
	case -4 : partName =  "anti-c"; break;
	case 5 : partName =  "b"; break;
	case -5 : partName =  "anti-b"; break;
	case 6 : partName =  "t"; break;
	case -6 : partName =  "anti-t"; break;
	case 7 : partName =  "b'"; break;
	case -7 : partName =  "anti-b'"; break;
	case 8 : partName =  "t'"; break;
	case -8 : partName =  "anti-t'"; break;
	case 11 : partName =  "e-"; break;
	case -11 : partName =  "e+"; break;
	case 12 : partName =  "nu_e"; break;
	case -12 : partName =  "anti-nu_e"; break;
	case 13 : partName =  "mu-"; break;
	case -13 : partName =  "mu+"; break;
	case 14 : partName =  "nu_mu"; break;
	case -14 : partName =  "anti-nu_mu"; break;
	case 15 : partName =  "tau-"; break;
	case -15 : partName =  "tau+"; break;
	case 16 : partName =  "nu_tau"; break;
	case -16 : partName =  "anti-nu_tau"; break;
	case 17 : partName =  "L-"; break;
	case -17 : partName =  "L+"; break;
	case 18 : partName =  "nu_L"; break;
	case -18 : partName =  "anti-nu_L"; break;
	case 21 : partName =  "g"; break;
	case 22 : partName =  "gamma"; break;
	case 10022 : partName =  "vpho"; break;
	case 20022 : partName =  "Cerenkov"; break;
	case 30022 : partName =  "radgam"; break;
	case 23 : partName =  "Z0"; break;
	case 24 : partName =  "W+"; break;
	case -24 : partName =  "W-"; break;
	case 25 : partName =  "Higgs0"; break;
	case 110 : partName =  "reggeon"; break;
	case 990 : partName =  "pomeron"; break;
	case 32 : partName =  "Z'0"; break;
	case 33 : partName =  "Z''0"; break;
	case 34 : partName =  "W'+"; break;
	case -34 : partName =  "W'-"; break;
	case 35 : partName =  "Higgs'0"; break;
	case 36 : partName =  "A0"; break;
	case 37 : partName =  "Higgs+"; break;
	case -37 : partName =  "Higgs-"; break;
	case 41 : partName =  "R0"; break;
	case -41 : partName =  "anti-R0"; break;
	case 61 : partName =  "Xu0"; break;
	case 62 : partName =  "Xu+"; break;
	case -62 : partName =  "Xu-"; break;
	case 81 : partName =  "specflav"; break;
	case 82 : partName =  "rndmflav"; break;
	case -82 : partName =  "anti-rndmflav"; break;
	case 83 : partName =  "phasespa"; break;
	case 84 : partName =  "c-hadron"; break;
	case -84 : partName =  "anti-c-hadron"; break;
	case 85 : partName =  "b-hadron"; break;
	case -85 : partName =  "anti-b-hadron"; break;
	case 86 : partName =  "t-hadron"; break;
	case -86 : partName =  "anti-t-hadron"; break;
	case 87 : partName =  "b'-hadron"; break;
	case -87 : partName =  "anti-b'-hadron"; break;
	case 88 : partName =  "t'-hadron"; break;
	case -88 : partName =  "anti-t'-hadron"; break;
	case 89 : partName =  "Wvirt+"; break;
	case -89 : partName =  "Wvirt-"; break;
	case 90 : partName =  "diquark"; break;
	case -90 : partName =  "anti-diquark"; break;
	case 91 : partName =  "cluster"; break;
	case 92 : partName =  "string"; break;
	case 93 : partName =  "indep"; break;
	case 94 : partName =  "CMshower"; break;
	case 95 : partName =  "SPHEaxis"; break;
	case 96 : partName =  "THRUaxis"; break;
	case 97 : partName =  "CLUSjet"; break;
	case 98 : partName =  "CELLjet"; break;
	case 99 : partName =  "table"; break;
	case 100 : partName =  "geantino"; break;
	case 1101 : partName =  "dd_0"; break;
	case -1101 : partName =  "anti-dd_0"; break;
	case 2101 : partName =  "ud_0"; break;
	case -2101 : partName =  "anti-ud_0"; break;
	case 2201 : partName =  "uu_0"; break;
	case -2201 : partName =  "anti-uu_0"; break;
	case 3101 : partName =  "sd_0"; break;
	case -3101 : partName =  "anti-sd_0"; break;
	case 3201 : partName =  "su_0"; break;
	case -3201 : partName =  "anti-su_0"; break;
	case 3301 : partName =  "ss_0"; break;
	case -3301 : partName =  "anti-ss_0"; break;
	case 4101 : partName =  "cd_0"; break;
	case -4101 : partName =  "anti-cd_0"; break;
	case 4201 : partName =  "cu_0"; break;
	case -4201 : partName =  "anti-cu_0"; break;
	case 4301 : partName =  "cs_0"; break;
	case -4301 : partName =  "anti-cs_0"; break;
	case 4401 : partName =  "cc_0"; break;
	case -4401 : partName =  "anti-cc_0"; break;
	case 5101 : partName =  "bd_0"; break;
	case -5101 : partName =  "anti-bd_0"; break;
	case 5201 : partName =  "bu_0"; break;
	case -5201 : partName =  "anti-bu_0"; break;
	case 5301 : partName =  "bs_0"; break;
	case -5301 : partName =  "anti-bs_0"; break;
	case 5401 : partName =  "bc_0"; break;
	case -5401 : partName =  "anti-bc_0"; break;
	case 5501 : partName =  "bb_0"; break;
	case -5501 : partName =  "anti-bb_0"; break;
	case 1103 : partName =  "dd_1"; break;
	case -1103 : partName =  "anti-dd_1"; break;
	case 2103 : partName =  "ud_1"; break;
	case -2103 : partName =  "anti-ud_1"; break;
	case 2203 : partName =  "uu_1"; break;
	case -2203 : partName =  "anti-uu_1"; break;
	case 3103 : partName =  "sd_1"; break;
	case -3103 : partName =  "anti-sd_1"; break;
	case 3203 : partName =  "su_1"; break;
	case -3203 : partName =  "anti-su_1"; break;
	case 3303 : partName =  "ss_1"; break;
	case -3303 : partName =  "anti-ss_1"; break;
	case 4103 : partName =  "cd_1"; break;
	case -4103 : partName =  "anti-cd_1"; break;
	case 4203 : partName =  "cu_1"; break;
	case -4203 : partName =  "anti-cu_1"; break;
	case 4303 : partName =  "cs_1"; break;
	case -4303 : partName =  "anti-cs_1"; break;
	case 4403 : partName =  "cc_1"; break;
	case -4403 : partName =  "anti-cc_1"; break;
	case 5103 : partName =  "bd_1"; break;
	case -5103 : partName =  "anti-bd_1"; break;
	case 5203 : partName =  "bu_1"; break;
	case -5203 : partName =  "anti-bu_1"; break;
	case 5303 : partName =  "bs_1"; break;
	case -5303 : partName =  "anti-bs_1"; break;
	case 5403 : partName =  "bc_1"; break;
	case -5403 : partName =  "anti-bc_1"; break;
	case 5503 : partName =  "bb_1"; break;
	case -5503 : partName =  "anti-bb_1"; break;
	case 9910113 : partName =  "rho_diff0"; break;
	case 9910211 : partName =  "pi_diff+"; break;
	case -9910211 : partName =  "pi_diff-"; break;
	case 9910223 : partName =  "omega_diff"; break;
	case 9910333 : partName =  "phi_diff"; break;
	case 9910443 : partName =  "psi_diff"; break;
	case 9912112 : partName =  "n_diffr"; break;
	case -9912112 : partName =  "anti-n_diffr"; break;
	case 9912212 : partName =  "p_diff+"; break;
	case -9912212 : partName =  "anti-p_diff-"; break;
	case 1011 : partName =  "deuteron"; break;
	case -1011 : partName =  "anti-deuteron"; break;
	case 1021 : partName =  "tritium"; break;
	case -1021 : partName =  "anti-tritium"; break;
	case 1012 : partName =  "He3"; break;
	case -1012 : partName =  "anti-He3"; break;
	case 1022 : partName =  "alpha"; break;
	case -1022 : partName =  "anti-alpha"; break;
	case 111 : partName =  "pi0"; break;
	case 211 : partName =  "pi+"; break;
	case -211 : partName =  "pi-"; break;
	case 10111 : partName =  "a_00"; break;
	case 10211 : partName =  "a_0+"; break;
	case -10211 : partName =  "a_0-"; break;
	case 100111 : partName =  "pi(2S)0"; break;
	case 100211 : partName =  "pi(2S)+"; break;
	case -100211 : partName =  "pi(2S)-"; break;
	case 113 : partName =  "rho0"; break;
	case 213 : partName =  "rho+"; break;
	case -213 : partName =  "rho-"; break;
	case 10113 : partName =  "b_10"; break;
	case 10213 : partName =  "b_1+"; break;
	case -10213 : partName =  "b_1-"; break;
	case 20113 : partName =  "a_10"; break;
	case 20213 : partName =  "a_1+"; break;
	case -20213 : partName =  "a_1-"; break;
	case 100113 : partName =  "rho(2S)0"; break;
	case 100213 : partName =  "rho(2S)+"; break;
	case -100213 : partName =  "rho(2S)-"; break;
	case 30113 : partName =  "rho(3S)0"; break;
	case 30213 : partName =  "rho(3S)+"; break;
	case -30213 : partName =  "rho(3S)-"; break;
	case 115 : partName =  "a_20"; break;
	case 215 : partName =  "a_2+"; break;
	case -215 : partName =  "a_2-"; break;
	case 9000221 : partName =  "f_0(600)"; break;
	case 221 : partName =  "eta"; break;
	case 331 : partName =  "eta'"; break;
	case 10221 : partName =  "f_0"; break;
	case 100221 : partName =  "eta(2S)"; break;
	case 10331 : partName =  "f'_0"; break;
	case 9020221 : partName =  "eta(1405)"; break;
	case 9030221 : partName =  "f_0(1500)"; break;
	case 223 : partName =  "omega"; break;
	case 333 : partName =  "phi"; break;
	case 10223 : partName =  "h_1"; break;
	case 20223 : partName =  "f_1"; break;
	case 10333 : partName =  "h'_1"; break;
	case 20333 : partName =  "f'_1"; break;
	case 100223 : partName =  "omega(2S)"; break;
	case 100333 : partName =  "phi(1680)"; break;
	case 225 : partName =  "f_2"; break;
	case 335 : partName =  "f'_2"; break;
	case 310 : partName =  "K_S0"; break;
	case 130 : partName =  "K_L0"; break;
	case 311 : partName =  "K0"; break;
	case -311 : partName =  "anti-K0"; break;
	case 321 : partName =  "K+"; break;
	case -321 : partName =  "K-"; break;
	case 90000311 : partName =  "K_0*(800)0"; break;
	case -90000311 : partName =  "anti-K_0*(800)0"; break;
	case 90000321 : partName =  "K_0*(800)+"; break;
	case -90000321 : partName =  "K_0*(800)-"; break;
	case 10311 : partName =  "K_0*0"; break;
	case -10311 : partName =  "anti-K_0*0"; break;
	case 10321 : partName =  "K_0*+"; break;
	case -10321 : partName =  "K_0*-"; break;
	case 30343 : partName =  "Xsd"; break;
	case -30343 : partName =  "anti-Xsd"; break;
	case 30353 : partName =  "Xsu"; break;
	case -30353 : partName =  "anti-Xsu"; break;
	case 30363 : partName =  "Xss"; break;
	case -30363 : partName =  "anti-Xss"; break;
	case 30643 : partName =  "Xdd"; break;
	case -30643 : partName =  "anti-Xdd"; break;
	case 30653 : partName =  "Xdu+"; break;
	case -30653 : partName =  "anti-Xdu-"; break;
	case 313 : partName =  "K*0"; break;
	case -313 : partName =  "anti-K*0"; break;
	case 323 : partName =  "K*+"; break;
	case -323 : partName =  "K*-"; break;
	case 10313 : partName =  "K_10"; break;
	case -10313 : partName =  "anti-K_10"; break;
	case 10323 : partName =  "K_1+"; break;
	case -10323 : partName =  "K_1-"; break;
	case 20313 : partName =  "K'_10"; break;
	case -20313 : partName =  "anti-K'_10"; break;
	case 20323 : partName =  "K'_1+"; break;
	case -20323 : partName =  "K'_1-"; break;
	case 100313 : partName =  "K'*0"; break;
	case -100313 : partName =  "anti-K'*0"; break;
	case 100323 : partName =  "K'*+"; break;
	case -100323 : partName =  "K'*-"; break;
	case 30313 : partName =  "K''*0"; break;
	case -30313 : partName =  "anti-K''*0"; break;
	case 30323 : partName =  "K''*+"; break;
	case -30323 : partName =  "K''*-"; break;
	case 315 : partName =  "K_2*0"; break;
	case -315 : partName =  "anti-K_2*0"; break;
	case 325 : partName =  "K_2*+"; break;
	case -325 : partName =  "K_2*-"; break;
	case 10315 : partName =  "K_2(1770)0"; break;
	case -10315 : partName =  "anti-K_2(1770)0"; break;
	case 10325 : partName =  "K_2(1770)+"; break;
	case -10325 : partName =  "K_2(1770)-"; break;
	case 20315 : partName =  "K_2(1820)0"; break;
	case -20315 : partName =  "anti-K_2(1820)0"; break;
	case 20325 : partName =  "K_2(1820)+"; break;
	case -20325 : partName =  "K_2(1820)-"; break;
	case 317 : partName =  "K_3*0"; break;
	case -317 : partName =  "anti-K_3*0"; break;
	case 327 : partName =  "K_3*+"; break;
	case -327 : partName =  "K_3*-"; break;
	case 319 : partName =  "K_4*0"; break;
	case -319 : partName =  "anti-K_4*0"; break;
	case 329 : partName =  "K_4*+"; break;
	case -329 : partName =  "K_4*-"; break;
	case 411 : partName =  "D+"; break;
	case -411 : partName =  "D-"; break;
	case 421 : partName =  "D0"; break;
	case -421 : partName =  "anti-D0"; break;
	case 10411 : partName =  "D_0*+"; break;
	case -10411 : partName =  "D_0*-"; break;
	case 10421 : partName =  "D_0*0"; break;
	case -10421 : partName =  "anti-D_0*0"; break;
	case 100411 : partName =  "D(2S)+"; break;
	case -100411 : partName =  "D(2S)-"; break;
	case 100421 : partName =  "D(2S)0"; break;
	case -100421 : partName =  "anti-D(2S)0"; break;
	case 413 : partName =  "D*+"; break;
	case -413 : partName =  "D*-"; break;
	case 423 : partName =  "D*0"; break;
	case -423 : partName =  "anti-D*0"; break;
	case 10413 : partName =  "D_1+"; break;
	case -10413 : partName =  "D_1-"; break;
	case 10423 : partName =  "D_10"; break;
	case -10423 : partName =  "anti-D_10"; break;
	case 20413 : partName =  "D'_1+"; break;
	case -20413 : partName =  "D'_1-"; break;
	case 20423 : partName =  "D'_10"; break;
	case -20423 : partName =  "anti-D'_10"; break;
	case 100413 : partName =  "D*(2S)+"; break;
	case -100413 : partName =  "D*(2S)-"; break;
	case 100423 : partName =  "D*(2S)0"; break;
	case -100423 : partName =  "anti-D*(2S)0"; break;
	case 415 : partName =  "D_2*+"; break;
	case -415 : partName =  "D_2*-"; break;
	case 425 : partName =  "D_2*0"; break;
	case -425 : partName =  "anti-D_2*0"; break;
	case 431 : partName =  "D_s+"; break;
	case -431 : partName =  "D_s-"; break;
	case 10431 : partName =  "D_s0*+"; break;
	case -10431 : partName =  "D_s0*-"; break;
	case 433 : partName =  "D_s*+"; break;
	case -433 : partName =  "D_s*-"; break;
	case 10433 : partName =  "D_s1+"; break;
	case -10433 : partName =  "D_s1-"; break;
	case 20433 : partName =  "D'_s1+"; break;
	case -20433 : partName =  "D'_s1-"; break;
	case 435 : partName =  "D_s2*+"; break;
	case -435 : partName =  "D_s2*-"; break;
	case 9000433 : partName =  "D_sj(2700)+"; break;
	case -9000433 : partName =  "D_sj(2700)-"; break;
	case 441 : partName =  "eta_c"; break;
	case 10441 : partName =  "chi_c0"; break;
	case 100441 : partName =  "eta_c(2S)"; break;
	case 443 : partName =  "J/psi"; break;
	case 10443 : partName =  "h_c"; break;
	case 20443 : partName =  "chi_c1"; break;
	case 100443 : partName =  "psi(2S)"; break;
	case 30443 : partName =  "psi(3770)"; break;
	case 9000443 : partName =  "psi(4040)"; break;
	case 9010443 : partName =  "psi(4160)"; break;
	case 9020443 : partName =  "psi(4415)"; break;
	case 445 : partName =  "chi_c2"; break;
	case 120443 : partName =  "X(3872)"; break;
	case 90000443 : partName =  "Y(3940)"; break;
	case 91000443 : partName =  "X(3940)"; break;
	case 511 : partName =  "B0"; break;
	case -511 : partName =  "anti-B0"; break;
	case 521 : partName =  "B+"; break;
	case -521 : partName =  "B-"; break;
	case 10511 : partName =  "B_0*0"; break;
	case -10511 : partName =  "anti-B_0*0"; break;
	case 10521 : partName =  "B_0*+"; break;
	case -10521 : partName =  "B_0*-"; break;
	case 513 : partName =  "B*0"; break;
	case -513 : partName =  "anti-B*0"; break;
	case 523 : partName =  "B*+"; break;
	case -523 : partName =  "B*-"; break;
	case 10513 : partName =  "B_10"; break;
	case -10513 : partName =  "anti-B_10"; break;
	case 10523 : partName =  "B_1+"; break;
	case -10523 : partName =  "B_1-"; break;
	case 20513 : partName =  "B'_10"; break;
	case -20513 : partName =  "anti-B'_10"; break;
	case 20523 : partName =  "B'_1+"; break;
	case -20523 : partName =  "B'_1-"; break;
	case 515 : partName =  "B_2*0"; break;
	case -515 : partName =  "anti-B_2*0"; break;
	case 525 : partName =  "B_2*+"; break;
	case -525 : partName =  "B_2*-"; break;
	case 531 : partName =  "B_s0"; break;
	case -531 : partName =  "anti-B_s0"; break;
	case 10531 : partName =  "B_s0*0"; break;
	case -10531 : partName =  "anti-B_s0*0"; break;
	case 533 : partName =  "B_s*0"; break;
	case -533 : partName =  "anti-B_s*0"; break;
	case 10533 : partName =  "B_s10"; break;
	case -10533 : partName =  "anti-B_s10"; break;
	case 20533 : partName =  "B'_s10"; break;
	case -20533 : partName =  "anti-B'_s10"; break;
	case 535 : partName =  "B_s2*0"; break;
	case -535 : partName =  "anti-B_s2*0"; break;
	case 541 : partName =  "B_c+"; break;
	case -541 : partName =  "B_c-"; break;
	case 10541 : partName =  "B_c0*+"; break;
	case -10541 : partName =  "B_c0*-"; break;
	case 543 : partName =  "B_c*+"; break;
	case -543 : partName =  "B_c*-"; break;
	case 10543 : partName =  "B_c1+"; break;
	case -10543 : partName =  "B_c1-"; break;
	case 20543 : partName =  "B'_c1+"; break;
	case -20543 : partName =  "B'_c1-"; break;
	case 545 : partName =  "B_c2*+"; break;
	case -545 : partName =  "B_c2*-"; break;
	case 551 : partName =  "eta_b"; break;
	case 10551 : partName =  "chi_b0"; break;
	case 100551 : partName =  "eta_b(2S)"; break;
	case 110551 : partName =  "chi_b0(2P)"; break;
	case 200551 : partName =  "eta_b(3S)"; break;
	case 210551 : partName =  "chi_b0(3P)"; break;
	case 553 : partName =  "Upsilon"; break;
	case 10553 : partName =  "h_b"; break;
	case 20553 : partName =  "chi_b1"; break;
	case 30553 : partName =  "Upsilon_1(1D)"; break;
	case 100553 : partName =  "Upsilon(2S)"; break;
	case 110553 : partName =  "h_b(2P)"; break;
	case 120553 : partName =  "chi_b1(2P)"; break;
	case 130553 : partName =  "Upsilon_1(2D)"; break;
	case 200553 : partName =  "Upsilon(3S)"; break;
	case 210553 : partName =  "h_b(3P)"; break;
	case 220553 : partName =  "chi_b1(3P)"; break;
	case 300553 : partName =  "Upsilon(4S)"; break;
	case 9000553 : partName =  "Upsilon(5S)"; break;
	case 555 : partName =  "chi_b2"; break;
	case 10555 : partName =  "eta_b2(1D)"; break;
	case 20555 : partName =  "Upsilon_2(1D)"; break;
	case 100555 : partName =  "chi_b2(2P)"; break;
	case 110555 : partName =  "eta_b2(2D)"; break;
	case 120555 : partName =  "Upsilon_2(2D)"; break;
	case 200555 : partName =  "chi_b2(3P)"; break;
	case 557 : partName =  "Upsilon_3(1D)"; break;
	case 100557 : partName =  "Upsilon_3(2D)"; break;
	case 2212 : partName =  "p+"; break;
	case -2212 : partName =  "anti-p-"; break;
	case 2112 : partName =  "n0"; break;
	case -2112 : partName =  "anti-n0"; break;
	case 12212 : partName =  "N(1440)+"; break;
	case -12212 : partName =  "anti-N(1440)-"; break;
	case 12112 : partName =  "N(1440)0"; break;
	case -12112 : partName =  "anti-N(1440)0"; break;
	case 2124 : partName =  "N(1520)+"; break;
	case -2124 : partName =  "anti-N(1520)-"; break;
	case 1214 : partName =  "N(1520)0"; break;
	case -1214 : partName =  "anti-N(1520)0"; break;
	case 22212 : partName =  "N(1535)+"; break;
	case -22212 : partName =  "anti-N(1535)-"; break;
	case 22112 : partName =  "N(1535)0"; break;
	case -22112 : partName =  "anti-N(1535)0"; break;
	case 2224 : partName =  "Delta++"; break;
	case -2224 : partName =  "anti-Delta--"; break;
	case 2214 : partName =  "Delta+"; break;
	case -2214 : partName =  "anti-Delta-"; break;
	case 2114 : partName =  "Delta0"; break;
	case -2114 : partName =  "anti-Delta0"; break;
	case 1114 : partName =  "Delta-"; break;
	case -1114 : partName =  "anti-Delta+"; break;
	case 32224 : partName =  "Delta(1600)++"; break;
	case -32224 : partName =  "anti-Delta(1600)--"; break;
	case 32214 : partName =  "Delta(1600)+"; break;
	case -32214 : partName =  "anti-Delta(1600)-"; break;
	case 32114 : partName =  "Delta(1600)0"; break;
	case -32114 : partName =  "anti-Delta(1600)0"; break;
	case 31114 : partName =  "Delta(1600)-"; break;
	case -31114 : partName =  "anti-Delta(1600)+"; break;
	case 2222 : partName =  "Delta(1620)++"; break;
	case -2222 : partName =  "anti-Delta(1620)--"; break;
	case 2122 : partName =  "Delta(1620)+"; break;
	case -2122 : partName =  "anti-Delta(1620)-"; break;
	case 1212 : partName =  "Delta(1620)0"; break;
	case -1212 : partName =  "anti-Delta(1620)0"; break;
	case 1112 : partName =  "Delta(1620)-"; break;
	case -1112 : partName =  "anti-Delta(1620)+"; break;
	case 3122 : partName =  "Lambda0"; break;
	case -3122 : partName =  "anti-Lambda0"; break;
	case 13122 : partName =  "Lambda(1405)0"; break;
	case -13122 : partName =  "anti-Lambda(1405)0"; break;
	case 3124 : partName =  "Lambda(1520)0"; break;
	case -3124 : partName =  "anti-Lambda(1520)0"; break;
	case 23122 : partName =  "Lambda(1600)0"; break;
	case -23122 : partName =  "anti-Lambda(1600)0"; break;
	case 33122 : partName =  "Lambda(1670)0"; break;
	case -33122 : partName =  "anti-Lambda(1670)0"; break;
	case 13124 : partName =  "Lambda(1690)0"; break;
	case -13124 : partName =  "anti-Lambda(1690)0"; break;
	case 43122 : partName =  "Lambda(1800)0"; break;
	case -43122 : partName =  "anti-Lambda(1800)0"; break;
	case 53122 : partName =  "Lambda(1810)0"; break;
	case -53122 : partName =  "anti-Lambda(1810)0"; break;
	case 3126 : partName =  "Lambda(1820)0"; break;
	case -3126 : partName =  "anti-Lambda(1820)0"; break;
	case 13126 : partName =  "Lambda(1830)0"; break;
	case -13126 : partName =  "anti-Lambda(1830)0"; break;
	case 3222 : partName =  "Sigma+"; break;
	case -3222 : partName =  "anti-Sigma-"; break;
	case 3212 : partName =  "Sigma0"; break;
	case -3212 : partName =  "anti-Sigma0"; break;
	case 3112 : partName =  "Sigma-"; break;
	case -3112 : partName =  "anti-Sigma+"; break;
	case 3224 : partName =  "Sigma*+"; break;
	case -3224 : partName =  "anti-Sigma*-"; break;
	case 3214 : partName =  "Sigma*0"; break;
	case -3214 : partName =  "anti-Sigma*0"; break;
	case 3114 : partName =  "Sigma*-"; break;
	case -3114 : partName =  "anti-Sigma*+"; break;
	case 13212 : partName =  "Sigma(1660)0"; break;
	case -13212 : partName =  "anti-Sigma(1660)0"; break;
	case 13214 : partName =  "Sigma(1670)0"; break;
	case -13214 : partName =  "anti-Sigma(1670)0"; break;
	case 23212 : partName =  "Sigma(1750)0"; break;
	case -23212 : partName =  "anti-Sigma(1750)0"; break;
	case 3216 : partName =  "Sigma(1775)0"; break;
	case -3216 : partName =  "anti-Sigma(1775)0"; break;
	case 3322 : partName =  "Xi0"; break;
	case -3322 : partName =  "anti-Xi0"; break;
	case 3312 : partName =  "Xi-"; break;
	case -3312 : partName =  "anti-Xi+"; break;
	case 3324 : partName =  "Xi*0"; break;
	case -3324 : partName =  "anti-Xi*0"; break;
	case 3314 : partName =  "Xi*-"; break;
	case -3314 : partName =  "anti-Xi*+"; break;
	case 13324 : partName =  "Xi(1820)0"; break;
	case -13324 : partName =  "anti-Xi(1820)0"; break;
	case 13314 : partName =  "Xi(1820)-"; break;
	case -13314 : partName =  "anti-Xi(1820)+"; break;
	case 3334 : partName =  "Omega-"; break;
	case -3334 : partName =  "anti-Omega+"; break;
	case 13334 : partName =  "Omega(2250)-"; break;
	case -13334 : partName =  "anti-Omega(2250)+"; break;
	case 4122 : partName =  "Lambda_c+"; break;
	case -4122 : partName =  "anti-Lambda_c-"; break;
	case 14122 : partName =  "Lambda_c(2593)+"; break;
	case -14122 : partName =  "anti-Lambda_c(2593)-"; break;
	case 4124 : partName =  "Lambda_c(2625)+"; break;
	case -4124 : partName =  "anti-Lambda_c(2625)-"; break;
	case 14124 : partName =  "Lambda_c(2765)+"; break;
	case -14124 : partName =  "anti-Lambda_c(2765)-"; break;
	case 24122 : partName =  "Lambda_c(2880)+"; break;
	case -24122 : partName =  "anti-Lambda_c(2880)-"; break;
	case 34122 : partName =  "Lambda_c(2940)+"; break;
	case -34122 : partName =  "anti-Lambda_c(2940)-"; break;
	case 4222 : partName =  "Sigma_c++"; break;
	case -4222 : partName =  "anti-Sigma_c--"; break;
	case 4212 : partName =  "Sigma_c+"; break;
	case -4212 : partName =  "anti-Sigma_c-"; break;
	case 4112 : partName =  "Sigma_c0"; break;
	case -4112 : partName =  "anti-Sigma_c0"; break;
	case 4224 : partName =  "Sigma_c*++"; break;
	case -4224 : partName =  "anti-Sigma_c*--"; break;
	case 4214 : partName =  "Sigma_c*+"; break;
	case -4214 : partName =  "anti-Sigma_c*-"; break;
	case 4114 : partName =  "Sigma_c*0"; break;
	case -4114 : partName =  "anti-Sigma_c*0"; break;
	case 14224 : partName =  "Sigma_c(2800)++"; break;
	case -14224 : partName =  "anti-Sigma_c(2800)--"; break;
	case 14214 : partName =  "Sigma_c(2800)+"; break;
	case -14214 : partName =  "anti-Sigma_c(2800)-"; break;
	case 14114 : partName =  "Sigma_c(2800)0"; break;
	case -14114 : partName =  "anti-Sigma_c(2800)0"; break;
	case 4232 : partName =  "Xi_c+"; break;
	case -4232 : partName =  "anti-Xi_c-"; break;
	case 4132 : partName =  "Xi_c0"; break;
	case -4132 : partName =  "anti-Xi_c0"; break;
	case 4322 : partName =  "Xi'_c+"; break;
	case -4322 : partName =  "anti-Xi'_c-"; break;
	case 4312 : partName =  "Xi'_c0"; break;
	case -4312 : partName =  "anti-Xi'_c0"; break;
	case 4324 : partName =  "Xi_c*+"; break;
	case -4324 : partName =  "anti-Xi_c*-"; break;
	case 4314 : partName =  "Xi_c*0"; break;
	case -4314 : partName =  "anti-Xi_c*0"; break;
	case 14232 : partName =  "Xi_c(2790)+"; break;
	case -14232 : partName =  "anti-Xi_c(2790)-"; break;
	//case 14232 : partName =  "Xi_c(2790)0"; break;
	//case -14232 : partName =  "anti-Xi_c(2790)0"; break;
	case 24332 : partName =  "Xi_c(2815)+"; break;
	case -24332 : partName =  "anti-Xi_c(2815)-"; break;
	case 24132 : partName =  "Xi_c(2815)0"; break;
	case -24132 : partName =  "anti-Xi_c(2815)0"; break;
	case 34232 : partName =  "Xi_c(2980)+"; break;
	case -34232 : partName =  "anti-Xi_c(2980)-"; break;
	case 34132 : partName =  "Xi_c(2980)0"; break;
	case -34132 : partName =  "anti-Xi_c(2980)0"; break;
	case 44232 : partName =  "Xi_c(3080)+"; break;
	case -44232 : partName =  "anti-Xi_c(3080)-"; break;
	case 44132 : partName =  "Xi_c(3080)0"; break;
	case -44132 : partName =  "anti-Xi_c(3080)0"; break;
	case 4332 : partName =  "Omega_c0"; break;
	case -4332 : partName =  "anti-Omega_c0"; break;
	case 4334 : partName =  "Omega_c*0"; break;
	case -4334 : partName =  "anti-Omega_c*0"; break;
	case 5122 : partName =  "Lambda_b0"; break;
	case -5122 : partName =  "anti-Lambda_b0"; break;
	case 5112 : partName =  "Sigma_b-"; break;
	case -5112 : partName =  "anti-Sigma_b+"; break;
	case 5114 : partName =  "Sigma_b*-"; break;
	case -5114 : partName =  "anti-Sigma_b*+"; break;
	case 5212 : partName =  "Sigma_b0"; break;
	case -5212 : partName =  "anti-Sigma_b0"; break;
	case 5214 : partName =  "Sigma_b*0"; break;
	case -5214 : partName =  "anti-Sigma_b*0"; break;
	case 5222 : partName =  "Sigma_b+"; break;
	case -5222 : partName =  "anti-Sigma_b-"; break;
	case 5224 : partName =  "Sigma_b*+"; break;
	case -5224 : partName =  "anti-Sigma_b*-"; break;
	case 5132 : partName =  "Xi_b-"; break;
	case -5132 : partName =  "anti-Xi_b+"; break;
	case 5232 : partName =  "Xi_b0"; break;
	case -5232 : partName =  "anti-Xi_b0"; break;
	case 5312 : partName =  "Xi'_b-"; break;
	case -5312 : partName =  "anti-Xi'_b+"; break;
	case 5314 : partName =  "Xi_b*-"; break;
	case -5314 : partName =  "anti-Xi_b*+"; break;
	case 5322 : partName =  "Xi'_b0"; break;
	case -5322 : partName =  "anti-Xi'_b0"; break;
	case 5324 : partName =  "Xi_b*0"; break;
	case -5324 : partName =  "anti-Xi_b*0"; break;
	case 5332 : partName =  "Omega_b-"; break;
	case -5332 : partName =  "anti-Omega_b+"; break;
	case 5334 : partName =  "Omega_b*-"; break;
	case -5334 : partName =  "anti-Omega_b*+"; break;

	default:
		char buf[10];
		sprintf(buf,"?(%d)",idhep);
		partName  =  buf;      // isajet or pdg number does not exist
		break;
	}

	return partName;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
