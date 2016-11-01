/*
 * constants.h
 *
 *  Created on: Aug 20, 2014
 *      Author: cervenkov
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

// The maximum lenght of a char array used for storing tuple names
static const int MAX_NAME_LENGHT = 100;

static const double PI = 3.14159265358979323846;

// The Y4S IDHEP code is redefined to a lower value, because
// the IDHEP plots are otherwise too "zoomed out"
static const int IDHEP_Y4S = 300553;
static const int IDHEP_Y4S_REDEFINED = 353;

// BASF identifies histograms via numbers
static const int BASF_HIST_ID_BASE = 1000; // The base for all the MCFlag separated histograms; BASE + MCFLAG
static const int BASF_HIST_ID_ALL = 1500; // All candidates together
static const int BASF_HIST_ID_BEST = 2000; // Best candidate

static const int NUM_MC_FLAGS = 19;
static const int MC_FLAGS_ASSIGNMENT[NUM_MC_FLAGS] = { -11, -10, -9, -5, -2, -1, 0, 1, 2, 3, 4, 5, 6, 10, 11, 20, 21,
		23, 24 }; // See geninfo.cc

// Signal MC flags; used for figure-of-merit calculation in DSRhoPeek
static const int MC_FLAG_SIGNAL = 1;
static const int MC_FLAG_SIGNAL_WITHOUT_FSR = 10;

// Nominal Mbc; used for Best Canditate Selection
static const double MBC = 5.28;

// Cuts
static const double MBC_CUT = 5.15;
static const double DELTA_E_CUT = 0.25;

// PID cuts
static const double KAON_ID = 0.2;
static const double PION_ID = 0.8;

// Impact parameter cuts
static const double IMPACT_PAR_DR = 0.5;
static const double IMPACT_PAR_DZ = 5;

// Cut on the mass difference between D* and its daughter D0
static const double DSD0_MASS_DIFF_LOW = 0.141;
static const double DSD0_MASS_DIFF_HIGH = 0.150;

// Mass window cuts; PDG mass +- window
static const double D0_MASS_WINDOW = 0.030;
static const double RHO_MASS_WINDOW = 0.500;

#endif /* CONSTANTS_H_ */
