/*
 * constants.cc
 *
 *  Created on: Jul 21, 2015
 *      Author: cervenkov
 */

#include "constants.h"

// The maximum lenght of a char array used for storing tuple names
extern const int MAX_NAME_LENGHT = 100;

extern const double PI = 3.14159265358979323846;

// The Y4S IDHEP code is redefined to a lower value, because
// the IDHEP plots are otherwise too "zoomed out"
extern const int IDHEP_Y4S = 300553;
extern const int IDHEP_Y4S_REDEFINED = 353;

extern const int IDHEP_B0 = 511;
extern const int IDHEP_DS = 413;
extern const int IDHEP_RHO = 213;
extern const int IDHEP_D0 = 421;
extern const int IDHEP_PI = 211;
extern const int IDHEP_PI0 = 111;
extern const int IDHEP_K = 321;
extern const int IDHEP_GAMMA = 22;

// BASF identifies histograms via numbers
extern const int BASF_HIST_ID_BEST = 2000; // Best candidate

// Nominal Mbc; used for Best Canditate Selection
// The sigma of Mbc is given mainly by beam energy measurement: 3 MeV.
extern const double MBC = 5.28;
extern const double MBCSIGMA = 0.003;

extern const double PI0MASS = 0.135;
extern const double PI0WIDTH = 0.0063;

extern const double FWD_ENDCAP_ANGLE = 0.55;
extern const double BKW_ENDCAP_ANGLE = 2.27;

//----------------- Cuts ------------------
// Values are in GeV, cm or dimensionless
extern const double MBC_CUT = 5.27;
extern const double MBC_SIDEBAND_CUT = 5.24;
extern const double DELTA_E_CUT_LOW = -0.15;
extern const double DELTA_E_CUT_HIGH = 0.10;

// PID cuts
extern const double KAON_ID = 0.2;
extern const double PION_ID = 0.8;

// Impact parameter cuts
extern const double IMPACT_PAR_DR = 0.5;
extern const double IMPACT_PAR_DZ = 5;

// Cuts on the minimal number of hits in the r-phi a z planes
extern const int N_HITS_SVD_R = 2;
extern const int N_HITS_SVD_Z = 2;

// Gamma energy cuts
extern const double GAMMA_E_CUT_BARREL = 0.05;
extern const double GAMMA_E_CUT_ENDCAP = 0.10;

// Cut on the mass difference between D* and its daughter D0
extern const double DSD0_MASS_DIFF_LOW = 0.143;
extern const double DSD0_MASS_DIFF_HIGH = 0.148;
extern const double DSD0_TIGHT_MASS_DIFF_LOW = 0.144;
extern const double DSD0_TIGHT_MASS_DIFF_HIGH = 0.147;

// Mass cuts; PDG mass +- window
extern const double D0_MASS_WINDOW = 0.020;
extern const double D0_PI0_MASS_WINDOW_LOW = 0.031;
extern const double D0_PI0_MASS_WINDOW_HIGH = 0.025;
extern const double D0_3PI_MASS_WINDOW = 0.014;
extern const double RHO_MASS_WINDOW_LOW = 0.340;
extern const double RHO_MASS_WINDOW_HIGH = 0.220;
extern const double PI0_MASS_CUT_LOW = 0.115;
extern const double PI0_MASS_CUT_HIGH = 0.150;

// Number of B candidates
extern const unsigned MAX_NUM_B_CANDIDATES = 5;

// Pi0s chi2 cuts
extern const int RHO_PI0_CHI2_CUT = 17;
extern const int D0_PI0_CHI2_CUT = 26;

// Continuum suppression cut
extern const double CS_CUT = -1000.6;

// Theta_b cut; see joural
extern const double THETA_B_CUT = 0.5;

// Track cuts
extern const int TRACK_NHITS_RPHI = 1;
extern const int TRACK_NHITS_Z = 2;

// Tag side track cuts
extern const double TAG_TRACK_DZ_CUT = -0.18;
extern const double TAG_TRACK_DZ_ERROR_CUT = 0.05;
extern const double TAG_TRACK_DR_CUT = 0.05; // Tokyo cut

