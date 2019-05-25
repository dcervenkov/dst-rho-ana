/**
 *  @file    constants.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-04-05
 *
 *  @brief Holds all the constants required by the program
 *
 */

#include "constants.h"

// Standard includes
#include <array>

// ROOT includes
#include "TMath.h"

namespace constants {

// Default transversity and cartesian parameters
const std::array<double, 16> par_input = {0.269, // ap
                                          0.56,  // apa
                                          0.941, // a0
                                          3.11,  // ata
                                          +0.00817, // xp
                                          +0.00533, // x0
                                          +0.00826, // xt
                                          -0.00661, // yp
                                          -0.00846, // y0
                                          +0.00372, // yt
                                          -0.01022, // xbp
                                          -0.00846, // xb0
                                          -0.00584, // xbt
                                          +0.00244, // ybp
                                          +0.00533, // yb0
                                          -0.00692};// ybt

// Transversity amplitudes
const double ap = 0.269;
const double a0 = 0.941;
const double at = 0.205;

// Format of the resulting images
const char format[] = ".pdf";

// B type: 0 == B0 and 1 == B+
const int btype = 0;

const double dm = 0.507;     // [ps^-1]
const double tau = 1.53439;  // [ps]
const double c = 0.0299792458;  // [cm/ps]

const double fraction_cr_of_crscf = 0.8585;
const double fraction_cr_of_crscfbkg = 0.7901;
const double fraction_scf_of_crscfbkg = 0.1302;

namespace cuts {

const double de_low = -0.14;  // [GeV]
const double de_high = 0.10;  // [GeV]

const double dt_low = -10;  // [ps]
const double dt_high = 10;  // [ps]

const double thetat_low = 0;            // [rad]
const double thetat_high = TMath::Pi(); // [rad]
const double thetab_low = 0.5;          // [rad]
const double thetab_high = 2.95;        // [rad]
const double phit_low = -TMath::Pi();   // [rad]
const double phit_high = TMath::Pi();   // [rad]

const double cs_bdtg = -0.6;

const int sig_vtx_h = 50;
const int tag_vtx_h = 50;

const double sig_vtx_multitrack_sigma_z = 0.02;   // [cm]
const double sig_vtx_singletrack_sigma_z = 0.05;  // [cm]
const double tag_vtx_multitrack_sigma_z = 0.02;   // [cm]
const double tag_vtx_singletrack_sigma_z = 0.05;  // [cm]

}  // namespace cuts

}  // namespace constants
