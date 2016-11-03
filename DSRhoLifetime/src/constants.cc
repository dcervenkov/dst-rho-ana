/**
 *  @file    constants.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-04-05
 *
 *  @brief Holds all the constants required by the program
 *
 */

#include "constants.h"

namespace constants {

const double kPi = 3.14159265358979323846;

const double kAp = 0.269;
const double kA0 = 0.941;
const double kAt = 0.205;

const char format[] = ".pdf";

//B type: 0 == B0 and 1 == B+
const int btype = 0;

const double dm = 0.507; // [ps^-1]
const double tau = 1.53439; // [ps]

const double cut_dt_low = -20; // [ps]
const double cut_dt_high = 20; // [ps]

const double fit_range_dt_low = -10; // [ps]
const double fit_range_dt_high = 10; // [ps]

const double c = 0.0299792458; // [cm/ps]

namespace cuts {

const int sig_vtx_h = 50;
const int tag_vtx_h = 50;

const double sig_vtx_multitrack_sigma_z = 0.02; // [cm]
const double sig_vtx_singletrack_sigma_z = 0.05; // [cm]
const double tag_vtx_multitrack_sigma_z = 0.02; // [cm]
const double tag_vtx_singletrack_sigma_z = 0.05; // [cm]

} // namespace cuts

} // namespace constants

