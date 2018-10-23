/**
 *  @file    constants.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-04-05
 *
 *  @brief Holds all the declarations of constants required by the program
 *
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

namespace constants {

extern const double pi;
extern const double ap;
extern const double at;
extern const double a0;
extern const char format[];
extern const int btype;
extern const double dm;
extern const double tau;
extern const double fit_range_dt_low;
extern const double fit_range_dt_high;
extern const double c;

namespace cuts {

extern const double dt_low;
extern const double dt_high;
extern const double thetab_low;
extern const double thetab_high;
extern const int sig_vtx_h;
extern const int tag_vtx_h;
extern const double sig_vtx_multitrack_sigma_z;
extern const double sig_vtx_singletrack_sigma_z;
extern const double tag_vtx_multitrack_sigma_z;
extern const double tag_vtx_singletrack_sigma_z;

} // namespace cuts

} // namespace constants

#endif /* CONSTANTS_H_ */
