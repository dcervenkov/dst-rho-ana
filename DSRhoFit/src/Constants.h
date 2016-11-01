#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED

#include "TString.h"

static const double PI = 3.141592;

static const int B0_IDHEP = 511;
static const int DS_IDHEP = 413;
static const int RHO_IDHEP = 213;
static const int D0_IDHEP = 421;
static const int PI0_IDHEP = 111;
static const int PI_IDHEP = 211;
static const int PHOTON_IDHEP = 22;

static const double Btau = 1.53439;// [ps]
static const double Bdm = 0.507;// [h-bar (ps)^-1]

static const double c = 0.29979; // [mm/ps]

static const double Bgamma = 1/Btau/c; // 2.17242; // [h-bar (mm/c)^-1]
static const double Bfreq = Bdm/c; //1.69117; //[h-bar/c^2 (mm/c)^-1] number taken directly from evtgen (in debug mode)

const double kThetaBMin = 0.5;
const double kThetaBMax = 2.95;

const TString format = ".pdf";
#endif
