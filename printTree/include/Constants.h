#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cstring>

#if defined(BELLE_NAMESPACE)
namespace Belle
{
#endif

// CONSTANTS - PARTICLES
const static double PIMASS         = 0.13957018; // GeV
const static double BMASS          = 5.27915;    // GeV
const static double BZEROMASS      = 5.27953;    // GeV
const static double ETACMASS       = 2.9803;     // GeV
const static double KZEROMASS      = 0.497614;   // GeV
const static double ETACWIDTH      = 0.0286;     // GeV

// CONSTANTS - EXPERIMENT
const static double ONRESLUM       = 710.526;
const static double ONRESENERGY    = 10.578;     // GeV
const static double OFFRESLUM      = 88.145;
const static double OFFRESENERGY   = 10.518;     // GeV
const static int    NBBAR          = 771581000;

// CONSTANTS - CUTS
const static double PIPIDCUT       = 0.9;
const static double KPIDCUT        = 0.6;
const static double EIDCUT         = 0.95;
const static double PPIDCUT        = 0.4;

const static double BRLPHOTONCUT   = 0.05;  // GeV
const static double FWDPHOTONCUT   = 0.10;  // GeV
const static double FWDPHOTONTCUT  = 0.20;  // GeV (Tighter cut)
const static double PI0CHI2CUT     = 50.0;

const static double IMPACTRPHICUT  = 0.5;   // cm
const static double IMPACTZCUT     = 3.0;   // cm

const static double KSHORTMASSCUTL = 0.482; // GeV
const static double KSHORTMASSCUTR = 0.514; // GeV

const static double PI0MASSCUTL    = 0.118; // GeV
const static double PI0MASSCUTR    = 0.150; // GeV

const static double ETACMASSCUTL   = 2.92;  // GeV
const static double ETACMASSCUTR   = 3.04;  // GeV

const static double MBCCUTL        = 5.20;  // GeV
const static double MBCCUTR        = 5.30;  // GeV
const static double MBCSIGCUTL     = 5.271; // GeV
const static double MBCSIGCUTR     = 5.29;  // GeV

const static double DECUT          = 0.25;  // GeV
const static double DESIGCUT       = 0.035; // GeV

// Random generator
const static int    RNDGENSEED     = 12345;

// Enums  

// 
// CHARGE
//   
enum CHARGE { ZERO=0, PLUS=1, MINUS=-1 };

//
// CWN - Number of candidates
//         
const static int BASFCWNCAN = 55;

enum CHTYPE { ETACK_KSKPI=0, ETACK_KKPI0=1, ETACK_PPBAR=2, ETACKS_KSKPI=10, ETACKS_KKPI0=11, ETACKS_PPBAR=12 };
   
//
// MCFLAG
//
enum MCFLAG { DATA=0, MC=1 };

//
// MC TRUTH INFORMATION
//
// -1 = not assigned yet
//  0 = one or more final states particles (FSP) have no link
//  1 = final state correctly reconstructed (final state might have scattered --> the same final state)
//  2 = final state correctly reconstructed (without final state radiation - gammas missing)
//  10= final state correctly reconstructed, but charge is opposite
//  20= final state correclty reconstructed (without final state radiation - gammas missing), but charge opposite
//  3 = final state assigned to one of it's decay products (created by Geant4 in material)
//  4 = one of more FSP misidentified, but final inv. mass OK
//  5 = one or more FSP misidentified
//  6 = just PID incorrectly assigned
//  7 = missing particle
//
enum MCTRUTH { NOASSIGN=-1, NOLINK=0, OK=1, OKWFR=2, OKCHRG=10, OKWFRCHRG=20, OKDEC=3, OKMISID=4, MISID=5, MISPID=6, MISSING=7};

//
// PARTICLE LUND NUMBER <-> name               
//                                                
// 0 = not assigned yet
//
enum PARTICLES { NOPART=0  , ELECTRON=11 , MUON   =13  , GAMMA  =22, 
                 PIZERO=111, PION    =211, RHO0   =113 , RHO    =213,
                 KSHORT=310, KLONG   =130, KZERO  =311 , KAON   =321,
                 ETAP  =331, ETAC    =441, DZERO  =421 , DSTAR  =413, 
		 BZERO =511, BMEZON  =521, PROTON =2212, NEUTRON=2112 };

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif // CONSTANTS_H
