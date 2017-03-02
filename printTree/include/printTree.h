//===================================================
//
//  Print tree
//
//  Routine to print decay tree
//
//  Created:  July 2001 - F.Ronga (IPHE)
//  Modified: J. Wicht
//            Z. Drasal --> compilable under gcc 4.4
//                      --> added parameters
//
//===================================================
//#ifndef __PRINT_TREE__
//#define __PRINT_TREE__

#include "belle.h"
#include "tuple/BelleTupleManager.h"
#include "event/BelleEvent.h"
#include "basf/module.h"

#include "tables/belletdf.h"
#include "tables/hepevt.h"
#include "tables/mdst.h"

#include "particle/Particle.h"

#include "belleutil/debugout.h"

#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

static const char PX='X';
static const char PY='Y';
static const char PZ='Z';
static const char VX='x';
static const char VY='y';
static const char VZ='z';
static const char E ='E';
static const char M ='M';
static const char ST='S';
static const char PI='P';
static const char ID='i';


// Module class definition
class printTreeModule : public Module
{  
 public:

// Constructor & destructor
   printTreeModule();
   ~printTreeModule();

// Initialization procedure - called once at the beginning
   void init( int* );

// Definition of histograms - called once at the beginning
   void hist_def( void );

// Process at the beginning of each run procedure - looped
   void begin_run( BelleEvent*, int* );

// Process each event procedure - looped
   void event( BelleEvent*, int* );

// Process at the end of each run - looped
   void end_run( BelleEvent*, int* ) {};

   void disp_stat( const char* ) {};
   void other( int*, BelleEvent*, int* ) {};

// Termination procedure - called once at the end
   void term( void );

// Parameters --> set from outside
   char outFileName[128]; // Output filename
   char strVars[11];      // String of variables to be plot

   int qqOnly;   // QQ event only
   int decDepth; // Decay length
   int whichB;   // Which B
   int runNo;    // Run number to print
   int evtNo;    // Event number to print
   int nEvts;    // + nEvts following events

 private:
  
//!Print the tree
   void print(HepAList<Gen_hepevt> list, int iExp, int iRun, int iEvt);

//!Get particle name based on idhep
   std::string getRealName(int idhep);
  
// Output file
   ofstream outFile;

   //DC: getEventInfo doesn't work for non-gsim data - returns only 0's; this is part of a workaround
   //int m_iEvt;

};

extern "C" Module_descr *mdcl_printtreemodule ()
{
    printTreeModule *module = new printTreeModule;
    Module_descr *dscr = new Module_descr ( "printtreemodule", module );

    dscr->define_param ("fileName", "Name of the file to save the output to.", 128, *(&module->outFileName));
    dscr->define_param ("qqOnly", "TODO: Write description.", &module->qqOnly);
    dscr->define_param ("strVars", "TODO: Write description.", 11, *(&module->strVars));
    dscr->define_param ("decDepth", "TODO: Write description.", &module->decDepth);
    dscr->define_param ("whichB", "TODO: Write description.", &module->whichB);
    dscr->define_param ("runNo", "Run number to process.", &module->runNo);
    dscr->define_param ("evtNo", "First to-be-processed event's number.", &module->evtNo);
    dscr->define_param ("nEvts", "Number of events to process.", &module->nEvts);

    basf_save_mdsc ( dscr );
    return dscr;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

//#endif // __PRINT_TREE__
