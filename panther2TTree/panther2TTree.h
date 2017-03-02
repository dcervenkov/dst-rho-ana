//===================================================
//
//  panther2TTree
//
//  Program to export panther table data from
//  EvtGen mdst file to ROOT TTree
//
//  Created:  July 2011 - D.Cervenkov
//
//===================================================
//#ifndef __PANTHER_2_TTREE__
//#define __PANTHER_2_TTREE__

//#include "belle.h"
//#include "tuple/BelleTupleManager.h"
//#include "event/BelleEvent.h"
//#include "basf/module.h"
//
//#include "tables/belletdf.h"
//#include "tables/hepevt.h"
//#include "tables/mdst.h"
//
//#include "particle/Particle.h"
//
//#include "belleutil/debugout.h"
//
//#include <iostream>
//#include <vector>
#include <fstream>

#include "belle.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "particle/Particle.h"
#include "tuple/BelleTupleManager.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

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

static const int maxParticles = 300;

/*
class MyParticle
{
   public:
      MyParticle()
      {	 itsId=0;itsId=0;itsIsthep=0;itsIdhep=0;itsMother=0;itsMo1=0;itsMo2=0;itsDa1=0;itsDa2=0;
	 itsP[0]=0;itsP[1]=0;itsP[2]=0;itsP[3]=0;itsP[4]=0;itsV[0]=0;itsV[1]=0;itsV[2]=0;itsV[3]=0;
      }
      ~MyParticle();

      void SetId(int id) {itsId = id;}
      void SetIsthep(int isthep) {itsIsthep = isthep;}
      void SetIdhep(int idhep) {itsIdhep = idhep;}
      void SetMother(int mother) {itsMother = mother;}
      void SetMo1(int mo1) {itsMo1 = mo1;}
      void SetMo2(int mo2) {itsMo2 = mo2;}
      void SetDa1(int da1) {itsDa1 = da1;}
      void SetDa2(int da2) {itsDa2 = da2;}
      void SetP(double p, int n) {itsP[n] = p;}
      void SetW(double v, int n) {itsV[n] = v;}

      int GetId() {return itsId;}
      int GetIsthep() {return itsIsthep;}
      int GetIdhep() {return itsIdhep;}
      int GetMother() {return itsMother;}
      int GetMo1() {return itsMo1;}
      int GetMo2() {return itsMo2;}
      int GetDa1() {return itsDa1;}
      int GetDa2() {return itsDa2;}
      double GetP(int n) {return itsP[n];}
      double GetW(int n) {return itsV[n];}

   private:
      int itsId;
      int itsIsthep;
      int itsIdhep;
      int itsMother;
      int itsMo1;
      int itsMo2;
      int itsDa1;
      int itsDa2;
      double itsP[5];
      double itsV[4];

   public:
      ClassDef(MyParticle,1)
};
*/

// Module class definition
class panther2TTreeModule : public Module
{  
 public:

// Constructor & destructor
   panther2TTreeModule();
   ~panther2TTreeModule();

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
   std::ofstream outFile;

//ROOT stuff
   TFile *rootFile;
   TTree *rootTree; 
   //TH1F *h1;

   //MyParticle particles[maxParticles];

   int r_numParticles;

   int r_id[maxParticles];
   int r_isthep[maxParticles];
   int r_idhep[maxParticles];
   int r_mother[maxParticles];
   int r_mo1[maxParticles];
   int r_mo2[maxParticles];
   int r_da1[maxParticles];
   int r_da2[maxParticles];
   double r_p[maxParticles][5];
   double r_v[maxParticles][4];

};

extern "C" Module_descr * mdcl_panther2TTree ()
{

	//   // Get instance of the module, wrap in a Module_descr
	panther2TTreeModule * module  = new panther2TTreeModule;
	Module_descr * moduleDscr = new Module_descr("panther2TTree", module);

	// Parameters:
	moduleDscr->define_param("fileName","Output file name" ,sizeof(module->outFileName)/sizeof(char), module->outFileName);
	moduleDscr->define_param("strVars" ,"Variables to plot",sizeof(module->strVars)/sizeof(char)    , module->strVars );
	moduleDscr->define_param("qqOnly"  ,"QQ event only"    ,&module->qqOnly     );
	moduleDscr->define_param("decDepth","Decay depth"      ,&module->decDepth   );
	moduleDscr->define_param("whichB"  ,"Which B to output",&module->whichB     );
	moduleDscr->define_param("runNo"   ,"Run number"       ,&module->runNo      );
	moduleDscr->define_param("evtNo"   ,"Event number"     ,&module->evtNo      );
	moduleDscr->define_param("nEvts"   ,"Number of events" ,&module->nEvts      );

	return moduleDscr;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

//#endif // __PANTHER_2_TTREE__
