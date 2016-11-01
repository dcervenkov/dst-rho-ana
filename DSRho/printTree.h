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
#ifndef PRINTTREE_H_
#define PRINTTREE_H_

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
static const char PD='P';
static const char ID='i';
static const char UD='u';


// Module class definition
class printTree {
 public:
   printTree(char* fileName);
   ~printTree();

   void print(HepAList<Gen_hepevt> list, int iExp, int iRun, int iEvt);
   void printParticleChildren(Particle& p);
   void printParticle(Particle& p);
   void printPiParticle(Particle& p);
   void printText(const char*);
   std::string getRealName(int idhep);

 private:
   ofstream outFile;

   // Parameters --> set from outside
   char outFileName[128]; // Output filename
   char strVars[11];      // String of variables to be plot

   int qqOnly;   // QQ event only
   int decDepth; // Decay length
   int whichB;   // Which B
   int runNo;    // Run number to print
   int evtNo;    // Event number to print
   int nEvts;    // + nEvts following events

   int level;


};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif // PRINTTREE_H_
