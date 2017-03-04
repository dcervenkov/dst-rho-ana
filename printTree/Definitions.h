#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <cstring>

// Include Belle
#include "belle.h"

// Include constants
#include "Constants.h"

// Include general MDST-tables 
#include "panther/panther.h"
#include "tables/belletdf.h"
#include "tables/mdst.h"
#include "mdst/mdst.h"

// Include general Particle ...
#include "particle/Particle.h"
#include "particle/utility.h"

#if defined(BELLE_NAMESPACE)
namespace Belle
{
#endif

// TypeDefs   
typedef std::vector<Particle> Particles;

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif // DEFINITIONS_H
