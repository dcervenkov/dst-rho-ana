/* vertex_info.h */


#ifndef VFIT_INFO_H
#define VFIT_INFO_H

#include "belle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

struct vfit_info
{
	bool         m_usable;
	int          m_status;
	HepPoint3D   m_pos;
	HepSymMatrix m_err;
	Hep3Vector   m_epar;
	Hep3Vector   m_eneg;
	Hep3Vector   m_epos;
	double       m_chi2;
	double		 m_effxi;
	int          m_ndf;
	int          m_effndf;
	int          m_ntrk;
	double       m_tmp[4];
	int          m_stat;

	vfit_info() :
		m_usable  ( false ),
		m_status  ( 0 ),
		m_pos     ( HepPoint3D(0,0,0) ),
		m_err     ( HepSymMatrix(3,0) ),
		m_epar    ( Hep3Vector(0.,0.,0.) ),
		m_eneg    ( Hep3Vector(0.,0.,0.) ),
		m_epos    ( Hep3Vector(0.,0.,0.) ),
		m_chi2    ( 0 ),
		m_effxi	  ( 0 ),
		m_ndf     ( 0 ),
		m_effndf  ( 0 ),
		m_ntrk    ( 0 ),
		m_stat    ( 0 )
		{ m_tmp[0] = m_tmp[1] = m_tmp[2] = m_tmp[3] = 0; }
};


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif /* VFIT_INFO_H */

