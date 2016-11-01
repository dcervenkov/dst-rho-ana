#include "UserInfo.h"

#include "belle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  
  // For Interface to Set UserInfo Class

  void setUserInfo(Particle &p, unsigned ch)
  {
    if(! &p.userInfo() )// return;
      p.userInfo(UserInfo(ch));
  }

  void setUserInfo(std::vector<Particle> &p, unsigned ch) {
    for(unsigned i=0;i<p.size();++i) setUserInfo(p[i],ch);
  }

  // UserInfo Class
  
  UserInfo::UserInfo()
    :     m_decayMode(0),
	  m_type(0),
    
	  m_dr(-999.),
	  m_dz(-999.),
	  m_eid(-1.),
	  m_pid(-1.),
	  m_muid(-1.),
	  m_protonID(-1.),

	  m_massbvf(0.0),
	  m_min_e9e25(0.0),
	  m_max_e9e25(0.0),
	  m_min_energy(0.0),
	  m_max_energy(0.0),
	  m_min_theta(-999.),
	  m_max_theta(-999.),
	  
	  m_decay_angle(-999.),

	  m_thetaT(-999.),
	  m_phiT(-999.),
	  m_thetaB(-999.),

	  m_mbc(-999.),
	  m_deltae(-999.)
  {
  }
  
  UserInfo::UserInfo(unsigned ch)
    :    m_decayMode(0),
	 m_type(0),
    
	 m_dr(-1.),
	 m_dz(-1.),
	 m_eid(-1.),
	 m_pid(-1.),
	 m_muid(-1.),
	 m_protonID(-1.),

	 m_massbvf(0.0),
	 m_min_e9e25(0.0),
	 m_max_e9e25(0.0),
	 m_min_energy(0.0),
	 m_max_energy(0.0),
	 m_min_theta(-999.),
	 m_max_theta(-999.),
    
	 m_decay_angle(-999.),

	 m_thetaT(-999.),
	 m_phiT(-999.),
	 m_thetaB(-999.),

	 m_mbc(-999.),
	 m_deltae(-999.)
  {
  }
  
  UserInfo::UserInfo(const UserInfo &x)
    :   	m_decayMode(x.m_decayMode),
		m_type(x.m_type),
	
		m_dr(x.m_dr),
		m_dz(x.m_dz),
		m_eid(x.m_eid),
		m_pid(x.m_pid),
		m_muid(x.m_muid),
		m_protonID(x.m_protonID),
	
		m_massbvf(x.m_massbvf),
		m_min_e9e25(x.m_min_e9e25),
		m_max_e9e25(x.m_max_e9e25),
		m_min_energy(x.m_min_energy),
		m_max_energy(x.m_max_energy),
		m_min_theta(x.m_min_theta),
		m_max_theta(x.m_max_theta),
	
		m_decay_angle(x.m_decay_angle),

		m_thetaT(x.m_thetaT),
		m_phiT(x.m_phiT),
		m_thetaB(x.m_thetaB),

		m_mbc(x.m_mbc),
		m_deltae(x.m_deltae)
  {
  }
  UserInfo::~UserInfo()
  {
  }
  
  UserInfo*
  UserInfo::clone(void) const
  {
    UserInfo *x = new UserInfo(*this);
    return x;
  }
  
  UserInfo &
  UserInfo::operator = (const UserInfo &x)
  {
    m_decayMode = x.m_decayMode;
    m_type   = x.m_type;
    
    m_dr     = x.m_dr;
    m_dz     = x.m_dz;
    m_eid = x.m_eid;
    m_pid = x.m_pid;
    m_muid= x.m_muid;
    m_protonID = x.m_protonID;

    m_massbvf   = x.m_massbvf;
    m_min_e9e25 = x.m_min_e9e25;
    m_max_e9e25 = x.m_max_e9e25;
    m_min_energy = x.m_min_energy;
    m_max_energy = x.m_max_energy;
    m_min_theta = x.m_min_theta;
    m_max_theta = x.m_max_theta;
    
    m_decay_angle = x.m_decay_angle;

    m_thetaT = x.m_thetaT;
	m_phiT = x.m_phiT;
	m_thetaB = x.m_thetaB;

    m_mbc         = x.m_mbc;
    m_deltae      = x.m_deltae;
    return *this;
  }
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
