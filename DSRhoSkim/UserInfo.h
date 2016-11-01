#ifndef USERINFO_H_
#define USERINFO_H_

#include "belle.h"
#include "particle/Particle.h"
#include "particle/ParticleUserInfo.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

  
  // For Interface to Set UserInfo Class
  
  void setUserInfo(Particle &p, unsigned ch=0);
  void setUserInfo(std::vector<Particle>& p, unsigned ch=0);
  
  // UserInfo Class

  class UserInfo : public ParticleUserInfo
  {
  public:
    /// Default constructor
    UserInfo();
    
    UserInfo(unsigned);
    
    /// Copy constructor
    UserInfo(const UserInfo &);
    
    /// Destructor
    virtual ~UserInfo();
    
    /// constructs self object.
    UserInfo * clone(void) const;
    
    /// Copy operator
    UserInfo & operator = (const UserInfo &);
    
public:
    const int & decayMode(void) const { return m_decayMode; }
    void decayMode(const int &v) { m_decayMode = v; }

    void setType(const int &v) { m_type = v; }
    const int & getType(void) const { return m_type; }

    void setDr(const double &v) { m_dr = v; }
    const double & getDr(void) const { return m_dr; }

    void setDz(const double &v) { m_dz = v; }
    const double & getDz(void) const { return m_dz; }

    void setEid(const double &v) { m_eid = v; }
    void setPid(const double &v) { m_pid = v; }
    void setMuid(const double &v) { m_muid = v; }
    void setProtonID(const double &v) { m_protonID = v; }
    const double & getEid(void) const { return m_eid; }
    const double & getPid(void) const { return m_pid; }
    const double & getMuid(void) const { return m_muid; }
    const double & getProtonID(void) const { return m_protonID; }

    void setMassBeforeVertexFit(const double &v) { m_massbvf = v; }
    const double & getMassBeforeVertexFit(void) const { return m_massbvf; }
    void setMinE9E25(const double &v) { m_min_e9e25 = v; }
    const double & getMinE9E25(void) const { return m_min_e9e25; }
    void setMaxE9E25(const double &v) { m_max_e9e25 = v; }
    const double & getMaxE9E25(void) const { return m_max_e9e25; }
    void setMinEnergy(const double &v) { m_min_energy = v; }
    const double & getMinEnergy(void) const { return m_min_energy; }
    void setMaxEnergy(const double &v) { m_max_energy = v; }
    const double & getMaxEnergy(void) const { return m_max_energy; }
    void setMinTheta(const double &v) { m_min_theta = v; }
    const double & getMinTheta(void) const { return m_min_theta; }
    void setMaxTheta(const double &v) { m_max_theta = v; }
    const double & getMaxTheta(void) const { return m_max_theta; }
    
    void setDecayAngle(const double &v) { m_decay_angle = v; }
    const double & getDecayAngle(void) const { return m_decay_angle; }

    void setMbc(const double &v) { m_mbc = v; }
    const double & getMbc(void) const { return m_mbc; }
    
    void setDeltaE(const double &v) { m_deltae = v; }
    const double & getDeltaE(void) const { return m_deltae; }

    void setTransversityAngles(const double& thetaT, const double& phiT, const double& thetaB) {
    	m_thetaT = thetaT;
    	m_phiT = phiT;
    	m_thetaB = thetaB;
    }
    const double& getThetaT(void) { return m_thetaT; }
    const double& getPhiT(void) { return m_phiT; }
    const double& getThetaB(void) { return m_thetaB; }

private:
	int m_decayMode;

	int      m_type;

	double   m_dr;
	double   m_dz;

	double   m_eid;
	double   m_pid;
	double   m_muid;
	double   m_protonID;

	double m_massbvf;
	double m_min_e9e25;
	double m_max_e9e25;
	double m_min_energy;
	double m_max_energy;
	double m_min_theta;
	double m_max_theta;

	double m_decay_angle;

	double m_thetaT;
	double m_phiT;
	double m_thetaB;

	double m_mbc;
	double m_deltae;
};
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif /* USERINFO_H_ */
