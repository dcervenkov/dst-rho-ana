#include "UserInfo.h"

#include "belle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


// For Interface to Set UserInfo Class

void setUserInfo(Particle &p, unsigned ch){
	if(! &p.userInfo() )// return;
		p.userInfo(UserInfo(ch));
}

void setUserInfo(std::vector<Particle> &p, unsigned ch) {
	for(unsigned i=0;i<p.size();++i) setUserInfo(p[i],ch);
}

// UserInfo Class

UserInfo::UserInfo() :
	m_decayMode(0),
	m_type(0),
	m_numCandidates(0),
	m_candidateSelection(-999),

	m_dr(-999.),
	m_dz(-999.),
	m_eid(-1.),
	m_pid(-1.),
	m_muid(-1.),
	m_protonID(-1.),

	m_massbvf(0.0),
	m_chisq(0.0),
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

	m_thetaTG(-999.),
	m_phiTG(-999.),
	m_thetaBG(-999.),

	m_mbc(-999.),
	m_deltae(-999.),

	m_cosThetaB(-999.),
	m_cosThetaTB(-999.),
	m_cosThetaT(-999.),

	m_contSuppressionBDTG(-3.),
	m_contSuppressionMLP(-3.),

	m_vertex_rec(),
	m_vertex_tag(),
	m_flavor(0),
	m_qr(0),
	m_wtag(0) {
}

UserInfo::UserInfo(unsigned ch) :
	m_decayMode(0),
	m_type(0),
	m_numCandidates(0),
	m_candidateSelection(-999),

	m_dr(-1.),
	m_dz(-1.),
	m_eid(-1.),
	m_pid(-1.),
	m_muid(-1.),
	m_protonID(-1.),

	m_massbvf(0.0),
	m_chisq(0.0),
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

	m_thetaTG(-999.),
	m_phiTG(-999.),
	m_thetaBG(-999.),

	m_mbc(-999.),
	m_deltae(-999.),

	m_cosThetaB(-999.),
	m_cosThetaTB(-999.),
	m_cosThetaT(-999.),

	m_contSuppressionBDTG(-3.),
	m_contSuppressionMLP(-3.),

	m_vertex_rec(),
	m_vertex_tag(),
	m_flavor(0),
	m_qr(0),
	m_wtag(0) {
}

UserInfo::UserInfo(const UserInfo &x) :
	m_decayMode(x.m_decayMode),
	m_type(x.m_type),
	m_numCandidates(x.m_numCandidates),
	m_candidateSelection(x.m_candidateSelection),

	m_dr(x.m_dr),
	m_dz(x.m_dz),
	m_eid(x.m_eid),
	m_pid(x.m_pid),
	m_muid(x.m_muid),
	m_protonID(x.m_protonID),

	m_massbvf(x.m_massbvf),
	m_chisq(x.m_chisq),
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

	m_thetaTG(x.m_thetaTG),
	m_phiTG(x.m_phiTG),
	m_thetaBG(x.m_thetaBG),

	m_mbc(x.m_mbc),
	m_deltae(x.m_deltae),

	m_ksfw(x.m_ksfw),

	m_cosThetaB(x.m_cosThetaB),
	m_cosThetaTB(x.m_cosThetaTB),
	m_cosThetaT(x.m_cosThetaT),

	m_contSuppressionBDTG(x.m_contSuppressionBDTG),
	m_contSuppressionMLP(x.m_contSuppressionMLP),

	m_vertex_rec(x.m_vertex_rec),
	m_vertex_tag(x.m_vertex_tag),
	m_flavor(x.m_flavor),
	m_qr(x.m_qr),
	m_wtag(x.m_wtag) {
}

UserInfo::~UserInfo() {
}

UserInfo* UserInfo::clone(void) const {
	UserInfo *x = new UserInfo(*this);
	return x;
}

UserInfo& UserInfo::operator = (const UserInfo &x) {
	m_decayMode = x.m_decayMode;
	m_type   = x.m_type;
	m_numCandidates = x.m_numCandidates;
	m_candidateSelection = x.m_candidateSelection;

	m_dr     = x.m_dr;
	m_dz     = x.m_dz;
	m_eid = x.m_eid;
	m_pid = x.m_pid;
	m_muid= x.m_muid;
	m_protonID = x.m_protonID;
	
	m_massbvf   = x.m_massbvf;
	m_chisq   = x.m_chisq;
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

	m_thetaTG = x.m_thetaTG;
	m_phiTG = x.m_phiTG;
	m_thetaBG = x.m_thetaBG;

	m_mbc         = x.m_mbc;
	m_deltae      = x.m_deltae;

	m_ksfw = x.m_ksfw;

	m_cosThetaB = x.m_cosThetaB;
	m_cosThetaTB = x.m_cosThetaTB;
	m_cosThetaT = x.m_cosThetaT;

	m_contSuppressionBDTG = x.m_contSuppressionBDTG;
	m_contSuppressionMLP = x.m_contSuppressionMLP;

	m_vertex_rec = x.m_vertex_rec;
	m_vertex_tag = x.m_vertex_tag;

	m_flavor = x.m_flavor;
	m_qr = x.m_qr;
	m_wtag = x.m_wtag;

	return *this;
}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
