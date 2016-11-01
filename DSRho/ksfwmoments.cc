
	      //////////////////////////////////////////
	      //                                      //
	      //  ksfwmoments.cc                      //
	      //                                      //
	      //  moment-calculation of the k_sfw     //
	      //  improved Super-Fox-Wolfram moments  //
	      //  to be used with rooksfw.{cc,h}      //
	      //                                      //
	      //  M. Nakao                            //
	      //                                      //
	      //////////////////////////////////////////

// Versions
// 2006111301  first version


#include <stdio.h>
#include <unistd.h>
#include <ip/IpProfile.h>
#include <particle/Particle.h>
#include <kfitter/kvertexfitter.h>
#include <mdst/findKs.h>
#include MDST_H

#include "ksfwmoments.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// ----------------------------------------------------------------------
// legendre
// ----------------------------------------------------------------------
inline double
legendre(const double z, const int i)
{
  switch (i) {
  case 0:  return 1.;
  case 1:  return z;
  case 2:  return 1.5*z*z - 0.5;
  case 3:  return z*(2.5*z*z - 1.5);
  case 4:  return (4.375*z*z*z*z - 3.75*z*z + 0.375);
  //(wrong) case 4:  return (4.735*z*z*z*z - 3.75*z*z + 0.375);
  default: return 0;
  }
}
// ----------------------------------------------------------------------
// constructor
// (copied from k_sfw.cc with minimum modification)
// ----------------------------------------------------------------------
ksfwmoments::ksfwmoments(const Particle &cand, double ebeam,
			 const Hep3Vector &cmboost)
{
  if (! IpProfile::usable()) {
    printf("======================================================\n");
    printf("rooksfw: IpProfile::begin_run() has not been called!!!\n");
    printf("rooksfw: IpProfile is not usuable!!!\n");
    printf("======================================================\n");
    exit(1);
  }
  
  m_uf = 0;

  // kinematically allowed maximum momentum for mbc>5.2
  // sqrt(5.29^2 - 5.2^2) ~ 1.0GeV
  // 10.58/4 + 1.0/2 = 3.145
  static const double P_MAX(3.2);
  static const double PIMASS(0.139570);

  //=====================
  // create particle list
  //=====================
  HepLorentzVector pb(cand.p());
  pb.boost(cmboost);
  const double Hso0_max(2*(2*ebeam-pb.e())); // 2 * (ebeam - delta_E)

  std::vector<int> Q_sigA, Q_sigB, Q_other;
  HepLorentzVector p_cms_missA(0, 0, 0, 2*ebeam);
  HepLorentzVector p_cms_missB(0, 0, 0, 2*ebeam);
  m_et[0] = m_et[1] = 0;

  Mdst_vee2_Manager& veeMgr(Mdst_vee2_Manager::get_manager());
  Mdst_charged_Manager& chMgr(Mdst_charged_Manager::get_manager());
  Mdst_gamma_Manager& gMgr(Mdst_gamma_Manager::get_manager());
  const std::vector<Particle *> & final 
    = cand.relation().finalStateParticles();

  // for signal (use_finalstate_for_sig == 0)
  for (unsigned i = 0; i < cand.nChildren(); i++){
    const Particle & dau(cand.child(i));
    HepLorentzVector p_cms(dau.p());
    p_cms.boost(cmboost);
    p_cms_sigA.push_back(p_cms.vect());
    Q_sigA.push_back((int)dau.charge());
    p_cms_missA -= p_cms;
    m_et[0] += p_cms.perp();
  }

  // ks
  FindKs fKs;
  std::vector<Mdst_charged *> ks_dau;
  for (std::vector<Mdst_vee2>::iterator i = veeMgr.begin();
       i != veeMgr.end(); i++){
    bool sig(false);
    for (std::vector<Particle *>::const_iterator j = final.begin();
	 j != final.end(); j++){
      if ((**j).mdstCharged() && 
	  ((**j).mdstCharged().get_ID() == (*i).chgd(0).get_ID() || 
	   (**j).mdstCharged().get_ID() == (*i).chgd(1).get_ID() )){
	sig = true;
	break;
      }
    }
    if (sig) continue;
    if ( (*i).kind() != 1 ) continue;
    fKs.candidates( *i, IpProfile::position() );
    if( !fKs.goodKs() ) continue;
    HepLorentzVector p_cms((*i).px(), (*i).py(), (*i).pz(), (*i).energy());
    p_cms.boost(cmboost);
    if (p_cms.rho() > P_MAX) continue;
    p_cms_other.push_back(p_cms.vect());
    Q_other.push_back(0);
    ks_dau.push_back(&(*i).chgd(0));
    ks_dau.push_back(&(*i).chgd(1));
    p_cms_missA -= p_cms;
    p_cms_missB -= p_cms;
    m_et[0] += p_cms.perp();
    m_et[1] += p_cms.perp();
  }

  // charged tracks
  HepVector a(5);
  HepSymMatrix Ea(5,0);
  for (std::vector<Mdst_charged>::const_iterator i = chMgr.begin();
       i != chMgr.end(); i++){
    bool ksdau(false);
    for (std::vector<Mdst_charged *>::const_iterator j = ks_dau.begin();
	 j != ks_dau.end(); j++){
      if ((**j).get_ID() == (*i).get_ID()){
	ksdau = true;
	break;
      }
    }
    if (ksdau) continue;
    bool sig(false);
    HepLorentzVector p_cms;
    for (std::vector<Particle *>::const_iterator j = final.begin();
	 j != final.end(); j++){
      if ((**j).mdstCharged()
	  && (**j).mdstCharged().get_ID() == (*i).get_ID()) {
	p_cms = (**j).p();
	sig = true;
      }
    }
    if (sig){
      // for signal (use_finalstate_for_sig == 1)
      p_cms.boost(cmboost);
      p_cms_sigB.push_back(p_cms.vect());
      Q_sigB.push_back((int)(*i).charge());
      p_cms_missB -= p_cms;
      m_et[1] += p_cms.perp();
    } else {
      HepPoint3D p((*i).trk().mhyp(2).pivot_x(),
		   (*i).trk().mhyp(2).pivot_y(),
		   (*i).trk().mhyp(2).pivot_z());
      a[0] = (*i).trk().mhyp(2).helix(0);
      a[1] = (*i).trk().mhyp(2).helix(1);
      a[2] = (*i).trk().mhyp(2).helix(2);
      a[3] = (*i).trk().mhyp(2).helix(3);
      a[4] = (*i).trk().mhyp(2).helix(4);
      Ea[0][0] = (*i).trk().mhyp(2).error(0);
      Ea[1][0] = (*i).trk().mhyp(2).error(1);
      Ea[1][1] = (*i).trk().mhyp(2).error(2);
      Ea[2][0] = (*i).trk().mhyp(2).error(3);
      Ea[2][1] = (*i).trk().mhyp(2).error(4);
      Ea[2][2] = (*i).trk().mhyp(2).error(5);
      Ea[3][0] = (*i).trk().mhyp(2).error(6);
      Ea[3][1] = (*i).trk().mhyp(2).error(7);
      Ea[3][2] = (*i).trk().mhyp(2).error(8);
      Ea[3][3] = (*i).trk().mhyp(2).error(9);
      Ea[4][0] = (*i).trk().mhyp(2).error(10);
      Ea[4][1] = (*i).trk().mhyp(2).error(11);
      Ea[4][2] = (*i).trk().mhyp(2).error(12);
      Ea[4][3] = (*i).trk().mhyp(2).error(13);
      Ea[4][4] = (*i).trk().mhyp(2).error(14);
      Helix hel(p, a, Ea);
      hel.pivot(IpProfile::position());
      if (std::fabs(hel.dr()) > 5.0 || std::fabs(hel.dz()) > 10.0) continue;

      HepSymMatrix Err(7, 0);
      HepPoint3D Pos;
      HepLorentzVector Mom(hel.momentum(0, PIMASS, Pos, Err));

      kvertexfitter kvf;
      kvf.initialVertex(IpProfile::position());
      kvf.beamProfile(IpProfile::position_err_b_life_smeared());
      kvf.addTrack(Mom, Pos, Err, (*i).charge(), PIMASS);
      
      p_cms = kvf.fit() ? Mom : kvf.momentum(0);
      //p_cms = Mom;
      p_cms.boost(cmboost);
      if (p_cms.rho() > P_MAX) continue;
      p_cms_other.push_back(p_cms.vect());
      Q_other.push_back((int)(*i).charge());
      p_cms_missA -= p_cms;
      p_cms_missB -= p_cms;
      m_et[0] += p_cms.perp();
      m_et[1] += p_cms.perp();
    }
  }

  //gammas
  for (std::vector<Mdst_gamma>::const_iterator i = gMgr.begin();
       i != gMgr.end(); i++){
    bool sig(false);
    HepLorentzVector p_cms;
    for (std::vector<Particle *>::const_iterator j = final.begin();
	 j != final.end(); j++){
      if ((**j).mdstGamma() && (**j).mdstGamma().get_ID() == (*i).get_ID()) {
	p_cms = (**j).p();
	sig = true;
      }
    }
    if (sig){
      // for signal (use_finalstate_for_sig == 1)
      p_cms.boost(cmboost);
      p_cms_sigB.push_back(p_cms.vect());
      Q_sigB.push_back(0);
      p_cms_missB -= p_cms;
      m_et[1] += p_cms.perp();
    } else {
      p_cms.setPx((*i).px());
      p_cms.setPy((*i).py());
      p_cms.setPz((*i).pz());
      p_cms.setE(p_cms.rho());
      if (p_cms.rho() < 0.05) continue;
      //if ((*i).ecl().match() || p_cms.rho() < 0.1) continue;
      p_cms.boost(cmboost);
      if (p_cms.rho() > P_MAX) continue;
      p_cms_other.push_back(p_cms.vect());
      Q_other.push_back(0);
      p_cms_missA -= p_cms;
      p_cms_missB -= p_cms;
      m_et[0] += p_cms.perp();
      m_et[1] += p_cms.perp();
    }
  }

  m_mm2[0] = p_cms_missA.e() > 0 
    ? p_cms_missA.mag2() 
    : -p_cms_missA.e()*p_cms_missA.e() - p_cms_missA.vect().mag2();
  m_mm2[1] = p_cms_missB.e() > 0 
    ? p_cms_missB.mag2() 
    : -p_cms_missB.e()*p_cms_missB.e() - p_cms_missB.vect().mag2();

  //========================
  // calculate discriminants
  //========================
  std::vector<Hep3Vector>::iterator pi, pj;
  std::vector<int>::iterator Qi, Qj;

  // calculate Hso components
  for (int i = 0; i < 3; i++){
    for (int k = 0; k < 5; k++){
      m_Hso[0][i][k] = m_Hso[1][i][k] = 0;
    }
  }

  // signal A (use_finalstate_for_sig == 0)
  for (pi = p_cms_sigA.begin(), Qi = Q_sigA.begin();
       pi != p_cms_sigA.end(); pi++, Qi++){
    const double pi_mag((*pi).mag());
    for (pj = p_cms_other.begin(), Qj = Q_other.begin();
	 pj != p_cms_other.end(); pj++, Qj++){
      const double pj_mag((*pj).mag());
      const double ij_cos((*pi)*(*pj)/pi_mag/pj_mag);
      const int c_or_n(0 == (*Qj) ? 1 : 0);  // 0: charged 1: neutral
      for (int k = 0; k < 5; k++){
	m_Hso[0][c_or_n][k] += ( k % 2 ) 
	  ? (*Qi)*(*Qj)*pj_mag*legendre(ij_cos, k)
	  : pj_mag*legendre(ij_cos, k);
      }
    }
    const double p_miss_mag(p_cms_missA.rho());
    const double i_miss_cos((*pi)*p_cms_missA.vect()/pi_mag/p_miss_mag);
    for (int k = 0; k < 5; k++){
      m_Hso[0][2][k] += ( k % 2 ) ? 0 : p_miss_mag*legendre(i_miss_cos, k);
    }
  }

  // signal B (use_finalstate_for_sig == 1)
  for (pi = p_cms_sigB.begin(), Qi = Q_sigB.begin();
       pi != p_cms_sigB.end(); pi++, Qi++){
    const double pi_mag((*pi).mag());
    for (pj = p_cms_other.begin(), Qj = Q_other.begin();
	 pj != p_cms_other.end(); pj++, Qj++){
      const double pj_mag((*pj).mag());
      const double ij_cos((*pi)*(*pj)/pi_mag/pj_mag);
      const int c_or_n(0 == (*Qj) ? 1 : 0);  // 0: charged 1: neutral
      for (int k = 0; k < 5; k++){
	m_Hso[1][c_or_n][k] += ( k % 2 ) 
	  ? (*Qi)*(*Qj)*pj_mag*legendre(ij_cos, k)
	  : pj_mag*legendre(ij_cos, k);
      }
    }
    const double p_miss_mag(p_cms_missB.rho());
    const double i_miss_cos((*pi)*p_cms_missB.vect()/pi_mag/p_miss_mag);
    for (int k = 0; k < 5; k++){
      m_Hso[1][2][k] += ( k % 2 ) ? 0 : p_miss_mag*legendre(i_miss_cos, k);
    }
  }

  // add missing to the lists
  std::vector<Hep3Vector> p_cms_otherA(p_cms_other), p_cms_otherB(p_cms_other);
  p_cms_otherA.push_back(p_cms_missA.vect());
  p_cms_otherB.push_back(p_cms_missB.vect());
  Q_other.push_back(0);

  // calculate Hoo components
  for (int k = 0; k < 5; k++){
    m_Hoo[0][k] = m_Hoo[1][k] = 0;
  }
  for (pi = p_cms_otherA.begin(), Qi = Q_other.begin();
       pi != p_cms_otherA.end(); pi++, Qi++){
    const double pi_mag((*pi).mag());
    for (pj = p_cms_otherA.begin(), Qj = Q_other.begin();
	 pj != pi; pj++, Qj++){
      const double pj_mag((*pj).mag());
      const double ij_cos((*pi)*(*pj)/pi_mag/pj_mag);
      for (int k = 0; k < 5; k++){
	m_Hoo[0][k] += ( k % 2 ) 
	  ? (*Qi)*(*Qj)*pi_mag*pj_mag*legendre(ij_cos, k)
	  : pi_mag*pj_mag*legendre(ij_cos, k);
      }
    }
  }
  for (pi = p_cms_otherB.begin(), Qi = Q_other.begin();
       pi != p_cms_otherB.end(); pi++, Qi++){
    const double pi_mag((*pi).mag());
    for (pj = p_cms_otherB.begin(), Qj = Q_other.begin();
	 pj != pi; pj++, Qj++){
      const double pj_mag((*pj).mag());
      const double ij_cos((*pi)*(*pj)/pi_mag/pj_mag);
      for (int k = 0; k < 5; k++){
	m_Hoo[1][k] += ( k % 2 ) 
	  ? (*Qi)*(*Qj)*pi_mag*pj_mag*legendre(ij_cos, k)
	  : pi_mag*pj_mag*legendre(ij_cos, k);
      }
    }
  }

  // nomalize so that it does not dependent on delta_e
  for (int k = 0; k < 5; k++){
    for (int j = 0; j < ((k%2) ? 1 : 3); j++){
      m_Hso[0][j][k] /= Hso0_max;
      m_Hso[1][j][k] /= Hso0_max;
    }
    m_Hoo[0][k] /= (Hso0_max*Hso0_max);
    m_Hoo[1][k] /= (Hso0_max*Hso0_max);
  }

  //printf("rooksfw::rooksfw: mm2=%f et=%f hoo2=%f hso02=%f\n",
  //	 m_mm2[0], m_et[0], m_Hoo[0][2], m_Hso[0][0][2]);
}
#if defined(BELLE_NAMESPACE)
}
#endif
