
	      //////////////////////////////////////////
	      //                                      //
	      //  ksfwmoments.h                       //
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

#ifndef KSFWMOMENTS_H_
#define KSFWMOMENTS_H_

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class Particle;
class Hep3Vector;

// ----------------------------------------------------------------------
// ksfwmoments
// ----------------------------------------------------------------------
class ksfwmoments {
public:
  ksfwmoments() : m_uf(0) {};
  ksfwmoments(const Particle &cand, double ebeam, const Hep3Vector &cmboost);
  ~ksfwmoments() {};

  int usefinal(int uf) { if (uf==0||uf==1) m_uf = uf; return m_uf; }

  int usefinal() const { return m_uf; }
  double mm2(int uf=-1) const { return m_mm2[uf<0?m_uf:uf]; }
  double et(int uf=-1) const { return m_et[uf<0?m_uf:uf]; }
  double Hoo(int i, int uf=-1) const {
    return m_Hoo[uf<0?m_uf:uf][i]; }
  double Hso(int i, int j, int uf=-1) const {
    return m_Hso[uf<0?m_uf:uf][i][j]; }
  double var(int i) { // the same order as the original k_sfw.cc
    if (i < 0)  return 0;
    if (i < 1)  return et();
    if (i < 4)  return Hso(i-1, 0);
    if (i < 5)  return Hso(0,   1);
    if (i < 8)  return Hso(i-5, 2);
    if (i < 9)  return Hso(0,   3);
    if (i < 12) return Hso(i-9, 4);
    if (i < 17) return Hoo(i-12);
    return 0; }
  std::vector<Hep3Vector> get_p_cms_sigA() const { return p_cms_sigA; }
  std::vector<Hep3Vector> get_p_cms_sigB() const { return p_cms_sigB; }
  std::vector<Hep3Vector> get_p_cms_other() const { return p_cms_other; }

private:
  int m_uf;
  double m_Hso[2][3][5];
  double m_Hoo[2][5];
  double m_et[2];
  double m_mm2[2];
  std::vector<Hep3Vector> p_cms_sigA, p_cms_sigB, p_cms_other;
};

#if defined(BELLE_NAMESPACE)
}
#endif

#endif /* KSFWMOMENTS_H_ */
