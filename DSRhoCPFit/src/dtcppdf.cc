/**
 *  @file    dtcppdf.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-06-25
 *
 *  @brief Class that calculates the full angular time-dependent efficiency corrected
 *  CP violating PDF and provides it's analytical integrals
 *
 */


#include "dtcppdf.h"

// ROOT includes
#include "TF1.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "Math/Functor.h"
#include "Math/IntegratorMultiDim.h"

// Local includes
#include "constants.h"

//ClassImp(DtCPPDF)

DtCPPDF::DtCPPDF(const char *name, const char *title, bool _CKM_favored, bool _perfect_tagging, int _efficiency_model,
		RooAbsReal& _tht,
		RooAbsReal& _thb,
		RooAbsReal& _phit,
		RooAbsReal& _ap,
		RooAbsReal& _apa,
		RooAbsReal& _a0,
		RooAbsReal& _ata,
		RooAbsReal& _xp,
		RooAbsReal& _x0,
		RooAbsReal& _xt,
		RooAbsReal& _yp,
		RooAbsReal& _y0,
		RooAbsReal& _yt,

		RooAbsReal& _wtag,
		RooAbsReal& _dt,
		RooAbsReal& _tau,
		RooAbsReal& _dm,
		RooAbsReal& _expmc,
		RooAbsReal& _expno,
		RooAbsReal& _costhetabz,
		RooAbsReal& _benergy,
		RooAbsReal& _mbc,
		RooAbsReal& _vrntrk,
		RooAbsReal& _vrzerr,
		RooAbsReal& _vrchi2,
		RooAbsReal& _vrndf,
		RooAbsReal& _vtntrk,
		RooAbsReal& _vtzerr,
		RooAbsReal& _vtchi2,
		RooAbsReal& _vtndf,
		RooAbsReal& _vtistagl) :
			RooAbsPdf(name,title),
			tht("tht","tht",this,_tht),
			thb("thb","thb",this,_thb),
			phit("phit","phit",this,_phit),
			ap("ap","ap",this,_ap),
			apa("apa","apa",this,_apa),
			a0("a0","a0",this,_a0),
			ata("ata","ata",this,_ata),
			xp("xp","xp",this,_xp),
			x0("x0","x0",this,_x0),
			xt("xt","xt",this,_xt),
			yp("yp","yp",this,_yp),
			y0("y0","y0",this,_y0),
			yt("yt","yt",this,_yt),

			wtag("wtag","wtag",this,_wtag),
			dt("dt","dt",this,_dt),
			tau("tau","tau",this,_tau),
			dm("dm","dm",this,_dm),
			expmc("expmc","expmc",this,_expmc),
			expno("expno","expno",this,_expno),
			costhetabz("costhetabz","costhetabz",this,_costhetabz),
			benergy("benergy","benergy",this,_benergy),
			mbc("mbc","mbc",this,_mbc),
			vrntrk("vrntrk","vrntrk",this,_vrntrk),
			vrzerr("vrzerr","vrzerr",this,_vrzerr),
			vrchi2("vrchi2","vrchi2",this,_vrchi2),
			vrndf("vrndf","vrndf",this,_vrndf),
			vtntrk("vtntrk","vtntrk",this,_vtntrk),
			vtzerr("vtzerr","vtzerr",this,_vtzerr),
			vtchi2("vtchi2","vtchi2",this,_vtchi2),
			vtndf("vtndf","vtndf",this,_vtndf),
			vtistagl("vtistagl","vtistagl",this,_vtistagl),

			efficiency_model(_efficiency_model),
			mixing(true),
			CKM_favored(_CKM_favored),
			perfect_tagging(_perfect_tagging)
{
	// The rest of this constructor computes angular integrals
	// of certain terms of the PDF. This is used to speed up
	// computation of normalization; see DtCPPDF::analyticalIntegral
	ROOT::Math::Functor wf1(this, &DtCPPDF::f1, 3);
	ROOT::Math::Functor wf2(this, &DtCPPDF::f2, 3);
	ROOT::Math::Functor wf3(this, &DtCPPDF::f3, 3);
	ROOT::Math::Functor wf4(this, &DtCPPDF::f4, 3);
	ROOT::Math::Functor wf5(this, &DtCPPDF::f5, 3);
	ROOT::Math::Functor wf6(this, &DtCPPDF::f6, 3);

	double a[] = {tht.min(), thb.min(), phit.min()};
	double b[] = {tht.max(), thb.max(), phit.max()};

	printf("INFO: Beginning precomputation of efficiency integrals... \n");
	ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
	ig.SetFunction(wf1);
	int_tht_thb_phit[0] = ig.Integral(a,b);
	ig.SetFunction(wf2);
	int_tht_thb_phit[1] = ig.Integral(a,b);
	ig.SetFunction(wf3);
	int_tht_thb_phit[2] = ig.Integral(a,b);
	ig.SetFunction(wf4);
	int_tht_thb_phit[3] = ig.Integral(a,b);
	ig.SetFunction(wf5);
	int_tht_thb_phit[4] = ig.Integral(a,b);
	ig.SetFunction(wf6);
	int_tht_thb_phit[5] = ig.Integral(a,b);
	printf("INFO: Efficiency integrals' precomputation finished. \n");
}

DtCPPDF::DtCPPDF(const char *name, const char *title,
		RooAbsReal& _dt,
		RooAbsReal& _tau,
		RooAbsReal& _expmc,
		RooAbsReal& _expno,
		RooAbsReal& _costhetabz,
		RooAbsReal& _benergy,
		RooAbsReal& _mbc,
		RooAbsReal& _vrntrk,
		RooAbsReal& _vrzerr,
		RooAbsReal& _vrchi2,
		RooAbsReal& _vrndf,
		RooAbsReal& _vtntrk,
		RooAbsReal& _vtzerr,
		RooAbsReal& _vtchi2,
		RooAbsReal& _vtndf,
		RooAbsReal& _vtistagl) :
			RooAbsPdf(name,title),
			dt("dt","dt",this,_dt),
			tau("tau","tau",this,_tau),
			expmc("expmc","expmc",this,_expmc),
			expno("expno","expno",this,_expno),
			costhetabz("costhetabz","costhetabz",this,_costhetabz),
			benergy("benergy","benergy",this,_benergy),
			mbc("mbc","mbc",this,_mbc),
			vrntrk("vrntrk","vrntrk",this,_vrntrk),
			vrzerr("vrzerr","vrzerr",this,_vrzerr),
			vrchi2("vrchi2","vrchi2",this,_vrchi2),
			vrndf("vrndf","vrndf",this,_vrndf),
			vtntrk("vtntrk","vtntrk",this,_vtntrk),
			vtzerr("vtzerr","vtzerr",this,_vtzerr),
			vtchi2("vtchi2","vtchi2",this,_vtchi2),
			vtndf("vtndf","vtndf",this,_vtndf),
			vtistagl("vtistagl","vtistagl",this,_vtistagl),

			mixing(false),
			CKM_favored(false),
			perfect_tagging(false)
{
}


DtCPPDF::DtCPPDF(const DtCPPDF& other, const char* name) :
    		RooAbsPdf(other,name),
			tht("tht",this,other.tht),
			thb("thb",this,other.thb),
			phit("phit",this,other.phit),
			ap("ap",this,other.ap),
			apa("apa",this,other.apa),
			a0("a0",this,other.a0),
			ata("ata",this,other.ata),
			xp("xp",this,other.xp),
			x0("x0",this,other.x0),
			xt("xt",this,other.xt),
			yp("yp",this,other.yp),
			y0("y0",this,other.y0),
			yt("yt",this,other.yt),

			wtag("wtag",this,other.wtag),
			dt("dt",this,other.dt),
			tau("tau",this,other.tau),
			dm("dm",this,other.dm),
			expmc("expmc",this,other.expmc),
			expno("expno",this,other.expno),
			costhetabz("costhetabz",this,other.costhetabz),
			benergy("benergy",this,other.benergy),
			mbc("mbc",this,other.mbc),
			vrntrk("vrntrk",this,other.vrntrk),
			vrzerr("vrzerr",this,other.vrzerr),
			vrchi2("vrchi2",this,other.vrchi2),
			vrndf("vrndf",this,other.vrndf),
			vtntrk("vtntrk",this,other.vtntrk),
			vtzerr("vtzerr",this,other.vtzerr),
			vtchi2("vtchi2",this,other.vtchi2),
			vtndf("vtndf",this,other.vtndf),
			vtistagl("vtistagl",this,other.vtistagl),

			efficiency_model(other.efficiency_model),
			mixing(other.mixing),
			CKM_favored(other.CKM_favored),
			perfect_tagging(other.perfect_tagging)
{
	for (int i = 0; i < 6; i++) {
        int_tht_thb_phit[i] = other.int_tht_thb_phit[i];
    }
}


Double_t DtCPPDF::evaluate() const {
	const bool mc = (expmc == 2) ? 1 : 0;
	const Belle::dtres_param_t* dtres_param = Belle::get_dtres_param((int)expno, mc);
	double ak, ck;
	const double pb_cms_sq = (benergy * benergy) - (mbc * mbc);
	const double Eb_cms = std::sqrt((Belle::dt_resol_global::mbzero * Belle::dt_resol_global::mbzero) + pb_cms_sq);
	Belle::CalcAkCk(costhetabz, Eb_cms, &ak, &ck, Belle::dt_resol_global::mbzero);

	double pdf = 0;
	double pdf_const = 0;
	double pdf_sin = 0;
	double pdf_cos = 0;

//	double norm = 0;
//	double nnorm = 0;
//	double norm_const = 0;
//	double norm_sin = 0;
//	double norm_cos = 0;

	if (mixing) {
		// AfRkRdetRnp_fullrec and MfRkRdetRnp_fullrec are supposedly (Sumisawa BAS 2010)
		// not normalized, hence the 0.5/tau factors
		pdf_const = EfRkRdetRnp_fullrec( dt, constants::btype,
				tau, ak, ck,
				vrntrk, vrzerr,	vrchi2, vrndf,
				vtntrk, vtzerr,	vtchi2, vtndf,
				vtistagl, dtres_param );
		pdf_sin = AfRkRdetRnp_fullrec( dt, constants::btype,
				tau, dm, ak, ck,
				vrntrk, vrzerr,	vrchi2, vrndf,
				vtntrk, vtzerr,	vtchi2, vtndf,
				vtistagl, dtres_param ) * 0.5/tau;
		pdf_cos = MfRkRdetRnp_fullrec( dt, constants::btype,
				tau, dm, ak, ck,
				vrntrk, vrzerr,	vrchi2, vrndf,
				vtntrk, vtzerr,	vtchi2, vtndf,
				vtistagl, dtres_param ) * 0.5/tau;

//		norm_const = norm_EfRkRdetRnp_fullrec(constants::cut_dt_low, constants::cut_dt_high, constants::btype,
//				tau, ak, ck,
//				vrntrk, vrzerr,	vrchi2, vrndf,
//				vtntrk, vtzerr,	vtchi2, vtndf,
//				vtistagl, dtres_param );
//		norm_sin = norm_AfRkRdetRnp_fullrec(constants::cut_dt_low, constants::cut_dt_high, constants::btype,
//				tau, dm, ak, ck,
//				vrntrk, vrzerr,	vrchi2, vrndf,
//				vtntrk, vtzerr,	vtchi2, vtndf,
//				vtistagl, dtres_param ) * 0.5/tau;
//		norm_cos = norm_MfRkRdetRnp_fullrec(constants::cut_dt_low, constants::cut_dt_high, constants::btype,
//				tau, dm, ak, ck,
//				vrntrk, vrzerr,	vrchi2, vrndf,
//				vtntrk, vtzerr,	vtchi2, vtndf,
//				vtistagl, dtres_param ) * 0.5/tau;

//		double alpha = 1;

//		pdf_const = Belle::AddOutlierWithBkg((int) expno, dt, 1, pdf_const, pdf_const, (int) vrntrk, (int) vtntrk, dtres_param,
//				norm_const / alpha, norm_const / alpha, constants::cut_dt_low, constants::cut_dt_high, alpha, 1);
//
//		pdf_sin = Belle::AddOutlierWithBkg((int) expno, dt, 1, pdf_sin, pdf_sin, (int) vrntrk, (int) vtntrk, dtres_param,
//				norm_sin / alpha, norm_sin / alpha, constants::cut_dt_low, constants::cut_dt_high, alpha, 1);
//
//		pdf_cos = Belle::AddOutlierWithBkg((int) expno, dt, 1, pdf_cos, pdf_cos, (int) vrntrk, (int) vtntrk, dtres_param,
//				norm_cos / alpha, norm_cos / alpha, constants::cut_dt_low, constants::cut_dt_high, alpha, 1);

//		double r = 1 - 2 * wtag;
//		int r_bin = GetRBin(r);
//		double wtag_binned = GetWTag(expno, r_bin, mc);
//		double delta_wtag_binned = GetDeltaWTag(expno, r_bin, mc);
//		double r_binned = 1 - 2 * wtag_binned;
		double sign = -1 + 2*CKM_favored;
//
//		if (perfect_tagging) {
//			delta_wtag_binned = 0;
//			r_binned = 1;
//		}

//		pdf = pdf_const * (1 + sign * delta_wtag_binned) - sign * r_binned * (S * pdf_sin + A * pdf_cos);
//		norm = norm_const * (1 + sign * delta_wtag_binned) - sign * r_binned * (S * norm_sin + A * norm_cos);


////////////////////////////////////////////////////////////////

		Double_t a0a = 0;
		Double_t at = sqrt(1-ap*ap-a0*a0);

		Double_t ap0r = ap*a0*cos(-apa+a0a);
		Double_t ap0i = ap*a0*sin(-apa+a0a);
		Double_t a0tr = a0*at*cos(-a0a+ata);
		Double_t a0ti = a0*at*sin(-a0a+ata);
		Double_t aptr = ap*at*cos(-apa+ata);
		Double_t apti = ap*at*sin(-apa+ata);

		Double_t At2 = 0;
		Double_t Ap2 = 0;
		Double_t A02 = 0;
		Double_t Ap0r = 0;
		Double_t A0ti = 0;
		Double_t Apti = 0;

		Ap2 = ap*ap*((1 + xp*xp + yp*yp)*pdf_const + (1 - xp*xp - yp*yp)*sign*pdf_cos + 2*yp*sign*pdf_sin);
		A02 = a0*a0*((1 + x0*x0 + y0*y0)*pdf_const + (1 - x0*x0 - y0*y0)*sign*pdf_cos + 2*y0*sign*pdf_sin);
		At2 = at*at*((1 + xt*xt + yt*yt)*pdf_const + (1 - xt*xt - yt*yt)*sign*pdf_cos + 2*yt*sign*pdf_sin);

		Ap0r = ap0r*((1 + xp*x0 + yp*y0)*pdf_const + (1 - xp*x0 - yp*y0)*sign*pdf_cos + (yp + y0)*sign*pdf_sin) -\
			   ap0i*((x0*yp - xp*y0)*(pdf_const - sign*pdf_cos) + (x0 - xp)*sign*pdf_sin);

		A0ti = a0ti*((1 + x0*xt + y0*yt)*pdf_const + (1 - x0*xt - y0*yt)*sign*pdf_cos + (y0 + yt)*sign*pdf_sin) +\
			   a0tr*((xt*y0 - x0*yt)*(pdf_const - sign*pdf_cos) + (xt - x0)*sign*pdf_sin);

		Apti = apti*((1 + xp*xt + yp*yt)*pdf_const + (1 - xp*xt - yp*yt)*sign*pdf_cos + (yp + yt)*sign*pdf_sin) +\
			   aptr*((xt*yp - xp*yt)*(pdf_const - sign*pdf_cos) + (xt - xp)*sign*pdf_sin);


		Double_t value =    (Ap2*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
							At2*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
							A02*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
							sqrt(2)*Ap0r*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
							sqrt(2)*A0ti*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
							2*Apti*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit));

		pdf = value * eff.GetEfficiency(tht, thb, phit, efficiency_model);

//		printf("DBG: value = %f, pdf = %f\n", value, pdf);

//		double A = ap*ap*(1 - xp*xp - yp*yp) + a0*a0*(1 - x0*x0 - y0*y0) + at*at*(1 - xt*xt - yt*yt);
//		double S = ap*ap*2*yp + a0*a0*2*y0 + at*at*2*yt;
//
//		printf("ap  = %f\n"
//				"apa = %f\n"
//				"a0  = %f\n"
//				"a0a = %f\n"
//				"at  = %f\n"
//				"ata = %f\n"
//				"xp  = %f\n"
//				"x0  = %f\n"
//				"xt  = %f\n"
//				"yp  = %f\n"
//				"y0  = %f\n"
//				"yt  = %f\n",
//				(double) ap,
//				(double) apa,
//				(double) a0,
//				(double) a0a,
//				(double) at,
//				(double) ata,
//				(double) xp,
//				(double) x0,
//				(double) xt,
//				(double) yp,
//				(double) y0,
//				(double) yt);
//
//		printf("\nA = %f\nS = %f\n\n", A, S);

//		double nAp2 = ap*ap*((1 + xp*xp + yp*yp)*norm_const + (1 - xp*xp - yp*yp)*sign*norm_cos + 2*yp*sign*norm_sin);
//		double nA02 = a0*a0*((1 + x0*x0 + y0*y0)*norm_const + (1 - x0*x0 - y0*y0)*sign*norm_cos + 2*y0*sign*norm_sin);
//		double nAt2 = at*at*((1 + xt*xt + yt*yt)*norm_const + (1 - xt*xt - yt*yt)*sign*norm_cos + 2*yt*sign*norm_sin);

//		norm = 64.0 * constants::kPi/9 * (Ap2 + A02 + At2);
//		nnorm = 64.0 * constants::kPi/9 * (nAp2 + nA02 + nAt2);

	} else {
		pdf = EfRkRdetRnp_fullrec( dt, constants::btype,
				tau, ak, ck,
				vrntrk, vrzerr,	vrchi2, vrndf,
				vtntrk, vtzerr,	vtchi2, vtndf,
				vtistagl, dtres_param );

//		norm = norm_EfRkRdetRnp_fullrec(constants::cut_dt_low, constants::cut_dt_high, constants::btype,
//					tau, ak, ck,
//					vrntrk, vrzerr,	vrchi2, vrndf,
//					vtntrk, vtzerr,	vtchi2, vtndf,
//					vtistagl, dtres_param );
	}


	// This is extremely dumb, as I calculate the normalization, normalize
	// the PDF (inside AddOutlier) and then de-normalize it, only to calculate
	// the normalization later in analyticalIntegral again. Don't have time to
	// fix the dumbness now.
	// TODO: Fix dumbness
	double alpha = 1;
	double nnorm = analyticalIntegral(12);
	pdf = Belle::AddOutlier(expno, dt, pdf, vrntrk, vtntrk, dtres_param,
			nnorm, constants::cuts::dt_low, constants::cuts::dt_high, alpha);
	return pdf*nnorm;

//	return Belle::AddOutlierWithBkg((int) expno, dt, 1, pdf, pdf, (int) vrntrk, (int) vtntrk, dtres_param,
//			nnorm / alpha, nnorm / alpha, constants::cut_dt_low, constants::cut_dt_high, alpha, 1);

	return pdf;

}

Int_t DtCPPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {

//	printf(" ******************************* DBG: Looking for analytical integral...\n");

	if(matchArgs(allVars,analVars,dt,tht,thb,phit)) return 1;

	if(matchArgs(allVars,analVars,thb,tht,phit)) return 2;

	// TODO: Before enabling the following analytical integrals, a proper
	// procedure to integrate out the angles (taking efficiency into account)
	// has to be implemented

//	if(matchArgs(allVars,analVars,dt,thb,tht)) return 3;
//	if(matchArgs(allVars,analVars,dt,tht,phit)) return 4;
//	if(matchArgs(allVars,analVars,dt,thb,phit)) return 5;
//
//	if(matchArgs(allVars,analVars,dt,tht)) return 6;
//	if(matchArgs(allVars,analVars,dt,thb)) return 7;
//	if(matchArgs(allVars,analVars,dt,phit)) return 8;
//
//	if(matchArgs(allVars,analVars,thb,tht)) return 9;
//	if(matchArgs(allVars,analVars,tht,phit)) return 10;
//	if(matchArgs(allVars,analVars,thb,phit)) return 11;

	if(matchArgs(allVars,analVars,dt)) return 12;

//	printf(" ******************************* DBG: Analytical integral not found!\n");

	return 0 ;
}

Double_t DtCPPDF::analyticalIntegral(Int_t code, const char* rangeName) const {

	const bool mc = (expmc == 2) ? 1 : 0;
	const Belle::dtres_param_t* dtres_param = Belle::get_dtres_param((int)expno, mc);
	double ak, ck;
	const double pb_cms_sq = (benergy * benergy) - (mbc * mbc);
	const double Eb_cms = std::sqrt((Belle::dt_resol_global::mbzero * Belle::dt_resol_global::mbzero) + pb_cms_sq);
	Belle::CalcAkCk(costhetabz, Eb_cms, &ak, &ck, Belle::dt_resol_global::mbzero);

	double pdf_const = 0;
	double pdf_sin = 0;
	double pdf_cos = 0;

	double norm_const = 0;
	double norm_sin = 0;
	double norm_cos = 0;

	if (mixing) {
		// There is no need to calculate both norm_{const,sin,cos} and pdf_{const,sin,cos}
		if (IsTimeIntegrated(code)) {
			norm_const = norm_EfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
					tau, ak, ck,
					vrntrk, vrzerr,	vrchi2, vrndf,
					vtntrk, vtzerr,	vtchi2, vtndf,
					vtistagl, dtres_param );
			// AfRkRdetRnp_fullrec and MfRkRdetRnp_fullrec are supposedly (Sumisawa BAS 2010)
			// not normalized, hence the 0.5/tau factors
			norm_sin = norm_AfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
					tau, dm, ak, ck,
					vrntrk, vrzerr,	vrchi2, vrndf,
					vtntrk, vtzerr,	vtchi2, vtndf,
					vtistagl, dtres_param ) * 0.5/tau;
			norm_cos = norm_MfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
					tau, dm, ak, ck,
					vrntrk, vrzerr,	vrchi2, vrndf,
					vtntrk, vtzerr,	vtchi2, vtndf,
					vtistagl, dtres_param ) * 0.5/tau;
		} else {
			pdf_const = EfRkRdetRnp_fullrec( dt, constants::btype,
					tau, ak, ck,
					vrntrk, vrzerr,	vrchi2, vrndf,
					vtntrk, vtzerr,	vtchi2, vtndf,
					vtistagl, dtres_param );
			pdf_sin = AfRkRdetRnp_fullrec( dt, constants::btype,
					tau, dm, ak, ck,
					vrntrk, vrzerr,	vrchi2, vrndf,
					vtntrk, vtzerr,	vtchi2, vtndf,
					vtistagl, dtres_param ) * 0.5/tau;
			pdf_cos = MfRkRdetRnp_fullrec( dt, constants::btype,
					tau, dm, ak, ck,
					vrntrk, vrzerr,	vrchi2, vrndf,
					vtntrk, vtzerr,	vtchi2, vtndf,
					vtistagl, dtres_param ) * 0.5/tau;
		}

//		double r = 1 - 2 * wtag;
//		int r_bin = GetRBin(r);
//		double wtag_binned = GetWTag(expno, r_bin, mc);
//		double delta_wtag_binned = GetDeltaWTag(expno, r_bin, mc);
//		double r_binned = 1 - 2 * wtag_binned;
		double sign = -1 + 2*CKM_favored;
//
//		if (perfect_tagging) {
//			delta_wtag_binned = 0;
//			r_binned = 1;
//		}

//		pdf = pdf_const * (1 + sign * delta_wtag_binned) - sign * r_binned * (S * pdf_sin + A * pdf_cos);
//		norm = norm_const * (1 + sign * delta_wtag_binned) - sign * r_binned * (S * norm_sin + A * norm_cos);


////////////////////////////////////////////////////////////////

		Double_t a0a = 0;
		Double_t at = sqrt(1-ap*ap-a0*a0);

		Double_t ap0r = ap*a0*cos(-apa+a0a);
		Double_t ap0i = ap*a0*sin(-apa+a0a);
		Double_t a0tr = a0*at*cos(-a0a+ata);
		Double_t a0ti = a0*at*sin(-a0a+ata);
		Double_t aptr = ap*at*cos(-apa+ata);
		Double_t apti = ap*at*sin(-apa+ata);

		Double_t At2 = 0;
		Double_t Ap2 = 0;
		Double_t A02 = 0;
		Double_t Ap0r = 0;
		Double_t A0ti = 0;
		Double_t Apti = 0;

		Ap2 = ap*ap*((1 + xp*xp + yp*yp)*pdf_const + (1 - xp*xp - yp*yp)*sign*pdf_cos + 2*yp*sign*pdf_sin);
		A02 = a0*a0*((1 + x0*x0 + y0*y0)*pdf_const + (1 - x0*x0 - y0*y0)*sign*pdf_cos + 2*y0*sign*pdf_sin);
		At2 = at*at*((1 + xt*xt + yt*yt)*pdf_const + (1 - xt*xt - yt*yt)*sign*pdf_cos + 2*yt*sign*pdf_sin);

		Ap0r = ap0r*((1 + xp*x0 + yp*y0)*pdf_const + (1 - xp*x0 - yp*y0)*sign*pdf_cos + (yp + y0)*sign*pdf_sin) -\
			   ap0i*((x0*yp - xp*y0)*(pdf_const - sign*pdf_cos) + (x0 - xp)*sign*pdf_sin);

		A0ti = a0ti*((1 + x0*xt + y0*yt)*pdf_const + (1 - x0*xt - y0*yt)*sign*pdf_cos + (y0 + yt)*sign*pdf_sin) +\
			   a0tr*((xt*y0 - x0*yt)*(pdf_const - sign*pdf_cos) + (xt - x0)*sign*pdf_sin);

		Apti = apti*((1 + xp*xt + yp*yt)*pdf_const + (1 - xp*xt - yp*yt)*sign*pdf_cos + (yp + yt)*sign*pdf_sin) +\
			   aptr*((xt*yp - xp*yt)*(pdf_const - sign*pdf_cos) + (xt - xp)*sign*pdf_sin);


		Double_t nAt2 = 0;
		Double_t nAp2 = 0;
		Double_t nA02 = 0;
		Double_t nAp0r = 0;
		Double_t nA0ti = 0;
		Double_t nApti = 0;

		nAp2 = ap*ap*((1 + xp*xp + yp*yp)*norm_const + (1 - xp*xp - yp*yp)*sign*norm_cos + 2*yp*sign*norm_sin);
		nA02 = a0*a0*((1 + x0*x0 + y0*y0)*norm_const + (1 - x0*x0 - y0*y0)*sign*norm_cos + 2*y0*sign*norm_sin);
		nAt2 = at*at*((1 + xt*xt + yt*yt)*norm_const + (1 - xt*xt - yt*yt)*sign*norm_cos + 2*yt*sign*norm_sin);

		nAp0r = ap0r*((1 + xp*x0 + yp*y0)*norm_const + (1 - xp*x0 - yp*y0)*sign*norm_cos + (yp + y0)*sign*norm_sin) -\
			   ap0i*((x0*yp - xp*y0)*(norm_const - sign*norm_cos) + (x0 - xp)*sign*norm_sin);

		nA0ti = a0ti*((1 + x0*xt + y0*yt)*norm_const + (1 - x0*xt - y0*yt)*sign*norm_cos + (y0 + yt)*sign*norm_sin) +\
			   a0tr*((xt*y0 - x0*yt)*(norm_const - sign*norm_cos) + (xt - x0)*sign*norm_sin);

		nApti = apti*((1 + xp*xt + yp*yt)*norm_const + (1 - xp*xt - yp*yt)*sign*norm_cos + (yp + yt)*sign*norm_sin) +\
			   aptr*((xt*yp - xp*yt)*(norm_const - sign*norm_cos) + (xt - xp)*sign*norm_sin);

		switch(code) {
		case 1: // Int[g,{dt,tht,thb,phit}]
			return  nAp2 * int_tht_thb_phit[0] +
					nAt2 * int_tht_thb_phit[1] +
					nA02 * int_tht_thb_phit[2] +
					nAp0r * int_tht_thb_phit[3] -
					nA0ti * int_tht_thb_phit[4] -
					nApti * int_tht_thb_phit[5];

		case 2: // Int[g,{tht,thb,phit}]
			return  Ap2 * int_tht_thb_phit[0] +
					At2 * int_tht_thb_phit[1] +
					A02 * int_tht_thb_phit[2] +
					Ap0r * int_tht_thb_phit[3] -
					A0ti * int_tht_thb_phit[4] -
					Apti * int_tht_thb_phit[5];

		// TODO: Coeficients are not right - efficiency not taken into account
		case 3: // Int[g,{dt,tht,thb}]
			return 16./9.*(nAt2 + 2*nA02*cos(phit)*cos(phit) + 2*nAp2*sin(phit)*sin(phit));

		case 4: // Int[g,{dt,tht,phit}]
			return 8.*constants::pi/3*((nAp2+nAt2)*sin(thb)*sin(thb) + 2*nA02*cos(thb)*cos(thb))*sin(thb);

		case 5: // Int[g,{dt,thb,phit}]
			return 8.*constants::pi/3*((nAp2+nA02)*sin(tht)*sin(tht) + 2*nAt2*cos(tht)*cos(tht))*sin(tht);

		case 6: // Int[g,{dt,tht}]
			return 4./3.*(4*nA02*cos(phit)*cos(phit)*cos(thb)*cos(thb) + \
					(nAt2 + 2*nAp2*sin(phit)*sin(phit))*sin(thb)*sin(thb) + sqrt(2)*ap0r*sin(2*phit)*sin(2*thb))*sin(thb);

		case 7: // Int[g,{dt,thb}]
			return 8./3.*(nAt2*cos(tht)*cos(tht) + nA02*cos(phit)*cos(phit)*sin(tht)*sin(tht) + \
					sin(phit)*(nAp2*sin(phit)*sin(tht)*sin(tht) - apti*sin(2*tht)))*sin(tht);

		case 8: // Int[g,{dt,phit}]
			return 2*constants::pi*(2*nAt2*cos(tht)*cos(tht)*sin(thb)*sin(thb) +
					(2*nA02*cos(thb)*cos(thb) + nAp2*sin(thb)*sin(thb))*sin(tht)*sin(tht))*(sin(tht)*sin(thb));

		case 9: // Int[g,{tht,thb}]
			return 16./9.*(At2 + 2*A02*cos(phit)*cos(phit) + 2*Ap2*sin(phit)*sin(phit));

		case 10: // Int[g,{tht,phit}]
			return 8.*constants::pi/3*((Ap2+At2)*sin(thb)*sin(thb) + 2*A02*cos(thb)*cos(thb))*sin(thb);

		case 11: // Int[g,{thb,phit}]
			return 8.*constants::pi/3*((Ap2+A02)*sin(tht)*sin(tht) + 2*At2*cos(tht)*cos(tht))*sin(tht);

		case 12: // Int[g,{dt}]
			return (nAp2*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
							nAt2*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
							nA02*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
							sqrt(2)*nAp0r*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
							sqrt(2)*nA0ti*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
							2*nApti*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)) * eff.GetEfficiency(tht, thb, phit, efficiency_model);

		default:
			return 0;
		}

	} else {
		return norm_EfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
					tau, ak, ck,
					vrntrk, vrzerr,	vrchi2, vrndf,
					vtntrk, vtzerr,	vtchi2, vtndf,
					vtistagl, dtres_param );
	}


//	return Belle::AddOutlier(expno, dt, pdf, vrntrk, vtntrk, dtres_param,
//			norm, constants::cut_dt_low, constants::cut_dt_high, alpha);

//	Belle::AddOutlierWithBkg((int) expno, dt, 1, pdf, pdf, (int) vrntrk, (int) vtntrk, dtres_param,
//			norm / alpha, norm / alpha, constants::cut_dt_low, constants::cut_dt_high, alpha, 1);



}

bool DtCPPDF::IsTimeIntegrated(int code) const {
    if (code == 1 || (code >= 3 && code <= 8) || code == 12) {
    	return true;
    } else {
    	return false;
    }
}

int DtCPPDF::GetRBin(double r) const {
    return (0. <= r && r <= 0.1 ? 0 :
    0.1 < r && r <= 0.25 ? 1 :
    0.25 < r && r <= 0.5 ? 2 :
    0.5 < r && r <= 0.625 ? 3 :
    0.625 < r && r <= 0.75 ? 4 :
    0.75 < r && r <= 0.875 ? 5 :
    0.875 < r && r <= 1.0 ? 6 : 7);
}

double DtCPPDF::GetWTag(int expno, int rbin, bool mc) const {
    double w_svd1_data[7] = {0.5, 0.418852, 0.329879, 0.233898, 0.170608, 0.099791, 0.0228501};

    double w_svd2_data[7] = {0.5, 0.418826, 0.319303, 0.222948, 0.163191, 0.104085, 0.0251454};

    double w_svd1_mc[7] = {0.5, 0.420827, 0.300296, 0.219317, 0.154636, 0.0916131, 0.0228891};

    double w_svd2_mc[7] = {0.5, 0.412222, 0.307838, 0.212765, 0.149933, 0.0913264, 0.0218754};

    if (mc) {
        if (expno < 30) {
            return w_svd1_mc[rbin];
        } else  {
            return w_svd2_mc[rbin];
        }
    } else {
        if (expno < 30) {
            return w_svd1_data[rbin];
        } else {
            return w_svd2_data[rbin];
        }
    }
}

double DtCPPDF::GetDeltaWTag(int expno, int rbin, bool mc) const {
    double dw_svd1_data[7] = {0., 0.0569661, 0.0126192, -0.0147724, -0.000550289, 0.00887704, 0.00465683};

    double dw_svd2_data[7] = {0., -0.00877001, 0.0103515, -0.0109253, -0.0186365, 0.00168037, -0.0036441};

    double dw_svd1_mc[7] = {0., 0.0583019, 0.00573998, -0.0392635, 0.00474508, -0.0118737, -0.00585326};

    double dw_svd2_mc[7] = {0., 0.00408778, 0.010326, -0.00479522, 0.00151989, 0.0143633, 0.00189979};
    if (mc) {
        if (expno < 30) {
            return dw_svd1_mc[rbin];
        } else {
            return dw_svd2_mc[rbin];
        }
    } else {
        if (expno < 30) {
            return dw_svd1_data[rbin];
        } else {
            return dw_svd2_data[rbin];
        }
    }
}

Double_t DtCPPDF::f1(const double * vars) {
	Double_t val = 2 * sin(vars[0]) * sin(vars[0]) * sin(vars[0]) * sin(vars[1]) * sin(vars[1]) * sin(vars[1]) * sin(vars[2]) * sin(vars[2]);
	return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t DtCPPDF::f2(const double * vars) {
	Double_t val = 2 * cos(vars[0]) * cos(vars[0]) * sin(vars[0]) * sin(vars[1]) * sin(vars[1]) * sin(vars[1]);
	return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t DtCPPDF::f3(const double * vars) {
	Double_t val = 4 * sin(vars[0]) * sin(vars[0]) * sin(vars[0]) * cos(vars[1]) * cos(vars[1]) * sin(vars[1]) * cos(vars[2]) * cos(vars[2]);
	return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t DtCPPDF::f4(const double * vars) {
	Double_t val = sqrt(2) * sin(vars[0]) * sin(vars[0]) * sin(vars[0]) * sin(2 * vars[1]) * sin(vars[1]) * sin(2 * vars[2]);
	return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t DtCPPDF::f5(const double * vars) {
	Double_t val = sqrt(2) * sin(2 * vars[0]) * sin(vars[0]) * sin(2 * vars[1]) * sin(vars[1]) * cos(vars[2]);
	return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}

Double_t DtCPPDF::f6(const double * vars) {
	Double_t val = 2 * sin(2 * vars[0]) * sin(vars[0]) * sin(vars[1]) * sin(vars[1]) * sin(vars[1]) * sin(vars[2]);
	return val * eff.GetEfficiency(vars[0], vars[1], vars[2], efficiency_model);
}
