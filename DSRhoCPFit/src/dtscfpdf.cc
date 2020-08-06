/**
 *  @file    dtscfpdf.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2017-12-25
 *
 *  @brief Class that holds the full angular time-dependent self-cross-feed PDF
 *
 */

#include "dtscfpdf.h"

// ROOT includes
#include "TF1.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "Math/Functor.h"
#include "Math/IntegratorMultiDim.h"

// Local includes
#include "constants.h"
#include "tools.h"

//ClassImp(DtSCFPDF)

DtSCFPDF::DtSCFPDF(const char *name, const char *title, bool _B_bar, bool _CKM_favored, bool _perfect_tagging,
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

            mixing(true),
            B_bar(_B_bar),
            CKM_favored(_CKM_favored),
            perfect_tagging(_perfect_tagging)
{
}

DtSCFPDF::DtSCFPDF(const char *name, const char *title,
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
            B_bar(false),
            CKM_favored(false),
            perfect_tagging(false)
{
}


DtSCFPDF::DtSCFPDF(const DtSCFPDF& other, const char* name) :
            RooAbsPdf(other,name),
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

            mixing(other.mixing),
            B_bar(other.B_bar),
            CKM_favored(other.CKM_favored),
            perfect_tagging(other.perfect_tagging)
{
}


Double_t DtSCFPDF::evaluate() const {
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

    if (mixing) {
        // AfRkRdetRnp_fullrec and MfRkRdetRnp_fullrec are supposedly (Sumisawa BAS 2010)
        // not normalized, hence the 0.5/tau factors
        pdf_const = EfRkRdetRnp_fullrec( dt, constants::btype,
                tau, ak, ck,
                vrntrk, vrzerr, vrchi2, vrndf,
                vtntrk, vtzerr, vtchi2, vtndf,
                vtistagl, dtres_param );
        pdf_sin = AfRkRdetRnp_fullrec( dt, constants::btype,
                tau, dm, ak, ck,
                vrntrk, vrzerr, vrchi2, vrndf,
                vtntrk, vtzerr, vtchi2, vtndf,
                vtistagl, dtres_param ) * 0.5/tau;
        pdf_cos = MfRkRdetRnp_fullrec( dt, constants::btype,
                tau, dm, ak, ck,
                vrntrk, vrzerr, vrchi2, vrndf,
                vtntrk, vtzerr, vtchi2, vtndf,
                vtistagl, dtres_param ) * 0.5/tau;

        Double_t At2 = 0;
        Double_t Ap2 = 0;
        Double_t A02 = 0;

        CalculateAmplitudeTerms(Ap2, A02, At2, pdf_const, pdf_cos, pdf_sin);

        pdf = (Ap2 + A02 + At2);

    } else {
        pdf = EfRkRdetRnp_fullrec( dt, constants::btype,
                tau, ak, ck,
                vrntrk, vrzerr, vrchi2, vrndf,
                vtntrk, vtzerr, vtchi2, vtndf,
                vtistagl, dtres_param );
    }

    return pdf;
}

Int_t DtSCFPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {

    if(matchArgs(allVars,analVars,dt)) return 1;

    return 0 ;
}

Double_t DtSCFPDF::analyticalIntegral(Int_t code, const char* rangeName) const {

    const bool mc = (expmc == 2) ? 1 : 0;
    const Belle::dtres_param_t* dtres_param = Belle::get_dtres_param((int)expno, mc);
    double ak, ck;
    const double pb_cms_sq = (benergy * benergy) - (mbc * mbc);
    const double Eb_cms = std::sqrt((Belle::dt_resol_global::mbzero * Belle::dt_resol_global::mbzero) + pb_cms_sq);
    Belle::CalcAkCk(costhetabz, Eb_cms, &ak, &ck, Belle::dt_resol_global::mbzero);

    double norm_const = 0;
    double norm_sin = 0;
    double norm_cos = 0;

    if (mixing) {
        norm_const = norm_EfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
                tau, ak, ck,
                vrntrk, vrzerr, vrchi2, vrndf,
                vtntrk, vtzerr, vtchi2, vtndf,
                vtistagl, dtres_param );
        // AfRkRdetRnp_fullrec and MfRkRdetRnp_fullrec are supposedly (Sumisawa BAS 2010)
        // not normalized, hence the 0.5/tau factors
        norm_sin = norm_AfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
                tau, dm, ak, ck,
                vrntrk, vrzerr, vrchi2, vrndf,
                vtntrk, vtzerr, vtchi2, vtndf,
                vtistagl, dtres_param ) * 0.5/tau;
        norm_cos = norm_MfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
                tau, dm, ak, ck,
                vrntrk, vrzerr, vrchi2, vrndf,
                vtntrk, vtzerr, vtchi2, vtndf,
                vtistagl, dtres_param ) * 0.5/tau;

        Double_t nAt2 = 0;
        Double_t nAp2 = 0;
        Double_t nA02 = 0;

        switch(code) {
        case 1: // Int[g,{dt}]
            CalculateAmplitudeTerms(nAp2, nA02, nAt2, norm_const, norm_cos, norm_sin);
            return nAp2 + nA02 + nAt2;

        default:
            return 0;
        }

    } else {
        return norm_EfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
                    tau, ak, ck,
                    vrntrk, vrzerr, vrchi2, vrndf,
                    vtntrk, vtzerr, vtchi2, vtndf,
                    vtistagl, dtres_param );
    }


//  return Belle::AddOutlier(expno, dt, pdf, vrntrk, vtntrk, dtres_param,
//          norm, constants::cut_dt_low, constants::cut_dt_high, alpha);

//  Belle::AddOutlierWithBkg((int) expno, dt, 1, pdf, pdf, (int) vrntrk, (int) vtntrk, dtres_param,
//          norm / alpha, norm / alpha, constants::cut_dt_low, constants::cut_dt_high, alpha, 1);



}

bool DtSCFPDF::IsTimeIntegrated(int code) const {
    if (code == 1 || (code >= 3 && code <= 8) || code == 12) {
        return true;
    } else {
        return false;
    }
}

void DtSCFPDF::CalculateAmplitudeTerms(double& Ap2, double& A02, double& At2,
                                      const double& constant, const double& cosine,
                                      const double& sine) const {

        const bool mc = (expmc == 2) ? 1 : 0;
		double r = 1 - 2 * wtag;
		int r_bin = tools::GetRBin(r);
		double wtag_binned = tools::GetWTag(expno, r_bin, mc);
		double delta_wtag_binned = tools::GetDeltaWTag(expno, r_bin, mc);
		double r_binned = 1 - 2 * wtag_binned;

		if (perfect_tagging) {
			delta_wtag_binned = 0;
			r_binned = 1;
        }

        Double_t at = sqrt(1-ap*ap-a0*a0);

        double sign = -1 + 2*CKM_favored;

		Ap2 = ap*ap*((1 + xp*xp + yp*yp)*constant*(1 - delta_wtag_binned) + ((1 - xp*xp - yp*yp)*sign*cosine + 2*yp*sign*sine)*r_binned);
		A02 = a0*a0*((1 + x0*x0 + y0*y0)*constant*(1 - delta_wtag_binned) + ((1 - x0*x0 - y0*y0)*sign*cosine + 2*y0*sign*sine)*r_binned);
		At2 = at*at*((1 + xt*xt + yt*yt)*constant*(1 - delta_wtag_binned) + ((1 - xt*xt - yt*yt)*sign*cosine + 2*yt*sign*sine)*r_binned);

}

