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

//ClassImp(DtSCFPDF)

DtSCFPDF::DtSCFPDF(const char *name, const char *title, bool _B_bar, bool _CKM_favored, bool _perfect_tagging,
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

            mixing(true),
            B_bar(_B_bar),
            CKM_favored(_CKM_favored),
            perfect_tagging(_perfect_tagging)
{
    // Self-cross-feed phit model
    scf_phit_poly_p2_ = new RooRealVar("scf_phit_poly_p2", "p_(2)", 0.856, -0.1, 2);
    scf_phit_f_ = new RooRealVar("scf_phit_f", "f_(poly)", 0.147, 0.1, 0.9);
    scf_phit_poly_ = new RooPolynomial("scf_phit_poly", "scf_phit_poly", _phit, *scf_phit_poly_p2_, 2);
    scf_phit_offset_ = new RooRealVar("scf_phit_offset", "#phi_(t)^(offset)", 0.056, -0.1, 0.1);
    scf_phit_phit_ = new RooFormulaVar("scf_phit_phit", "scf_phit_phit", "phit - scf_phit_offset",
                                 RooArgList(_phit, *scf_phit_offset_));
    scf_phit_cos_ = new RooGenericPdf("scf_phit_cos", "scf_phit_cos", "cos(scf_phit_phit)^2",
                                RooArgList(*scf_phit_phit_));

    // Self-cross-feed thetat model
    scf_thetat_f_ = new RooRealVar("scf_thetat_f", "#theta_(t)^(w)", -0.051, -0.1, 0.1);
    scf_thetat_thetat_ = new RooFormulaVar("scf_thetat_thetat", "scf_thetat_thetat",
                                     "(thetat - 1.5708)*(1+scf_thetat_f) + 1.5708",
                                     RooArgList(_tht, *scf_thetat_f_));

    // Self-cross-feed thetab model
    scf_thetab_gaus_mu_ = new RooRealVar("scf_thetab_gaus_mu", "#mu", 2.885, 1.5, 3);
    scf_thetab_gaus_sigma_l_ = new RooRealVar("scf_thetab_gaus_sigma_l", "#sigma_(L)", 0.411, 0, 3);
    scf_thetab_gaus_sigma_r_ = new RooRealVar("scf_thetab_gaus_sigma_r", "#sigma_(R)", 0.094, 0, 3);
    scf_thetab_gaus_ = new RooBifurGauss(
        "scf_thetab_gaus",  "scf_thetab_gaus",       _thb,
        *scf_thetab_gaus_mu_, *scf_thetab_gaus_sigma_l_, *scf_thetab_gaus_sigma_r_);
    scf_thetab_exp_alpha_ = new RooRealVar("scf_thetab_exp_alpha", "#alpha", -4.63, -10, 0.0);
    scf_thetab_exp_ = new RooExponential("scf_thetab_exp", "scf_thetab_exp", _thb,
                                   *scf_thetab_exp_alpha_);
    scf_thetab_f_ = new RooRealVar("scf_thetab_f", "f_(exp)", 0.625, 0, 1);

    scf_phit_model_ = new RooAddPdf("scf_phit_model", "scf_phit_model",
                              RooArgList(*scf_phit_poly_, *scf_phit_cos_), RooArgList(*scf_phit_f_));

    scf_thetat_model_ = new RooGenericPdf("scf_thetat_model", "scf_thetat_model",
                                    "sin(scf_thetat_thetat)^3", RooArgList(*scf_thetat_thetat_));

    scf_thetab_model_ = new RooAddPdf("scf_thetab_model", "scf_thetab_model",
                              RooArgList(*scf_thetab_exp_, *scf_thetab_gaus_), RooArgList(*scf_thetab_f_));
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
            
            // Self-cross-feed phit model
            scf_phit_poly_p2_(other.scf_phit_poly_p2_),
            scf_phit_f_(other.scf_phit_f_),
            scf_phit_poly_(other.scf_phit_poly_),
            scf_phit_offset_(other.scf_phit_offset_),
            scf_phit_phit_(other.scf_phit_phit_),
            scf_phit_cos_(other.scf_phit_cos_),
            scf_phit_model_(other.scf_phit_model_),

            // Self-cross-feed thetat model
            scf_thetat_f_(other.scf_thetat_f_),
            scf_thetat_thetat_(other.scf_thetat_thetat_),
            scf_thetat_model_(other.scf_thetat_model_),

            // Self-cross-feed thetab model
            scf_thetab_gaus_mu_(other.scf_thetab_gaus_mu_),
            scf_thetab_gaus_sigma_l_(other.scf_thetab_gaus_sigma_l_),
            scf_thetab_gaus_sigma_r_(other.scf_thetab_gaus_sigma_r_),
            scf_thetab_gaus_(other.scf_thetab_gaus_),
            scf_thetab_exp_alpha_(other.scf_thetab_exp_alpha_),
            scf_thetab_exp_(other.scf_thetab_exp_),
            scf_thetab_f_(other.scf_thetab_f_),
            scf_thetab_model_(other.scf_thetab_model_),

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

        pdf = (Ap2 + A02 + At2) * scf_thetat_model_->getVal() * scf_thetab_model_->getVal() * scf_phit_model_->getVal();

        double dtht = tht;
        RooRealVar rootht("thetat", "thetat", dtht);
        printf("tht = %f, tht_model = %f\n", dtht, scf_thetat_model_->getVal(RooArgSet(rootht)));

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

//  printf(" ******************************* DBG: Looking for analytical integral...\n");

    if(matchArgs(allVars,analVars,dt,tht,thb,phit)) return 1;

    if(matchArgs(allVars,analVars,thb,tht,phit)) return 2;

    // TODO: Before enabling the following analytical integrals, a proper
    // procedure to integrate out the angles (taking efficiency into account)
    // has to be implemented

    if(matchArgs(allVars,analVars,dt,thb,tht)) return 3;
    if(matchArgs(allVars,analVars,dt,tht,phit)) return 4;
    if(matchArgs(allVars,analVars,dt,thb,phit)) return 5;
//
//  if(matchArgs(allVars,analVars,dt,tht)) return 6;
//  if(matchArgs(allVars,analVars,dt,thb)) return 7;
//  if(matchArgs(allVars,analVars,dt,phit)) return 8;
//
//  if(matchArgs(allVars,analVars,thb,tht)) return 9;
//  if(matchArgs(allVars,analVars,tht,phit)) return 10;
//  if(matchArgs(allVars,analVars,thb,phit)) return 11;

    if(matchArgs(allVars,analVars,dt)) return 12;

//  printf(" ******************************* DBG: Analytical integral not found!\n");

    return 0 ;
}

Double_t DtSCFPDF::analyticalIntegral(Int_t code, const char* rangeName) const {

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
        } else {
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
        }

        Double_t At2 = 0;
        Double_t Ap2 = 0;
        Double_t A02 = 0;

        Double_t nAt2 = 0;
        Double_t nAp2 = 0;
        Double_t nA02 = 0;

        switch(code) {
        case 1: // Int[g,{dt,tht,thb,phit}]
            CalculateAmplitudeTerms(nAp2, nA02, nAt2, norm_const, norm_cos, norm_sin);
            return nAp2 + nA02 + nAt2;

        case 2: // Int[g,{tht,thb,phit}]
            CalculateAmplitudeTerms(Ap2, A02, At2, pdf_const, pdf_cos, pdf_sin);
            return Ap2 + A02 + At2;

        case 3: // Int[g,{dt,tht,thb}]
            CalculateAmplitudeTerms(nAp2, nA02, nAt2, norm_const, norm_cos, norm_sin);
            return (nAp2 + nA02 + nAt2) * scf_phit_model_->getVal();

        case 4: // Int[g,{dt,tht,phit}]
            CalculateAmplitudeTerms(nAp2, nA02, nAt2, norm_const, norm_cos, norm_sin);
            return (nAp2 + nA02 + nAt2) * scf_thetab_model_->getVal();

        case 5: // Int[g,{dt,thb,phit}]
            CalculateAmplitudeTerms(nAp2, nA02, nAt2, norm_const, norm_cos, norm_sin);
            return (nAp2 + nA02 + nAt2) * scf_thetat_model_->getVal();

        case 6: // Int[g,{dt,tht}]

        case 7: // Int[g,{dt,thb}]

        case 8: // Int[g,{dt,phit}]

        case 9: // Int[g,{tht,thb}]

        case 10: // Int[g,{tht,phit}]

        case 11: // Int[g,{thb,phit}]

        case 12: // Int[g,{dt}]
            CalculateAmplitudeTerms(nAp2, nA02, nAt2, norm_const, norm_cos, norm_sin);
            return (nAp2 + nA02 + nAt2) * scf_thetat_model_->getVal() * scf_thetab_model_->getVal() * scf_phit_model_->getVal();

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

int DtSCFPDF::GetRBin(double r) const {
    return (0. <= r && r <= 0.1 ? 0 :
    0.1 < r && r <= 0.25 ? 1 :
    0.25 < r && r <= 0.5 ? 2 :
    0.5 < r && r <= 0.625 ? 3 :
    0.625 < r && r <= 0.75 ? 4 :
    0.75 < r && r <= 0.875 ? 5 :
    0.875 < r && r <= 1.0 ? 6 : 7);
}

double DtSCFPDF::GetWTag(int expno, int rbin, bool mc) const {
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

double DtSCFPDF::GetDeltaWTag(int expno, int rbin, bool mc) const {
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

void DtSCFPDF::CalculateAmplitudeTerms(double& Ap2, double& A02, double& At2,
                                      const double& constant, const double& cosine, 
                                      const double& sine) const {

        const bool mc = (expmc == 2) ? 1 : 0;
		double r = 1 - 2 * wtag;
		int r_bin = GetRBin(r);
		double wtag_binned = GetWTag(expno, r_bin, mc);
		double delta_wtag_binned = GetDeltaWTag(expno, r_bin, mc);
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
    
