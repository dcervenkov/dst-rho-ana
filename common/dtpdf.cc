/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "dtpdf.h"

// ROOT includes
#include "TF1.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

// Local includes
#include "constants.h"
#include "tools.h"

//ClassImp(DSRhoPDF)

DtPDF::DtPDF(const char *name, const char *title, bool _CKM_favored, bool _perfect_tagging,
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
        CKM_favored(_CKM_favored),
        perfect_tagging(_perfect_tagging)
{
}

DtPDF::DtPDF(const char *name, const char *title,
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


DtPDF::DtPDF(const DtPDF& other, const char* name) :
            RooAbsPdf(other,name),
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
            CKM_favored(other.CKM_favored),
            perfect_tagging(other.perfect_tagging)
{
}


Double_t DtPDF::evaluate() const {
    const bool mc = (expmc == 2) ? 1 : 0;
    const Belle::dtres_param_t* dtres_param = Belle::get_dtres_param((int)expno, mc);
    double ak, ck;
    const double pb_cms_sq = (benergy * benergy) - (mbc * mbc);
    const double Eb_cms = std::sqrt((Belle::dt_resol_global::mbzero * Belle::dt_resol_global::mbzero) + pb_cms_sq);
    Belle::CalcAkCk(costhetabz, Eb_cms, &ak, &ck, Belle::dt_resol_global::mbzero);

    double pdf = 0;
    double pdf_const = 0;
    double pdf_cos = 0;

    double norm = 0;
    double norm_const = 0;
    double norm_cos = 0;

    if (mixing) {
        // AfRkRdetRnp_fullrec and MfRkRdetRnp_fullrec are supposedly (Sumisawa BAS 2010)
        // not normalized, hence the 0.5/tau factors
        pdf_const = EfRkRdetRnp_fullrec( dt, constants::btype,
                tau, ak, ck,
                vrntrk, sqrt(vrzerr),	vrchi2, vrndf,
                vtntrk, sqrt(vtzerr),	vtchi2, vtndf,
                vtistagl, dtres_param );
        pdf_cos = MfRkRdetRnp_fullrec( dt, constants::btype,
                tau, dm, ak, ck,
                vrntrk, sqrt(vrzerr),	vrchi2, vrndf,
                vtntrk, sqrt(vtzerr),	vtchi2, vtndf,
                vtistagl, dtres_param ) * 0.5/tau;

        norm_const = norm_EfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
                tau, ak, ck,
                vrntrk, sqrt(vrzerr),	vrchi2, vrndf,
                vtntrk, sqrt(vtzerr),	vtchi2, vtndf,
                vtistagl, dtres_param );
        norm_cos = norm_MfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
                tau, dm, ak, ck,
                vrntrk, sqrt(vrzerr),	vrchi2, vrndf,
                vtntrk, sqrt(vtzerr),	vtchi2, vtndf,
                vtistagl, dtres_param ) * 0.5/tau;

        double r = 1 - 2 * wtag;
        int r_bin = tools::GetRBin(r);
        double wtag_binned = tools::GetWTag(expno, r_bin, mc);
        double r_binned = 1 - 2 * wtag_binned;
        double sign = -1 + 2*CKM_favored;

        if (perfect_tagging) {
            r_binned = 1;
        }

        pdf  = pdf_const  + sign * r_binned * pdf_cos;
        norm = norm_const + sign * r_binned * norm_cos;

    } else {
        pdf = EfRkRdetRnp_fullrec( dt, constants::btype,
                tau, ak, ck,
                vrntrk, sqrt(vrzerr),	vrchi2, vrndf,
                vtntrk, sqrt(vtzerr),	vtchi2, vtndf,
                vtistagl, dtres_param );

        norm = norm_EfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
                    tau, ak, ck,
                    vrntrk, sqrt(vrzerr),	vrchi2, vrndf,
                    vtntrk, sqrt(vtzerr),	vtchi2, vtndf,
                    vtistagl, dtres_param );
    }

    double alpha = 1;

    return Belle::AddOutlier(expno, dt, pdf, vrntrk, vtntrk, dtres_param, norm,
                             constants::cuts::dt_low, constants::cuts::dt_high, alpha);

    // return Belle::AddOutlierWithBkg((int) expno, dt, 1, pdf, pdf, (int) vrntrk, (int) vtntrk, dtres_param,
    //         norm / alpha, norm / alpha, constants::cuts::dt_low, constants::cuts::dt_high, alpha, 1);

}

Int_t DtPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {
    if(matchArgs(allVars,analVars,dt)) return 1;

    return 0 ;
}

Double_t DtPDF::analyticalIntegral(Int_t code, const char* rangeName) const {
    switch(code) {
    case 1: // Int[g,{dt}]
        return 1;
    default:
        return 0;
    }

}
