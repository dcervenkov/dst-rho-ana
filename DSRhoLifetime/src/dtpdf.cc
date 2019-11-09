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

//ClassImp(DSRhoPDF)

DtPDF::DtPDF(const char *name, const char *title, bool _CKM_favored, bool _perfect_tagging,
        RooAbsReal& _S,
        RooAbsReal& _A,
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
        S("S","S",this,_S),
        A("A","A",this,_A),
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
            S("S",this,other.S),
            A("A",this,other.A),
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
    double pdf_sin = 0;
    double pdf_cos = 0;

    double norm = 0;
    double norm_const = 0;
    double norm_sin = 0;
    double norm_cos = 0;

    if (mixing) {
        // AfRkRdetRnp_fullrec and MfRkRdetRnp_fullrec are supposedly (Sumisawa BAS 2010)
        // not normalized, hence the 0.5/tau factors
        pdf_const = EfRkRdetRnp_fullrec( dt, constants::btype,
                tau, ak, ck,
                vrntrk, sqrt(vrzerr),	vrchi2, vrndf,
                vtntrk, sqrt(vtzerr),	vtchi2, vtndf,
                vtistagl, dtres_param );
        pdf_sin = AfRkRdetRnp_fullrec( dt, constants::btype,
                tau, dm, ak, ck,
                vrntrk, sqrt(vrzerr),	vrchi2, vrndf,
                vtntrk, sqrt(vtzerr),	vtchi2, vtndf,
                vtistagl, dtres_param ) * 0.5/tau;
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
        norm_sin = norm_AfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
                tau, dm, ak, ck,
                vrntrk, sqrt(vrzerr),	vrchi2, vrndf,
                vtntrk, sqrt(vtzerr),	vtchi2, vtndf,
                vtistagl, dtres_param ) * 0.5/tau;
        norm_cos = norm_MfRkRdetRnp_fullrec(constants::cuts::dt_low, constants::cuts::dt_high, constants::btype,
                tau, dm, ak, ck,
                vrntrk, sqrt(vrzerr),	vrchi2, vrndf,
                vtntrk, sqrt(vtzerr),	vtchi2, vtndf,
                vtistagl, dtres_param ) * 0.5/tau;

        double r = 1 - 2 * wtag;
        int r_bin = GetRBin(r);
        double wtag_binned = GetWTag(expno, r_bin, mc);
        double delta_wtag_binned = GetDeltaWTag(expno, r_bin, mc);
        double r_binned = 1 - 2 * wtag_binned;
        double sign = -1 + 2*CKM_favored;

        if (perfect_tagging) {
            delta_wtag_binned = 0;
            r_binned = 1;
        }

        pdf = pdf_const * (1 + sign * delta_wtag_binned) - sign * r_binned * (S * pdf_sin + A * pdf_cos);
        norm = norm_const * (1 + sign * delta_wtag_binned) - sign * r_binned * (S * norm_sin + A * norm_cos);

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

	// return Belle::AddOutlier(expno, dt, pdf, vrntrk, vtntrk, dtres_param,
	// 		norm, constants::cuts::dt_low, constants::cuts::dt_high, alpha);

    return Belle::AddOutlierWithBkg((int) expno, dt, 1, pdf, pdf, (int) vrntrk, (int) vtntrk, dtres_param,
            norm / alpha, norm / alpha, constants::cuts::dt_low, constants::cuts::dt_high, alpha, 1);

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

int DtPDF::GetRBin(double r) const {
    return (0. <= r && r <= 0.1 ? 0 :
    0.1 < r && r <= 0.25 ? 1 :
    0.25 < r && r <= 0.5 ? 2 :
    0.5 < r && r <= 0.625 ? 3 :
    0.625 < r && r <= 0.75 ? 4 :
    0.75 < r && r <= 0.875 ? 5 :
    0.875 < r && r <= 1.0 ? 6 : 7);
}

double DtPDF::GetWTag(int expno, int rbin, bool mc) const {
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

double DtPDF::GetDeltaWTag(int expno, int rbin, bool mc) const {
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


