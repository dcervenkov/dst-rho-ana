#include "RooTatamiHelper.h"
#include "tatami/tatami.h"
#include <string>
#include <regex.h>
#include <vector>
#include <utility>
#include <boost/foreach.hpp>

RooDataSet ascii2DataSet(const std::string& asciifilename,
		bool mc,
		RooRealVar& dt,
		RooRealVar& mbc,
		RooRealVar& deltaE,
		RooCategory& expno,
		RooRealVar& CosThetaB,
		RooRealVar& ntrk_sig,
		RooRealVar& ntrk_tag,
		RooRealVar& z_err_sig,
		RooRealVar& z_err_tag,
		RooRealVar& chisq_tracks_sig,
		RooRealVar& chisq_tracks_tag,
		RooRealVar& dgf_tracks_sig,
		RooRealVar& dgf_tracks_tag,
		RooRealVar& keeptagl_tag,
		RooRealVar& ecms,
		RooCategory& rbin,
		RooRealVar& wtag,
		RooRealVar& delta_wtag,
		RooCategory& flavor_tag,
		RooRealVar& ebeam) {

	FILE *fp = fopen(asciifilename.c_str(), "r");
	if (!fp) {
		std::cout << "Failed to open " << asciifilename << std::endl;
		perror("fopen");
		exit(1);
	} else {
		std::cout << "Creating new RooDataSet data set with " << asciifilename << " as input." << std::endl;
	}

	RooArgSet m_varset(dt, mbc, deltaE, expno, CosThetaB, ntrk_sig, ntrk_tag, z_err_sig, z_err_tag);
	m_varset.add(chisq_tracks_sig);
	m_varset.add(chisq_tracks_tag);
	m_varset.add(dgf_tracks_sig);
	m_varset.add(dgf_tracks_tag);
	m_varset.add(keeptagl_tag);
	m_varset.add(ecms);
	m_varset.add(rbin);
	m_varset.add(wtag);
	m_varset.add(delta_wtag);
	m_varset.add(flavor_tag);
	m_varset.add(ebeam);
	RooDataSet r_set("r_set", "r_set", m_varset);

	int m_icpv_exp, m_icpv_run, m_icpv_evt, m_icpv_frm, m_icpv_hepevt, m_icpv_mc_flv;
	double m_cos_theta_tr, m_cos_theta_1, m_icpv_mc_delta_t, m_icpv_mc_pbs, m_icpv_mc_costh;
	int m_icpv_decay, m_icpv_flv, m_icpv_rec_vt_ntrk, m_icpv_asc_vt_ntrk, m_icpv_rec_vt_ndf, m_icpv_asc_vt_ndf, m_icpv_keeptagl;
	double m_icpv_rec_vt_pos, m_icpv_rec_vt_err, m_icpv_rec_vt_chi2, m_icpv_rec_vt_cl;
	double m_icpv_asc_vt_pos, m_icpv_asc_vt_err, m_icpv_asc_vt_chi2, m_icpv_asc_vt_cl;
	double m_icpv_wtag_mc, m_icpv_deltae, m_icpv_mbc, m_icpv_pbs, m_icpv_costh, m_icpv_ebeam, m_icpv_costhtr;
	double m_icpv_dt, m_icpv_r, m_icpv_wtag, m_icpv_delta_wtag, m_icpv_ecms;
	int m_icpv_rbin;
	float m_dummy_f;

	char buffer[1024];

	while (fgets(buffer, sizeof (buffer), fp))
	{
		int ntoken = 0;

		if (!(pattern_match(buffer)) == true) continue;

		ntoken = sscanf(
				buffer,
				"%d %d %d %d "
				"%d "
				"%d "
				"%le %le %le "
				"%le %le "
				"%d "
				"%le %le %d %le %d %le "
				"%le %le %d %le %d %le "
				"%f "
				"%d  %d %le "
				"%le %le %le %le %le %le \n",

				&m_icpv_exp, &m_icpv_run, &m_icpv_evt, &m_icpv_frm,
				&m_icpv_hepevt,
				&m_icpv_mc_flv,
				&m_cos_theta_tr, &m_cos_theta_1,
				&m_icpv_mc_delta_t,
				&m_icpv_mc_pbs, &m_icpv_mc_costh,
				&m_icpv_decay,
				&m_icpv_rec_vt_pos, &m_icpv_rec_vt_err, &m_icpv_rec_vt_ntrk, &m_icpv_rec_vt_chi2, &m_icpv_rec_vt_ndf, &m_icpv_rec_vt_cl,
				&m_icpv_asc_vt_pos, &m_icpv_asc_vt_err, &m_icpv_asc_vt_ntrk, &m_icpv_asc_vt_chi2, &m_icpv_asc_vt_ndf, &m_icpv_asc_vt_cl,
				&m_dummy_f,
				&m_icpv_flv, &m_icpv_keeptagl, &m_icpv_wtag_mc,
				&m_icpv_deltae, &m_icpv_mbc, &m_icpv_pbs, &m_icpv_costh,
				&m_icpv_ebeam, &m_icpv_costhtr
		);

		if (ntoken != 34) continue; // have read ascii record for event, now calculate needed quantities (e.g. determine r bin, wtag, sigfrac, ...)

		//calculate variables
		m_icpv_dt = (m_icpv_rec_vt_pos - m_icpv_asc_vt_pos) * Belle::dt_resol_global::inv_bgc;

		m_icpv_r = 1 - 2 * m_icpv_wtag_mc;
		m_icpv_rbin = get_rbin(m_icpv_r);
		m_icpv_wtag = get_wtag(m_icpv_exp, m_icpv_rbin, mc);
		m_icpv_delta_wtag = get_delta_wtag(m_icpv_exp, m_icpv_rbin, mc);

		m_icpv_ecms = sqrt(m_icpv_pbs * m_icpv_pbs + Belle::dt_resol_global::mbzero * Belle::dt_resol_global::mbzero);

		//perform ICPV cuts
		if (m_icpv_rbin == 7) continue;
		if (m_icpv_rec_vt_err <= 0) continue;
		if (m_icpv_asc_vt_err <= 0) continue;

		if (m_icpv_rec_vt_ndf > 0 && m_icpv_rec_vt_chi2 / m_icpv_rec_vt_ndf >= 50) continue;
		if (m_icpv_asc_vt_ndf > 0 && m_icpv_asc_vt_chi2 / m_icpv_asc_vt_ndf >= 50) continue;

		if (m_icpv_rec_vt_ntrk == 1 && m_icpv_rec_vt_err >= 0.05) continue;
		if (m_icpv_asc_vt_ntrk == 1 && m_icpv_asc_vt_err >= 0.05) continue;

		if (m_icpv_rec_vt_ntrk > 1 && m_icpv_rec_vt_err >= 0.02) continue;
		if (m_icpv_asc_vt_ntrk > 1 && m_icpv_asc_vt_err >= 0.02) continue;

		typedef std::pair<double, RooRealVar*> valuevar;
		std::vector<valuevar> vars;
		vars.push_back(std::make_pair(m_icpv_dt, &dt));
		vars.push_back(std::make_pair(m_icpv_deltae, &deltaE));
		vars.push_back(std::make_pair(m_icpv_mbc, &mbc));
		//vars.push_back(std::make_pair(m_icpv_exp, &expno));
		vars.push_back(std::make_pair(m_icpv_costh, &CosThetaB));
		vars.push_back(std::make_pair(m_icpv_rec_vt_ntrk, &ntrk_sig));
		vars.push_back(std::make_pair(m_icpv_asc_vt_ntrk, &ntrk_tag));
		vars.push_back(std::make_pair(m_icpv_rec_vt_err, &z_err_sig));
		vars.push_back(std::make_pair(m_icpv_asc_vt_err, &z_err_tag));
		vars.push_back(std::make_pair(m_icpv_rec_vt_chi2, &chisq_tracks_sig));
		vars.push_back(std::make_pair(m_icpv_asc_vt_chi2, &chisq_tracks_tag));
		vars.push_back(std::make_pair(m_icpv_rec_vt_ndf, &dgf_tracks_sig));
		vars.push_back(std::make_pair(m_icpv_asc_vt_ndf, &dgf_tracks_tag));
		vars.push_back(std::make_pair(m_icpv_keeptagl, &keeptagl_tag));
		vars.push_back(std::make_pair(m_icpv_ecms, &ecms));
		vars.push_back(std::make_pair(m_icpv_wtag, &wtag));
		vars.push_back(std::make_pair(m_icpv_delta_wtag, &delta_wtag));
		vars.push_back(std::make_pair(m_icpv_ebeam, &ebeam));

		bool in_range = true;


		BOOST_FOREACH(valuevar &p, vars) {
			if (p.first < p.second->getMin() || p.first > p.second->getMax()) {
				in_range = false;
				break;
			} else
				p.second->setVal(p.first);
		}

		if (flavor_tag.setIndex(m_icpv_flv, false)) in_range = false;
		if (rbin.setIndex(m_icpv_rbin, false)) in_range = false;

		if (m_icpv_exp > 30)
			expno.setIndex(31, false);
		else
			expno.setIndex(7, false);

		//if (ntrk_sig.setIndex(m_icpv_rec_vt_ntrk, false)) in_range = false;
		//if (ntrk_tag.setIndex(m_icpv_asc_vt_ntrk, false)) in_range = false;
		//if (dgf_tracks_sig.setIndex(m_icpv_rec_vt_ndf, false)) in_range = false;
		//if (dgf_tracks_tag.setIndex(m_icpv_asc_vt_ndf, false)) in_range = false;
		//if (keeptagl_tag.setIndex(m_icpv_keeptagl, false)) in_range = false;

		if (in_range) r_set.add(m_varset);

	}
	std::cout << "number of events for fit: " << r_set.numEntries() << std::endl;
	fclose(fp);
	return r_set;
}

RooDataSet ascii2DataSet(const std::string& asciifilename,
		bool mc,
		RooRealVar& dt,
		RooRealVar& mbc,
		RooRealVar& deltaE,
		RooCategory& expno,
		RooRealVar& CosThetaB,
		RooRealVar& ntrk_sig,
		RooRealVar& ntrk_tag,
		RooRealVar& z_err_sig,
		RooRealVar& z_err_tag,
		RooRealVar& chisq_tracks_sig,
		RooRealVar& chisq_tracks_tag,
		RooRealVar& dgf_tracks_sig,
		RooRealVar& dgf_tracks_tag,
		RooRealVar& keeptagl_tag,
		RooRealVar& ecms,
		RooCategory& rbin,
		RooRealVar& wtag,
		RooRealVar& delta_wtag,
		RooCategory& flavor_tag,
		RooRealVar& ebeam,
		std::vector<read_in_variable> add_var) {

	FILE *fp = fopen(asciifilename.c_str(), "r");
	if (!fp) {
		std::cout << "Failed to open " << asciifilename << std::endl;
		perror("fopen");
		exit(1);
	} else {
		std::cout << "Creating new RooDataSet data set with " << asciifilename << " as input." << std::endl;
	}

	RooArgSet m_varset(dt, mbc, deltaE, expno, CosThetaB, ntrk_sig, ntrk_tag, z_err_sig, z_err_tag);
	m_varset.add(chisq_tracks_sig);
	m_varset.add(chisq_tracks_tag);
	m_varset.add(dgf_tracks_sig);
	m_varset.add(dgf_tracks_tag);
	m_varset.add(keeptagl_tag);
	m_varset.add(ecms);
	m_varset.add(rbin);
	m_varset.add(wtag);
	m_varset.add(delta_wtag);
	m_varset.add(flavor_tag);
	m_varset.add(ebeam);

	BOOST_FOREACH(read_in_variable& ri_var, add_var) {
		m_varset.add(*(ri_var.first));
	}

	RooDataSet r_set("r_set", "r_set", m_varset);

	double d_readin[34];
	int i_readin[34];
	//int m_icpv_exp, m_icpv_run, m_icpv_evt, m_icpv_frm, m_icpv_hepevt, m_icpv_mc_flv;
	//double m_cos_theta_tr, m_cos_theta_1, m_icpv_mc_delta_t, m_icpv_mc_pbs, m_icpv_mc_costh;
	//int m_icpv_decay, m_icpv_flv, m_icpv_rec_vt_ntrk, m_icpv_asc_vt_ntrk, m_icpv_rec_vt_ndf, m_icpv_asc_vt_ndf, m_icpv_keeptagl;
	//double m_icpv_rec_vt_pos, m_icpv_rec_vt_err, m_icpv_rec_vt_chi2, m_icpv_rec_vt_cl;
	//double m_icpv_asc_vt_pos, m_icpv_asc_vt_err, m_icpv_asc_vt_chi2, m_icpv_asc_vt_cl;
	//double m_icpv_wtag_mc, m_icpv_deltae, m_icpv_mbc, m_icpv_pbs, m_icpv_costh, m_icpv_ebeam, m_icpv_costhtr;
	double m_icpv_dt, m_icpv_r, m_icpv_wtag, m_icpv_delta_wtag, m_icpv_ecms;
	int m_icpv_rbin;
	//float m_dummy_f;

	char buffer[1024];

	while (fgets(buffer, sizeof (buffer), fp))
	{
		int ntoken = 0;

		if (!(pattern_match(buffer)) == true) continue;

		std::string scan_string = "%d %d %d %d %d %d %le %le %le %le %le %d %le %le %d %le %d %le %le %le %d %le %d %le %f %d %d %le %le %le %le %le %le %le \n";

		ntoken = sscanf(
				buffer,
				scan_string.c_str(),

				&i_readin[0], &d_readin[1], &d_readin[2], &d_readin[3], &d_readin[4], &d_readin[5], &d_readin[6], &d_readin[7], &d_readin[8], &d_readin[9],
				&d_readin[10], &d_readin[11], &d_readin[12], &d_readin[13], &i_readin[14], &d_readin[15], &i_readin[16], &d_readin[17], &d_readin[18], &d_readin[19],
				&i_readin[20], &d_readin[21], &i_readin[22], &d_readin[23], &d_readin[24], &i_readin[25], &i_readin[26], &d_readin[27], &d_readin[28], &d_readin[29],
				&d_readin[30], &d_readin[31], &d_readin[32], &d_readin[33]
				                                                       //&m_icpv_exp, &m_icpv_run, &m_icpv_evt, &m_icpv_frm,
				                                                       //&m_icpv_hepevt,
				                                                       //&m_icpv_mc_flv,
				                                                       //&m_cos_theta_tr, &m_cos_theta_1,
				                                                       //&m_icpv_mc_delta_t,
				                                                       //&m_icpv_mc_pbs, &m_icpv_mc_costh,
				                                                       //&m_icpv_decay,
				                                                       //&m_icpv_rec_vt_pos, &m_icpv_rec_vt_err, &m_icpv_rec_vt_ntrk, &m_icpv_rec_vt_chi2, &m_icpv_rec_vt_ndf, &m_icpv_rec_vt_cl,
				                                                       //&m_icpv_asc_vt_pos, &m_icpv_asc_vt_err, &m_icpv_asc_vt_ntrk, &m_icpv_asc_vt_chi2, &m_icpv_asc_vt_ndf, &m_icpv_asc_vt_cl,
				                                                       //&m_dummy_f,
				                                                       //&m_icpv_flv, &m_icpv_keeptagl, &m_icpv_wtag_mc,
				                                                       //&m_icpv_deltae, &m_icpv_mbc, &m_icpv_pbs, &m_icpv_costh,
				                                                       //&m_icpv_ebeam, &m_icpv_costhtr
		);

		if (ntoken != 34) continue; // have read ascii record for event, now calculate needed quantities (e.g. determine r bin, wtag, sigfrac, ...)

		//calculate variables
		double m_icpv_rec_vt_pos = d_readin[12];
		double m_icpv_asc_vt_pos = d_readin[18];
		double m_icpv_wtag_mc = d_readin[27];
		int m_icpv_exp = i_readin[0];
		double m_icpv_pbs = d_readin[30];
		double m_icpv_rec_vt_err = d_readin[13];
		double m_icpv_asc_vt_err = d_readin[19];
		int m_icpv_rec_vt_ndf = i_readin[16];
		int m_icpv_asc_vt_ndf = i_readin[22];
		double m_icpv_rec_vt_chi2 = d_readin[15];
		double m_icpv_asc_vt_chi2 = d_readin[21];
		int m_icpv_rec_vt_ntrk = i_readin[14];
		int m_icpv_asc_vt_ntrk = i_readin[20];
		double m_icpv_deltae = d_readin[28];
		double m_icpv_mbc = d_readin[29];
		double m_icpv_costh = d_readin[31];
		int m_icpv_keeptagl = i_readin[26];
		double m_icpv_ebeam = d_readin[32];
		int m_icpv_flv = i_readin[25];

		m_icpv_dt = (m_icpv_rec_vt_pos - m_icpv_asc_vt_pos) * Belle::dt_resol_global::inv_bgc;

		m_icpv_r = 1 - 2 * m_icpv_wtag_mc;
		m_icpv_rbin = get_rbin(m_icpv_r);
		m_icpv_wtag = get_wtag(m_icpv_exp, m_icpv_rbin, mc);
		m_icpv_delta_wtag = get_delta_wtag(m_icpv_exp, m_icpv_rbin, mc);

		m_icpv_ecms = sqrt(m_icpv_pbs * m_icpv_pbs + Belle::dt_resol_global::mbzero * Belle::dt_resol_global::mbzero);

		//perform ICPV cuts
		if (m_icpv_rbin == 7) continue;
		if (m_icpv_rec_vt_err <= 0) continue;
		if (m_icpv_asc_vt_err <= 0) continue;

		if (m_icpv_rec_vt_ndf > 0 && m_icpv_rec_vt_chi2 / m_icpv_rec_vt_ndf >= 50) continue;
		if (m_icpv_asc_vt_ndf > 0 && m_icpv_asc_vt_chi2 / m_icpv_asc_vt_ndf >= 50) continue;

		if (m_icpv_rec_vt_ntrk == 1 && m_icpv_rec_vt_err >= 0.05) continue;
		if (m_icpv_asc_vt_ntrk == 1 && m_icpv_asc_vt_err >= 0.05) continue;

		if (m_icpv_rec_vt_ntrk > 1 && m_icpv_rec_vt_err >= 0.02) continue;
		if (m_icpv_asc_vt_ntrk > 1 && m_icpv_asc_vt_err >= 0.02) continue;

		typedef std::pair<double, RooRealVar*> valuevar;
		std::vector<valuevar> vars;
		vars.push_back(std::make_pair(m_icpv_dt, &dt));
		vars.push_back(std::make_pair(m_icpv_deltae, &deltaE));
		vars.push_back(std::make_pair(m_icpv_mbc, &mbc));
		//vars.push_back(std::make_pair(m_icpv_exp, &expno));
		vars.push_back(std::make_pair(m_icpv_costh, &CosThetaB));
		vars.push_back(std::make_pair(m_icpv_rec_vt_ntrk, &ntrk_sig));
		vars.push_back(std::make_pair(m_icpv_asc_vt_ntrk, &ntrk_tag));
		vars.push_back(std::make_pair(m_icpv_rec_vt_err, &z_err_sig));
		vars.push_back(std::make_pair(m_icpv_asc_vt_err, &z_err_tag));
		vars.push_back(std::make_pair(m_icpv_rec_vt_chi2, &chisq_tracks_sig));
		vars.push_back(std::make_pair(m_icpv_asc_vt_chi2, &chisq_tracks_tag));
		vars.push_back(std::make_pair(m_icpv_rec_vt_ndf, &dgf_tracks_sig));
		vars.push_back(std::make_pair(m_icpv_asc_vt_ndf, &dgf_tracks_tag));
		vars.push_back(std::make_pair(m_icpv_keeptagl, &keeptagl_tag));
		vars.push_back(std::make_pair(m_icpv_ecms, &ecms));
		vars.push_back(std::make_pair(m_icpv_wtag, &wtag));
		vars.push_back(std::make_pair(m_icpv_delta_wtag, &delta_wtag));
		vars.push_back(std::make_pair(m_icpv_ebeam, &ebeam));

		BOOST_FOREACH(read_in_variable& ri_var, add_var) {
			vars.push_back(std::make_pair(d_readin[ri_var.second - 1], ri_var.first));
		}

		bool in_range = true;

		BOOST_FOREACH(valuevar &p, vars) {
			if (p.first < p.second->getMin() || p.first > p.second->getMax()) {
				in_range = false;
				break;
			} else
				p.second->setVal(p.first);
		}

		if (flavor_tag.setIndex(m_icpv_flv, false)) in_range = false;
		if (rbin.setIndex(m_icpv_rbin, false)) in_range = false;

		if (m_icpv_exp > 30)
			expno.setIndex(31, false);
		else
			expno.setIndex(1, false);

		//if (ntrk_sig.setIndex(m_icpv_rec_vt_ntrk, false)) in_range = false;
		//if (ntrk_tag.setIndex(m_icpv_asc_vt_ntrk, false)) in_range = false;
		//if (dgf_tracks_sig.setIndex(m_icpv_rec_vt_ndf, false)) in_range = false;
		//if (dgf_tracks_tag.setIndex(m_icpv_asc_vt_ndf, false)) in_range = false;
		//if (keeptagl_tag.setIndex(m_icpv_keeptagl, false)) in_range = false;

		if (in_range) r_set.add(m_varset);

	}
	std::cout << "number of events for fit: " << r_set.numEntries() << std::endl;
	fclose(fp);
	return r_set;
}

bool pattern_match(const std::string& temp_buf) {
	static regex_t preg;
	const char *buf = temp_buf.c_str();

	static bool first = true;
	if (first) {
		regcomp(&preg, "^[-+.0-9 \tEe]*\n$", 0);
		first = false;
	}

	if (regexec(&preg, buf, 0, NULL, 0) == REG_NOMATCH) {
		fprintf(stderr, "skipped bad event record (regexec): %s", buf);
		return false;
	}

	return true;
}

inline int get_rbin(double r) {
	return (0. <= r && r <= 0.1 ? 0 :
	0.1 < r && r <= 0.25 ? 1 :
	0.25 < r && r <= 0.5 ? 2 :
	0.5 < r && r <= 0.625 ? 3 :
	0.625 < r && r <= 0.75 ? 4 :
	0.75 < r && r <= 0.875 ? 5 :
	0.875 < r && r <= 1.0 ? 6 : 7);
}

inline double get_wtag(int expno, int rbin, bool mc) {
	double w_svd1_data[7] = {0.5, 0.418852, 0.329879, 0.233898, 0.170608, 0.099791, 0.0228501};

	double w_svd2_data[7] = {0.5, 0.418826, 0.319303, 0.222948, 0.163191, 0.104085, 0.0251454};

	double w_svd1_mc[7] = {0.5, 0.420827, 0.300296, 0.219317, 0.154636, 0.0916131, 0.0228891};

	double w_svd2_mc[7] = {0.5, 0.412222, 0.307838, 0.212765, 0.149933, 0.0913264, 0.0218754};

	if (mc) {
		if (expno < 30) {
			return w_svd1_mc[rbin];
		} else 	{
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

inline double get_delta_wtag(int expno, int rbin, bool mc) {
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
