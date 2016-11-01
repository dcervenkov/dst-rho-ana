/* 
 * File:   RooTatamiGeneral.h
 * Author: kronenbitter
 *
 * Created on October 28, 2011, 1:28 PM
 */

#ifndef ROOTATAMIHELPER_H
#define ROOTATAMIHELPER_H

#include "RooDataSet.h"
#include <string>
#include <vector>
#include <utility> 
#include <fstream>
#include "RooRealVar.h"
#include "RooCategory.h"

typedef std::pair<RooRealVar*,int> read_in_variable;

RooDataSet ascii2DataSet(const std::string& asciifilename,
		bool mc,
		RooRealVar& dt,
		RooRealVar& dmbc,
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
		RooRealVar& ebeam);

RooDataSet ascii2DataSet(const std::string& asciifilename,
		bool mc,
		RooRealVar& dt,
		RooRealVar& dmbc,
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
		std::vector<read_in_variable> add_var);

bool pattern_match(const std::string& temp_buf);

inline int get_rbin(double r);

inline double get_wtag(int expno, int rbin, bool mc);

inline double get_delta_wtag(int expno, int rbin, bool mc);

#endif	/* ROOTATAMIHELPER_H */
