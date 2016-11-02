/*
 * Constants.h
 *
 *  Created on: Oct 9, 2015
 *      Author: cervenkov
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const char* OUTPUT_FORMAT = ".gif";
const int MIN_ENTRIES = 0;

const int NUM_MC_FLAGS = 9;
const int MC_FLAGS_ASSIGNMENT[NUM_MC_FLAGS] = { -1, 1, 2, 3,  4, 5,  6, 7,  8};
const int COLOR_MAP[NUM_MC_FLAGS] = 		  { 13, 4, 7, 6, 20, 2, 30, 8, 31};

// Signal MC flags; used for figure-of-merit calculation
const int MC_FLAG_SIGNAL = 1;
const int MC_FLAG_SIGNAL_WITHOUT_FSR = 10;

// BASF identifies histograms via numbers
const int BASF_HIST_ID_BEST = 2000; // Best candidate

#endif /* CONSTANTS_H_ */
