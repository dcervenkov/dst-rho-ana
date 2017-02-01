/**
 *  @file    colors.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2014-08-27
 *
 *  @brief File defining perceptually equidistant colors
 *
 */

#ifndef COLORS_H_
#define COLORS_H_

#include "TROOT.h"
#include "TColor.h"
#include "TStyle.h"

namespace colors {

const double pastelColors[8][3] = {
		{114,147,203},
		{225,151, 76},
		{132,186, 91},
		{211, 94, 96},
		{128,133,133},
		{144,103,167},
		{171,104, 87},
		{204,194, 16},
};

const double saturatedColors[8][3] = {
		{ 57,106,177},
		{218,124, 48},
		{ 62,150, 81},
		{204, 37, 41},
		{ 83, 81, 84},
		{107, 76,154},
		{146, 36, 40},
		{148,139, 61},
};

/**
 * Replace the first 8 default ROOT colors with 8 perceptually equidistant
 * colors suitable for plotting. Replace the default 'colz' palette by a 
 * perceptually uniform one.
 */
void setColors(){
	TColor *color;
	for (int i = 0; i < 8; ++i) {
		color = gROOT->GetColor(i+2);
		color->SetRGB(pastelColors[i][0]/255., pastelColors[i][1]/255., pastelColors[i][2]/255.);
	}

	gStyle->SetPalette(kViridis);
}

} // namespace colors

#endif /* COLORS_H_ */
