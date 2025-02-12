/*
 * SPDX-FileCopyrightText: 2025 Alexandra Mielke <alexandra.mielke@smail.emt.h-brs.de>
 *
 * SPDX-License-Identifier: MIT
 */

/*
2021_Jribi.cpp: Application entry point.

This file is part of the implementation of Jribi et al. 2021: 
An SE(3) invariant description for 3D face recognition.

Developed for the institute of security research at the Bonn-Rhein-Sieg university of applied
sciences. See the COPYRIGHT file at the top-level directory of this distribution
for details on code ownership.
*/

//project imports
#include "../inc/2021Jribi.h"

int main()
{
	const char filename[] = "/home/nforce/Downloads/2024-10-14-synthetic-test.stl";
	//const char filename[] = "/home/nforce/Downloads/20240724141021698_13_Unknown_zuschnitt.stl";
	//const char filename[] = "/home/nforce/Projekte/5_Messdaten/Testdaten/NikolausHausASCII.stl";
	const auto [vertices, faces] = import_mesh(filename);

	//------------------------------------------- 1. Landmark detection
	//manually pick landmarks on ITIGKKIF
	std::map<std::string, jribi::Landmark> landmarks;	
	landmarks["Pronasale"] = { 0.517606, 0.067231, 1.632528, 0 }; //nose-tip
	landmarks["Endocanthion_left"] = { 0.497978, 0.109879, 1.660569, 0 };
	landmarks["Exocanthion_left"] = { 0.472765, 0.114076, 1.662901, 0 };
	landmarks["Endocanthion_right"] = { 0.537204, 0.109998, 1.660317, 0 };
	landmarks["Exocanthion_right"] = { 0.562480, 0.113942, 1.662759, 0 };

	//--------------------------- Create level curves
	math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> sampled_curve_nose;
	math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> sampled_curve_lefteye;
	math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> sampled_curve_righteye;	

	jribi::FaceMesh* facemesh = new jribi::FaceMesh(vertices, faces);
	facemesh->createLevelCurves(landmarks, sampled_curve_nose, 
										sampled_curve_lefteye, 
										sampled_curve_righteye);	
	delete facemesh;

	//--------------------------- Compare
	double mean_sampled = (math::sum(sampled_curve_nose,1,2) + math::sum(sampled_curve_lefteye,1,2) + math::sum(sampled_curve_righteye,1,2))/
							(jribi::NR_LEVEL_CURVES_NOSE*jribi::NR_SUBSAMPLES_NOSE + 2*jribi::NR_LEVEL_CURVES_EYE*jribi::NR_SUBSAMPLES_EYE);

	math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> read_curve_nose;
	math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> read_curve_lefteye;
	math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> read_curve_righteye;	
	const auto [mean_read, flag] = jribi::read_curves("./out/levelcurves.csv", read_curve_nose, read_curve_lefteye, read_curve_righteye);
	if(!flag) return 0;
	
	double sum_top = 0;
	double sum_sampled = 0;
	double sum_read = 0;
	int nr_level_curves = std::max(jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_LEVEL_CURVES_EYE);
	int nr_subsamples = std::max(jribi::NR_SUBSAMPLES_NOSE, jribi::NR_SUBSAMPLES_EYE);
	for(int i=0; i < nr_level_curves; i++)
	{
		for(int j=0; j < nr_subsamples; j++)
		{
			if((nr_level_curves < jribi::NR_LEVEL_CURVES_NOSE) && (nr_subsamples < jribi::NR_SUBSAMPLES_NOSE))
			{
				sum_top += (sampled_curve_nose[i][j][1] - mean_sampled) * (read_curve_nose[i][j][1] - mean_read);
				sum_sampled += std::pow((sampled_curve_nose[i][j][1] - mean_sampled),2);
				sum_read += std::pow((read_curve_nose[i][j][1] - mean_read),2);
			}
			if((nr_level_curves < jribi::NR_LEVEL_CURVES_EYE) && (nr_subsamples < jribi::NR_SUBSAMPLES_EYE))
			{
				sum_top += (sampled_curve_lefteye[i][j][1] - mean_sampled) * (read_curve_lefteye[i][j][1] - mean_read);
				sum_top += (sampled_curve_righteye[i][j][1] - mean_sampled) * (read_curve_righteye[i][j][1] - mean_read);
				sum_sampled += std::pow((sampled_curve_lefteye[i][j][1] - mean_sampled),2)
							+ std::pow((sampled_curve_righteye[i][j][1] - mean_sampled),2);
				sum_read += std::pow((read_curve_lefteye[i][j][1] - mean_read),2)
							+ std::pow((read_curve_righteye[i][j][1] - mean_read),2);
			}
		}
	}
	double correlation = sum_top/std::sqrt(sum_read*sum_sampled);
	std::cout << "Correlation: " << correlation << std::endl;
	return 1;
}

