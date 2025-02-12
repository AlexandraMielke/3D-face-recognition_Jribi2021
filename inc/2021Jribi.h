/*
 * SPDX-FileCopyrightText: 2025 Alexandra Mielke <alexandra.mielke@smail.emt.h-brs.de>
 *
 * SPDX-License-Identifier: MIT
 */

#ifndef _2021_JRIBI_H_
#define _2021_JRIBI_H_

//============================================== DEPENDENCIES ==============================================
//General dependencies
#include <map>
#include <vector>

//Library includes
#include <boost/filesystem.hpp>
#include "../lib/openstl/openstl.h"

//Project includes
#include "../inc/LevelCurves.h"
#include "../inc/Math_helper.h"

//============================================== FUNCTION DECLARATIONS =====================================
/**
 * @brief Read mesh from file and create geodesic::Mesh with it.
 * 
 * @param filename Path to mesh file
 * @return Imported mesh with std::vectors of vertices and faces.
 */
inline const auto import_mesh(const char *filename)
{
	std::cout << "Reading STL file" << std::endl;
	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open()) {
		std::cerr << "Error: Unable to open file '" << filename << "'" << std::endl;
	}

	// Deserialize the triangles in either binary or ASCII format
	std::vector<openstl::Triangle> triangles = openstl::deserializeStl(file);
	file.close();

	if (TEST_PRINT) std::cout << "Imported " << triangles.size() << " triangles" << std::endl;

	//convert triangles to vertices and faces
	return openstl::convertToVerticesAndFaces(triangles);
}

/**
 * @brief Calculates and prints the 3D invariant description of the face located at the filename.
 * 
 * @param filename Path to STL file.
 * @param landmarks Landmarks for the STL file.
 * @return int Returns 0 if successful and 1 otherwise.
 */
int jribi_description(std::string filename, std::map<std::string, jribi::Landmark> landmarks);

#endif // _2021_JRIBI_H_

