/*
 * SPDX-FileCopyrightText: 2025 Alexandra Mielke <alexandra.mielke@smail.emt.h-brs.de>
 *
 * SPDX-License-Identifier: MIT
 */

// LevelCurves.h: Function definition for everything level curve related including Geodist struct and FaceMesh class.

#ifndef _LEVELCURVES_H_
#define _LEVELCURVES_H_

//============================================== DEPENDENCIES ==============================================
//General dependencies
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <map>

//Library includes
#include "../lib/geodesic/geodesic_algorithm_base.h"
#include "../lib/geodesic/geodesic_algorithm_dijkstra.h"
#include "../lib/geodesic/geodesic_algorithm_subdivision.h"
#include "../lib/geodesic/geodesic_algorithm_exact.h"
#include "../lib/openstl/openstl.h"
#include "../inc/Math_helper.h"
#include "../lib/SplineLib/Sources/Splines/b_spline.hpp"

//================================================ DEFINES ==================================================
#ifndef TEST_PRINT
#define TEST_PRINT true
#endif

#ifndef SUBDIVISION
#define SUBDIVISION
#endif

namespace jribi
{
#ifdef DIJKSTRA
typedef geodesic::GeodesicAlgorithmDijkstra GeodesicAlgorithm;
#endif

#ifdef EXACT
	typedef geodesic::GeodesicAlgorithmExact GeodesicAlgorithm;
#endif

#ifdef SUBDIVISION
	typedef geodesic::GeodesicAlgorithmSubdivision GeodesicAlgorithm;
#endif

typedef std::vector<openstl::Vec3> LevelCurve;


const int NR_LEVEL_CURVES_EYE = 9;
const int NR_LEVEL_CURVES_NOSE = 15;
const int NR_SUBSAMPLES_EYE = 20;
const int NR_SUBSAMPLES_NOSE = 40;
const int SUBDIVISION_LEVEL = 4;

enum FACEREGION {NOSE, LEFTEYE, RIGHTEYE};

struct Landmark{
	float x, y, z;
	unsigned int vertex_id;
	double dist_to_vertex = 1000;
};
struct Geodist{
	openstl::Vec3 xyz;
	float dist_lm1, dist_lm2, dist_lm3; //just for sorting purposes
	float sum,angle;
	int level_curve_nr;

	bool operator<(const Geodist &a)
	{
		if (level_curve_nr == a.level_curve_nr)
		{
			if(angle < a.angle) return true;
			else return false;
		}
		else
		{
			if(level_curve_nr < a.level_curve_nr) return true;
			else return false;
		}
		if(level_curve_nr < a.level_curve_nr) return true;
		else return false;
	}

	bool operator<=(const Geodist &a)
	{
		if (level_curve_nr == a.level_curve_nr)
		{
			if(angle <= a.angle) return true;
			else return false;
		}
		else
		{
			if(level_curve_nr < a.level_curve_nr) return true;
			else return false;
		}
	}

	bool operator>(const Geodist &a)
	{
		if (level_curve_nr == a.level_curve_nr)
		{
			if(angle > a.angle) return true;
			else return false;
		}
		else
		{
			if(level_curve_nr > a.level_curve_nr) return true;
			else return false;
		}
	}

	bool operator>=(const Geodist &a)
	{
		if (level_curve_nr == a.level_curve_nr)
		{
			if(angle >= a.angle) return true;
			else return false;
		}
		else
		{
			if(level_curve_nr > a.level_curve_nr) return true;
			else return false;
		}
	}

	bool operator==(const Geodist &a)
	{
		if (level_curve_nr == a.level_curve_nr)
		{
			if(angle == a.angle) return true;
			else return false;
		}
		else
		{
			return false;
		}
	}
};

//============================================== FUNCTION DECLARATIONS =====================================
class FaceMesh
{
public:
	/**
	 * @brief Construct a new Face Mesh object. 
	 * 
	 * @param vertices 
	 * @param faces 
	 */
	FaceMesh(const std::vector<openstl::Vec3> &vertices,
		const std::vector<openstl::Face> &faces);
	~FaceMesh();

	/**
	 * @brief Create a Level Curves object
	 * 
	 * @param landmarks Dictionary with landmarks of face. Should contain nose, eye center, inner and outer corner.
	 * @param sampled_curve_nose Output array for the nose region that stores level curves and their respective curvature.
	 * @param sampled_curve_lefteye Output array for the left-eye region that stores level curves and their respective curvature.
	 * @param sampled_curve_righteye Output array for the right-eye region that stores level curves and their respective curvature.
	 */
	void createLevelCurves(std::map<std::string, Landmark> &landmarks,
		math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> &sampled_curve_nose,
		math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_lefteye,
		math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_righteye);

	/**
	 * @brief Prints equal level curves from file. File is an ASCII csv with the following format:
	 * Curve_nr,shape_index,flag_valid\n
	 * nose,,\n
	 * ...data...
	 * lefteye,,\n
	 * ...data...
	 * righteye,,\n
	 * ...data...
	 * With data stored as comma separated values.
	 * 
	 * @param filename Path to file.
	 * @param sampled_curve_nose Array with curvature data of the equal level curves for the nose region.
	 * @param sampled_curve_lefteye Array with curvature data of the equal level curves for the lefteye region.
	 * @param sampled_curve_righteye Array with curvature data of the equal level curves for the righteye region.
	 */
	void print_curves(const char *filename, 
		math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> &sampled_curve_nose,
		math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_lefteye,
		math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_righteye);

private:
	geodesic::Mesh mesh;
	openstl::Vec3 face_sep;
	openstl::Vec3 eye2eye_vec;
	openstl::Vec3 facenormal_vec;
	openstl::Vec3 mid_nose;
	openstl::Vec3 mid_lefteye;
	openstl::Vec3 mid_righteye;

	//=================================================== GEODESIC FUNCTIONS
	/**
	 * @brief Propagate geodesic distance from a landmark.
	 * 
	 * @param vertex_id ID of the landmark's vertex within the mesh.
	 * @param mesh 
	 * @param propagation_limit Maximum propagation distance in the mesh's unit. Defaults to 1.0E100
	 * @return GeodesicAlgorithm 
	 */
	GeodesicAlgorithm propagate_from_landmark(unsigned int vertex_id, 
												geodesic::Mesh &mesh,
												double propagation_limit = geodesic::GEODESIC_INF);
	/**
	 * @brief Get ID of closest vertex to a landmark. The landmark's coordinates do not have to be vertices 
	 * of the mesh.
	 * 
	 * @param x X coordinate of landmark.
	 * @param y Y coordinate of landmark.
	 * @param z Z coordinate of landmark.
	 */
	int closest_vertexID(float x, float y, float z);

	/**
	 * @brief Get the ID of the vertex nearest to the point in the path at path_index.
	 * 
	 * @param path_index Index of point in path.
	 * @param path Geodesic path.
	 * @return unsigned int ID of nearest vertex.
	 */
	unsigned int closest_vertexID_inMesh(unsigned int path_index, std::vector<geodesic::SurfacePoint> path);

	/**
	 * @brief Calculate cartesian distance between two vertices.
	 * 
	 * @param vertexA First vertex.
	 * @param vertexB Second vertex. Order not important.
	 * @return double Cartesian distance.
	 */
	double cartesian_dist(geodesic::Vertex* vertexA, geodesic::Vertex* vertexB);

	//=================================================== BIN & SORT CURVE FUNCTIONS
	/**
	 * @brief Bins the multipolar parametrization distances into nr_level_curves+1 level curves.
	 * The respective level curve number is stored in the Geodist struct under level_curve_nr. If the
	 * point is not part of a level curve, level_curve_nr is set to 0.
	 * 
	 * @param geodist Geodist struct with multipolar distance.
	 * @param level_curves Vector with the level curves. Each level curve is a vector of openstl::Vec3.
	 * @param faceregion Type of region in face. Can be any of jribi::FACEREGION.
	 */
	std::vector<LevelCurve> bin_level_curves(std::vector<Geodist> &geodist, 
			std::map<std::string, Landmark> &landmarks,
			const int faceregion);

	/**
	 * @brief Calculates and sums up the geodesic distance between target and two landmarks.
	 * 
	 * @param geomesh_landmark1 Propagated mesh from landmark 1.
	 * @param geomesh_landmark2 Propagated mesh from landmark 2.
	 * @param target Target vertice of the mesh.
	 * @param geodist_buff Buffer to store temporary geodesic distances.
	 * @param geodist Vector with the vertices and their respective mutlipolar distances. 
	 */
	void bipolar_param(	GeodesicAlgorithm &geomesh_landmark1, 
					GeodesicAlgorithm &geomesh_landmark2,
					geodesic::SurfacePoint &target,
					Geodist &geodist_buff,
					std::vector<Geodist> &geodist);

	/**
	 * @brief Calculates and sums up the geodesic distance between target and two landmarks.
	 * 
	 * @param geomesh_landmark1 Propagated mesh from landmark 1.
	 * @param geomesh_landmark2 Propagated mesh from landmark 2.
	 * @param geomesh_landmark3 Propahated mesh from landmark 3.
	 * @param target Target vertice of the mesh.
	 * @param geodist_buff Buffer to store temporary geodesic distances.
	 * @param geodist Vector with the vertices and their respective mutlipolar distances.
	 */
	void threepolar_param(GeodesicAlgorithm &geomesh_landmark1, 
					GeodesicAlgorithm &geomesh_landmark2,
					GeodesicAlgorithm &geomesh_landmark3,
					geodesic::SurfacePoint &target,
					Geodist &geodist_buff,
					std::vector<Geodist> &geodist);

	/**
	 * @brief Calculates the angle.
	 * 
	 * @param geodist Geodist struct with multipolar distance.
	 * @param mid Midpoint of the face region.
	 */
	void calc_angle(std::vector<Geodist> &geodist, openstl::Vec3 mid);

	//=================================================== INTERPOLATION FUNCTIONS
	/**
	 * @brief Returns B-Spline interpolation of the coordinates in curve_data. Samples of the B-Spline
	 * are stored in sampled_curve.
	 * 
	 * @param curve_data Storage vector for a single level curves' coordinates.
	 * @param sampled_curve Storage vector for the B-spline samples. Returned are the closest vertex_id
	 * and this point's respective principal curvatures (min and max).
	 * of 
	 * @param faceregion Face region determines the number of samples. Can be either of jribri::FACEREGION.
	 */
	void bspline_interpolation(LevelCurve &curve_data, double sampled_curve[][3], int faceregion);

	/**
	 * @brief Handles the B-Spline interpolation and calls bspline_interpolation() for each level curve in a
	 * face region.
	 * 
	 * @param level_curves Storage vector for all level curves of a face region.
	 * @param faceregion Face region determines the number of levelcurves and samples. Can be either
	 * of jribri::FACEREGION.
	 */
	template<std::size_t NR_LEVEL_CURVES, std::size_t NR_SUBSAMPLES>
	void interpolate_curves(std::vector<LevelCurve> &level_curves, int faceregion,
							math::Array3d<double, NR_LEVEL_CURVES, NR_SUBSAMPLES, 3> &sampled_curve)
	{
		double interpolated_curve_i[NR_SUBSAMPLES][3]; //temporary storage for interpolated curve
		LevelCurve level_curves_i; //temporary storage for single level curve to shrink the curve
		math::Array2d<double, NR_LEVEL_CURVES, NR_SUBSAMPLES> buff;
		for(int k=0; k<NR_LEVEL_CURVES; k++) 
		{
			level_curves_i = level_curves[k];
			//Shrink curve 
			size_t i = 0;
			while(i<level_curves_i.size())
			{
				std::vector<size_t> toRemove;
				float count = 0.0; //to count the number of points within the distance threshold

				// Check subsequent points for distance and accumulate weight
				for (size_t j = i + 1; j < level_curves_i.size(); ++j) 
				{
					if (math::cartesian_dist(level_curves_i[i],level_curves_i[j]) <= 2*this->mesh.average_edge) 
					{
						level_curves_i[i].x += level_curves_i[j].x;
						level_curves_i[i].y += level_curves_i[j].y;
						level_curves_i[i].z += level_curves_i[j].z;
						count += 1.0;
						toRemove.push_back(j);
					}
				}

				// Update the weight of the current point
				level_curves_i[i].x = level_curves_i[i].x/(count+1);
				level_curves_i[i].y = level_curves_i[i].y/(count+1);
				level_curves_i[i].z = level_curves_i[i].z/(count+1);

				// Remove all points within the distance threshold
				for (auto it = toRemove.rbegin(); it != toRemove.rend(); ++it) {
					level_curves_i.erase(level_curves_i.begin() + *it);
				}
				i++;
			}

			bspline_interpolation(level_curves_i, interpolated_curve_i, faceregion);
			for(int j=0; j<NR_SUBSAMPLES; j++)
			{
				sampled_curve[k][j][0] = interpolated_curve_i[j][0]; //index
				buff[k][j] = interpolated_curve_i[j][1]; //shape index
				sampled_curve[k][j][2] = interpolated_curve_i[j][2]; //valid flag
			}
		}
		float w0 = 100.0;
		float w1 = 10.0;
		float w2 = 1.0;
		float sum1, sum2, count1, count2;
		//To Do: ignore invalid points
		for(int i=0; i<NR_LEVEL_CURVES; i++)
		{
			for(int j=0; j<NR_SUBSAMPLES; j++)
			{
				int j_p1, j_m1, j_p2, j_m2; //indices for the ring, needed for overflow
				j_m1 = j-1;
				j_p1 = j+1;
				j_m2 = j-2;
				j_p2 = j+2;
				if(j-1 < 0) j_m1 = NR_SUBSAMPLES-j-1;
				if(j+1 > NR_SUBSAMPLES-1) j_p1 = j+1-NR_SUBSAMPLES;
				if(j-2 < 0) j_m2 = NR_SUBSAMPLES-j-2;
				if(j+2 > NR_SUBSAMPLES-1) j_p2 = j+2-NR_SUBSAMPLES;
				
				if((i > 0) && (i < NR_LEVEL_CURVES-1)) //default for ring1
				{
					sum1 = buff[i-1][j] + buff[i+1][j] + buff[i][j_m1] + buff[i][j_p1]; //normal case ring1
					count1 = 4;
					if((i > 1) && (i < NR_LEVEL_CURVES-2)) //default for ring2
					{
						sum2 = buff[i][j_m2] + buff[i+1][j_m1] + buff[i+2][j] + buff[i+1][j_p1] +
							buff[i][j_p2] + buff[i-1][j_p1] + buff[i-2][j] + buff[i-1][j_m1];
						count2 = 8;
					}
					else //exception for ring2 for level curves 2 and N-2 
					{
						if(i == 1) sum2 = buff[i+2][j] + buff[i][j_m2] + buff[i][j_p2] + buff[i-1][j_m1]
										+ buff[i-1][j_p1] + buff[i+1][j_m1] + buff[i+1][j_p1];
						if(i == NR_LEVEL_CURVES-2) sum2 = buff[i-2][j] + buff[i][j_m2] + buff[i][j_p2]
										+ buff[i-1][j_m1] + buff[i-1][j_p1] + buff[i+1][j_m1] + buff[i+1][j_p1];
						count2 = 7;
					}
				}
				else 
				{						
					if(i == 0) 
					{
						sum1 = buff[i+1][j] + buff[i][j_m1] + buff[i][j_p1];
						sum2 = buff[i+2][j] + buff[i][j_m2] + buff[i][j_p2] + buff[i+1][j_m1] + buff[i+1][j_p1];
					}
					if(i == NR_LEVEL_CURVES-1) 
					{
						sum1 = buff[i-1][j] + buff[i][j_m1] + buff[i][j_p1];
						sum2 = buff[i-2][j] + buff[i-1][j_p1] + buff[i][j_p2] + buff[i][j_m2] + buff[i-1][j_m1];
					}
					count1 = 3;
					count2 = 5;
				}
				sampled_curve[i][j][1] = (w0*buff[i][j] + w1*sum1 + w2*sum2)/(w0 + count1*w1 + count2*w2); //copy to output array
			}
		}
	}

	/**
	 * @brief Determines the curves starting point. This point lies on the plane with eye2eye_vector as normal and
	 * the mid of the respective face region as point. This ensures that all samples are started at the same position.
	 * The point is determined by choosing the closest point to the plane, either first or last data of the curve, and
	 * then calculating the footpoint of this point on the plane.
	 * 
	 * @param first_curve_data First coordinate of the respective curve.
	 * @param last_curve_data Last coordinate of the respective curve.
	 * @param faceregion Face region determines the number of levelcurves and samples. Can be either
	 * of jribri::FACEREGION.
	 * @return splinelib::Array<splinelib::Coordinate, 3UL> 
	 */
	splinelib::Array<splinelib::Coordinate, 3UL> curve_startingPoint(openstl::Vec3 first_curve_data, openstl::Vec3 last_curve_data, int faceregion);

	/**
	 * @brief Determines Gaussian and mean curvature of a point.
	 * 
	 * @param vertex_id Respective point of the mesh.
	 * @return std::tuple<double, double> Gaussian curvature and mean curvature as float.
	 */
	std::tuple<double, double> principalCurvature(int vertex_id);
	//=================================================== PRINT FUNCTIONS
	/**
	 * @brief Prints geodesic path information to ASCII csv. Meant for debugging purposes.
	 * 
	 * @param filename Path to output file.
	 * @param path Path of gedodesic. 
	 */
	void print_geodesicpath(const char *filename, std::vector<geodesic::SurfacePoint> path);

	/**
	 * @brief Prints sum of multipolar distances and respective vertex coordinates to ASCII csv. Meant for debugging purposes.
	 * 
	 * @param filename Path to output file.
	 * @param geodist Vector with multipolar parametrization distance.
	 */
	void print_multipolar_dist(const char *filename, std::vector<Geodist> &geodist);

	/**
	 * @brief Prints the sorted level curves of a specific face region to an ASCII csv. Meant for debugging purposes.
	 * 
	 * @param filename Path to output file.
	 * @param level_curves Vector with the coordinates of the level curve.
	 * @param nr_level_curve Number of level curves for the respective face region.
	 */
	void print_sortedcurve(const char *filename, std::vector<LevelCurve> level_curves, int nr_level_curve);

	/**
	 * @brief Prints the interpolated level curves of all face regions to a single ASCII csv for each face region. 
	 * Meant for debugging purposes.
	 * 
	 * @param sampled_curve_nose 
	 * @param sampled_curve_lefteye 
	 * @param sampled_curve_righteye 
	 */
	void print_interpolatedCurve(
	math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> &sampled_curve_nose,
	math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_lefteye,
	math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_righteye);
};	

//=================================================== MAIN
/**
 * @brief Read equal level curves from file. File should be an ASCII csv with the following format:
 * Curve_nr,shape_index,flag_valid\n
 * nose,,\n
 * ...data...
 * lefteye,,\n
 * ...data...
 * righteye,,\n
 * ...data...
 * With data stored as comma separated values.
 * 
 * @param filename Path to file including filename.
 * @param sampled_curve_nose Array to which the equal level curves for the nose region are stored. 
 * @param sampled_curve_lefteye Array to which the equal level curves for the lefteye region are stored. 
 * @param sampled_curve_righteye Array to which the equal level curves for the righteye region are stored. 
 * @return std::tuple<float, int> If second value is 1, the first returns the mean curvature. If the second is 0, an error occured.
 */
std::tuple<double, int> read_curves(const char *filename, 
			math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> &sampled_curve_nose,
			math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_lefteye,
			math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_righteye);

} //namespace jribi

#endif // _LEVELCURVES_H_

