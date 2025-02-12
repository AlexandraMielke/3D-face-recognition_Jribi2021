/*
 * SPDX-FileCopyrightText: 2025 Alexandra Mielke <alexandra.mielke@smail.emt.h-brs.de>
 *
 * SPDX-License-Identifier: MIT
 */

//project imports
#include "../inc/LevelCurves.h"

namespace jribi
{
//============================================== FACEMESH CLASS =====================================
//=================================================== PUBLIC FUNCTIONS
FaceMesh::FaceMesh(const std::vector<openstl::Vec3> &vertices,
	            const std::vector<openstl::Face> &faces)
{
	std::cout << "Initialising mesh data..." << std::endl;
	this->mesh.initialize_mesh_data(vertices, faces);		//create internal mesh data structure including edges
	std::cout << "Done.\n" << std::endl;
}
FaceMesh::~FaceMesh(){}

void FaceMesh::createLevelCurves(std::map<std::string, Landmark> &landmarks,
			math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> &sampled_curve_nose,
			math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_lefteye,
			math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_righteye)	
{
	//Determine eye centers as middle between Endocanthion and Exocanthion
	landmarks["leftEye_center"] = { (landmarks["Endocanthion_left"].x + landmarks["Exocanthion_left"].x)/2, 
									(landmarks["Endocanthion_left"].y + landmarks["Exocanthion_left"].y)/2, 
									(landmarks["Endocanthion_left"].z + landmarks["Exocanthion_left"].z)/2, 0 };
	landmarks["rightEye_center"] = { (landmarks["Endocanthion_right"].x + landmarks["Exocanthion_right"].x)/2,
									(landmarks["Endocanthion_right"].y + landmarks["Exocanthion_right"].y)/2,
									(landmarks["Endocanthion_right"].z + landmarks["Exocanthion_right"].z)/2, 0 };
	this->face_sep = { //get vector that seperates the face in the middle. We use nose to mid of eyes as an estimate here.
			std::abs(landmarks["Pronasale"].x - (landmarks["Endocanthion_right"].x - landmarks["Endocanthion_left"].x)/2), 
			std::abs(landmarks["Pronasale"].y - (landmarks["Endocanthion_right"].y - landmarks["Endocanthion_left"].y)/2), 
			std::abs(landmarks["Pronasale"].z - (landmarks["Endocanthion_right"].z - landmarks["Endocanthion_left"].z)/2)};
	this->eye2eye_vec = { //get normal to seperating vector to determine if angle is positive or negative
			landmarks["Endocanthion_left"].x - landmarks["Endocanthion_right"].x, 
			landmarks["Endocanthion_left"].y - landmarks["Endocanthion_right"].y, 
			landmarks["Endocanthion_left"].z - landmarks["Endocanthion_right"].z};
	this->facenormal_vec = math::cross_product(this->eye2eye_vec, this->face_sep);
	//Get closest vertex to each landmark
	for(std::map<std::string, Landmark>::iterator it = landmarks.begin(); it != landmarks.end(); ++it)
	{
		it->second.vertex_id = closest_vertexID(it->second.x, it->second.y, it->second.z);
	}
	std::cout << "Calculating geodesic distances..." << std::endl;
	//Level curves
	//Set propagation maximum to maximum distance between eye centers and nose
	double const distance_limit_nose = 1.3*std::max(
		cartesian_dist(&mesh.vertices()[landmarks["Pronasale"].vertex_id], &mesh.vertices()[landmarks["leftEye_center"].vertex_id]), 
		cartesian_dist(&mesh.vertices()[landmarks["Pronasale"].vertex_id], &mesh.vertices()[landmarks["rightEye_center"].vertex_id])); 

	double const distance_limit_eye = 1.5*std::max(
		cartesian_dist(&mesh.vertices()[landmarks["Endocanthion_left"].vertex_id], &mesh.vertices()[landmarks["Exocanthion_left"].vertex_id]), 
		cartesian_dist(&mesh.vertices()[landmarks["Endocanthion_right"].vertex_id], &mesh.vertices()[landmarks["Exocanthion_right"].vertex_id])); 

	GeodesicAlgorithm geomesh_nose = propagate_from_landmark(landmarks["Pronasale"].vertex_id, mesh, distance_limit_nose);
	GeodesicAlgorithm geomesh_leftEye_inner = propagate_from_landmark(landmarks["Endocanthion_left"].vertex_id, mesh, distance_limit_nose);
	GeodesicAlgorithm geomesh_rightEye_inner = propagate_from_landmark(landmarks["Endocanthion_right"].vertex_id, mesh, distance_limit_nose);
	GeodesicAlgorithm geomesh_leftEye_outer = propagate_from_landmark(landmarks["Exocanthion_left"].vertex_id, mesh, distance_limit_eye);
	GeodesicAlgorithm geomesh_rightEye_outer = propagate_from_landmark(landmarks["Exocanthion_right"].vertex_id, mesh, distance_limit_eye);
	GeodesicAlgorithm geomesh_rightEye_center = propagate_from_landmark(landmarks["rightEye_center"].vertex_id, mesh, distance_limit_eye);
	GeodesicAlgorithm geomesh_leftEye_center = propagate_from_landmark(landmarks["leftEye_center"].vertex_id, mesh, distance_limit_eye);	

	std::vector<Geodist> geodist_nose;
	std::vector<Geodist> geodist_lefteye;
	std::vector<Geodist> geodist_righteye;
	for(int i=0; i<mesh.vertices().size(); i++) 
	{
		geodesic::SurfacePoint target(&mesh.vertices()[i]);
		//nose
		Geodist geodist_buffer;
		geodist_buffer.sum = 0;
		threepolar_param(geomesh_nose, geomesh_leftEye_inner, geomesh_rightEye_inner, target, geodist_buffer, geodist_nose);	
		//left eye
		geodist_buffer.sum = 0;
		bipolar_param(geomesh_leftEye_outer, geomesh_leftEye_center, target, geodist_buffer, geodist_lefteye);
		//right eye
		geodist_buffer.sum = 0;
		bipolar_param(geomesh_rightEye_outer, geomesh_rightEye_center, target, geodist_buffer, geodist_righteye);
	}
	std::cout << "Done." << std::endl;
	
	std::cout << "Binning level curves..." << std::endl;
	std::vector<jribi::LevelCurve> levelcurves_righteye(jribi::NR_LEVEL_CURVES_EYE);
	std::vector<jribi::LevelCurve> levelcurves_lefteye(jribi::NR_LEVEL_CURVES_EYE);
	std::vector<jribi::LevelCurve> levelcurves_nose(jribi::NR_LEVEL_CURVES_NOSE);
	levelcurves_nose = bin_level_curves(geodist_nose, landmarks, FACEREGION::NOSE);
	levelcurves_lefteye = bin_level_curves(geodist_lefteye, landmarks, FACEREGION::LEFTEYE);
	levelcurves_righteye = bin_level_curves(geodist_righteye, landmarks , FACEREGION::RIGHTEYE);
	std::cout << "Done." << std::endl;

	if(TEST_PRINT)
	{
		print_multipolar_dist("./out/geodist2lefteye_binned.csv", geodist_lefteye);
		print_multipolar_dist("./out/geodist2righteye_binned.csv", geodist_righteye);
		print_multipolar_dist("./out/geodist2nose_binned.csv", geodist_nose);

		print_sortedcurve("./out/levelcurve_sorted_nose.csv", levelcurves_nose, NR_LEVEL_CURVES_NOSE);
		print_sortedcurve("./out/levelcurve_sorted_lefteye.csv", levelcurves_lefteye, NR_LEVEL_CURVES_EYE);
		print_sortedcurve("./out/levelcurve_sorted_righteye.csv", levelcurves_righteye, NR_LEVEL_CURVES_EYE);
	}

	std::cout << "Interpolating level curves..." << std::endl;
	interpolate_curves(levelcurves_nose, jribi::FACEREGION::NOSE, sampled_curve_nose);
	interpolate_curves(levelcurves_lefteye, jribi::FACEREGION::LEFTEYE, sampled_curve_lefteye);
	interpolate_curves(levelcurves_righteye, jribi::FACEREGION::RIGHTEYE, sampled_curve_righteye);
	if (TEST_PRINT)
	{
		print_interpolatedCurve(sampled_curve_nose, sampled_curve_lefteye, sampled_curve_righteye);
	}
	std::cout << "Done." << std::endl;
}

//=================================================== GEODESIC FUNCTIONS
GeodesicAlgorithm FaceMesh::propagate_from_landmark(unsigned int vertex_id, 
											geodesic::Mesh &mesh,
											double propagation_limit)
{
	unsigned const subdivision_level = SUBDIVISION_LEVEL;
	std::vector<geodesic::SurfacePoint> sources;
	GeodesicAlgorithm geomesh(&mesh, subdivision_level);
	sources.push_back(geodesic::SurfacePoint(&mesh.vertices()[vertex_id]));
	geomesh.propagate(sources, propagation_limit);		//cover the whole mesh
	return geomesh;
}

int FaceMesh::closest_vertexID(float x, float y, float z)
{
	double error_point;
	double max_error = std::numeric_limits<double>::max();
	int vertex_id;

	//Iterate through mesh and find closest vertex
	for(unsigned i=0; i<this->mesh.vertices().size(); ++i)
	{
		error_point = abs(this->mesh.vertices()[i].x()-x)
							+ abs(this->mesh.vertices()[i].y()-y)
							+ abs(this->mesh.vertices()[i].z()-z);
		if(error_point < max_error)
		{
			vertex_id = i;
			max_error = error_point;
		}		
	}
	return vertex_id;
}

unsigned int FaceMesh::closest_vertexID_inMesh(unsigned int path_index, std::vector<geodesic::SurfacePoint> path)
{
	unsigned int pointID;					//convert to actual vertex ID
	float err = 10000;
	if(path[path_index].base_element()->type() == geodesic::PointType::VERTEX)
	{
		pointID = path[path_index].base_element()->id(); //id of vertex
	}
	else //if base is edge or face, find closest vertex
	{
		geodesic::base_pointer p = path[path_index].base_element();
		err = 10000;
		for(int k=0; k<p->adjacent_vertices().size(); ++k)
		{
			double error_point = abs(p->adjacent_vertices()[k]->x()-path[path_index].x())
					+ abs(p->adjacent_vertices()[k]->y()-path[path_index].y())
					+ abs(p->adjacent_vertices()[k]->z()-path[path_index].z());
			
			if(error_point < err)
			{
				pointID = p->adjacent_vertices()[k]->id();
				err = error_point;
			}
		}
	}
	return pointID;
}

double FaceMesh::cartesian_dist(geodesic::Vertex* vertexA, geodesic::Vertex* vertexB)
{
	double determinant = std::pow(vertexA->x() - vertexB->x(),2) + 
						std::pow(vertexA->y() - vertexB->y(),2) + 
						std::pow(vertexA->z() - vertexB->z(),2);
	return std::sqrt(determinant);
}

//=================================================== BIN & SORT CURVE FUNCTIONS
std::vector<LevelCurve> FaceMesh::bin_level_curves(std::vector<Geodist> &geodist, 
				std::map<std::string, Landmark> &landmarks,
				const int faceregion)
{
	//Step 1: Calculate angle
	//Use smallest element as midpoint. Smallest and largest element are also used to normalize the sums before binning.
	openstl::Vec3 mid; 
	double min_sum = std::numeric_limits<double>::max();
	double max_sum = 0;
	for(std::vector<Geodist>::iterator it = geodist.begin(); it != geodist.end(); ++it)
	{
		if(it->sum < min_sum) 
		{
			min_sum = it->sum;
			mid = it->xyz;
		}
		if(it->sum > max_sum) max_sum = it->sum;
	}
	int nr_level_curves;
	switch (faceregion)
	{
		case FACEREGION::NOSE:
			this->mid_nose = mid;
			nr_level_curves = NR_LEVEL_CURVES_NOSE;
			break;
		case FACEREGION::LEFTEYE:
			this->mid_lefteye = mid;
			nr_level_curves = NR_LEVEL_CURVES_EYE;
			break;
		case FACEREGION::RIGHTEYE:
			this->mid_righteye = mid;
			nr_level_curves = NR_LEVEL_CURVES_EYE;
			break;
	}

	calc_angle(geodist, mid);

	//Step 2: Assign level curves to data, but keep them in original vector as they are not sorted yet.
	int nr_bins = 10 * nr_level_curves;
	int skip = int(nr_bins/1.5);
	std::vector<double> target_levels = math::linspace(1, nr_bins-skip, nr_level_curves);

	int count[nr_level_curves] = {0};
	for(std::vector<Geodist>::iterator it = geodist.begin(); it != geodist.end(); ++it)
	{
		//loop through all points i Geodist
		double binned_sum = std::floor((it->sum - min_sum) / (max_sum - min_sum) * nr_bins);
		it->level_curve_nr = 0; //stays 0 if not in level curve
		for(std::vector<double>::size_type i = 0; i != target_levels.size(); i++)
		{
			if(binned_sum == int(target_levels[i]))
			{
				it->level_curve_nr = i+1;
				count[i]++;
				break;
			}
		}
	}

	//Step 3: Sort by level curve number and angle
	std::sort(geodist.begin(), geodist.end());

	//Step 4: Bin into new curve arrays. Arrays are sorted by angle.
	std::vector<LevelCurve> new_level_curves(nr_level_curves);
	for(int i=0; i<nr_level_curves; i++)
	{
		new_level_curves[i].reserve(count[i]);
	}
	for(std::vector<Geodist>::iterator it = geodist.begin(); it != geodist.end(); ++it)
	{
		if(it->level_curve_nr > 0)
		{
			new_level_curves[it->level_curve_nr-1].push_back(it->xyz);
		}
	}

	//Step 5: Check if angle sorting was correct and fix if necessary
	for(int i=0; i<nr_level_curves; i++)
	{
		std::size_t curve_size = new_level_curves[i].size();
		this->facenormal_vec = math::cross_product(this->face_sep, this->eye2eye_vec);
		for(int n=curve_size; n>1; n--)
		{
			for(int j=0; j<n-1; j++)
			{
				openstl::Vec3 a = {new_level_curves[i][j].x - mid.x, new_level_curves[i][j].y - mid.y, new_level_curves[i][j].z - mid.z};
				openstl::Vec3 b = {new_level_curves[i][j+1].x - mid.x, new_level_curves[i][j+1].y - mid.y, new_level_curves[i][j+1].z - mid.z};
				if(math::determinant3x3(a, b, facenormal_vec) < 0)
				{
					openstl::Vec3 temp = new_level_curves[i][j];
					new_level_curves[i][j] = new_level_curves[i][j+1];
					new_level_curves[i][j+1] = temp;
				}
			}
		}
	}	
	return new_level_curves;
}

void FaceMesh::bipolar_param(GeodesicAlgorithm &geomesh_landmark1, 
						GeodesicAlgorithm &geomesh_landmark2,
						geodesic::SurfacePoint &target,
						Geodist &geodist_buff,
						std::vector<Geodist> &geodist)
{
	std::vector<geodesic::SurfacePoint> path;
	//copy all data as mesh will be destroyed after level curves are generated
	geodist_buff.xyz.x = target.x();
	geodist_buff.xyz.y = target.y();
	geodist_buff.xyz.z = target.z();
	geomesh_landmark1.trace_back(target, path);
	double length1 = geodesic::length(path);
	path.clear();
	//vertices that are further away than the propagation limit result in zeros and should be skipped
	if (length1 != 0)
	{
		geomesh_landmark2.trace_back(target, path);
		double length2 = geodesic::length(path);
		path.clear();
		if(length2 != 0)
		{
			geodist_buff.sum = length1 + length2 + geodist_buff.sum;
			geodist.push_back(geodist_buff);
		}
	}
}

void FaceMesh::threepolar_param(GeodesicAlgorithm &geomesh_landmark1, 
						GeodesicAlgorithm &geomesh_landmark2,
						GeodesicAlgorithm &geomesh_landmark3,
						geodesic::SurfacePoint &target,
						Geodist &geodist_buff,
						std::vector<Geodist> &geodist)
{
	std::vector<geodesic::SurfacePoint> path;
	geomesh_landmark1.trace_back(target, path);
	geodist_buff.sum = geodesic::length(path); 
	path.clear();

	//vertices that are further away than the propagation limit result in zeros and should be skipped
	if (geodist_buff.sum != 0) bipolar_param(geomesh_landmark2, geomesh_landmark3, target, geodist_buff, geodist);
}

void FaceMesh::calc_angle(std::vector<Geodist> &geodist, openstl::Vec3 mid)
{
	float len_face_sep = std::sqrt(this->face_sep.x*this->face_sep.x + this->face_sep.y*this->face_sep.y + this->face_sep.z*this->face_sep.z);

	for(auto it = geodist.begin(); it != geodist.end(); it++)
	{
		if (it->level_curve_nr != 0)
		{
			float x = it->xyz.x - mid.x;
			float y = it->xyz.y - mid.y;
			float z = it->xyz.z - mid.z;
			float dot = x*this->face_sep.x + y*this->face_sep.y + z*this->face_sep.z;
			float len_a = std::sqrt(x*x + y*y + z*z);
			float angle = std::acos(dot/(len_a*len_face_sep));

			//check if vector a is in the same direction as n
			if (x*this->eye2eye_vec.x + y*this->eye2eye_vec.y + z*this->eye2eye_vec.z > 0) it->angle = angle;
			else it->angle = -angle; 
		}
		else it->angle = 0;
	}
}

//=================================================== INTERPOLATION FUNCTIONS
void FaceMesh::bspline_interpolation(LevelCurve &curve_data, double sampled_curve[][3], int faceregion)
{
	using BSpline = splinelib::sources::splines::BSpline<1, 3>;
	using ParameterSpace = BSpline::ParameterSpace_;
	using VectorSpace = BSpline::VectorSpace_;
	using Coordinates = VectorSpace::Coordinates_; //array of coordinates --> vector space in here?
	using Degrees = ParameterSpace::Degrees_;
	using KnotVectors = ParameterSpace::KnotVectors_;
	using Coordinate = Coordinates::value_type; //standard coordinate?
	using Degree = Degrees::value_type;
	using KnotVector = KnotVectors::value_type::element_type;
	using Knots = KnotVector::Knots_;
	using ScalarCoordinate = Coordinate::value_type;
	using Knot = Knots::value_type;

	int nr_subsamples;
	switch (faceregion)
	{
		case FACEREGION::NOSE:
			nr_subsamples = NR_SUBSAMPLES_NOSE;
			break;
		case FACEREGION::LEFTEYE:
			nr_subsamples = NR_SUBSAMPLES_EYE;
			break;
		case FACEREGION::RIGHTEYE:
			nr_subsamples = NR_SUBSAMPLES_EYE;
			break;
	}

	splinelib::Vector<splinelib::Array<splinelib::Coordinate, 3UL>> coordinates_new;
	coordinates_new.reserve(curve_data.size()+2);
	//add point that lies on the face seperator plane to start and end of curve data for nose
	Coordinate startpoint = curve_startingPoint(curve_data[0], curve_data[curve_data.size()-1], faceregion);
	coordinates_new.push_back(startpoint);		
	for(auto it = curve_data.begin(); it != curve_data.end(); it++)
	{
		Coordinate kCoordinate{ScalarCoordinate{it->x}, ScalarCoordinate{it->y}, ScalarCoordinate{it->z}};
		coordinates_new.push_back(kCoordinate);
	}
	coordinates_new.push_back(startpoint);

	const int degree = 2;
	constexpr Knot const kKnot0_0{0.0}, kKnot1_0{1.0};
	std::vector<Knot> knots_new;
	for(int i = 0; i < degree; i++)
	{
		knots_new.push_back(kKnot0_0);
	}
	int len_data = curve_data.size()+2-(degree-1); //-3 bei degree 5?
	for(int i = 0; i < len_data; i++)
	{
		double expr = double(double(i)/double(len_data-1));
		knots_new.push_back(Knot{expr});
	}
	for(int i = 0; i < degree; i++)
	{
		knots_new.push_back(kKnot1_0);
	}
	constexpr ScalarCoordinate const kCoordinate0_0{0.0}, kCoordinate1_0{1.0};
	constexpr Degree const kDegree{degree};
	splinelib::SharedPointer<KnotVector> const knot_vector{std::make_shared<KnotVector>(knots_new)};
	constexpr Degrees const kDegrees{kDegree};
	KnotVectors const knot_vectors{knot_vector};
	splinelib::SharedPointer<ParameterSpace> const parameter_space{std::make_shared<ParameterSpace>(knot_vectors, kDegrees)}; //knot_vectors and kDegrees here
	splinelib::SharedPointer<VectorSpace> const vector_space{std::make_shared<VectorSpace>(coordinates_new)}; //coordinates should be in here
	BSpline const b_spline{parameter_space, vector_space}; //This seems to be the actual b-spline Calculation

	// Evaluate the planar B-spline curve.
	using ParametricCoordinate = BSpline::ParametricCoordinate_;
	using ScalarParametricCoordinate = ParametricCoordinate::value_type;

	constexpr ScalarParametricCoordinate const kParametricCoordinate0_25{0.25};
	constexpr ParametricCoordinate const kParametricCoordinate{kParametricCoordinate0_25};
	//Coordinate const &evaluated = b_spline(kParametricCoordinate);
	int vertex_id;
	for(int i= 0; i < nr_subsamples; i++)
	{
		Coordinate const &evaluated = b_spline(ParametricCoordinate{ScalarParametricCoordinate{double(i)/double(nr_subsamples-1)}});
		//std::cout << "Evaluated: " << evaluated[0].Get() << " " << evaluated[1].Get() << " " << evaluated[2].Get() << std::endl;
		vertex_id = closest_vertexID(float(evaluated[0].Get()), float(evaluated[1].Get()), float(evaluated[2].Get()));
		if (mesh.vertices()[vertex_id].is_boundary())
		{
			sampled_curve[i][0] = -1.0; //set index to -1 to indicate that vertex is a boundary vertex
			sampled_curve[i][1] = 0.0;
			std::cout << "WARNING: Vertex " << i << " is a boundary vertex. Cannot determine curvature." << std::endl;
		}
		else
		{
			const auto [k1, k2]  = principalCurvature(vertex_id);
			if(vertex_id<0) std::runtime_error("Error: Vertex ID of interpolated vertix is negative.");
			sampled_curve[i][0] = double(vertex_id);
			if(k1 == k2)
				sampled_curve[i][1] = 0.0; //results in illegal division by zero
				sampled_curve[i][2] = 1.0;
			if(k1>k2)
				sampled_curve[i][1] =  double(0.5 - (atan2((k1 + k2) , (k1 - k2))/M_PI));
				sampled_curve[i][2] = 1.0;
			if(k1<k2)
				sampled_curve[i][1] =  double(0.5 - (atan2((k2 + k1) , (k2 - k1))/M_PI));
				sampled_curve[i][2] = 1.0;
			if(std::isnan(sampled_curve[i][1])) 
			{
				sampled_curve[i][1] = 0.0;
				sampled_curve[i][2] = 0.0;
			}
		}
	}
}

splinelib::Array<splinelib::Coordinate, 3UL> FaceMesh::curve_startingPoint(openstl::Vec3 first_curve_data, openstl::Vec3 last_curve_data, int faceregion)
{
	float d;
	switch (faceregion)
	{
		case FACEREGION::NOSE:
			d = math::dot_product(this->eye2eye_vec, this->mid_nose);
			break;
		case FACEREGION::LEFTEYE:
			d = math::dot_product(this->eye2eye_vec, this->mid_lefteye);
			break;
		case FACEREGION::RIGHTEYE:
			d = math::dot_product(this->eye2eye_vec, this->mid_righteye);
			break;
	}
	float k = (std::pow(this->eye2eye_vec.x,2) + std::pow(this->eye2eye_vec.y,2) + std::pow(this->eye2eye_vec.z,2));
  	float t_first = (d - math::dot_product(this->eye2eye_vec, first_curve_data)) / k;
	float t_last = (d - math::dot_product(this->eye2eye_vec, last_curve_data)) / k;
	if(std::abs(t_last) < std::abs(t_first))
	{
		splinelib::Coordinate x{last_curve_data.x + (t_last * this->eye2eye_vec.x)};
		splinelib::Coordinate y{last_curve_data.y + (t_last * this->eye2eye_vec.y)};
		splinelib::Coordinate z{last_curve_data.z + (t_last * this->eye2eye_vec.z)};
		return {x,y,z};
	}
	else
	{
		splinelib::Coordinate x{first_curve_data.x + (t_first * this->eye2eye_vec.x)};
		splinelib::Coordinate y{first_curve_data.y + (t_first * this->eye2eye_vec.y)};
		splinelib::Coordinate z{first_curve_data.z + (t_first * this->eye2eye_vec.z)};
		return {x,y,z};
	}
} 

std::tuple<double, double> FaceMesh::principalCurvature(int vertex_id)
{
	
	//Calculate Gaussian Curvature by iterating through all faces of a vertex to get sum of angles
	double sum_angle = 0;
	double sum_area = 0;
	geodesic::SimpleVector<geodesic::face_pointer> m_faces = mesh.vertices()[vertex_id].adjacent_faces();
	for(unsigned i=0; i<m_faces.size(); ++i)
	{
		geodesic::Face* f = m_faces[i];
		double angle = 0;
		std::vector<geodesic::Vertex*> vertexes;
		for(unsigned j=0; j<3; ++j) //Iterate through all vertices of face
		{
			if (f->adjacent_vertices()[j]->id() == mesh.vertices()[vertex_id].id())
			{ //Save angle, that belongs to the current vertex
				angle = f->corner_angles()[j];
			}
			else
			{ //save both other vertexes to calculate the area and the sign of the mean curvature
				vertexes.push_back(f->adjacent_vertices()[j]);
			}
		}
		double dists = cartesian_dist(vertexes[0], &mesh.vertices()[vertex_id]) * cartesian_dist(vertexes[1], &mesh.vertices()[vertex_id]);
		double tri_area = (1.0/6.0)* double(std::sin(angle)) * dists;                        
		sum_area = sum_area + tri_area;
		sum_angle += angle;
	}
	double gauss_curvature = ((2*M_PI - sum_angle)/sum_area);

	//Determine mean curvature
	geodesic::SimpleVector<geodesic::edge_pointer> m_edges = mesh.vertices()[vertex_id].adjacent_edges();
	double sum_x = 0;
	double sum_y = 0;
	double sum_z = 0;
	for(unsigned i=0; i<m_edges.size(); ++i)
	{
		geodesic::Edge* e = m_edges[i];
		geodesic::Vertex* opposite_vertex = e->opposite_vertex(&mesh.vertices()[vertex_id]);
		std::vector<double> vertex_angles;
		//get angle of vertex that is neither current vertex, nor main vertex
		for(unsigned j=0; j<2; j++)
		{
			geodesic::Face* f = e->adjacent_faces()[j];
			for(unsigned l=0; l<3; l++)
			{
				//get angle of vertex that is neither current vertex, nor main vertex
				if((f->adjacent_vertices()[l]->id() != mesh.vertices()[vertex_id].id()) && (f->adjacent_vertices()[l]->id() != opposite_vertex->id()))
				{
					vertex_angles.push_back(f->corner_angles()[l]);                                
				}
			}
		}
		double factor = (std::cos(vertex_angles[0])/std::sin(vertex_angles[0])) + (std::cos(vertex_angles[1])/std::sin(vertex_angles[1]));
		sum_x += factor * (opposite_vertex->x() - mesh.vertices()[vertex_id].x());
		sum_y += factor * (opposite_vertex->y() - mesh.vertices()[vertex_id].y());
		sum_z += factor * (opposite_vertex->z() - mesh.vertices()[vertex_id].z());
	}
	//Calculate norm of absolute mean curvature
	openstl::Vec3 p = openstl::Vec3{float(sum_x/(2*sum_area)), float(sum_y/(2*sum_area)), float(sum_z/(2*sum_area))};
	double mean_curvature = (std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z)/2.0);
	//Calculate sign of mean curvature by averaging face's normals around a vertex
	if (math::dot_product(this->facenormal_vec, p) < 0) mean_curvature = -mean_curvature;
	double min_curvature = mean_curvature - std::sqrt(mean_curvature*mean_curvature - gauss_curvature);
	double max_curvature = mean_curvature + std::sqrt(mean_curvature*mean_curvature - gauss_curvature);
	return std::make_tuple(min_curvature, max_curvature);	
}      

//=================================================== PRINT FUNCTIONS
void FaceMesh::print_geodesicpath(const char *filename, std::vector<geodesic::SurfacePoint> path)
{
	std::ofstream outputfile;
	outputfile.open(filename);
	outputfile << "x,y,z\n";
	for(int j=0; j < path.size(); j++)
	{
		outputfile << path[j].x() << "," << path[j].y() << "," << path[j].z() << "\n";
	}
	outputfile.close();
}

void FaceMesh::print_multipolar_dist(const char *filename, std::vector<Geodist> &geodist)
{
	std::ofstream outputfile;
	outputfile.open(filename);
	outputfile << "x,y,z,dist, level_curve_nr\n";
	for(int j = 0; j < geodist.size(); j++)
	{
		outputfile 	<< geodist[j].xyz.x << "," 
					<< geodist[j].xyz.y << "," 
					<< geodist[j].xyz.z << ","
					<< geodist[j].sum << ","
					<< geodist[j].level_curve_nr << "\n";
	}
	outputfile.close();
}

void FaceMesh::print_sortedcurve(const char *filename, std::vector<LevelCurve> level_curves, int nr_level_curve)
{
	std::ofstream outputfile;
	outputfile.open(filename);
	outputfile << "x,y,z,levelcurve_nr\n";
	for(int i=0; i<nr_level_curve; i++) 
	{
		for(int j = 0; j < level_curves[i].size(); j++)
		{
			outputfile 	<< level_curves[i][j].x << "," 
						<< level_curves[i][j].y << "," 
						<< level_curves[i][j].z << ","
						<< i << "\n";
		}
	}
	outputfile.close();
}

void FaceMesh::print_interpolatedCurve(
			math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> &sampled_curve_nose,
			math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_lefteye,
			math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_righteye)
{
	for(int k=0; k<3; k++)
	{
		std::ofstream outputfile;
		int samples, nr_curves, vertex_id;
		if(k==FACEREGION::NOSE) 
		{
			outputfile.open("./out/interpolated_curve_nose.csv");
			samples = NR_SUBSAMPLES_NOSE;
			nr_curves = NR_LEVEL_CURVES_NOSE;
		}
		else if(k==FACEREGION::LEFTEYE) 
		{
			outputfile.open("./out/interpolated_curve_lefteye.csv");
			samples = NR_SUBSAMPLES_EYE;
			nr_curves = NR_LEVEL_CURVES_EYE;
		}
		else 
		{
			outputfile.open("./out/interpolated_curve_righteye.csv");
			samples = NR_SUBSAMPLES_EYE;
			nr_curves = NR_LEVEL_CURVES_EYE;
		}
		outputfile << "curve_nr,x,y,z,shape_index\n";
		for (int i = 0; i < nr_curves; i++)
		{
			for(int j = 0; j < samples; j++)
			{
				if(k==FACEREGION::NOSE) vertex_id = sampled_curve_nose[i][j][0];
				else if(k==FACEREGION::LEFTEYE) vertex_id = sampled_curve_lefteye[i][j][0];
				else vertex_id = sampled_curve_righteye[i][j][0];
				if (vertex_id < 0)
				{
					outputfile 	<< i << ","
							<< 0.0 << "," 
							<< 0.0 << "," 
							<< 0.0 << ","
							<< 0.0 << "\n";
				}
				else
				{
					if(k==FACEREGION::NOSE) 
					{
						outputfile 	<< i << ","
							<< mesh.vertices()[vertex_id].x() << "," 
							<< mesh.vertices()[vertex_id].y() << "," 
							<< mesh.vertices()[vertex_id].z() << ","
							<< sampled_curve_nose[i][j][1] << "\n";
					}
					else if(k==FACEREGION::LEFTEYE) 
					{
						outputfile 	<< i << ","
							<< mesh.vertices()[vertex_id].x() << "," 
							<< mesh.vertices()[vertex_id].y() << "," 
							<< mesh.vertices()[vertex_id].z() << ","
							<< sampled_curve_lefteye[i][j][1] << "\n";
					}
					else 
					{
						outputfile 	<< i << ","
							<< mesh.vertices()[vertex_id].x() << "," 
							<< mesh.vertices()[vertex_id].y() << "," 
							<< mesh.vertices()[vertex_id].z() << ","
							<< sampled_curve_righteye[i][j][1] << "\n";
					}
				}
			}
		}
		outputfile.close();
	}
}

void FaceMesh::print_curves(const char *filename, 
			math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> &sampled_curve_nose,
			math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_lefteye,
			math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_righteye)
{
	std::ofstream outputfile;
	outputfile.open(filename);
	outputfile << "Curve_nr,shape_index,flag_valid\n";
	for(int k = 0; k<3; k++)
	{
		int samples = 0;
		int nr_curves = 0;
		int vertex_id = 0;
		if(k==0) 
		{
			outputfile << "nose,,\n";
			samples = jribi::NR_SUBSAMPLES_NOSE;
			nr_curves = jribi::NR_LEVEL_CURVES_NOSE;
		}
		else
		{
			if(k==1) outputfile << "lefteye,,\n";
			else outputfile << "righteye,,\n";
			samples = jribi::NR_SUBSAMPLES_EYE;
			nr_curves = jribi::NR_LEVEL_CURVES_EYE;
		}
		for(int i = 0; i < nr_curves; i++)
		{
			for(int j = 0; j < samples; j++)
			{
				if(k==0) vertex_id = sampled_curve_nose[i][j][0];
				else if(k==1) vertex_id = sampled_curve_lefteye[i][j][0];
				else vertex_id = sampled_curve_righteye[i][j][0];

				if (vertex_id < 0)
				{
					outputfile 	<< double(i) << ","
							<< 0.0 << "," 
							<< 1.0 << "\n";
				}
				else
				{
					if(k==0) outputfile << i << ","
								<< sampled_curve_nose[i][j][1] << "," 
								<< sampled_curve_nose[i][j][2] << "\n";
					else if(k==1) outputfile 	<< i << ","
								<< sampled_curve_lefteye[i][j][1] << "," 
								<< sampled_curve_lefteye[i][j][2] << "\n";
					else outputfile 	<< i << ","
								<< sampled_curve_righteye[i][j][1] << "," 
								<< sampled_curve_righteye[i][j][2] << "\n";
				}
			}
		}
	}
	outputfile << "mean" << "," 
		<< double((math::sum(sampled_curve_nose,1,2) + math::sum(sampled_curve_lefteye,1,2) + math::sum(sampled_curve_righteye,1,2))/
				(jribi::NR_LEVEL_CURVES_NOSE*jribi::NR_SUBSAMPLES_NOSE + 2*jribi::NR_LEVEL_CURVES_EYE*jribi::NR_SUBSAMPLES_EYE)) << ","
		<< 1.0 << "\n"; //Print mean value to reduce calculation time at comparison
	outputfile.close();
}

std::tuple<double, int> read_curves(const char *filename, 
				math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> &sampled_curve_nose,
				math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_lefteye,
				math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> &sampled_curve_righteye)
{
	std::ifstream inputfile(filename);
	if(!inputfile.is_open())
	{
		std::cerr << "Error: Could not open file " << filename << std::endl;
		return std::make_tuple(-1.0, 0);
	}
	std::string curve_nr,shape_index,flag_valid;
	std::string newline;
	std::getline(inputfile, newline);
	if(newline != "Curve_nr,shape_index,flag_valid")
	{
		std::cerr << "Error: Wrong file format. First line should be 'Curve_nr,shape_index,flag_valid', but found: " << newline << std::endl;
		return std::make_tuple(-1.0, 0);
	}
	std::getline(inputfile, newline);
	if(newline != "nose,,")
	{
		std::cerr << "Error: Wrong file format. Second line should be 'nose,nose,nose', but found: " << newline << std::endl;
		return std::make_tuple(-1.0, 0);
	}
	for(int i=0; i<NR_LEVEL_CURVES_NOSE; i++)
	{
		for(int j=0; j<NR_SUBSAMPLES_NOSE; j++)
		{
			std::getline(inputfile, curve_nr, ',');
			std::getline(inputfile, shape_index, ',');
			std::getline(inputfile, flag_valid, '\n');
			sampled_curve_nose[i][j][0] = std::stod(curve_nr);
			sampled_curve_nose[i][j][1] = std::stod(shape_index);
			sampled_curve_nose[i][j][2] = std::stod(flag_valid);
		}
	}
	std::getline(inputfile, newline);
	if(newline != "lefteye,,")
	{
		std::cerr << "Error: Wrong file format. Line should be 'lefteye,,', but found: " << newline << std::endl;
		return std::make_tuple(-1.0, 0);
	}
	for(int i=0; i<NR_LEVEL_CURVES_EYE; i++)
	{
		for(int j=0; j<NR_SUBSAMPLES_EYE; j++)
		{
			std::getline(inputfile, curve_nr, ',');
			std::getline(inputfile, shape_index, ',');
			std::getline(inputfile, flag_valid, '\n');
			sampled_curve_lefteye[i][j][0] = std::stod(curve_nr);
			sampled_curve_lefteye[i][j][1] = std::stod(shape_index);
			sampled_curve_lefteye[i][j][2] = std::stod(flag_valid);
		}
	}
	std::getline(inputfile, newline);
	if(newline != "righteye,,")
	{
		std::cerr << "Error: Wrong file format. Line should be 'righteye,,', but found: " << newline << std::endl;
		return std::make_tuple(-1.0, 0);
	}
	for(int i=0; i<NR_LEVEL_CURVES_EYE; i++)
	{
		for(int j=0; j<NR_SUBSAMPLES_EYE; j++)
		{
			std::getline(inputfile, curve_nr, ',');
			std::getline(inputfile, shape_index, ',');
			std::getline(inputfile, flag_valid, '\n');
			sampled_curve_righteye[i][j][0] = std::stod(curve_nr);
			sampled_curve_righteye[i][j][1] = std::stod(shape_index);
			sampled_curve_righteye[i][j][2] = std::stod(flag_valid);
		}
	}	
	std::getline(inputfile, newline, ',');
	if(newline != "mean")
	{
		std::cerr << "Error: Wrong file format. Line should be 'mean', but found: " << newline << std::endl;
		return std::make_tuple(-1.0, 0);
	}
	std::getline(inputfile, shape_index, ',');
	std::getline(inputfile, flag_valid, '\n');
	return std::make_tuple(std::stof(shape_index), 1);
}
} //namespace jribi

