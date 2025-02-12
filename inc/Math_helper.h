/*
 * SPDX-FileCopyrightText: 2025 Alexandra Mielke <alexandra.mielke@smail.emt.h-brs.de>
 *
 * SPDX-License-Identifier: MIT
 */
/*
 * SPDX-FileCopyrightText: 2017 Akavil
 *
 * SPDX-License-Identifier: CC-BY-SA-3.0
 */

#ifndef _MATH_HELPER_H_
#define _MATH_HELPER_H_

#include <vector>
#include <math.h>

namespace math
{
template <typename T, std::size_t Row, std::size_t Col, std::size_t Dep>
using Array3d = std::array<std::array<std::array<T, Dep>, Col>, Row>;

template <typename T, std::size_t Row, std::size_t Col>
using Array2d = std::array<std::array<T, Col>, Row>;

/**
 * @author Akavil 2014, edited 2017 CC BY-SA https://stackoverflow.com/questions/27028226/python-linspace-in-c
 * @brief Returns linspace between two values.
 * 
 * @tparam T Either float or int derivative.
 * @param start_in Start value of linspace.
 * @param end_in End value of linspace.
 * @param num_in Number of elements in linspace.
 * @return std::vector<double> linspaced elements.
 */
template<typename T>
inline std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // Ensure that start and end are exactly the same as the input
  return linspaced;
}

/**
 * @author Alexandra Mielke
 * @brief Returns determinant of x,y and z vector. If a,b,c are sorted in a right-hand system, the determinant is positive.
 * 
 * @tparam T Vector struct with x,y,z members accessed via . operator.
 * @param a Vector a.
 * @param b Vector b.
 * @param c Vector c. 
 * @return float 
 */
template<typename T>
inline float determinant3x3(T a, T b, T c)
{
  return a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) + a.z * (b.x * c.y - b.y * c.x);
}

/**
 * @author Alexandra Mielke
 * @brief Returns cross product of two vectors, e.g. normal vector to them.
 * 
 * @tparam T Vector struct with x,y,z members accessed via . operator. 
 * @param a Vector a.
 * @param b Vector b.
 * @return T Vector normal to a and b.
 */
template<typename T>
inline T cross_product(T a, T b)
{
  T c;
  c.x = a.y*b.z - a.z*b.y;
  c.y = a.z*b.x - a.x*b.z;
  c.z = a.x*b.y - a.y*b.x;
  return c;
}

/**
 * @author Alexandra Mielke
 * @brief Returns dot product of two vectors.
 * 
 * @tparam T Vector struct with x,y,z members accessed via . operator.
 * @param a Vector a.
 * @param b Vector b.
 * @return float Dot product, e.g. scalar value.
 */
template<typename T>
inline float dot_product(T &a, T &b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

/**
 * @author Alexandra Mielke
 * @brief Determines the scalar product of a vector and a scalar.
 * 
 * @tparam T Vector struct with x,y,z members accessed via . operator.
 * @tparam U Type of scalar. Should be float or int.
 * @param a 3D Vector.
 * @param b Scalar.
 * @return T Vector multiplied by scalar b.
 */
template<typename T, typename U>
inline T scalar_product(T &a, U b)
{
  T c;
  c.x = a.x*b;
  c.y = a.y*b;
  c.z = a.z*b;
  return c;
}

/**
 * @author Alexandra Mielke
 * @brief Returns the signed angle between two vectors. Sign is determined via the normal vector.
 * 
 * @tparam T Vector struct with x, y, z members accessed via . operator.
 * @param a Vector a. 
 * @param b Vector b.
 * @param n Normal vector.
 * @return float 
 */
template<typename T>
inline float signed_angle(T &a, T &b, T &n)
{
	float dot = a.x*b.x + a.y*b.y + a.z*b.z;
	float len_a = std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
	float len_b = std::sqrt(b.x*b.x + b.y*b.y + b.z*b.z);
	float angle = std::acos(dot/(len_a*len_b));

	//check if vector a is in the same direction as n
	if (a.x*n.x + a.y*n.y + a.z*n.z > 0) return angle;
	else return -angle; 
}

/**
 * @author Alexandra Mielke
 * @brief Determine angle between two vectors.
 * 
 * @tparam T Vector struct with x, y, z members accessed via . operator.
 * @param a Vector a.
 * @param b Vector b. 
 * @return float 
 */
template<typename T>
inline float standard_angle(T &a, T &b)
{
	float dot = a.x*b.x + a.y*b.y + a.z*b.z;
	float len_a = std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
	float len_b = std::sqrt(b.x*b.x + b.y*b.y + b.z*b.z);
	return std::acos(dot/(len_a*len_b));
}

template<typename T>
inline float cartesian_dist(T &a, T &b)
{
  double determinant = std::pow(a.x - b.x,2) + 
                            std::pow(a.y - b.y,2) + 
                            std::pow(a.z - b.z,2);
        return std::sqrt(determinant);
}

/**
 * @author Alexandra Mielke
 * @brief Calculates the mean of a 3D array.
 * 
 * @tparam T type of array. Should be float, double or int.
 * @tparam ROW Number of rows in 3D array.
 * @tparam COL Number of columns in 3D array.
 * @tparam DEP Depth of 3D array.
 * @param arr array.
 * @return double Mean value for the slice indexed by depth_index.
 */
template<typename T, std::size_t ROW, std::size_t COL, std::size_t DEP>
double mean(Array3d<T, ROW, COL, DEP> &arr)
{
  double sum = 0;
  for(int i=0; i<ROW; i++)
  {
    for(int j=0; j<COL; j++)
    {
      for(int k=0; k<DEP; k++)
      {
        sum += double(arr[i][j][k]);
      }
    }
  }
  return sum/double(ROW*COL*DEP);
}

/**
 * @author Alexandra Mielke
 * @brief Calculates the sum of all members of a sliced 3D array at a specific depth index.
 * 
 * @tparam T type of array. Should be float, double or int.
 * @tparam ROW Number of rows in 3D array.
 * @tparam COL Number of columns in 3D array.
 * @tparam DEP Depth of 3D array. This is the dimension that is sliced.
 * @param arr array.
 * @param depth_index Depth index that is used to calculate the mean.
 * @param flag_index Flag index that is used to indicate if value at depth index is legit.
 * @return double Mean value for the slice indexed by depth_index.
 */
template<typename T, std::size_t ROW, std::size_t COL, std::size_t DEP>
double sum(Array3d<T, ROW, COL, DEP> &arr, int depth_index, int flag_index)
{
  double sum = 0;
  for(int i=0; i<ROW; i++)
  {
    for(int j=0; j<COL; j++)
    {
      if (int(arr[i][j][flag_index]) == 1) sum += double(arr[i][j][depth_index]);
    }
  }
  return sum;
}

} //namespace math

#endif // _MATH_HELPER_H_