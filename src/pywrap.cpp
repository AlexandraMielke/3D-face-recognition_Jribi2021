/*
 * SPDX-FileCopyrightText: 2025 Alexandra Mielke <alexandra.mielke@smail.emt.h-brs.de>
 *
 * SPDX-License-Identifier: MIT
 */
/*
 * SPDX-FileCopyrightText: 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

/*PyBind11:

Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>, All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Please also refer to the file .github/CONTRIBUTING.md, which clarifies licensing of
external contributions to this project including patches, pull requests, etc.*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "../inc/2021Jribi.h"

namespace py = pybind11;

PYBIND11_MODULE(JribiDescription, m) 
{
    m.doc() = "Jribi module for invariant description of 3D face"; // optional module docstring
    m.def("jribi_description", 
    [](std::string filename, std::string filenameout, const py::dict landmarks_dict) -> int
    {
        // Convert the Python dictionary to a C++ map
        std::map<std::string, jribi::Landmark> landmarks;
        for (auto item : landmarks_dict)
        {
            auto array = item.second.cast<py::array_t<float>>();
            auto buf = array.request();
            float* ptr = (float*)buf.ptr;
            jribi::Landmark landmark_buff;
            landmark_buff.x = ptr[0];
            landmark_buff.y = ptr[1];
            landmark_buff.z = ptr[2];
            landmark_buff.vertex_id = 0;
            landmarks[std::string(py::str(item.first))] = landmark_buff;
        }
        //Copy to landmarks
        std::cout << "Reading STL file" << std::endl;
        const auto [vertices, faces] = import_mesh(filename.c_str());
        std::cout << "Done." << std::endl;

        //parse json file
        std::string filename_stl = filename + ".JSON";

        // Open the file
        std::ifstream file(filename);
        if (!file) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return 1;
        }

        // Read the entire file into a string
        std::stringstream buffer;
        buffer << file.rdbuf();
        std::string jsonString = buffer.str();

        //create level curves
        math::Array3d<double, jribi::NR_LEVEL_CURVES_NOSE, jribi::NR_SUBSAMPLES_NOSE, 3> sampled_curve_nose;
        math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> sampled_curve_lefteye;
        math::Array3d<double, jribi::NR_LEVEL_CURVES_EYE, jribi::NR_SUBSAMPLES_EYE, 3> sampled_curve_righteye;	

        jribi::FaceMesh* facemesh = new jribi::FaceMesh(vertices, faces);
        facemesh->createLevelCurves(landmarks, sampled_curve_nose, 
                                            sampled_curve_lefteye, 
                                            sampled_curve_righteye);	
        facemesh->print_curves(filenameout.c_str(), sampled_curve_nose, sampled_curve_lefteye, sampled_curve_righteye);
        return 0;
    }, 
        "A function that calculates geodesic distances for a given mesh and landmarks and prints the invariant description of the face.");
}