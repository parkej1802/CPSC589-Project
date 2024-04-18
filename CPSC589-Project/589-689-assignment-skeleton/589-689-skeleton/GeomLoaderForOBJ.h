#pragma once

//------------------------------------------------------------------------------
// This file contains a single-function namespace for loading .obj files
// into your program via the "tinyobjloader" library.
//------------------------------------------------------------------------------

#include <string>

#include "Geometry3D.h"

namespace GeomLoaderForOBJ {
	//CPU_Geometry loadIntoCPUGeometry(std::string filename);
	CPU_Geometry3D loadIntoCPUGeometry(std::string filename);
};
