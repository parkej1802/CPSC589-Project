#pragma once

//------------------------------------------------------------------------------
// This file contains a single-function namespace for loading .obj files
// into your program via the "tinyobjloader" library.
//------------------------------------------------------------------------------

#include <string>

#include "Geometry.h"

namespace GeomLoaderForOBJ {
	CPU_Geometry loadIntoCPUGeometry(std::string filename);
};
