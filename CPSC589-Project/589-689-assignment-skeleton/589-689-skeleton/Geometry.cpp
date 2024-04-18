

#include "Geometry.h"

#include <utility>


/*
GPU_Geometry::GPU_Geometry()
	: vao()
	, vertBuffer(0, 3, GL_FLOAT)

	, uvBuffer(1, 2, GL_FLOAT)
	, normalsBuffer(2, 3, GL_FLOAT)
	, colBuffer(3, 3, GL_FLOAT)



{}
*/
GPU_Geometry::GPU_Geometry()
	: vao()
	, vertBuffer(0, 3, GL_FLOAT)
	, colBuffer(1, 3, GL_FLOAT)

{}


void GPU_Geometry::setVerts(const std::vector<glm::vec3>& verts) {
	vertBuffer.uploadData(sizeof(glm::vec3) * verts.size(), verts.data(), GL_STATIC_DRAW);
}


void GPU_Geometry::setCols(const std::vector<glm::vec3>& cols) {
	colBuffer.uploadData(sizeof(glm::vec3) * cols.size(), cols.data(), GL_STATIC_DRAW);
}



/*
#include "Geometry.h"

#include <utility>


GPU_Geometry::GPU_Geometry()
	: vao()
	, vertBuffer(0, 3, GL_FLOAT)
	, colBuffer(1, 3, GL_FLOAT)
{}


void GPU_Geometry::setVerts(const std::vector<std::vector<glm::vec3>>& verts) {
	for (const auto& v : verts) {
		size_t vSize = v.size() * sizeof(glm::vec3);
		vertBuffer.uploadData(vSize, v.data(), GL_DYNAMIC_DRAW);
	}
}


void GPU_Geometry::setCols(const std::vector<std::vector<glm::vec3>>& cols) {
	for (const auto& c : cols) {
		size_t cSize = c.size() * sizeof(glm::vec3);
		colBuffer.uploadData(cSize, c.data(), GL_DYNAMIC_DRAW);
	}
}
*/
