#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "GeomLoaderForOBJ.h"

#include "Log.h";

#include "tiny_obj_loader.h"

// Most of this function is just boilerplate from tinyobjloader's GitHub README.
CPU_Geometry GeomLoaderForOBJ::loadIntoCPUGeometry(std::string filename) {
	CPU_Geometry geom;

	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string warn, err;

	if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename.c_str())) {
		Log::error("Error loading OBJ file {}: {}", filename, err);
		Log::warn("Warnings for OBJ file {}: {}", filename, warn);
		throw std::runtime_error(warn + "\n" + err);
	}

	// Loop over "shapes" (described further down).
	for (const auto& shape : shapes) {
		size_t vertIndexOffset = 0;

		// Each shape has multiple faces.
		// These should be triangulated by default by the loader, but we might as
		// well still read in the size from the list. You may want to change this
		// to check if the number is 3 and throw an exception if not. Up to you.
		for (size_t face = 0; face < shape.mesh.num_face_vertices.size(); face++) {
			size_t numFaceVerts = size_t(shape.mesh.num_face_vertices[face]);

			// Loop over vertices in the face. Again, unless there's a bug with
			// the obj loader library or something, this *should* be 3 verts.
			for (size_t v = 0; v < numFaceVerts; v++) {
				// Get the vertex by its index. Get the x, y, & z components.
				tinyobj::index_t idx = shape.mesh.indices[vertIndexOffset + v];
				tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
				tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
				tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];

				geom.verts.push_back(glm::vec3(vx, vy, vz));

				// You could have a .obj file without normal data.
				// In this case, `normal_index` would be negative.
				if (idx.normal_index >= 0) {
					// If there is normal data, we add it to our CPU_Geometry object.
					tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
					tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
					tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
					geom.normals.push_back(glm::vec3(nx, ny, nz));
				}
				else {
					// The skeleton is set up so that normals are required.
					// So right now, we throw an exception if the .obj file
					// does not have any normals.
					// You could change this if you want. E.g., you could have
					// a shader not relying on normals, or you could
					// auto-generate the normals here.
					Log::error("Missing normals in OBJ file {}!", filename);
					throw std::runtime_error("No normal data supplied!");
				}

				// The skeleton can tolerate missing texture coordinates.
				// So we simply check if they exist and add them if so.
				if (idx.texcoord_index >= 0) {
					tinyobj::real_t tx = attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
					tinyobj::real_t ty = attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
					geom.uvs.push_back(glm::vec2(tx, ty));
				}

			}
			vertIndexOffset += numFaceVerts;

			// If your .obj comes with a material, here is where you'd access
			// it and do something with it. To see how to use this
			// tinyobj::material_t, you can see the source code at:
			// https://github.com/tinyobjloader/tinyobjloader/blob/v2.0.0rc10/tiny_obj_loader.h#L181
			// Currently, the skeleton does not use them, but you may change that.
			//
			// shape.mesh.material_ids[face]; // Materials are per-face. 
		}
	}
	return geom;
}
