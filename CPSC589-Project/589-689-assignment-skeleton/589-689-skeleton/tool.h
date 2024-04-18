#pragma once

float findFrontMaxY(const std::vector<glm::vec3>& verts) {
	float maxY = std::numeric_limits<float>::lowest();

	for (const auto& vert : verts) {
		if (vert.y > maxY) {
			maxY = vert.y;
		}
	}

	return maxY;
}

float findFrontMaxX(const std::vector<glm::vec3>& verts) {
	float maxX = std::numeric_limits<float>::lowest();

	for (const auto& vert : verts) {
		if (vert.x > maxX) {
			maxX = vert.x;
		}
	}

	return maxX;
}

float findFrontMinY(const std::vector<glm::vec3>& verts) {
	float minY = std::numeric_limits<float>::max();

	for (const auto& vert : verts) {
		if (vert.y < minY) {
			minY = vert.y;
		}
	}

	return minY;
}

float findFrontMinX(const std::vector<glm::vec3>& verts) {
	float minX = std::numeric_limits<float>::max();

	for (const auto& vert : verts) {
		if (vert.x < minX) {
			minX = vert.x;
		}
	}

	return minX;
}

float findSideMaxY(const std::vector<glm::vec3>& verts) {
	float maxY = std::numeric_limits<float>::lowest();

	for (const auto& vert : verts) {
		if (vert.y > maxY) {
			maxY = vert.y;
		}
	}

	return maxY;
}

float findSideMaxZ(const std::vector<glm::vec3>& verts) {
	float maxZ = std::numeric_limits<float>::lowest();

	for (const auto& vert : verts) {
		if (vert.z > maxZ) {
			maxZ = vert.z;
		}
	}

	return maxZ;
}
