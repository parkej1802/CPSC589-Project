#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>
#include <unordered_map>
#include <map>
#include <fstream>
#include <set>
#include <random>


// Window.h `#include`s ImGui, GLFW, and glad in correct order.
#include "Window.h"

#include "CDT.h"

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Camera.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "tool.h"
#include "CDTUtils.h"


float PI = 3.14159265359;
bool clear = false;
bool drew = false;
bool showDraw = false;

CPU_Geometry grid;
glm::vec3 red = glm::vec3(1, 0, 0);
glm::vec3 green = glm::vec3(0, 1, 0);
glm::vec3 blue = glm::vec3(0, 0, 1);
glm::vec3 white = glm::vec3(1, 1, 1);

//***************************************************************************|
//							STRUCTS
//***************************************************************************|

struct CustomPoint2D
{
	double data[2];
};

struct CustomEdge
{
	std::pair<std::size_t, std::size_t> vertices;
};

struct Mesh;
struct HalfEdge;
struct Vertex;
struct Face;
struct Edge;

struct Mesh
{
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<std::vector<int>> triangles;
};

struct HalfEdgeMesh
{
	std::list<HalfEdge> halfEdges;
	std::list<Vertex> vertices;
	std::list<Face> faces;
	std::list<Edge> edges;
};

struct HalfEdge {
	Vertex* vertex;
	Face* face;
	HalfEdge* next;
	HalfEdge* pair;

	Edge* edge;
};

struct Vertex {
	HalfEdge* halfEdge;
	glm::vec3 position;
	glm::vec3 normal;

	bool isNew;
	glm::vec3 newPosition;
};

struct Face {
	HalfEdge* halfEdge;
};

struct Edge {
	HalfEdge* halfEdge;
	bool isNew;
	glm::vec3 newPosition;
};

struct EdgeKeyHash {
	std::size_t operator()(const std::pair<int, int>& key) const {
		return std::hash<int>()(key.first) ^ (std::hash<int>()(key.second) << 1);
	}
};

// CALLBACKS
class MyCallbacks : public CallbackInterface {

public:
	// Constructor. We use values of -1 for attributes that, at the start of
	// the program, have no meaningful/"true" value.
	MyCallbacks(ShaderProgram& shader, int screenWidth, int screenHeight)
		: shader(shader)
		, currentFrame(0)
		, leftMouseActiveVal(false)
		, lastLeftPressedFrame(-1)
		, lastRightPressedFrame(-1)
		, screenMouseX(-1.0)
		, screenMouseY(-1.0)
		, screenWidth(screenWidth)
		, screenHeight(screenHeight)
		, x_angle(0)
		, y_angle(0)
		, camera(glm::radians(45.f), glm::radians(45.f), 3.0)
		, rightMouseDown(false)
		, aspect(1.0f)

	{
		updateUniformLocations();
	}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (key == GLFW_KEY_R && action == GLFW_PRESS) {
			x_angle = 0.0f;
			y_angle = 0.0f;
			shader.recompile();
			updateUniformLocations();
		}
	}

	virtual void mouseButtonCallback(int button, int action, int mods) {
		// If we click the mouse on the ImGui window, we don't want to log that
		// here. But if we RELEASE the mouse over the window, we do want to
		// know that!

		auto& io = ImGui::GetIO();
		if (io.WantCaptureMouse && action == GLFW_PRESS) return;


		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
			leftMouseActiveVal = true;
			lastLeftPressedFrame = currentFrame;
		}

		if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
			lastRightPressedFrame = currentFrame;
		}


		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
			leftMouseActiveVal = false;
		}


		if (button == GLFW_MOUSE_BUTTON_RIGHT) {
			if (action == GLFW_PRESS)			rightMouseDown = true;
			else if (action == GLFW_RELEASE)	rightMouseDown = false;
		}
	}

	// Updates the screen width and height, in screen coordinates
	// (not necessarily the same as pixels)
	virtual void windowSizeCallback(int width, int height) {
		screenWidth = width;
		screenHeight = height;
		aspect = float(width) / float(height);
	}

	// Sets the new cursor position, in screen coordinates
	virtual void cursorPosCallback(double xpos, double ypos) {
		if (rightMouseDown) {
			float xoffset = xpos - screenMouseX;
			float yoffset = ypos - screenMouseY;

			const float sensitivity = 0.003f;

			x_angle += xoffset * sensitivity;
			y_angle += yoffset * sensitivity;

			screenMouseX = xpos;
			screenMouseY = ypos;
		}
		else {
			screenMouseX = xpos;
			screenMouseY = ypos;
		}
	}

	virtual void scrollCallback(double xoffset, double yoffset) {
		camera.incrementR(yoffset);
	}

	void viewPipeline() {
		glm::mat4 M = glm::mat4(1.0);
		glm::mat4 V = camera.getView();
		glm::mat4 P = glm::perspective(glm::radians(45.0f), aspect, 0.01f, 1000.f);
		glUniformMatrix4fv(mLoc, 1, GL_FALSE, glm::value_ptr(M));
		glUniformMatrix4fv(vLoc, 1, GL_FALSE, glm::value_ptr(V));
		glUniformMatrix4fv(pLoc, 1, GL_FALSE, glm::value_ptr(P));
	}

	// Whether the left mouse was pressed down this frame.
	bool leftMouseJustPressed() {
		return lastLeftPressedFrame == currentFrame;
	}

	// Whether the left mouse button is being pressed down at all.
	bool leftMouseActive() {
		return leftMouseActiveVal;
	}

	// Whether the right mouse button was pressed down this frame.
	bool rightMouseJustPressed() {
		return lastRightPressedFrame == currentFrame;
	}

	// Tell the callbacks object a new frame has begun.
	void incrementFrameCount() {
		currentFrame++;
	}

	float getXAngle() {
		return x_angle;
	}

	float getYAngle() {
		return y_angle;
	}

	// Converts the cursor position from screen coordinates to GL coordinates
	// and returns the result.
	glm::vec2 getCursorPosGL() {
		glm::vec2 screenPos(screenMouseX, screenMouseY);
		// Interpret click as at centre of pixel.
		glm::vec2 centredPos = screenPos + glm::vec2(0.5f, 0.5f);
		// Scale cursor position to [0, 1] range.
		glm::vec2 scaledToZeroOne = centredPos / glm::vec2(screenWidth, screenHeight);

		glm::vec2 flippedY = glm::vec2(scaledToZeroOne.x, 1.0f - scaledToZeroOne.y);

		// Go from [0, 1] range to [-1, 1] range.
		return 2.f * flippedY - glm::vec2(1.f, 1.f);
	}

	// Takes in a list of points, given in GL's coordinate system,
	// and a threshold (in screen coordinates) 
	// and then returns the index of the first point within that distance from
	// the cursor.
	// Returns -1 if no such point is found.
	int indexOfPointAtCursorPos(const std::vector<glm::vec3>& glCoordsOfPointsToSearch, float screenCoordThreshold) {
		// First, we conver thte points from GL to screen coordinates.
		std::vector<glm::vec3> screenCoordVerts;
		for (const auto& v : glCoordsOfPointsToSearch) {
			screenCoordVerts.push_back(glm::vec3(glPosToScreenCoords(v), 0.f));
		}

		// We make sure we interpret the cursor position as at the centre of
		// the relevant pixel, for consistency with getCursorPosGL().
		glm::vec3 cursorPosScreen(screenMouseX + 0.5f, screenMouseY + 0.5f, 0.f);


		for (size_t i = 0; i < screenCoordVerts.size(); i++) {
			// Return i if length of difference vector within threshold.
			glm::vec3 diff = screenCoordVerts[i] - cursorPosScreen;
			if (glm::length(diff) < screenCoordThreshold) {
				return i;
			}
		}
		return -1; // No point within threshold found.
	}

	Camera camera;

private:

	void updateUniformLocations() {
		mLoc = glGetUniformLocation(shader, "M");
		vLoc = glGetUniformLocation(shader, "V");
		pLoc = glGetUniformLocation(shader, "P");;
		/*
		lightPosLoc = glGetUniformLocation(shader, "lightPos");;
		lightColLoc = glGetUniformLocation(shader, "lightCol");;
		diffuseColLoc = glGetUniformLocation(shader, "diffuseCol");;
		ambientStrengthLoc = glGetUniformLocation(shader, "ambientStrength");;
		texExistenceLoc = glGetUniformLocation(shader, "texExistence");;
		*/
	}

	int screenWidth;
	int screenHeight;

	double screenMouseX;
	double screenMouseY;

	int currentFrame;

	bool leftMouseActiveVal;

	int lastLeftPressedFrame;
	int lastRightPressedFrame;

	float x_angle;
	float y_angle;

	ShaderProgram& shader;

	// Converts GL coordinates to screen coordinates.
	glm::vec2 glPosToScreenCoords(glm::vec2 glPos) {
		// Convert the [-1, 1] range to [0, 1]
		glm::vec2 scaledZeroOne = 0.5f * (glPos + glm::vec2(1.f, 1.f));

		glm::vec2 flippedY = glm::vec2(scaledZeroOne.x, 1.0f - scaledZeroOne.y);
		glm::vec2 screenPos = flippedY * glm::vec2(screenWidth, screenHeight);
		return screenPos;
	}

	bool rightMouseDown;

	float aspect;

	GLint mLoc;
	GLint vLoc;
	GLint pLoc;
};

bool operator==(const glm::vec3& a, const glm::vec3& b) {
	return a.x == b.x && a.y == b.y && a.z == b.z;

}

// get control points of a line user provided
float minDist = 0.1;
std::vector<glm::vec3> get_control_points(std::vector<glm::vec3>& line, int smoothness) {
	std::vector<glm::vec3> result;
	int step_size = line.size() / smoothness;
	int index = 0;

	if (!line.empty()) {
		result.push_back(line[0]);
	}

	for (int index = step_size; index < line.size(); index += step_size) {
		const glm::vec3& lastPoint = result.back();
		const glm::vec3& newPoint = line[index];

		if (glm::distance(lastPoint, newPoint) >= minDist) {
			result.push_back(newPoint);
		}
	}

	return result;
}


// flatten all line vectors to one vector
std::vector<glm::vec3> flattenLineVerts(std::vector<std::vector<glm::vec3>>& lineVerts) {
	std::vector<glm::vec3> flattenedVerts;

	for (const auto& vertGroup : lineVerts) {
		flattenedVerts.insert(flattenedVerts.end(), vertGroup.begin(), vertGroup.end());
	}

	return flattenedVerts;
}

// let user provide cross sections
void draw_cross_sections(
	std::shared_ptr<MyCallbacks>& cb,
	std::vector<std::vector<glm::vec3>>& lineVerts,
	CPU_Geometry& cpuGeom,
	int& cross_section) {

	if (cb->leftMouseActive()) {
		lineVerts[cross_section].push_back(glm::vec3(cb->getCursorPosGL(), 0.f));
		cpuGeom.verts = flattenLineVerts(lineVerts);

		cpuGeom.cols.clear();
		for (int i = 0; i < lineVerts.size(); i++) {
			for (int j = 0; j < lineVerts[i].size(); j++) {
				if (i == 0) {
					cpuGeom.cols.push_back(red);
				}
				else if (i == 1) {
					cpuGeom.cols.push_back(green);
				}
				else {
					cpuGeom.cols.push_back(blue);
				}
			}
		}
	}
}


void transform(
	std::vector<std::vector<glm::vec3>>& lineVerts,
	std::vector<std::vector<glm::vec3>>& transformedVerts) {

	//transformedVerts.clear();
	transformedVerts[0] = lineVerts[0];

	for (int i = 1; i < 3; i++) {
		for (int j = 0; j < lineVerts[i].size(); j++) {
			float x = lineVerts[i][j].x;
			float y = lineVerts[i][j].y;
			if (i == 1) {
				glm::vec3 side = glm::vec3(0.f, y, x);
				transformedVerts[i].push_back(side);
			}
			else {
				glm::vec3 top = glm::vec3(x, 0.f, y);
				transformedVerts[i].push_back(top);
			}
		}
	}
}

// show the user cross sections in 3D
void combine(
	std::shared_ptr<MyCallbacks>& cb,
	std::vector<std::vector<glm::vec3>>& lineVerts,
	CPU_Geometry& cpuGeom) {

	cpuGeom.verts.clear();
	cpuGeom.cols.clear();

	float x_angle = cb->getXAngle();
	glm::mat3 X_R = glm::mat3(cos(x_angle), 0, sin(x_angle),
		0, 1, 0,
		-sin(x_angle), 0, cos(x_angle));

	float y_angle = cb->getYAngle();
	glm::mat3 Y_R = glm::mat3(1, 0, 0,
		0, cos(y_angle), -sin(y_angle),
		0, sin(y_angle), cos(y_angle));

	cpuGeom.verts = flattenLineVerts(lineVerts);
	for (auto& vert : cpuGeom.verts) {
		vert = vert * X_R;
		vert = vert * Y_R;
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < lineVerts[i].size(); j++) {
			if (i == 0) {
				cpuGeom.cols.push_back(red);
			}
			else if (i == 1) {
				cpuGeom.cols.push_back(green);
			}
			else {
				cpuGeom.cols.push_back(blue);
			}
		}
	}
	return;
}

std::vector<glm::vec3> calculateVertexNormals(const Mesh& mesh) {
	std::vector<glm::vec3> normals(mesh.vertices.size(), glm::vec3(0.0f, 0.0f, 0.0f));

	for (const auto& triangle : mesh.triangles) {
		auto v1_index = triangle[0] - 1;
		auto v2_index = triangle[1] - 1;
		auto v3_index = triangle[2] - 1;

		const glm::vec3& v1 = mesh.vertices[v1_index];
		const glm::vec3& v2 = mesh.vertices[v2_index];
		const glm::vec3& v3 = mesh.vertices[v3_index];

		glm::vec3 normal = glm::normalize(glm::cross(v2 - v1, v3 - v1));

		normals[v1_index] += normal;
		normals[v2_index] += normal;
		normals[v3_index] += normal;
	}

	for (auto& normal : normals) {
		normal = glm::normalize(normal);
	}

	return normals;
}

void saveMeshToOBJ(const Mesh& mesh, const std::string& filename) {
	std::ofstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Failed to open file for writing: " << filename << std::endl;
		return;
	}

	for (const auto& vertex : mesh.vertices) {
		file << "v " << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
	}

	for (const auto& normal : mesh.normals) {
		file << "vn " << normal.x << " " << normal.y << " " << normal.z << std::endl;
	}

	file << "s 1" << std::endl;

	for (const auto& triangle : mesh.triangles) {
		file << "f";
		for (int i = 0; i < 3; ++i) {
			int vertexIndex = triangle[i];
			file << " " << vertexIndex << "//" << vertexIndex;
		}
		file << std::endl;
	}

	file.close();
}

Mesh get_front_mesh(const CDT::Triangulation<double>& cdt) {
	Mesh mesh;

	for (const auto& vertex : cdt.vertices) {
		glm::vec3 temp = glm::vec3(vertex.x, vertex.y, 0.0f);
		mesh.vertices.push_back(temp);
	}

	for (const auto& vertex : cdt.vertices) {
		glm::vec3 temp = glm::vec3(0.0, 0.0, 1.0f);
		mesh.normals.push_back(temp);
	}

	for (const auto& triangle : cdt.triangles) {
		std::vector<int> temp;
		for (int i = 0; i < 3; ++i) {
			int vertexIndex = triangle.vertices[i] + 1;
			temp.push_back(vertexIndex);
		}
		mesh.triangles.push_back(temp);
	}

	return mesh;
}

Mesh get_back_mesh(const CDT::Triangulation<double>& cdt) {
	Mesh mesh;

	for (const auto& vertex : cdt.vertices) {
		glm::vec3 temp = glm::vec3(vertex.x, vertex.y, 0.0f);
		mesh.vertices.push_back(temp);
	}

	for (const auto& vertex : cdt.vertices) {
		glm::vec3 temp = glm::vec3(0.0, 0.0, 1.0f);
		mesh.normals.push_back(temp);
	}

	for (const auto& triangle : cdt.triangles) {
		std::vector<int> temp;
		for (int i = 2; i >= 0; --i) {
			int vertexIndex = triangle.vertices[i] + 1;
			temp.push_back(vertexIndex);
		}
		mesh.triangles.push_back(temp);
	}

	return mesh;
}


void front_inflation_side(Mesh& mesh, const std::vector<glm::vec3>& controlPoints) {
	std::vector<glm::vec3> filteredControlPoints;
	std::copy_if(controlPoints.begin(), controlPoints.end(), std::back_inserter(filteredControlPoints),
		[](const glm::vec3& cp) { return cp.z > 0; });

	if (filteredControlPoints.empty()) return;

	std::sort(filteredControlPoints.begin(), filteredControlPoints.end(),
		[](const glm::vec3& a, const glm::vec3& b) { return a.y < b.y; });

	for (auto& vert : mesh.vertices) {
		auto it = std::lower_bound(filteredControlPoints.begin(), filteredControlPoints.end(), vert,
			[](const glm::vec3& cp, const glm::vec3& vert) { return cp.y < vert.y; });

		if (it == filteredControlPoints.begin()) {
			vert.z = it->z;
		}
		else if (it == filteredControlPoints.end()) {
			vert.z = (it - 1)->z;
		}
		else {
			const glm::vec3& low = *(it - 1);
			const glm::vec3& high = *it;
			float ratio = (vert.y - low.y) / (high.y - low.y);
			vert.z = low.z + ratio * (high.z - low.z);
		}
	}
}

void front_inflation_top(Mesh& mesh, const std::vector<glm::vec3>& controlPoints) {
	std::vector<glm::vec3> filteredControlPoints;
	std::copy_if(controlPoints.begin(), controlPoints.end(), std::back_inserter(filteredControlPoints),
		[](const glm::vec3& cp) { return cp.z > 0; });

	if (filteredControlPoints.empty()) return;

	std::sort(filteredControlPoints.begin(), filteredControlPoints.end(),
		[](const glm::vec3& a, const glm::vec3& b) { return a.x < b.x; });

	for (auto& vert : mesh.vertices) {
		auto it = std::lower_bound(filteredControlPoints.begin(), filteredControlPoints.end(), vert,
			[](const glm::vec3& cp, const glm::vec3& vert) { return cp.x < vert.x; });

		if (it == filteredControlPoints.begin()) {
			if (vert.z > it->z) {
				vert.z = it->z;
			}
		}
		else if (it == filteredControlPoints.end()) {
			if (vert.z > (it - 1)->z) {
				vert.z = (it - 1)->z;
			}
		}
		else {
			const glm::vec3& low = *(it - 1);
			const glm::vec3& high = *it;
			float ratio = (vert.x - low.x) / (high.x - low.x);
			float newZ = low.z + ratio * (high.z - low.z);
			if (vert.z > newZ) {
				vert.z = newZ;
			}
		}
	}
}


void back_inflation_side(Mesh& mesh, const std::vector<glm::vec3>& controlPoints) {
	std::vector<glm::vec3> filteredControlPoints;
	std::copy_if(controlPoints.begin(), controlPoints.end(), std::back_inserter(filteredControlPoints),
		[](const glm::vec3& cp) { return cp.z < 0; });

	if (filteredControlPoints.empty()) return;

	std::sort(filteredControlPoints.begin(), filteredControlPoints.end(),
		[](const glm::vec3& a, const glm::vec3& b) { return a.y < b.y; });

	for (auto& vert : mesh.vertices) {
		auto it = std::lower_bound(filteredControlPoints.begin(), filteredControlPoints.end(), vert,
			[](const glm::vec3& cp, const glm::vec3& vert) { return cp.y < vert.y; });

		if (it == filteredControlPoints.begin()) {
			vert.z = it->z;
		}
		else if (it == filteredControlPoints.end()) {
			vert.z = (it - 1)->z;
		}
		else {
			const glm::vec3& low = *(it - 1);
			const glm::vec3& high = *it;
			float ratio = (vert.y - low.y) / (high.y - low.y);
			vert.z = low.z + ratio * (high.z - low.z);
		}
	}
}

void back_inflation_top(Mesh& mesh, const std::vector<glm::vec3>& controlPoints) {
	std::vector<glm::vec3> filteredControlPoints;
	std::copy_if(controlPoints.begin(), controlPoints.end(), std::back_inserter(filteredControlPoints),
		[](const glm::vec3& cp) { return cp.z < 0; });

	if (filteredControlPoints.empty()) return;

	std::sort(filteredControlPoints.begin(), filteredControlPoints.end(),
		[](const glm::vec3& a, const glm::vec3& b) { return a.x < b.x; });

	for (auto& vert : mesh.vertices) {
		auto it = std::lower_bound(filteredControlPoints.begin(), filteredControlPoints.end(), vert,
			[](const glm::vec3& cp, const glm::vec3& vert) { return cp.x < vert.x; });

		if (it == filteredControlPoints.begin()) {
			if (vert.z < it->z) {
				vert.z = it->z;
			}
		}
		else if (it == filteredControlPoints.end()) {
			if (vert.z < (it - 1)->z) {
				vert.z = (it - 1)->z;
			}
		}
		else {
			const glm::vec3& low = *(it - 1);
			const glm::vec3& high = *it;
			float ratio = (vert.x - low.x) / (high.x - low.x);
			float newZ = low.z + ratio * (high.z - low.z);
			if (vert.z < newZ) {
				vert.z = newZ;
			}
		}
	}
}

bool isCounterClockwise(const std::vector<glm::vec3>& controlPoints) {
	if (controlPoints.size() < 3) return false;  // 최소 3점이 필요

	float sum = 0;
	for (size_t i = 0; i < controlPoints.size(); ++i) {
		const glm::vec3& A = controlPoints[i];
		const glm::vec3& B = controlPoints[(i + 1) % controlPoints.size()];
		const glm::vec3& C = controlPoints[(i + 2) % controlPoints.size()];

		float direction = (B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y);
		sum += direction;
	}

	return sum > 0;
}

/*
void connectPlanarMeshes(Mesh& frontMesh, Mesh& backMesh, const std::vector<glm::vec3>& controlPoints) {
	int frontVerticesCount = frontMesh.vertices.size();

	for (auto& vert : frontMesh.vertices) {
		vert.z = 0.01f;
	}

	for (auto& vert : backMesh.vertices) {
		vert.z = -0.01f;
	}

	// backMesh의 정점과 삼각형을 frontMesh에 추가
	frontMesh.vertices.insert(frontMesh.vertices.end(), backMesh.vertices.begin(), backMesh.vertices.end());
	frontMesh.normals.insert(frontMesh.normals.end(), backMesh.normals.begin(), backMesh.normals.end());

	for (const auto& triangle : backMesh.triangles) {
		std::vector<int> newTriangle = { triangle[0] + frontVerticesCount, triangle[1] + frontVerticesCount, triangle[2] + frontVerticesCount };
		frontMesh.triangles.push_back(newTriangle);
	}


	for (int i = 0; i < controlPoints.size(); i++) {
		// Find vertices in frontMesh and backMesh that match controlPoints[i] and controlPoints[i+1]
		int frontIndexA = -1, backIndexA = -1;
		int frontIndexB = -1, backIndexB = -1;

		for (int j = 0; j < frontVerticesCount; j++) {
			if (glm::vec2(frontMesh.vertices[j].x, frontMesh.vertices[j].y) == glm::vec2(controlPoints[i].x, controlPoints[i].y)) {
				frontIndexA = j;
				break; // Found the matching vertex, no need to continue searching
			}
		}


		if (i < controlPoints.size() - 1) {
			for (int j = 0; j < frontVerticesCount; j++) {
				if (glm::vec2(frontMesh.vertices[j].x, frontMesh.vertices[j].y) == glm::vec2(controlPoints[i + 1].x, controlPoints[i + 1].y)) {
					frontIndexB = j;
					break; // Found the matching vertex, no need to continue searching
				}
			}
		}
		else {
			frontIndexB = 0;
		}

		for (int j = 0; j < backMesh.vertices.size(); j++) {
			if (glm::vec2(backMesh.vertices[j].x, backMesh.vertices[j].y) == glm::vec2(controlPoints[i].x, controlPoints[i].y)) {
				backIndexA = j + frontVerticesCount;
				break; // Found the matching vertex, no need to continue searching
			}
		}

		if (i < controlPoints.size() - 1) {
			for (int j = 0; j < backMesh.vertices.size(); j++) {
				if (glm::vec2(backMesh.vertices[j].x, backMesh.vertices[j].y) == glm::vec2(controlPoints[i + 1].x, controlPoints[i + 1].y)) {
					backIndexB = j + frontVerticesCount;
					break; // Found the matching vertex, no need to continue searching
				}
			}
		}
		else {
			backIndexB = frontVerticesCount;
		}

		// Connect the vertices with triangles
		if (isCounterClockwise(controlPoints)) {
			frontMesh.triangles.push_back({ frontIndexA + 1, backIndexA + 1, frontIndexB + 1 });
			frontMesh.triangles.push_back({ frontIndexB + 1, backIndexA + 1, backIndexB + 1 });
		}
		else {
			frontMesh.triangles.push_back({ frontIndexA + 1, frontIndexB + 1, backIndexA + 1 });
			frontMesh.triangles.push_back({ frontIndexB + 1, backIndexB + 1, backIndexA + 1 });
		}
	}

	frontMesh.normals = calculateVertexNormals(frontMesh);
}
*/
void generate3dMesh(Mesh& frontMesh, const std::vector<int> controlPoints) {
	int frontVerticesCount = frontMesh.vertices.size();
	std::vector<glm::vec3> points;

	for (auto& p : controlPoints) {
		glm::vec3 temp = frontMesh.vertices[p - 1];
		points.push_back(temp);
	}

	// generate back mesh
	Mesh backMesh = frontMesh;
	for (auto& triangle : backMesh.triangles) {
		int v1 = triangle[0];
		int v2 = triangle[1];
		int v3 = triangle[2];

		std::vector<int> newTriangle = { v3, v2, v1 };
		triangle = newTriangle;
	}

	for (auto& vert : frontMesh.vertices) {
		vert.z = 0.1;
	}

	for (auto& vert : backMesh.vertices) {
		vert.z = -0.1;
	}

	frontMesh.vertices.insert(frontMesh.vertices.end(), backMesh.vertices.begin(), backMesh.vertices.end());
	frontMesh.normals.insert(frontMesh.normals.end(), backMesh.normals.begin(), backMesh.normals.end());

	for (const auto& triangle : backMesh.triangles) {
		std::vector<int> newTriangle = { triangle[0] + frontVerticesCount, triangle[1] + frontVerticesCount, triangle[2] + frontVerticesCount };
		frontMesh.triangles.push_back(newTriangle);
	}

	int frontA;
	int frontB;
	int backA;
	int backB;
	for (int i = 0; i < controlPoints.size(); i++) {
		if (i < controlPoints.size()-1) {
			frontA = controlPoints[i];
			frontB = controlPoints[i + 1];
			backA = frontA + frontVerticesCount;
			backB = frontB + frontVerticesCount;
		}
		else {
			frontA = controlPoints[i];
			frontB = controlPoints[0];
			backA = frontA + frontVerticesCount;
			backB = frontB + frontVerticesCount;
		}

		if (isCounterClockwise(points)) {
			frontMesh.triangles.push_back({ frontA, backA, frontB });
			frontMesh.triangles.push_back({ frontB, backA, backB });
		}
		else {
			frontMesh.triangles.push_back({ frontA, frontB, backA });
			frontMesh.triangles.push_back({ frontB, backB, backA });
		}
		
	}

	frontMesh.normals = calculateVertexNormals(frontMesh);
}

bool isPointInsidePolygon(const glm::vec3& point, const std::vector<glm::vec3>& polygon) {
	int intersections = 0;
	for (size_t i = 0; i < polygon.size(); ++i) {
		const glm::vec3& start = polygon[i];
		const glm::vec3& end = polygon[(i + 1) % polygon.size()];

		// 포인트와 폴리곤 에지가 y-좌표에서 겹치는지 확인
		if ((start.y <= point.y && point.y < end.y) || (end.y <= point.y && point.y < start.y)) {
			// 수평선이 에지와 교차하는 x-좌표 계산
			float x = start.x + (point.y - start.y) * (end.x - start.x) / (end.y - start.y);
			// 교차점이 포인트의 오른쪽에 있으면 교차 횟수 증가
			if (point.x < x) {
				++intersections;
			}
		}
	}
	return (intersections % 2) != 0;
}


void randomDart(CDT::Triangulation<double>& cdt, const std::vector<glm::vec3>& lineVert, float minDist, int newPointsCount) {
	std::mt19937 rng(std::random_device{}());
	std::uniform_real_distribution<float> distX, distY;

	glm::vec3 minPos = lineVert[0];
	glm::vec3 maxPos = lineVert[0];
	for (const auto& point : lineVert) {
		minPos.x = std::min(minPos.x, point.x);
		minPos.y = std::min(minPos.y, point.y);
		maxPos.x = std::max(maxPos.x, point.x);
		maxPos.y = std::max(maxPos.y, point.y);
	}

	distX = std::uniform_real_distribution<float>(minPos.x, maxPos.x);
	distY = std::uniform_real_distribution<float>(minPos.y, maxPos.y);

	std::vector<glm::vec3> placedPoints(lineVert.begin(), lineVert.end());

	for (int i = 0; i < newPointsCount; ++i) {
		bool pointAccepted = false;
		glm::vec3 newPos;
		int maxAttempts = 100;
		int attempts = 0;

		while (!pointAccepted && attempts < maxAttempts) {
			double x = distX(rng);
			double y = distY(rng);

			newPos = glm::vec3(x, y, 0.0f);
			pointAccepted = true;

			if (isPointInsidePolygon(newPos, lineVert)) {
				for (const auto& placedPoint : placedPoints) {
					if (glm::distance(newPos, placedPoint) < minDist) {
						pointAccepted = false;
						break;
					}
				}
			}
			else {
				pointAccepted = false;
			}

			attempts++;
		}

		if (pointAccepted) {
			placedPoints.push_back(newPos);
			cdt.insertVertices({ CDT::V2d<double>::make(newPos.x, newPos.y) });
		}
	}
}

void inflation_side(Mesh& mesh, const std::vector<glm::vec3>& controlPoints) {
	const float threshold = 0.01f; // 정점 중복성 검사를 위한 임계값

	std::vector<glm::vec3> positiveVertices;
	std::vector<glm::vec3> negativeVertices;

	std::copy_if(controlPoints.begin(), controlPoints.end(), std::back_inserter(positiveVertices),
		[](const glm::vec3& cp) { return cp.z > 0; });
	std::copy_if(controlPoints.begin(), controlPoints.end(), std::back_inserter(negativeVertices),
		[](const glm::vec3& cp) { return cp.z < 0; });

	if (positiveVertices.empty() || negativeVertices.empty()) return;

	std::sort(positiveVertices.begin(), positiveVertices.end(),
		[](const glm::vec3& a, const glm::vec3& b) { return a.y < b.y; });
	std::sort(negativeVertices.begin(), negativeVertices.end(),
		[](const glm::vec3& a, const glm::vec3& b) { return a.y < b.y; });

	for (auto& vert : mesh.vertices) {
		glm::vec3 proposedChange = vert; // 변경을 위한 임시 변수

		// 정점 위치 변경 결정 로직
		if (vert.z > 0) {
			// 기존 로직에 따라 proposedChange.z를 적절히 설정
			auto it = std::lower_bound(positiveVertices.begin(), positiveVertices.end(), vert,
				[](const glm::vec3& cp, const glm::vec3& vert) { return cp.y < vert.y; });

			if (it == positiveVertices.begin()) {
				proposedChange.z = it->z;
			}
			else if (it == positiveVertices.end()) {
				proposedChange.z = (it - 1)->z;
			}
			else {
				const glm::vec3& low = *(it - 1);
				const glm::vec3& high = *it;
				float ratio = (vert.y - low.y) / (high.y - low.y);
				proposedChange.z = low.z + ratio * (high.z - low.z);
			}
		}
		else if (vert.z < 0) {
			// 기존 로직에 따라 proposedChange.z를 적절히 설정
			auto it = std::lower_bound(negativeVertices.begin(), negativeVertices.end(), vert,
				[](const glm::vec3& cp, const glm::vec3& vert) { return cp.y < vert.y; });

			if (it == negativeVertices.begin()) {
				proposedChange.z = it->z;
			}
			else if (it == negativeVertices.end()) {
				proposedChange.z = (it - 1)->z;
			}
			else {
				const glm::vec3& low = *(it - 1);
				const glm::vec3& high = *it;
				float ratio = (vert.y - low.y) / (high.y - low.y);
				proposedChange.z = low.z + ratio * (high.z - low.z);
			}
		}

		// 기존 메쉬 정점들 중에서 proposedChange와 "충분히 가까운" 정점이 있는지 검사
		bool isCloseVertexExist = std::any_of(mesh.vertices.begin(), mesh.vertices.end(),
			[&](const glm::vec3& existingVert) {
				return glm::distance(existingVert, proposedChange) < threshold;
			});

		// "충분히 가까운" 정점이 없는 경우에만 변경 적용
		if (!isCloseVertexExist) {
			vert.z = proposedChange.z;
		}

	}
}

void inflation_top(Mesh& mesh, const std::vector<glm::vec3>& controlPoints) {
	const float threshold = 0.01f;

	std::vector<glm::vec3> positiveVertices;
	std::vector<glm::vec3> negativeVertices;

	std::copy_if(controlPoints.begin(), controlPoints.end(), std::back_inserter(positiveVertices),
		[](const glm::vec3& cp) { return cp.z > 0; });
	std::copy_if(controlPoints.begin(), controlPoints.end(), std::back_inserter(negativeVertices),
		[](const glm::vec3& cp) { return cp.z < 0; });

	if (positiveVertices.empty()) return;
	if (negativeVertices.empty()) return;

	std::sort(positiveVertices.begin(), positiveVertices.end(),
		[](const glm::vec3& a, const glm::vec3& b) { return a.x < b.x; });
	std::sort(negativeVertices.begin(), negativeVertices.end(),
		[](const glm::vec3& a, const glm::vec3& b) { return a.x < b.x; });

	for (auto& vert : mesh.vertices) {
		glm::vec3 proposedChange = vert;

		if (vert.z > 0) {
			auto it = std::lower_bound(positiveVertices.begin(), positiveVertices.end(), vert,
				[](const glm::vec3& cp, const glm::vec3& vert) { return cp.x < vert.x; });

			if (it == positiveVertices.begin()) {
				if (proposedChange.z > it->z) {
					proposedChange.z = it->z;
				}
			}
			else if (it == positiveVertices.end()) {
				if (proposedChange.z > (it - 1)->z) {
					proposedChange.z = (it - 1)->z;
				}
			}
			else {
				const glm::vec3& low = *(it - 1);
				const glm::vec3& high = *it;
				float ratio = (vert.x - low.x) / (high.x - low.x);
				float newZ = low.z + ratio * (high.z - low.z);
				if (proposedChange.z > newZ) {
					proposedChange.z = newZ;
				}
			}
		}
		else if (vert.z < 0) {
			auto it = std::lower_bound(negativeVertices.begin(), negativeVertices.end(), vert,
				[](const glm::vec3& cp, const glm::vec3& vert) { return cp.x < vert.x; });

			if (it == negativeVertices.begin()) {
				if (proposedChange.z < it->z) {
					proposedChange.z = it->z;
				}
			}
			else if (it == negativeVertices.end()) {
				if (proposedChange.z < (it - 1)->z) {
					proposedChange.z = (it - 1)->z;
				}
			}
			else {
				const glm::vec3& low = *(it - 1);
				const glm::vec3& high = *it;
				float ratio = (vert.x - low.x) / (high.x - low.x);
				float newZ = low.z + ratio * (high.z - low.z);
				if (proposedChange.z < newZ) {
					proposedChange.z = newZ;
				}
			}
		}

		bool isCloseVertexExist = std::any_of(mesh.vertices.begin(), mesh.vertices.end(),
			[&](const glm::vec3& existingVert) {
				return glm::distance(existingVert, proposedChange) < threshold;
			});

		// "충분히 가까운" 정점이 없는 경우에만 변경 적용
		if (!isCloseVertexExist) {
			vert.z = proposedChange.z;
		}
		else {
			if (proposedChange.z > 0) {
				vert.z = proposedChange.z + 0.015;
			}
			else {
				vert.z = proposedChange.z - 0.015;
			}
		}
	}
}

void inflate_front_back(
	Mesh& mesh,
	const std::vector<glm::vec3>& side_line,
	const std::vector<glm::vec3>& top_line) {
	inflation_side(mesh, side_line);
	inflation_top(mesh, top_line);

	mesh.normals = calculateVertexNormals(mesh);
}

// EdgeKey 타입 정의 및 makeEdgeKey 함수
typedef std::pair<int, int> EdgeKey;
EdgeKey makeEdgeKey(int v1, int v2) {
	return { std::min(v1, v2), std::max(v1, v2) };
}

HalfEdgeMesh MeshToHalfEdge(Mesh& mesh) {
	HalfEdgeMesh heMesh;
	std::unordered_map<EdgeKey, HalfEdge*, EdgeKeyHash> edgeMap;
	std::vector<Vertex*> vertexPtrs;
	std::vector<Face*> facePtrs;
	std::vector<HalfEdge*> halfEdgePtrs;

	// setting the vertices
	for (int i = 0; i < mesh.vertices.size(); i++) {
		Vertex v;
		v.halfEdge = nullptr;
		v.position = mesh.vertices[i];
		v.normal = mesh.normals[i];
		v.isNew = false;
		heMesh.vertices.push_back(v);

		Vertex* vertexPtr = &heMesh.vertices.back();
		vertexPtrs.push_back(vertexPtr);
	}

	// setting the faces and halfedges
	for (int i = 0; i < mesh.triangles.size(); i++) {
		Face newFace;
		heMesh.faces.push_back(newFace);
		Face* newFacePtr = &heMesh.faces.back();
		facePtrs.push_back(newFacePtr);

		newFacePtr->halfEdge = nullptr;


		for (int j = 0; j < 3; j++) {
			HalfEdge newHE;
			heMesh.halfEdges.push_back(newHE);
			HalfEdge* newHEPtr = &heMesh.halfEdges.back();
			halfEdgePtrs.push_back(newHEPtr);

			int vertexIndex = mesh.triangles[i][j] - 1;
			int nextVertexIndex = mesh.triangles[i][(j + 1) % 3] - 1;

			newHEPtr->vertex = vertexPtrs[vertexIndex];
			newHEPtr->face = newFacePtr;
			newHEPtr->next = nullptr;
			newHEPtr->edge = nullptr;

			if (vertexPtrs[vertexIndex]->halfEdge == nullptr) vertexPtrs[vertexIndex]->halfEdge = newHEPtr;

			// setting up pairs
			EdgeKey key = makeEdgeKey(vertexIndex, nextVertexIndex);

			// Before adding the new half-edge to heMesh, check if its pair already exists
			auto it = edgeMap.find(key);
			if (it != edgeMap.end()) {
				HalfEdge* existingHE = it->second;
				newHEPtr->pair = existingHE;
				existingHE->pair = newHEPtr;

				// setting edges
				Edge edge;
				heMesh.edges.push_back(edge);
				Edge* edgePtr = &heMesh.edges.back();

				edgePtr->halfEdge = newHEPtr;
				edgePtr->isNew = false;

				newHEPtr->edge = edgePtr;
			}
			else {
				// Pair not found, add the new half-edge to the map for future pairing
				edgeMap[key] = newHEPtr;
			}

			if (newFacePtr->halfEdge == nullptr) newFacePtr->halfEdge = newHEPtr;
			//if (vertexPtrs[vertexIndex]->halfEdge == nullptr) vertexPtrs[vertexIndex]->halfEdge = newHEPtr;

			// update halfedge.next
			if (j > 0) {
				// set current half edge to the next half edge of the previous one
				halfEdgePtrs[halfEdgePtrs.size() - 2]->next = newHEPtr;
			}
		}
		halfEdgePtrs[halfEdgePtrs.size() - 1]->next = halfEdgePtrs[halfEdgePtrs.size() - 3];
	}


	return heMesh;
}

std::vector<Vertex*> getNeighbours(Vertex* vertex) {
	std::vector<Vertex*> neighbours;

	HalfEdge* startEdge = vertex->halfEdge->pair;
	HalfEdge* currentEdge = startEdge;

	do {
		neighbours.push_back(currentEdge->vertex);
		currentEdge = currentEdge->next->pair;

	} while (currentEdge != startEdge);

	return neighbours;
}

Vertex* splitEdge(Edge* e, HalfEdgeMesh& he) {
	HalfEdge* StN = e->halfEdge;
	HalfEdge* NtS = StN->pair;

	HalfEdge* NtW = StN->next;
	HalfEdge* WtS = NtW->next;

	HalfEdge* StE = NtS->next;
	HalfEdge* EtN = StE->next;

	Vertex* N = NtS->vertex;
	Vertex* S = StN->vertex;
	Vertex* E = EtN->vertex;
	Vertex* W = WtS->vertex;

	Face* L = StN->face;
	Face* R = NtS->face;

	// Renaming
	HalfEdge* CtN = StN;
	HalfEdge* NtC = NtS;

	Face* UL = L;
	Face* UR = R;

	// New Elements
	Vertex newVertex;
	he.vertices.push_back(newVertex);
	Vertex* C = &he.vertices.back();
	C->isNew = true;

	Face newFace1;
	he.faces.push_back(newFace1);
	Face* BL = &he.faces.back();
	Face newFace2;
	he.faces.push_back(newFace2);
	Face* BR = &he.faces.back();

	HalfEdge newHalfEdge1;
	he.halfEdges.push_back(newHalfEdge1);
	HalfEdge* StC = &he.halfEdges.back();
	HalfEdge newHalfEdge2;
	he.halfEdges.push_back(newHalfEdge2);
	HalfEdge* CtS = &he.halfEdges.back();
	HalfEdge newHalfEdge3;
	he.halfEdges.push_back(newHalfEdge3);
	HalfEdge* WtC = &he.halfEdges.back();
	HalfEdge newHalfEdge4;
	he.halfEdges.push_back(newHalfEdge4);
	HalfEdge* CtW = &he.halfEdges.back();
	HalfEdge newHalfEdge5;
	he.halfEdges.push_back(newHalfEdge5);
	HalfEdge* EtC = &he.halfEdges.back();
	HalfEdge newHalfEdge6;
	he.halfEdges.push_back(newHalfEdge6);
	HalfEdge* CtE = &he.halfEdges.back();

	// Also, add three Edges using the above halfedges.
	Edge newEdge1;
	he.edges.push_back(newEdge1);
	Edge* edgeS = &he.edges.back();
	Edge newEdge2;
	he.edges.push_back(newEdge2);
	Edge* edgeW = &he.edges.back();
	Edge newEdge3;
	he.edges.push_back(newEdge3);
	Edge* edgeE = &he.edges.back();

	edgeS->halfEdge = StC;
	edgeW->halfEdge = WtC;
	edgeE->halfEdge = EtC;

	edgeS->isNew = false;
	edgeW->isNew = true;
	edgeE->isNew = true;

	// Now, the switching
	CtN->vertex = C;
	CtE->vertex = C;
	CtW->vertex = C;
	CtS->vertex = C;

	WtC->vertex = W;
	StC->vertex = S;
	EtC->vertex = E;

	C->halfEdge = CtN;
	S->halfEdge = StE;
	//W->halfEdge = WtS;
	//E->halfEdge = EtN;

	UL->halfEdge = CtN;
	UR->halfEdge = NtC;

	BL->halfEdge = StC;
	BR->halfEdge = CtS;

	WtC->face = UL;
	CtE->face = UR;

	StC->face = CtW->face = WtS->face = BL;
	CtS->face = StE->face = EtC->face = BR;

	WtC->next = CtN;
	NtW->next = WtC;

	NtC->next = CtE;
	CtE->next = EtN;

	// Then all 3 next's for each of the 2 bottom faces...
	StC->next = CtW;
	CtW->next = WtS;
	WtS->next = StC;

	CtS->next = StE;
	StE->next = EtC;
	EtC->next = CtS;

	// Then you set .pair for all 6 of the new halfEdges...
	StC->pair = CtS;
	CtS->pair = StC;
	WtC->pair = CtW;
	CtW->pair = WtC;
	EtC->pair = CtE;
	CtE->pair = EtC;

	return C;
}

void flipEdge(Edge* e) {
	HalfEdge* StN = e->halfEdge;
	HalfEdge* NtS = StN->pair;

	HalfEdge* NtW = StN->next;
	HalfEdge* WtS = NtW->next;

	HalfEdge* StE = NtS->next;
	HalfEdge* EtN = StE->next;

	// added
	HalfEdge* NtE = EtN->pair;
	HalfEdge* StW = WtS->pair;

	Vertex* N = NtS->vertex;
	Vertex* S = StN->vertex;
	Vertex* E = EtN->vertex;
	Vertex* W = WtS->vertex;

	Face* f0 = StN->face;
	Face* f1 = NtS->face;

	// Renaming
	HalfEdge* WtE = StN;
	HalfEdge* EtW = NtS;

	N->halfEdge = NtW;
	S->halfEdge = StE;

	WtE->vertex = W;
	EtW->vertex = E;

	f0->halfEdge = WtE;
	f1->halfEdge = EtW;

	WtE->face = f0;
	EtN->face = f0;
	NtW->face = f0;

	EtW->face = f1;
	WtS->face = f1;
	StE->face = f1;

	// HE connectivity changes.
	WtE->next = EtN;
	EtN->next = NtW;
	NtW->next = WtE;

	EtW->next = WtS;
	WtS->next = StE;
	StE->next = EtW;
}

struct Vec3Hash {
	size_t operator()(const glm::vec3& v) const {
		// 간단한 해시 함수 예시; 실제 사용 시에는 더 나은 충돌 회피 방법을 고려해야 합니다.
		return std::hash<float>()(v.x) ^ std::hash<float>()(v.y) ^ std::hash<float>()(v.z);
	}
};
Mesh HalfEdgeToMesh(HalfEdgeMesh& he) {
	Mesh mesh;
	std::unordered_map<glm::vec3, int, Vec3Hash> vertexIndex;

	// adding vertices
	int index = 0;
	for (auto& vert : he.vertices) {
		mesh.vertices.push_back(vert.position);
		vertexIndex[vert.position] = index;
		index++;
	}

	// adding triangles
	for (auto& face : he.faces) {
		glm::vec3 v1 = face.halfEdge->vertex->position;
		glm::vec3 v2 = face.halfEdge->next->vertex->position;
		glm::vec3 v3 = face.halfEdge->next->next->vertex->position;

		auto it = vertexIndex.find(v1);
		int v1Index = it->second + 1;
		it = vertexIndex.find(v2);
		int v2Index = it->second + 1;
		it = vertexIndex.find(v3);
		int v3Index = it->second + 1;

		std::vector<int> triangle = { v1Index , v2Index , v3Index };
		mesh.triangles.push_back(triangle);
	}

	mesh.normals = calculateVertexNormals(mesh);

	return mesh;
}

void fake_loopSubdivision(Mesh& mesh) {
	std::vector<glm::vec3> new_vertices = mesh.vertices;
	std::vector<std::vector<int>> new_triangles;
	std::map<std::pair<int, int>, int> edgeMidpointIndices;

	// 에지의 중점을 계산하고 새로운 정점으로 추가
	for (const auto& triangle : mesh.triangles) {
		for (int i = 0; i < 3; ++i) {
			int v1 = triangle[i];
			int v2 = triangle[(i + 1) % 3];
			std::pair<int, int> edge = std::make_pair(std::min(v1, v2), std::max(v1, v2));

			if (edgeMidpointIndices.find(edge) == edgeMidpointIndices.end()) {
				glm::vec3 midpoint = (mesh.vertices[v1 - 1] + mesh.vertices[v2 - 1]) * 0.5f;
				new_vertices.push_back(midpoint);
				edgeMidpointIndices[edge] = new_vertices.size();
			}
		}
	}

	// 새로운 삼각형을 생성
	for (const auto& triangle : mesh.triangles) {
		int midpoints[3];
		for (int i = 0; i < 3; ++i) {
			int v1 = triangle[i];
			int v2 = triangle[(i + 1) % 3];
			std::pair<int, int> edge = std::make_pair(std::min(v1, v2), std::max(v1, v2));
			midpoints[i] = edgeMidpointIndices[edge];
		}

		new_triangles.push_back({ triangle[0], midpoints[0], midpoints[2] });
		new_triangles.push_back({ triangle[1], midpoints[1], midpoints[0] });
		new_triangles.push_back({ triangle[2], midpoints[2], midpoints[1] });
		new_triangles.push_back({ midpoints[0], midpoints[1], midpoints[2] });
	}

	// 메시 업데이트
	mesh.vertices = new_vertices;
	mesh.triangles = new_triangles;
	mesh.normals = calculateVertexNormals(mesh);
}

void loopSubdivision(Mesh& mesh, bool faceSplit) {
	HalfEdgeMesh heMesh = MeshToHalfEdge(mesh);

	// 1. Mark all vertices as belonging to the original mesh
	for (Vertex& v : heMesh.vertices) {
		v.isNew = false;
	}

	// 2. Compute updated positions for all vertex - vertex positions and *store *
	for (Vertex& vert : heMesh.vertices) {
		auto neighbours = getNeighbours(&vert);

		glm::vec3 neighbourSum(0.f);
		for (Vertex* neighbour : neighbours) {
			neighbourSum += neighbour->position;
		}

		float n = (float)(neighbours.size());
		float termToSquare = (3.f / 8.f) + 0.25f * cos((2.f * PI) / n);
		float alpha = (1.f / n) * ((5.f / 8.f) - termToSquare * termToSquare);

		vert.newPosition = (1.f - n * alpha) * vert.position + alpha * neighbourSum;
	}

	if (faceSplit) {
		// 3. Compute and *store* edge-vertex positions
		for (Edge& edge : heMesh.edges) {
			HalfEdge* h = edge.halfEdge;
			glm::vec3 v1 = h->vertex->position;
			glm::vec3 v2 = h->pair->vertex->position;
			glm::vec3 v3 = h->pair->next->pair->vertex->position;
			glm::vec3 v4 = h->next->pair->vertex->position;

			edge.newPosition = (3.f / 8.f) * (v1 + v2) + (1.f / 8.f) * (v3 + v4);
		}


		// 4. Split every edge in the mesh.
		std::list<Edge> originalEdges = heMesh.edges;
		for (Edge& e : originalEdges) {
			glm::vec3 newPositionCopy = e.newPosition;

			Vertex* newVert = splitEdge(&e, heMesh);
			newVert->newPosition = newPositionCopy;
		}

		// 5. Flip any *new* edge that connects an old and new vertex.
		for (Edge& e : heMesh.edges) {
			if (e.isNew) {
				bool firstVertNew = e.halfEdge->vertex->isNew;
				bool secondVertNew = e.halfEdge->pair->vertex->isNew;
				if (firstVertNew != secondVertNew) flipEdge(&e);
			}
		}
	}

	// 6. Finally, copy the new vertex positions (Vertex::newPosition) into the
	// usual vertex positions (Vertex::position).
	for (Vertex& v : heMesh.vertices) {
		v.position = v.newPosition;
	}

	mesh = HalfEdgeToMesh(heMesh);


}

void laplacian(Mesh& mesh) {
	HalfEdgeMesh heMesh = MeshToHalfEdge(mesh);

	for (auto& vert : heMesh.vertices) {

		auto neighbours = getNeighbours(&vert);
		int n = neighbours.size();

		glm::vec3 m = glm::vec3(0.f, 0.f, 0.f);
		for (auto& neighbour : neighbours) {
			m.x = m.x + neighbour->position.x;
			m.y = m.y + neighbour->position.y;
			m.z = m.z + neighbour->position.z;
		}

		m.x = m.x / n;
		m.y = m.y / n;
		m.z = m.z / n;

		vert.position.x = m.x;
		vert.position.y = m.y;
		vert.position.z = m.z;
	}

	mesh = HalfEdgeToMesh(heMesh);
}

std::pair<int, int> make_ordered_edge(int v1, int v2) {
	return v1 < v2 ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
}
std::vector<int> boundary_vertices(Mesh& mesh) {
	std::unordered_map<std::pair<int, int>, int, EdgeKeyHash> edgeCount;
	std::unordered_map<int, std::vector<int>> adjacencyList;
	std::vector<int> result;

	// count edge
	for (auto& triangle : mesh.triangles) {
		for (int i = 0; i < triangle.size(); ++i) {
			auto edge = make_ordered_edge(triangle[i], triangle[(i + 1) % triangle.size()]);
			edgeCount[edge]++;
		}
	}

	for (const auto& item : edgeCount) {
		if (item.second == 1) { // 바운더리 엣지
			adjacencyList[item.first.first].push_back(item.first.second);
			adjacencyList[item.first.second].push_back(item.first.first);
		}
	}

	// 바운더리 정점 순서대로 추출
	if (!adjacencyList.empty()) {
		std::unordered_set<int> visited;
		int start = adjacencyList.begin()->first;
		int current = start;
		int previous = -1;
		do {
			result.push_back(current);
			visited.insert(current);
			bool foundNext = false;
			for (int next : adjacencyList[current]) {
				if (next != previous && visited.find(next) == visited.end()) {
					previous = current;
					current = next;
					foundNext = true;
					break;
				}
			}
			if (!foundNext) break; // 중단 조건 추가
		} while (current != start);
	}

	return result;
}

void drawCommonElements(
	GPU_Geometry& gpuGeom,
	CPU_Geometry& lineCpu,
	const std::vector<glm::vec3>& gridVerts,
	const std::vector<glm::vec3>& gridCols
) {

	gpuGeom.setVerts(gridVerts);
	gpuGeom.setCols(gridCols);
	glDrawArrays(GL_LINE_STRIP, 0, GLsizei(2));
	glDrawArrays(GL_LINE_STRIP, 2, GLsizei(2));

}

int Eulerian_Trail = 0;

void phaseFront(std::shared_ptr<MyCallbacks>& cb,
	std::vector<std::vector<glm::vec3>>& lineVerts, CPU_Geometry& lineCpu,
	std::vector<std::vector<glm::vec3>>& controlPointVerts, CPU_Geometry& controlPointCpu,
	std::vector<std::vector<glm::vec3>>& transformedVerts,
	GPU_Geometry& gpuGeom,
	int& cross_section) {

	if (clear) {
		printf("front cleared");
		lineVerts[0].clear();
		clear = false;
		Eulerian_Trail = 0;
	}

	if (Eulerian_Trail == 1) {
		draw_cross_sections(cb, lineVerts, lineCpu, cross_section);
	}

	gpuGeom.setVerts(lineCpu.verts);
	gpuGeom.setCols(lineCpu.cols);
	glDrawArrays(GL_LINE_STRIP, 0, GLsizei(lineVerts[0].size()));
}

void phaseSide(std::shared_ptr<MyCallbacks>& cb,
	std::vector<std::vector<glm::vec3>>& lineVerts, CPU_Geometry& lineCpu,
	std::vector<std::vector<glm::vec3>>& controlPointVerts, CPU_Geometry& controlPointCpu,
	std::vector<std::vector<glm::vec3>>& transformedVerts,
	GPU_Geometry& gpuGeom,
	int& cross_section) {


	if (clear) {
		printf("side cleared");
		lineVerts[1].clear();
		clear = false;
		Eulerian_Trail = 2;
	}

	if (Eulerian_Trail == 3) {
		draw_cross_sections(cb, lineVerts, lineCpu, cross_section);
	}


	gpuGeom.setVerts(lineCpu.verts);
	gpuGeom.setCols(lineCpu.cols);

	glDrawArrays(GL_LINE_STRIP, GLsizei(lineVerts[0].size()), GLsizei(lineVerts[1].size()));

}

void phaseTop(std::shared_ptr<MyCallbacks>& cb,
	std::vector<std::vector<glm::vec3>>& lineVerts, CPU_Geometry& lineCpu,
	std::vector<std::vector<glm::vec3>>& controlPointVerts, CPU_Geometry& controlPointCpu,
	std::vector<std::vector<glm::vec3>>& transformedVerts,
	GPU_Geometry& gpuGeom,
	int& cross_section) {



	if (clear) {
		printf("side cleared");
		lineVerts[2].clear();
		clear = false;
		Eulerian_Trail = 4;
	}

	if (Eulerian_Trail == 5) {
		draw_cross_sections(cb, lineVerts, lineCpu, cross_section);
	}


	gpuGeom.setVerts(lineCpu.verts);
	gpuGeom.setCols(lineCpu.cols);

	glDrawArrays(GL_LINE_STRIP, GLsizei(lineVerts[0].size() + lineVerts[1].size()), GLsizei(lineVerts[2].size()));


}

void showDrew(std::shared_ptr<MyCallbacks>& cb,
	std::vector<std::vector<glm::vec3>>& lineVerts, CPU_Geometry& lineCpu,
	std::vector<std::vector<glm::vec3>>& controlPointVerts, CPU_Geometry& controlPointCpu,
	std::vector<std::vector<glm::vec3>>& transformedVerts,
	GPU_Geometry& gpuGeom,
	int& cross_section) {

	glDrawArrays(GL_LINE_STRIP, 0, GLsizei(lineVerts[0].size()));
	glDrawArrays(GL_LINE_STRIP, GLsizei(lineVerts[0].size()), GLsizei(lineVerts[1].size()));
}

void phaseCreateMesh(std::shared_ptr<MyCallbacks>& cb,
	std::vector<std::vector<glm::vec3>>& lineVerts, CPU_Geometry& lineCpu,
	std::vector<std::vector<glm::vec3>>& controlPointVerts, CPU_Geometry& controlPointCpu,
	std::vector<std::vector<glm::vec3>>& transformedVerts,
	GPU_Geometry& gpuGeom,
	int& cross_section) {

	for (int i = 0; i < 3; i++) {
		controlPointVerts[i] = get_control_points(lineVerts[i], 50);
	}
	transform(lineVerts, transformedVerts);
	cross_section++;


	// CDT
	auto cdt = CDT::Triangulation<double>(
		CDT::VertexInsertionOrder::Auto,
		CDT::IntersectingConstraintEdges::TryResolve,
		0.);
	cdt.insertVertices(
		controlPointVerts[0].begin(),
		controlPointVerts[0].end(),
		[](const glm::vec3& p) { return p.x; },
		[](const glm::vec3& p) { return p.y; }
	);

	//insert_Vertices(cdt, lineVerts[0]);
	randomDart(cdt, lineVerts[0], 0.15, 100);

	struct CustomEdge
	{
		std::pair<std::size_t, std::size_t> vertices;
	};

	std::vector<CustomEdge> edges;
	for (int i = 0; i < controlPointVerts[0].size() - 1; i++) {
		CustomEdge edge;
		edge.vertices = std::make_pair(i, i + 1);
		edges.push_back(edge);
	}

	CustomEdge edge;
	edge.vertices = std::make_pair(controlPointVerts[0].size() - 1, 0);
	edges.push_back(edge);
	cdt.insertEdges(
		edges.begin(),
		edges.end(),
		[](const CustomEdge& e) { return e.vertices.first; },
		[](const CustomEdge& e) { return e.vertices.second; }
	);

	cdt.eraseOuterTrianglesAndHoles();


	// generating 3D model starts here
	Mesh front_mesh = get_front_mesh(cdt);
	fake_loopSubdivision(front_mesh);
	fake_loopSubdivision(front_mesh);
	//front_inflation_side(front_mesh, transformedVerts[1]);
	//front_inflation_top(front_mesh, transformedVerts[2]);
	//front_mesh.normals = calculateVertexNormals(front_mesh);

	auto boundary = boundary_vertices(front_mesh);
	generate3dMesh(front_mesh, boundary);

	//Mesh back_mesh = get_back_mesh(cdt);
	//fake_loopSubdivision(back_mesh);
	//back_inflation_side(back_mesh, transformedVerts[1]);
	//back_inflation_top(back_mesh, transformedVerts[2]);
	//back_mesh.normals = calculateVertexNormals(back_mesh);

	//connectPlanarMeshes(front_mesh, back_mesh, controlPointVerts[0]);
	//fake_loopSubdivision(front_mesh);

	//inflate_front_back(front_mesh, transformedVerts[1], transformedVerts[2]);
	//laplacian(front_mesh);
	//loopSubdivision(front_mesh, false);
	
	for (int i = 0; i < 5; i++) {
		inflate_front_back(front_mesh, transformedVerts[1], transformedVerts[2]);
		loopSubdivision(front_mesh, false);
		//laplacian(front_mesh);
	}

	//fake_loopSubdivision(front_mesh);
	//loopSubdivision(front_mesh, true);
	//inflate_front_back(front_mesh, transformedVerts[1], transformedVerts[2]);

	//loopSubdivision(front_mesh, true);
	//loopSubdivision(front_mesh, true);
	//loopSubdivision(front_mesh, true);
	//loopSubdivision(front_mesh, true);
	//loopSubdivision(front_mesh, true);
	//fake_loopSubdivision(front_mesh);
	//front_mesh.normals = calculateVertexNormals(front_mesh);
	//loopSubdivision(front_mesh, false);
	//loopSubdivision(front_mesh, false);
	//loopSubdivision(front_mesh, false);
	//loopSubdivision(front_mesh, false);

	//fake_loopSubdivision(front_mesh);
	//inflation_side(front_mesh, transformedVerts[1]);
	//inflation_top(front_mesh, transformedVerts[2]);
	//front_mesh.normals = calculateVertexNormals(front_mesh);
	//loopSubdivision(front_mesh);
	//inflation_side(front_mesh, transformedVerts[1]);
	//inflation_top(front_mesh, transformedVerts[2]);
	saveMeshToOBJ(front_mesh, "C:/Users/dhktj/OneDrive/Desktop/after.obj");
	//saveMeshToOBJ(front_mesh, "C:/Users/U/Documents/ImaginationModeling/589-689-3D-skeleton/models/merged.obj");
	

	CDT::extractEdgesFromTriangles(cdt.triangles);
}

void showCombined(std::shared_ptr<MyCallbacks>& cb,
	std::vector<std::vector<glm::vec3>>& lineVerts, CPU_Geometry& lineCpu,
	std::vector<std::vector<glm::vec3>>& controlPointVerts, CPU_Geometry& controlPointCpu,
	std::vector<std::vector<glm::vec3>>& transformedVerts,
	GPU_Geometry& gpuGeom,
	int& cross_section) {
	combine(cb, transformedVerts, lineCpu);
	gpuGeom.setVerts(lineCpu.verts);
	gpuGeom.setCols(lineCpu.cols);
	int start_index = 0;
	for (int i = 0; i < 3; i++) {
		glDrawArrays(GL_LINE_STRIP, start_index, GLsizei(transformedVerts[i].size()));
		start_index += transformedVerts[i].size();
	}
}


int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit();
	Window window(800, 800, "CPSC 589/689"); // could set callbacks at construction if desired

	GLDebug::enable();

	// SHADERS
	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	auto cb = std::make_shared<MyCallbacks>(shader, window.getWidth(), window.getHeight());
	// CALLBACKS
	window.setCallbacks(cb);

	window.setupImGui(); // Make sure this call comes AFTER GLFW callbacks set.

	// Variables that ImGui will alter.
	float pointSize = 10.0f; // Diameter of drawn points
	float color[3] = { 1.f, 0.f, 0.f }; // Color of new points
	bool drawLines = true; // Whether to draw connecting lines
	int selectedPointIndex = -1; // Used for point dragging & deletion


	// GEOMETRY
	int cross_section = 0;
	CPU_Geometry lineCpu;
	CPU_Geometry controlPointCpu;
	CPU_Geometry bsplineCurveCpu;

	std::vector<std::vector<glm::vec3>> lineVerts(3);
	std::vector<std::vector<glm::vec3>> transformedVerts(3);
	std::vector<std::vector<glm::vec3>> controlPointVerts(3);
	std::vector<std::vector<glm::vec3>> bsplineCurveVerts(3);

	grid.verts = { glm::vec3(-1.0f, .0f, .0f), glm::vec3(1.0f, .0f, .0f), glm::vec3(.0f, 1.0f, .0f), glm::vec3(.0f, -1.0f, .0f) };
	grid.cols = { glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(1.0f, 1.0f, 1.0f), glm::vec3(1.0f, 1.0f, 1.0f) };

	GPU_Geometry gpuGeom;

	//================================================================================================================//

	GLfloat vertices[] = {

		 0.8f, -0.8f, 0.0f,  1.0f, 1.0f, 1.0f,
		 0.95f, -0.8f, 0.0f,  1.0f, 1.0f, 1.0f,
		 0.95f, -0.95f, 0.0f,  1.0f, 1.0f, 1.0f,
		 0.8f, -0.95f, 0.0f,  1.0f, 1.0f, 1.0f,
	};


	GLuint indices[] = {
		0, 1, 2,
		0, 2, 3
	};

	GLfloat clearVertices[] = {

		 0.8f, 0.8f, 0.0f,  1.0f, 1.0f, 1.0f,
		 0.95f, 0.8f, 0.0f,  1.0f, 1.0f, 1.0f,
		 0.95f, 0.95f, 0.0f,  1.0f, 1.0f, 1.0f,
		 0.8f, 0.95f, 0.0f,  1.0f, 1.0f, 1.0f,
	};

	GLuint clearIndices[] = {
		0, 1, 2,
		0, 2, 3
	};

	GLfloat showVertices[] = {
		 0.8f, 0.6f, 0.0f,  1.0f, 1.0f, 1.0f,
		 0.95f, 0.6f, 0.0f,  1.0f, 1.0f, 1.0f,
		 0.95f, 0.75f, 0.0f,  1.0f, 1.0f, 1.0f,
		 0.8f, 0.75f, 0.0f,  1.0f, 1.0f, 1.0f,
	};

	GLuint showIndices[] = {
		0, 1, 2,
		0, 2, 3
	};

	// VAO, VBO, EBO 설정
	GLuint VAO, VBO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	// position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// color attribute
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	//=========================================================================================================================//

	GLuint VAO2, VBO2, EBO2;
	glGenVertexArrays(1, &VAO2);
	glGenBuffers(1, &VBO2);
	glGenBuffers(1, &EBO2);

	glBindVertexArray(VAO2);

	glBindBuffer(GL_ARRAY_BUFFER, VBO2);
	glBufferData(GL_ARRAY_BUFFER, sizeof(clearVertices), clearVertices, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO2);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(clearIndices), clearIndices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	bool button = false;

	//=========================================================================================================================//

	GLuint VAO3, VBO3, EBO3;
	glGenVertexArrays(1, &VAO3);
	glGenBuffers(1, &VBO3);
	glGenBuffers(1, &EBO3);

	glBindVertexArray(VAO3);

	glBindBuffer(GL_ARRAY_BUFFER, VBO3);
	glBufferData(GL_ARRAY_BUFFER, sizeof(showVertices), showVertices, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO3);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(showIndices), showIndices, GL_STATIC_DRAW);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);


	//=========================================================================================================================//
	// RENDER LOOP
	while (!window.shouldClose()) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Tell callbacks object a new frame's begun BEFORE polling events!
		cb->incrementFrameCount();
		glfwPollEvents();

		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);

		glBindVertexArray(VAO2);
		glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
		/*

		// If mouse just went down, see if it was on a point.
		if (cb->leftMouseJustPressed() || cb->rightMouseJustPressed()) {
			// We use the point DIAMETER as the threshold, meaning the user
			// can click anywhere within 2x radius to select.
			// You may want to change that.
			float threshold = pointSize;

			selectedPointIndex = cb->indexOfPointAtCursorPos(cpuGeom.verts, threshold);
		}

		*/

		// when the left button gets pressed increase the cross_section
		if (cb->leftMouseJustPressed()) {
			std::cout << "Position: " << glm::vec3(cb->getCursorPosGL(), 0.f) << std::endl;
			if (cross_section == 0 && Eulerian_Trail < 2) {
				Eulerian_Trail++;
			}

			if (cross_section == 1 && Eulerian_Trail < 4) {
				Eulerian_Trail++;
			}

			if (cross_section == 2 && Eulerian_Trail < 6) {
				Eulerian_Trail++;
			}

			if (cb->getCursorPosGL().x >= 0.8f && cb->getCursorPosGL().x <= 0.95f && cb->getCursorPosGL().y <= -0.8f && cb->getCursorPosGL().y >= -0.95f) {
				button = true;
				if (cross_section < 5) cross_section++;



			}

			if (cb->getCursorPosGL().x >= 0.8f && cb->getCursorPosGL().x <= 0.95f && cb->getCursorPosGL().y >= 0.8f && cb->getCursorPosGL().y <= 0.95f) {
				printf("clear");
				clear = true;


			}


		}


		if (cb->rightMouseJustPressed()) {

			showDraw = !showDraw;
			//std::cout << "Show draw: " << showDraw << std::endl;
		}

		/*
		else if (cb->rightMouseJustPressed()) {
			if (selectedPointIndex >= 0) {
				// If we right-clicked on a vertex, erase it.
				controlPointcpu.verts.erase(controlPointcpu.verts.begin() + selectedPointIndex);
				controlPointcpu.cols.erase(controlPointcpu.cols.begin() + selectedPointIndex);
				selectedPointIndex = -1; // So that we don't drag in next frame.

				controlPointgpu.setVerts(controlPointcpu.verts);
				controlPointgpu.setCols(controlPointcpu.cols);
			}
		}

		else if (cb->leftMouseActive() && selectedPointIndex >= 0) {
			// Drag selected point.
			cpuGeom.verts[selectedPointIndex] = glm::vec3(cb->getCursorPosGL(), 0.f);
			gpuGeom.setVerts(cpuGeom.verts);

		}
		*/


		bool change = false; // Whether any ImGui variable's changed.

		// Three functions that must be called each new frame.
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		ImGui::Begin("Sample window.");

		ImGui::Text("Sample text.");

		change |= ImGui::SliderFloat("Point size", &pointSize, 1.f, 20.f);

		change |= ImGui::ColorEdit3("New pt color", (float*)&color);

		change |= ImGui::Checkbox("Draw lines", &drawLines);

		change |= ImGui::Checkbox("Draw lines", &drawLines);

		if (ImGui::Button("clear pts")) {
			change = true;
			lineVerts[0].clear();
		}

		ImGui::Text("Average %.1f ms/frame (%.1f fps)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

		ImGui::End();
		ImGui::Render();

		shader.use();
		gpuGeom.bind();

		if (change)
		{
			//cb->updateShadingUniforms(lightPos, lightCol, diffuseCol, ambientStrength, texExistence);
		}
		cb->viewPipeline();

		glPointSize(pointSize);

		//glEnable(GL_FRAMEBUFFER_SRGB);
		//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// actual drawing happening here
		/*
		draw(
			cb,
			lineVerts, lineCpu,
			controlPointVerts, controlPointCpu,
			bsplineCurveVerts, bsplineCurveCpu,
			transformedVerts,
			gpuGeom,
			cross_section);

		*/
		drawCommonElements(gpuGeom, lineCpu, grid.verts, grid.cols);
		if (cross_section == 0) {

			phaseFront(
				cb,
				lineVerts, lineCpu,
				controlPointVerts, controlPointCpu,
				transformedVerts,
				gpuGeom,
				cross_section);
		}
		else if (cross_section == 1) {
			phaseSide(
				cb,
				lineVerts, lineCpu,
				controlPointVerts, controlPointCpu,
				transformedVerts,
				gpuGeom,
				cross_section);
		}
		else if (cross_section == 2) {
			phaseTop(
				cb,
				lineVerts, lineCpu,
				controlPointVerts, controlPointCpu,
				transformedVerts,
				gpuGeom,
				cross_section);
		}
		else if (cross_section == 3) {
			phaseCreateMesh(
				cb,
				lineVerts, lineCpu,
				controlPointVerts, controlPointCpu,
				transformedVerts,
				gpuGeom,
				cross_section);
		}
		else {
			showCombined(cb,
				lineVerts, lineCpu,
				controlPointVerts, controlPointCpu,
				transformedVerts,
				gpuGeom,
				cross_section);
		}

		if (showDraw) {
			showDrew(cb,
				lineVerts, lineCpu,
				controlPointVerts, controlPointCpu,
				transformedVerts,
				gpuGeom,
				cross_section);
		}

		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		window.swapBuffers();
	}

	// Cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	glfwTerminate();
	return 0;
}
