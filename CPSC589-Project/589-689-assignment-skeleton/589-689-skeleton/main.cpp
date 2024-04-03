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
#include <random>

CPU_Geometry grid;
glm::vec3 red = glm::vec3(1, 0, 0);
glm::vec3 green = glm::vec3(0, 1, 0);
glm::vec3 blue = glm::vec3(0, 0, 1);
glm::vec3 white = glm::vec3(1, 1, 1);

int k = 4;
int m;
float ui = 0.01f;

int delta(float u, const std::vector<float>& knots, int m, int k) {
	for (int i = 0; i < m + k - 1; i++) {
		if (u >= knots[i] && u < knots[i + 1]) {
			return i;
		}
	}
	return -1;
}

std::vector<float> knotSequence(int m, int k) {
	int n = k + m;
	std::vector<float> knotSequence(n + 1);

	for (int i = 0; i < k - 1; i++) {
		knotSequence[i] = 0.0f;
	}

	for (int i = 1; i < m - 1; i++) {
		knotSequence[k - 1 + i] = (float)i / (m - 1);
	}

	for (int i = n - k + 1; i <= n; i++) {
		knotSequence[i] = 1.0f;
	}

	return knotSequence;
}


void Bspline(std::vector<glm::vec3>& E, std::vector<glm::vec3>& result, int k, float ui) {
	int m = E.size() - 1;
	result.clear();
	if (m < 1) return;
	if (k > m + 1) return;


	std::vector<float> ks = knotSequence(m, k);
	std::vector<glm::vec3> c(k);

	for (float u = ks[k - 1]; u <= ks[m + 1]; u += ui) {

		int d = delta(u, ks, m, k);
		if (d >= E.size()) {
			return;
		}

		for (int i = 0; i <= k - 1; i++) {
			c[i] = E[d - i];
		}

		for (int r = k; r >= 2; r--) {
			int i = d;
			for (int s = 0; s <= r - 2; s++) {
				float omega = (u - ks[i]) / (ks[i + r - 1] - ks[i]);
				c[s] = omega * c[s] + (1 - omega) * c[s + 1];
				i--;
			}
		}
		result.push_back(c[0]);
	}

	if (!result.empty()) {
		result.push_back(E.back());
	}

	return;
}

struct CustomPoint2D
{
	double data[2];
};

struct CustomEdge
{
	std::pair<std::size_t, std::size_t> vertices;
};

struct Mesh
{
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	std::vector<std::vector<int>> triangles;
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

	float getXAngle(){
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
std::vector<glm::vec3> flattenLineVerts(std::vector<std::vector<glm::vec3>> &lineVerts) {
	std::vector<glm::vec3> flattenedVerts;

	for (const auto& vertGroup : lineVerts) {
		flattenedVerts.insert(flattenedVerts.end(), vertGroup.begin(), vertGroup.end());
	}

	return flattenedVerts;
}

// let user provide cross sections
void draw_cross_sections(
	std::shared_ptr<MyCallbacks> &cb,
	std::vector<std::vector<glm::vec3>> &lineVerts,
	CPU_Geometry &cpuGeom,
	int &cross_section) {

	if (cb->leftMouseActive()) {
		// load cpu geometry
		lineVerts[cross_section].push_back(glm::vec3(cb->getCursorPosGL(), 0.f));
		cpuGeom.verts = flattenLineVerts(lineVerts);

		// load cpu geometry colors
		if (cross_section == 0) {
			cpuGeom.cols.push_back(red);
		}
		else if (cross_section == 1) {
			cpuGeom.cols.push_back(green);
		}
		else {
			cpuGeom.cols.push_back(blue);
		}
	}

	return;
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
	CPU_Geometry& cpuGeom){

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


void inflation_side(Mesh& mesh, const std::vector<glm::vec3>& controlPoints) {
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

void inflation_top(Mesh& mesh, const std::vector<glm::vec3>& controlPoints) {
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

void connectPlanarMeshes(Mesh& frontMesh, const Mesh& backMesh, const std::vector<glm::vec3>& controlPoints) {
	int frontVerticesCount = frontMesh.vertices.size();

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

		for (int j = 0; j < frontMesh.vertices.size(); j++) {
			if (glm::vec2(frontMesh.vertices[j].x, frontMesh.vertices[j].y) == glm::vec2(controlPoints[i].x, controlPoints[i].y)) {
				frontIndexA = j;
				break; // Found the matching vertex, no need to continue searching
			}
		}

		for (int j = 0; j < frontMesh.vertices.size(); j++) {
			if (i < controlPoints.size() - 1) {
				if (glm::vec2(frontMesh.vertices[j].x, frontMesh.vertices[j].y) == glm::vec2(controlPoints[i + 1].x, controlPoints[i + 1].y)) {
					frontIndexB = j;
					break; // Found the matching vertex, no need to continue searching
				}
			}
			else {
				frontIndexB = 0;
			}
			
		}

		for (int j = 0; j < backMesh.vertices.size(); j++) {
			if (glm::vec2(backMesh.vertices[j].x, backMesh.vertices[j].y) == glm::vec2(controlPoints[i].x, controlPoints[i].y)) {
				backIndexA = j + frontVerticesCount;
				break; // Found the matching vertex, no need to continue searching
			}
		}

		for (int j = 0; j < backMesh.vertices.size(); j++) {
			if (i < controlPoints.size() - 1) {
				if (glm::vec2(backMesh.vertices[j].x, backMesh.vertices[j].y) == glm::vec2(controlPoints[i + 1].x, controlPoints[i + 1].y)) {
					backIndexB = j + frontVerticesCount;
					break; // Found the matching vertex, no need to continue searching
				}
			}
			else {
				backIndexB = frontVerticesCount;
			}
		}

		// Connect the vertices with triangles
		if (frontIndexA != -1 && frontIndexB != -1 && backIndexA != -1 && backIndexB != -1) {
			frontMesh.triangles.push_back({ frontIndexA + 1, frontIndexB + 1, backIndexA + 1 });
			frontMesh.triangles.push_back({ frontIndexB + 1, backIndexB + 1, backIndexA + 1 });
		}
	}
	
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

void loopSubdivision(Mesh& mesh) {
	std::vector<glm::vec3> new_vertices = mesh.vertices;
	std::vector<std::vector<int>> new_triangles;
	std::map<std::pair<int, int>, int> edgeMidpointIndices;

	// create edge-vertices
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

	// create new triangles
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

	// update the mesh
	mesh.vertices = new_vertices;
	mesh.triangles = new_triangles;
}


// std::pair<int, int>에 대한 사용자 정의 해시 함수
struct pair_hash {
	template <class T1, class T2>
	std::size_t operator () (const std::pair<T1, T2>& pair) const {
		auto hash1 = std::hash<T1>{}(pair.first);
		auto hash2 = std::hash<T2>{}(pair.second);
		return hash1 ^ (hash2 << 1); // 오른쪽으로 비트를 1번 시프트하여 더 좋은 해시 분포를 얻음
	}
};

std::vector<glm::vec3> findBoundaryVertices(const Mesh& mesh) {
	std::unordered_set<int> boundaryVerticesSet;
	// 사용자 정의 해시 함수를 사용하여 unordered_map 선언
	std::unordered_map<std::pair<int, int>, int, pair_hash> edgeCount;

	for (const auto& triangle : mesh.triangles) {
		for (int i = 0; i < 3; ++i) {
			int v1 = std::min(triangle[i], triangle[(i + 1) % 3]);
			int v2 = std::max(triangle[i], triangle[(i + 1) % 3]);
			std::pair<int, int> edge(v1, v2);

			if (edgeCount.find(edge) == edgeCount.end()) {
				edgeCount[edge] = 1;
			}
			else {
				edgeCount[edge]++;
			}
		}
	}

	for (const auto& pair : edgeCount) {
		if (pair.second == 1) { // 에지가 한 번만 나타나면 경계에 있는 것임
			boundaryVerticesSet.insert(pair.first.first);
			boundaryVerticesSet.insert(pair.first.second);
		}
	}

	std::vector<int> boundaryVertices(boundaryVerticesSet.begin(), boundaryVerticesSet.end());
	std::vector<glm::vec3> result;
	for (int index : boundaryVertices) {
		result.push_back(mesh.vertices[index]);
	}
	return result;
}

void draw(
	std::shared_ptr<MyCallbacks>& cb,
	std::vector<std::vector<glm::vec3>>& lineVerts, CPU_Geometry& lineCpu,
	std::vector<std::vector<glm::vec3>>& controlPointVerts, CPU_Geometry& controlPointCpu,
	std::vector<std::vector<glm::vec3>>& bsplineVerts, CPU_Geometry& bsplineCurveCpu,
	std::vector<std::vector<glm::vec3>>& transformedVerts,
	GPU_Geometry& gpuGeom,
	int& cross_section) {

	if (cross_section < 3) {
		// draw grid
		gpuGeom.setVerts(grid.verts);
		gpuGeom.setCols(grid.cols);
		glDrawArrays(GL_LINE_STRIP, 0, GLsizei(2));
		glDrawArrays(GL_LINE_STRIP, 2, GLsizei(2));

		// draw user input lines 
		draw_cross_sections(cb, lineVerts, lineCpu, cross_section);
		gpuGeom.setVerts(lineCpu.verts);
		gpuGeom.setCols(lineCpu.cols);
		int start_index = 0;
		for (int i = 0; i < 3; i++) {
			glDrawArrays(GL_LINE_STRIP, start_index, GLsizei(lineVerts[i].size()));
			start_index += lineVerts[i].size();
		}
	}
	else if(cross_section == 3){
		// get control points
		for (int i = 0; i < 3; i++) {
			controlPointVerts[i] = get_control_points(lineVerts[i], 50);
			Bspline(controlPointVerts[i], bsplineVerts[i], k, ui);
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
		randomDart(cdt, lineVerts[0], 0.15, 50);
		
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

		Mesh front_mesh = get_front_mesh(cdt);

		std::vector<glm::vec3> boundary = findBoundaryVertices(front_mesh);

		inflation_side(front_mesh, transformedVerts[1]);
		inflation_top(front_mesh, transformedVerts[2]);
		front_mesh.normals = calculateVertexNormals(front_mesh);
		saveMeshToOBJ(front_mesh, "C:/Users/dhktj/OneDrive/Desktop/front.obj");


		Mesh back_mesh = get_back_mesh(cdt);
		back_inflation_side(back_mesh, transformedVerts[1]);
		back_inflation_top(back_mesh, transformedVerts[2]);
		back_mesh.normals = calculateVertexNormals(back_mesh);
		saveMeshToOBJ(back_mesh, "C:/Users/dhktj/OneDrive/Desktop/back.obj");


		connectPlanarMeshes(front_mesh, back_mesh, controlPointVerts[0]);
		loopSubdivision(front_mesh);
		loopSubdivision(front_mesh);
		front_mesh.normals = calculateVertexNormals(front_mesh);
		saveMeshToOBJ(front_mesh, "C:/Users/dhktj/OneDrive/Desktop/output.obj");


		cdt.triangles;
		cdt.vertices;
		cdt.fixedEdges;
		CDT::extractEdgesFromTriangles(cdt.triangles);

	}
	else {
		/*
		combine(cb, bsplineVerts, bsplineCurveCpu);
		gpuGeom.setVerts(bsplineCurveCpu.verts);
		gpuGeom.setCols(bsplineCurveCpu.cols);
		int start_index = 0;
		for (int i = 0; i < 3; i++) {
			glDrawArrays(GL_LINE_STRIP, start_index, GLsizei(bsplineVerts[i].size()));
			start_index += bsplineVerts[i].size();
		}
		*/
		
		// show cross sections in 3D
		combine(cb, transformedVerts, lineCpu);
		gpuGeom.setVerts(lineCpu.verts);
		gpuGeom.setCols(lineCpu.cols);
		int start_index = 0;
		for (int i = 0; i < 3; i++) {
			glDrawArrays(GL_LINE_STRIP, start_index, GLsizei(transformedVerts[i].size()));
			start_index += transformedVerts[i].size();
		}
		
	}

	return;
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
	int cross_section = -1;
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

	// RENDER LOOP
	while (!window.shouldClose()) {

		// Tell callbacks object a new frame's begun BEFORE polling events!
		cb->incrementFrameCount();
		glfwPollEvents();

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
			if (cross_section < 5) cross_section++;
		}

		if (cb->rightMouseJustPressed()) {

			float maxFY = findFrontMaxY(transformedVerts[0]);
			float maxFX = findFrontMaxX(transformedVerts[0]);
			float minFY = findFrontMinY(transformedVerts[0]);
			float minFX = findFrontMinX(transformedVerts[0]);
			std::cout << "Max Front X : " << maxFX << std::endl;
			std::cout << "Max Front Y : " << maxFY << std::endl;
			std::cout << "Min Front X : " << minFX << std::endl;
			std::cout << "Min Front Y : " << minFY << std::endl;

			float maxSY = findSideMaxY(transformedVerts[1]);
			float maxSZ = findSideMaxZ(transformedVerts[1]);
			std::cout << "Max Side Y : " << maxSY << std::endl;
			std::cout << "Max Side Z : " << maxSZ << std::endl;


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

		change |= ImGui::SliderInt("Curve's Order", &k, 1, 10);

		change |= ImGui::SliderFloat("u", &ui, 0.01f, 0.99f);

		change |= ImGui::Checkbox("Draw lines", &drawLines);

		if (ImGui::Button("clear pts")) {
			change = true;
			lineVerts.clear();
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

		glEnable(GL_FRAMEBUFFER_SRGB);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// actual drawing happening here
		draw(
			cb,
			lineVerts, lineCpu,
			controlPointVerts, controlPointCpu,
			bsplineCurveVerts, bsplineCurveCpu,
			transformedVerts,
			gpuGeom,
			cross_section);


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
