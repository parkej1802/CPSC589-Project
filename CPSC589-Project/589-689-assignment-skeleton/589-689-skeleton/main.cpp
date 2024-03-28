#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>
#include <unordered_map>
#include <fstream>

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
		/*
		if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
			lastRightPressedFrame = currentFrame;
		}
		*/

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
std::vector<glm::vec3> get_control_points(std::vector<glm::vec3>& line, int smoothness) {
	std::vector<glm::vec3> result;
	int step_size = line.size() / smoothness;
	int index = 0;

	while (index < line.size()) {
		result.push_back(line[index]);
		index += step_size;
	}

	auto it = std::unique(result.begin(), result.end());
	result.erase(it, result.end());

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

std::vector<glm::vec3> calculateVertexNormals(const CDT::Triangulation<double>& cdt) {
	std::vector<glm::vec3> normals(cdt.vertices.size(), glm::vec3(0.0f, 0.0f, 0.0f));

	for (const auto& triangle : cdt.triangles) {
		auto v1_index = triangle.vertices[0];
		auto v2_index = triangle.vertices[1];
		auto v3_index = triangle.vertices[2];

		glm::vec3 v1(cdt.vertices[v1_index].x, cdt.vertices[v1_index].y, 0.0f);
		glm::vec3 v2(cdt.vertices[v2_index].x, cdt.vertices[v2_index].y, 0.0f);
		glm::vec3 v3(cdt.vertices[v3_index].x, cdt.vertices[v3_index].y, 0.0f);

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

Mesh get_mesh(const CDT::Triangulation<double>& cdt) {
	Mesh mesh;

	for (const auto& vertex : cdt.vertices) {
		glm::vec3 temp = glm::vec3(vertex.x, vertex.y, 0.0f);
		mesh.vertices.push_back(temp);
	}

	auto normals = calculateVertexNormals(cdt);
	for (const auto& normal : normals) {
		glm::vec3 temp = glm::vec3(normal.x, normal.y, normal.z);
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

void inflation(Mesh& mesh, std::vector<glm::vec3> Verts) {
	std::unordered_map<float, float> yToZMap;
	for (const auto& vert : Verts) {
		if (vert.z > 0) {
			yToZMap[vert.y] = vert.z;
		}	
	}

	for (auto& vert : mesh.vertices) {
		auto it = yToZMap.find(vert.y);
		if (it != yToZMap.end()) {
			vert.z = it->second;
		}
		
	}
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
		//CDT::Triangulation<double> cdt;
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

		for (auto vertex : cdt.vertices) {
			std::cout << vertex.x << ", " << vertex.y << std::endl;
		}

		for (auto triangle : cdt.triangles) {
			std::cout << triangle.vertices[0] << ", " << triangle.vertices[1] << ", " << triangle.vertices[2] << std::endl;
		}

		for (auto edge : cdt.fixedEdges) {
			std::cout << edge.v1() << ", " << edge.v2() << std::endl;
		}


		Mesh mesh = get_mesh(cdt);
		inflation(mesh, transformedVerts[1]);
		saveMeshToOBJ(mesh, "C:/Users/dhktj/OneDrive/Desktop/output3.obj");
		
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


struct Face {
	int v1, v2, v3;
};



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
			if (cross_section < 5) cross_section++;
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
