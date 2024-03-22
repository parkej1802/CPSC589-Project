#include <iostream>

// Window.h `#include`s ImGui, GLFW, and glad in correct order.
#include "Window.h"

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include <vector>


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


void Bspline(CPU_Geometry& controlPointcpu, GPU_Geometry& controlPointgpu, std::vector<glm::vec3> E, int m, int k, float ui) {
	controlPointcpu.verts.clear();
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
		controlPointcpu.verts.push_back(c[0]);
		controlPointcpu.cols.push_back(glm::vec3(1.f, 1.f, 1.f));
	}

	if (!controlPointcpu.verts.empty()) {
		controlPointcpu.verts.push_back(E.back());
		controlPointcpu.cols.push_back(glm::vec3(1.f, 1.f, 1.f));

	}

	controlPointgpu.setVerts(controlPointcpu.verts);
	controlPointgpu.setCols(controlPointcpu.cols);
}


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
	{}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (key == GLFW_KEY_R && action == GLFW_PRESS) {
			shader.recompile();
		}

		if (key == GLFW_KEY_RIGHT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			x_angle += 0.01;
		}

		if (key == GLFW_KEY_LEFT && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			x_angle -= 0.01;
		}

		if (key == GLFW_KEY_UP && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			y_angle += 0.01;
		}

		if (key == GLFW_KEY_DOWN && (action == GLFW_PRESS || action == GLFW_REPEAT)) {
			y_angle -= 0.01;
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
	}

	// Updates the screen width and height, in screen coordinates
	// (not necessarily the same as pixels)
	virtual void windowSizeCallback(int width, int height) {
		screenWidth = width;
		screenHeight = height;
	}

	// Sets the new cursor position, in screen coordinates
	virtual void cursorPosCallback(double xpos, double ypos) {
		screenMouseX = xpos;
		screenMouseY = ypos;
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

private:
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
};


std::vector<glm::vec3> get_control_points(std::vector<glm::vec3> &line,int smoothness) {
	std::vector<glm::vec3> result;

	int step_size = line.size() / smoothness;
	int index = 0;

	while (index < line.size()) {
		result.push_back(line[index]);
		index += step_size;
	}
	return result; 
}

std::vector<glm::vec3> flattenLineVerts(std::vector<std::vector<glm::vec3>> &lineVerts) {
	std::vector<glm::vec3> flattenedVerts;

	for (const auto& vertGroup : lineVerts) {
		flattenedVerts.insert(flattenedVerts.end(), vertGroup.begin(), vertGroup.end());
	}

	return flattenedVerts;
}


void draw_cross_sections(std::shared_ptr<MyCallbacks> &cb,
	std::vector<std::vector<glm::vec3>> &lineVerts,
	CPU_Geometry &cpuGeom,
	GPU_Geometry &gpuGeom,
	int &cross_section) {

	cpuGeom.verts.clear();
	cpuGeom.cols.clear();

	if (cb->leftMouseActive()) {
		lineVerts[cross_section].push_back(glm::vec3(cb->getCursorPosGL(), 0.f));
		cpuGeom.verts = flattenLineVerts(lineVerts);

		for (int i = 0; i < 3; i++) {
			glm::vec3 color;
			if (i == 0) {
				color = glm::vec3(1, 0, 0);
			}
			else if (i == 1) {
				color = glm::vec3(0, 1, 0);
			}
			else {
				color = glm::vec3(0, 0, 1);
			}
			for (int j = 0; j < lineVerts[i].size(); j++) {
				cpuGeom.cols.push_back(color);
			}
		}
		gpuGeom.setVerts(cpuGeom.verts);
		gpuGeom.setCols(cpuGeom.cols);
	}
}

void draw_lines(std::vector<std::vector<glm::vec3>>& lineVerts) {
	int start_index = 0;

	for (int i = 0; i < 3; i++) {
		glDrawArrays(GL_LINE_STRIP, start_index, GLsizei(lineVerts[i].size()));
		start_index += lineVerts[i].size();
	}
}

void combine(std::vector<std::vector<glm::vec3>> &lineVerts,
	CPU_Geometry &cpuGeom,
	GPU_Geometry &gpuGeom,
	std::shared_ptr<MyCallbacks> &cb){
		cpuGeom.verts.clear();
		cpuGeom.cols.clear();

		float x_angle = cb -> getXAngle(); 
		glm::mat3 X_R = glm::mat3(cos(-x_angle),0, sin(-x_angle),
								0,1,0,
								-sin(-x_angle),0,cos(-x_angle));

		float y_angle = cb->getYAngle();
		glm::mat3 Y_R = glm::mat3(1, 0,0,
								0, cos(y_angle), -sin(y_angle),
								0, sin(y_angle), cos(y_angle));

		std::vector<glm::vec3> flattenedVerts;
		for (int i=0; i <3; i++) {
			for(int j=0;j<lineVerts[i].size();j++){
				if(i == 0){
					glm::vec3 coord = lineVerts[i][j] * X_R;
					coord = coord * Y_R;
					flattenedVerts.push_back(coord);
				}else if(i == 1){
					float x = lineVerts[i][j].x;
					float y = lineVerts[i][j].y;
					glm::vec3 side = glm::vec3(0.f,y,x);
					side = side * X_R;
					side = side * Y_R;
					flattenedVerts.push_back(side);
				}else{
					float x = lineVerts[i][j].x;
					float y = lineVerts[i][j].y;
					glm::vec3 top = glm::vec3(x,0.f,y);
					top = top * X_R;
					top = top * Y_R;
					flattenedVerts.push_back(top);
				}
			}
		}

		cpuGeom.verts = flattenedVerts;

		for (int i = 0; i < 3; i++) {
			glm::vec3 color;
			if (i == 0) {
				color = glm::vec3(1, 0, 0);
			}
			else if (i == 1) {
				color = glm::vec3(0, 1, 0);
			}
			else {
				color = glm::vec3(0, 0, 1);
			}
			for (int j = 0; j < lineVerts[i].size(); j++) {
				cpuGeom.cols.push_back(color);
			}
		}
		gpuGeom.setVerts(cpuGeom.verts);
		gpuGeom.setCols(cpuGeom.cols);

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

	// GEOMETRY
	CPU_Geometry cpuGeom;
	GPU_Geometry gpuGeom;

	CPU_Geometry controlPointcpu;
	GPU_Geometry controlPointgpu;

	// Variables that ImGui will alter.
	float pointSize = 10.0f; // Diameter of drawn points
	float color[3] = { 1.f, 0.f, 0.f }; // Color of new points
	bool drawLines = true; // Whether to draw connecting lines
	int selectedPointIndex = -1; // Used for point dragging & deletion

	int cross_section = -1; 

	std::vector<std::vector<glm::vec3>> lineVerts(3);
	// RENDER LOOP
	while (!window.shouldClose()) {

		// Tell callbacks object a new frame's begun BEFORE polling events!
		cb->incrementFrameCount();
		glfwPollEvents();

		// If mouse just went down, see if it was on a point.
		if (cb->leftMouseJustPressed() || cb->rightMouseJustPressed()) {
			// We use the point DIAMETER as the threshold, meaning the user
			// can click anywhere within 2x radius to select.
			// You may want to change that.
			float threshold = pointSize;

			selectedPointIndex = cb->indexOfPointAtCursorPos(cpuGeom.verts, threshold);
		}

		/* CROSS SECTION */
		// when the left button gets pressed increase the cross_section
		if (cb->leftMouseJustPressed()) {
			if (cross_section < 3) cross_section++;
		}

		if(cross_section < 3){
			draw_cross_sections(cb, lineVerts, cpuGeom, gpuGeom, cross_section);
		}else{
			combine(lineVerts, cpuGeom, gpuGeom, cb);
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

		glPointSize(pointSize);

		glEnable(GL_FRAMEBUFFER_SRGB);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//if (drawLines) glDrawArrays(GL_LINE_STRIP, 0, GLsizei(cpuGeom.verts.size()));

		//std::vector<glm::vec3> temp;

		/*
		if (!lineVerts.empty()) {
			controlPointcpu.verts = get_control_points(lineVerts[0], 30);

			 if (drawLines) {
				 if (!controlPointcpu.verts.empty()) {
					 gpuGeom.setVerts(controlPointcpu.verts);
					 glDrawArrays(GL_POINTS, 0, GLsizei(controlPointcpu.verts.size()));

				 }
			 }
		}
		*/

		// if (drawLines) {
		// 	for (const auto& line : lineVerts) {
		// 		gpuGeom.setVerts(line);
		// 		glDrawArrays(GL_LINE_STRIP, 0, GLsizei(line.size()));
		// 	}
		// }

		// if (drawLines) {
		// 	if (!cpuGeom.verts.empty()) {
		// 		gpuGeom.setVerts(cpuGeom.verts);
		// 		glDrawArrays(GL_LINE_STRIP, 0, GLsizei(cpuGeom.verts.size()));
		// 	}
		// }
		
		//m = controlPointcpu.verts.size() - 1;

		//Bspline(controlPointcpu, controlPointgpu, controlPointcpu.verts, m, k, ui);

		// controlPointgpu.bind();
		// glDrawArrays(GL_LINE_STRIP, 0, GLsizei(controlPointcpu.verts.size()));

		// glDrawArrays(GL_POINTS, 0, GLsizei(cpuGeom.verts.size()));

		draw_lines(lineVerts);

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
