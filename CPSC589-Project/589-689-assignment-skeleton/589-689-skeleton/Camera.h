#pragma once

//------------------------------------------------------------------------------
// This file contains an implementation of a spherical camera
//------------------------------------------------------------------------------

//#include <GL/glew.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

class Camera {
public:

	Camera(float t, float p, float r);

	glm::mat4 getView();
	glm::vec3 getPos();
	glm::mat4 getViewMatrix();
	glm::mat4 getProjectionMatrix();
	void incrementTheta(float dt);
	void incrementPhi(float dp);
	void incrementR(float dr);

private:

	float theta;
	float phi;
	float radius;
};
