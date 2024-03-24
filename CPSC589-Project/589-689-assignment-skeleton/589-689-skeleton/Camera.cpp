#include "Camera.h"

#define _USE_MATH_DEFINES
#include <math.h>

#include <iostream>

#include "glm/gtc/matrix_transform.hpp"

Camera::Camera(float t, float p, float r) : theta(t), phi(p), radius(r) {
}

glm::mat4 Camera::getView() {
	glm::vec3 eye = radius * glm::vec3(std::cos(theta) * std::sin(phi), std::sin(theta), std::cos(theta) * std::cos(phi));
	glm::vec3 at = glm::vec3(0.0f, 0.0f, 0.0f);
	glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

	return glm::lookAt(eye, at, up);
}

glm::vec3 Camera::getPos() {
	return radius * glm::vec3(std::cos(theta) * std::sin(phi), std::sin(theta), std::cos(theta) * std::cos(phi));
}

void Camera::incrementTheta(float dt) {
	if (theta + (dt / 100.0f) < M_PI_2 && theta + (dt / 100.0f) > -M_PI_2) {
		theta += dt / 100.0f;
	}
}

void Camera::incrementPhi(float dp) {
	phi -= dp / 100.0f;
	if (phi > 2.0 * M_PI) {
		phi -= 2.0 * M_PI;
	}
	else if (phi < 0.0f) {
		phi += 2.0 * M_PI;
	}
}

void Camera::incrementR(float dr) {
	radius -= dr;
}

glm::mat4 Camera::getViewMatrix() {
	// Calculate and return the view matrix based on the camera's current state
	glm::vec3 cameraPosition; // should be set somewhere in the Camera class
	glm::vec3 cameraTarget;   // the point the camera is looking at
	glm::vec3 upVector;       // usually (0, 1, 0) for y-up coordinate systems

	return glm::lookAt(cameraPosition, cameraTarget, upVector);
}

glm::mat4 Camera::getProjectionMatrix() {
	// Calculate and return the projection matrix
	float fov = glm::radians(45.0f); // Example field of view
	float aspectRatio = 800.0f / 600.0f; // Example aspect ratio
	float nearPlane = 0.1f; // Example near clipping plane
	float farPlane = 100.0f; // Example far clipping plane

	return glm::perspective(fov, aspectRatio, nearPlane, farPlane);
}
