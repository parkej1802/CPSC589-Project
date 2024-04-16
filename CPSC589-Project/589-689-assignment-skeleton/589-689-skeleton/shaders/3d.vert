#version 330 core
layout (location = 0) in vec3 pos;
layout (location = 1) in vec2 uv;
layout (location = 2) in vec3 normal;

uniform mat4 M;
uniform mat4 V;
uniform mat4 P;

out vec3 fragPos;
out vec2 fragUV;
out vec3 n;

void main() {
	fragUV = uv;

	// If you need extra efficiency, you may want to calculate the mat3 on
	// the CPU instead and then upload it as a uniform.
	n = mat3(transpose(inverse(M))) * normal;  

	fragPos = vec3(M * vec4(pos, 1.0));

	gl_Position = P * V * M * vec4(pos, 1.0);
}
