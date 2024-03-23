#pragma once

#include "GLHandles.h"

#include <glad/glad.h>


class VertexBuffer {

public:
	VertexBuffer(GLuint index, GLint size, GLenum dataType);

	// Because we're using the VertexBufferHandle to do RAII for the buffer for us
	// and our other types are trivial or provide their own RAII
	// we don't have to provide any specialized functions here. Rule of zero
	//
	// https://en.cppreference.com/w/cpp/language/rule_of_three
	// https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#Rc-zero

	// Public interface
	void bind() const { glBindBuffer(GL_ARRAY_BUFFER, bufferID); }
	void uploadData(GLsizeiptr size, const void* data, GLenum usage);

private:
	VertexBufferHandle bufferID;

	// The use of this bool is to only enable or disable "AttribArray" when
	// absolutely required. It assumes that no other areas of the code will
	// be touching the respective "AttribArray" at this index for the VAO
	// that this buffer is part of. If your code breaks that assumption,
	// then you may want to remove this bool and work something else out.
	bool attribArrayEnabled;

	GLuint attribIndex;
	GLint attribSize;
	GLenum attribDataType;
};

