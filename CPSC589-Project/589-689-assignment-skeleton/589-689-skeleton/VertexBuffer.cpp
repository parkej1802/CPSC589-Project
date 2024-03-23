#include "VertexBuffer.h"

#include <utility>


VertexBuffer::VertexBuffer(GLuint index, GLint size, GLenum dataType)
	: bufferID{}
	, attribArrayEnabled(false)
	, attribIndex(index)
	, attribSize(size)
	, attribDataType(dataType)
{}


void VertexBuffer::uploadData(GLsizeiptr size, const void* data, GLenum usage) {
	if (size > 0) {
		bind();
		glBufferData(GL_ARRAY_BUFFER, size, data, usage);

		// If we have data and did not yet set up and enable the AttribArray
		// with this vertex buffer, then we do so now.
		if (!attribArrayEnabled) {
			glVertexAttribPointer(attribIndex, attribSize, attribDataType, GL_FALSE, 0, (void*)0);
			glEnableVertexAttribArray(attribIndex);
			attribArrayEnabled = true;
		}
	}
	else if (attribArrayEnabled) {
		// If there's no data, we disable use of an array for this attrib
		// and instead set a constant default value.
		// OpenGL sets the default to (0, 0, 0, 1), or a portion of this
		// if the size is less than 4. To set your own default values,
		// look at:
		// https://registry.khronos.org/OpenGL-Refpages/gl4/html/glVertexAttrib.xhtml
		// A good, longer description of default behaviour is given at:
		// https://stackoverflow.com/questions/51721058/confusion-about-glvertexattrib-functions
		glDisableVertexAttribArray(attribIndex);
		attribArrayEnabled = false;
	}

}
