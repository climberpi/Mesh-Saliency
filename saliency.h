/*
 * author: Yupan Liu
 * date: Dec 26, 2015
 * brief: normal, mean curvature and mesh saliency
 */

#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library
#include <assimp/cimport.h> // C importer
#include <assimp/scene.h> // collects data
#include <assimp/postprocess.h> // various extra operations
#include "lib/maths_funcs.h"
#include "QSlim.h"
#include <queue>
#include <iostream>
using namespace std;

#define _USE_MATH_DEFINES

extern GLfloat* normals; // array of vertex normals
extern GLfloat* meanCurvature;
extern GLfloat* smoothSaliency;
extern GLfloat* points;

// Load a mesh using the assimp library, and calculate its mesh saliency
bool load_mesh (const char* file_name, GLuint* vao, int* point_count);

extern int vertexCnt;

extern float fMin(float x, float y);
extern float fMax(float x, float y);
extern void updateDisplayType(int type);
