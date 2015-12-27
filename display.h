/*
 * author: Yupan Liu
 * date: Dec 27, 2015
 * brief: normal, mean curvature and mesh saliency
 * comment: modified by "Anton's OpenGL 4 Tutorials"
 */

#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library
#include <math.h>
#include <assert.h>
#include "lib/maths_funcs.h"
#include "lib/gl_utils.h"

#define _USE_MATH_DEFINES
#define GL_LOG_FILE "gl.log"
#define VERTEX_SHADER_FILE "shader/test_vs.glsl"
#define FRAGMENT_SHADER_FILE "shader/test_fs.glsl"
#define MESH_FILE "object/cow.obj"
#define oo 88888888.0f

extern float fMin(float x, float y);
extern float fMax(float x, float y);

extern GLfloat* normals; // array of vertex normals
extern GLfloat* meanCurvature;
extern GLfloat* smoothSaliency;

extern int vertexCnt;

extern GLuint objVAO;
extern int ObjPointCount;
	
extern float cam_speed;
extern float cam_yaw_speed;
extern float cam_pos[3];
extern float cam_xaw;
extern float cam_yaw;
extern float cam_zaw;

extern int view_mat_location;
extern int proj_mat_location;

extern double xLoc, yLoc;
extern double dx, dy;

// Update display mode depends on keyboard state
void updateDisplayType(int type);

// Check the pixel's color whether is legal or not
inline void checkLegal(float &x);

// Get RGB value from the YUV value
void getYUVtoRGB(bool flag, float Y, float U, float V, GLfloat& R, GLfloat& G, GLfloat& B);

void initializeOpenGL();

void createShaders(GLuint shader_programme);

inline void updateView();

void keyBoardEvent(double elapsed_seconds);

void glfw_mouse_pos_callback(GLFWwindow* window, double xNow, double yNow);
