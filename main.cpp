/*
 * author: Yupan Liu
 * date: Dec 26, 2015
 * brief: normal, mean curvature and mesh saliency
 * comment: modified by "Anton's OpenGL 4 Tutorials"
 */
#include <stdio.h>
#include <stdlib.h>
#include "saliency.h"
#include "display.h"
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
using namespace std;

GLfloat* normals = NULL; // array of vertex normals
GLfloat* meanCurvature = NULL;
GLfloat* smoothSaliency = NULL;

// keep track of window size for things like the viewport and the mouse cursor
int g_gl_width = 640;
int g_gl_height = 480;
GLFWwindow* g_window = NULL;

float cam_speed = 1.0f; // 1 unit per second
float cam_yaw_speed = 10.0f; // 10 degrees per second
float cam_pos[3] = {0.0f, 0.0f, 5.0f}; // don't start at zero, or we will be too close
float cam_xaw = 0.0f; // x-rotation in degrees
float cam_yaw = 0.0f; // y-rotation in degrees
float cam_zaw = 0.0f; // z-rotation in degrees

int view_mat_location, proj_mat_location;

double xLoc, yLoc;
double dx, dy;

int vertexCnt = 0;

GLuint objVAO;
int ObjPointCount = 0;

float fMin(float x, float y) {return x < y ? x : y;}
float fMax(float x, float y) {return x > y ? x : y;}

int main (int argc, char* argv[]) {
	initializeOpenGL();
    glfwSetCursorPosCallback(g_window, glfw_mouse_pos_callback);
	// load the mesh using assimp

    assert (argc >= 2);
	assert (load_mesh (argv[1], &objVAO, &ObjPointCount));
	GLuint shader_programme = create_programme_from_files (
		VERTEX_SHADER_FILE, FRAGMENT_SHADER_FILE
	);
	createShaders(shader_programme);	

    string filename(argv[1]);
    
    // Output the saliency file
    filename += string(".saliency");
    ofstream ofs;
    ofs.open(filename.c_str());
    vector<float> saliencys;
    for(int i = 0; i < ObjPointCount; i++)
        saliencys.push_back(smoothSaliency[i]);
    sort(saliencys.begin(), saliencys.end());
    float threshold = saliencys[(int)(ObjPointCount * 0.7f)];
    for(int i = 0; i < ObjPointCount; i++)
        if(smoothSaliency[i] >= threshold)
            ofs << 1000000.0f * smoothSaliency[i] << endl;
        else 
            ofs << 10000.0f * smoothSaliency[i] << endl;
    ofs.close();

	while (!glfwWindowShouldClose (g_window)) {
		static double previous_seconds = glfwGetTime ();
		double current_seconds = glfwGetTime ();
		double elapsed_seconds = current_seconds - previous_seconds;
		previous_seconds = current_seconds;
	
		_update_fps_counter (g_window);
		// wipe the drawing surface clear
		glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	    glViewport (0, 0, g_gl_width, g_gl_height);

		glUseProgram (shader_programme);
		glBindVertexArray (objVAO);
		glDrawArrays (GL_TRIANGLES, 0, ObjPointCount);
		// update other events like input handling 
		glfwPollEvents ();
		// control keys
		keyBoardEvent(elapsed_seconds);
		// put the stuff we've been drawing onto the display
		glfwSwapBuffers (g_window);
	}
	
	// close GL context and any other GLFW resources
	glfwTerminate();
	return 0;
}
