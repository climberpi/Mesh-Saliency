/*
 * author: Yupan Liu
 * date: Dec 27, 2015
 * brief: normal, mean curvature and mesh saliency
 * comment: modified by "Anton's OpenGL 4 Tutorials"
 */
#include "display.h"

inline void checkLegal(float &x) {
    if(x > 0.9999f) x = 0.9999f;
    if(x < -0.9999f) x = -0.9999f;
}

// Change the YUV to RGB, and guarantee all R,G,B values in [-1, 1].
void getYUVtoRGB(bool flag, float Y, float U, float V, GLfloat& R, GLfloat& G, GLfloat& B) {
    R = (GLfloat)(Y + 1.4075 * (V - 128.0));
    G = (GLfloat)(Y - 0.3455 * (U - 128.0) - (0.7169 * (V - 128.0f)));
    B = (GLfloat)(Y + 1.7790 * (U - 128.0));
    if(!flag)
        R = 255.0f - R,
        B *= 0.8f;
    R = (R-127.5f)/255.0f;
    G = (G-127.5f)/255.0f;
    B = (B-127.5f)/255.0f;
    checkLegal(R), checkLegal(G), checkLegal(B);   
    //printf("$(%f %f %f)\n", R, G, B);
}

// Update display mode depends on keyboard state
void updateDisplayType(int type) {
    if(type != 1 && type != 2 && type != 3 && type != 4 && type != 5) {
        fprintf(stderr, "Wrong Display Type!\n");
        return;
    }

    int vertexCntHere = (type >= 4 ? simplifiedVertexCnt : vertexCnt);
    // Copy all vertecies in mesh data into VBOs
    {
        GLuint vbo;
        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData (
            GL_ARRAY_BUFFER,
            3 * vertexCntHere * sizeof (GLfloat),
            type >= 4 ? simplifiedPoints : points,
            GL_STATIC_DRAW
        );
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
        glEnableVertexAttribArray(0);
    }
    // Copy all vertices' color values in mesh data into VBOs
    {
        GLuint vbo;
        glGenBuffers (1, &vbo);
        glBindBuffer (GL_ARRAY_BUFFER, vbo);
        float* colors = NULL;
        colors = (float*)malloc(vertexCntHere * 3 * sizeof(float));
        switch(type) {
            // Normals mode
            case 1: case 4:{
                        for(int i = 0; i < vertexCntHere; i++)
                            colors[3*i+0] = normals[3*i+0],
                            colors[3*i+1] = normals[3*i+1],
                            colors[3*i+2] = normals[3*i+2];
                        break;
                    }
            // Mean curvature mode
            case 2: {
                        float xMin = oo, xMax = -oo;
                        for(int i = 0; i < vertexCntHere; i++)
                            xMin = fMin(xMin, fabs(meanCurvature[i])),
                            xMax = fMax(xMax, fabs(meanCurvature[i]));
                        // Normalize the Y value
                        for(int i = 0; i < vertexCntHere; i++) {
                            float Y = 255.0f*((log(1e-8+fabs(meanCurvature[i]))-log(1e-8+xMin))
                                    / (log(1e-8+xMax)-log(1e-8+xMin)));
                            getYUVtoRGB(meanCurvature[i]>0.0, Y, 255.0f, 255.0f, colors[3*i], colors[3*i+1], colors[3*i+2]);
                        }
                        break;
                    }
            // Mesh saliency mode
            case 3: {
                        float xMin = oo, xMax = -oo;
                        for(int i = 0; i < vertexCntHere; i++)
                            xMin = fMin(xMin, fabs(smoothSaliency[i])),
                            xMax = fMax(xMax, fabs(smoothSaliency[i]));
                        // Normalize the Y value
                        for(int i = 0; i < vertexCntHere; i++) {
                            float Y = 255.0f*(log(1e-8+smoothSaliency[i])-log(1e-8+xMin))
                                    / (log(1e-8+xMax)-log(1e-8+xMin));
                            getYUVtoRGB(smoothSaliency[i]>0, Y, 255.0f, 255.0f, colors[3*i], colors[3*i+1], colors[3*i+2]);
                        }
                        break;
                    }
        }
        glBufferData (
            GL_ARRAY_BUFFER,
            3 * vertexCntHere * sizeof (GLfloat),
            colors,
            GL_STATIC_DRAW
        );
        glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
        glEnableVertexAttribArray (1);
    }
}

void initializeOpenGL() {
	assert (restart_gl_log ());
	assert (start_gl ());
	glEnable (GL_DEPTH_TEST); // enable depth-testing
	glDepthFunc (GL_LESS); // depth-testing interprets a smaller value as "closer"
	glEnable (GL_CULL_FACE); // cull face
	glCullFace (GL_BACK); // cull back face
	glFrontFace (GL_CCW); // set counter-clock-wise vertex order to mean the front
	glClearColor (0.2, 0.2, 0.2, 1.0); // grey background to help spot mistakes
	glViewport (0, 0, g_gl_width, g_gl_height);
}

void createShaders(GLuint shader_programme) {	
	#define ONE_DEG_IN_RAD (2.0 * M_PI) / 360.0 // 0.017444444
	// input variables
	float near = 0.1f; // clipping plane
	float far = 100.0f; // clipping plane
	float fov = 67.0f * ONE_DEG_IN_RAD; // convert 67 degrees to radians
	float aspect = (float)g_gl_width / (float)g_gl_height; // aspect ratio
	// matrix components
	float range = tan (fov * 0.5f) * near;
	float Sx = (2.0f * near) / (range * aspect + range * aspect);
	float Sy = near / range;
	float Sz = -(far + near) / (far - near);
    float Pz = -(2.0f * far * near) / (far - near);
	GLfloat proj_mat[] = {
		Sx, 0.0f, 0.0f, 0.0f,
		0.0f, Sy, 0.0f, 0.0f,
		0.0f, 0.0f, Sz, -1.0f,
		0.0f, 0.0f, Pz, 0.0f
	};
	
    mat4 T = translate (identity_mat4 (), vec3 (-cam_pos[0], -cam_pos[1], -cam_pos[2]));
    mat4 R = rotate_z_deg (identity_mat4 (), -cam_zaw);
    mat4 view_mat = T * R;
	view_mat_location = glGetUniformLocation (shader_programme, "view");
	glUseProgram (shader_programme);
	glUniformMatrix4fv (view_mat_location, 1, GL_FALSE, view_mat.m);
	proj_mat_location = glGetUniformLocation (shader_programme, "proj");
	glUseProgram (shader_programme);
	glUniformMatrix4fv (proj_mat_location, 1, GL_FALSE, proj_mat);
}

// Maintain translations and rotations.
inline void updateView() {
    mat4 T = translate (identity_mat4 (), vec3 (-cam_pos[0], -cam_pos[1], -cam_pos[2])); // cam translation
    mat4 R = rotate_z_deg(identity_mat4(), -cam_zaw)
            * rotate_y_deg(identity_mat4(), -cam_yaw)
            * rotate_x_deg(identity_mat4(), -cam_xaw);
    mat4 view_mat = T * R;
    glUniformMatrix4fv (view_mat_location, 1, GL_FALSE, view_mat.m);
}

void keyBoardEvent(double elapsed_seconds) {
    // Update the display mode.
    if (glfwGetKey(g_window, GLFW_KEY_1))
        updateDisplayType(1);
    else if (glfwGetKey(g_window, GLFW_KEY_2))
        updateDisplayType(2);
    else if (glfwGetKey(g_window, GLFW_KEY_3))
        updateDisplayType(3);
    else if (glfwGetKey(g_window, GLFW_KEY_4))
        updateDisplayType(4);
    // Update the object's location and angle.
    bool cam_moved = false;
    if (glfwGetKey (g_window, GLFW_KEY_A)) {
        cam_pos[0] -= cam_speed * elapsed_seconds;
        cam_moved = true;
    }
    if (glfwGetKey (g_window, GLFW_KEY_D)) {
        cam_pos[0] += cam_speed * elapsed_seconds;
        cam_moved = true;
    }
    if (glfwGetKey (g_window, GLFW_KEY_Z)){
        cam_pos[1] += cam_speed * elapsed_seconds;
        cam_moved = true;
    }
    if (glfwGetKey (g_window, GLFW_KEY_C)) {
        cam_pos[1] -= cam_speed * elapsed_seconds;
        cam_moved = true;
    }
    if (glfwGetKey (g_window, GLFW_KEY_W)) {
        cam_pos[2] -= cam_speed * elapsed_seconds;
        cam_moved = true;
    }
    if (glfwGetKey (g_window, GLFW_KEY_S)) {
        cam_pos[2] += cam_speed * elapsed_seconds;
        cam_moved = true;
    }
    if (glfwGetKey (g_window, GLFW_KEY_LEFT)) {
        cam_zaw += cam_yaw_speed * elapsed_seconds;
        cam_moved = true;
    }
    if (glfwGetKey (g_window, GLFW_KEY_RIGHT)) {
        cam_zaw -= cam_yaw_speed * elapsed_seconds;
        cam_moved = true;
    }
    // update view matrix
    if (cam_moved) updateView(); 
    
    if (GLFW_PRESS == glfwGetKey (g_window, GLFW_KEY_ESCAPE)) {
        glfwSetWindowShouldClose (g_window, 1);
    }
}

void mouseEventHandler(GLFWwindow* window, double xNow, double yNow) {
    float EPS = 1e-6;
    float distRatio = 1.0f;
    float angleRatio = 0.3f;
    dx = angleRatio*(xNow-xLoc), dy = angleRatio*(yNow-yLoc);
    //printf("[%f,%f] pos:(%f %f %f) angle:(%f %f %f)\n", dx, dy, cam_pos[0], cam_pos[1], cam_pos[2], cam_xaw, cam_yaw, cam_zaw);
    bool isMoved = (fabs(dx) > EPS) || (fabs(dy) > EPS);
    if(isMoved) {
        cam_xaw += distRatio * dy;// * cos(dx);//distRatio * cos(dy) * cos(dx);// * cos(dy);
        //cam_yaw += dy * sin(dx);//distRatio * cos(dy) * sin(dx);
        cam_zaw += distRatio * dx;//distRatio * sin(dy);// * cos(dy);
        updateView();
    }
    xLoc = xNow, yLoc = yNow;
    dx = dy = 0.0;
}

