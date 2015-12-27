/*
 * author: Yupan Liu
 * date: Dec 26, 2015
 * brief: normal, mean curvature and mesh saliency
 * comment: modified by "Anton's OpenGL 4 Tutorials"
 */
#include "lib/maths_funcs.h"
#include "lib/gl_utils.h"
#include <assimp/cimport.h> // C importer
#include <assimp/scene.h> // collects data
#include <assimp/postprocess.h> // various extra operations
#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define GL_LOG_FILE "gl.log"
#define VERTEX_SHADER_FILE "shader/test_vs.glsl"
#define FRAGMENT_SHADER_FILE "shader/test_fs.glsl"
#define MESH_FILE "object/bunny.obj"

// keep track of window size for things like the viewport and the mouse cursor
int g_gl_width = 640;
int g_gl_height = 480;
GLFWwindow* g_window = NULL;

/* load a mesh using the assimp library */
bool load_mesh (const char* file_name, GLuint* vao, int* point_count) {
	const aiScene* scene = aiImportFile (file_name, aiProcess_Triangulate);
	if (!scene) {
		fprintf (stderr, "ERROR: reading mesh %s\n", file_name);
		return false;
	}
   
	/* get first mesh in file only */
	const aiMesh* mesh = scene->mMeshes[0];
	printf ("    %i vertices in mesh[0]\n", mesh->mNumVertices);
	printf ("    %i faces in mesh[0]\n", mesh->mNumFaces);
	
	/* pass back number of vertex points in mesh */
	*point_count = mesh->mNumVertices;
	
	/* generate a VAO, using the pass-by-reference parameter that we give to the
	function */
	glGenVertexArrays (1, vao);
	glBindVertexArray (*vao);
	
	/* we really need to copy out all the data from AssImp's funny little data
	structures into pure contiguous arrays before we copy it into data buffers
	because assimp's texture coordinates are not really contiguous in memory.
	i allocate some dynamic memory to do this. */
	GLfloat* points = NULL; // array of vertex points
	GLfloat* normals = NULL; // array of vertex normals
    if (!mesh->HasPositions()) {
        fprintf(stderr, "ERROR: mesh %s don't have vertex data!\n", file_name);
        return false;
    }

    // Get the mesh's verteices
    points = (GLfloat*)malloc (*point_count * 3 * sizeof (GLfloat));
    for (int i = 0; i < *point_count; i++) {
        const aiVector3D* vp = &(mesh->mVertices[i]);
        points[i * 3] = (GLfloat)vp->x;
        points[i * 3 + 1] = (GLfloat)vp->y;
        points[i * 3 + 2] = (GLfloat)vp->z;
    }

    // Calculate the mesh's normal
    normals = (GLfloat*)malloc(*point_count * 3 * sizeof(GLfloat));
    for(int i = 0; i < mesh->mNumFaces; i++) {
        int idx[3];
        for(int k = 0; k < 3; k++)
            idx[k] = mesh->mFaces[i].mIndices[k];
        // get all vertecies' location
        const aiVector3D* v1 = &(mesh->mVertices[idx[0]]);
        const aiVector3D* v2 = &(mesh->mVertices[idx[1]]);
        const aiVector3D* v3 = &(mesh->mVertices[idx[2]]);
        // vector1
        float x1 = v2->x - v1->x;
        float y1 = v2->y - v1->y;
        float z1 = v2->z - v1->z;
        // vector2
        float x2 = v3->x - v2->x;
        float y2 = v3->y - v2->y;
        float z2 = v3->z - v2->z;
        for(int k = 0; k < 3; k++) {
            normals[idx[k]*3+0] += (GLfloat)(y1*z2-y2*z1),
            normals[idx[k]*3+1] += (GLfloat)(x2*z1-x1*z2),
            normals[idx[k]*3+2] += (GLfloat)(x1*y2-x2*y1);
        }
    }

    for(int i = 0; i < *point_count; i++) {
        float norm = 0.0f;
        for(int k = 0; k < 3; k++)
            norm += normals[i*3+k]*normals[i*3+k];
        for(int k = 0; k < 3; k++)
            normals[i*3+k] /= sqrt(norm);
    }
	
	// Copy all vertecies in mesh data into VBOs
    {
        GLuint vbo;
        glGenBuffers (1, &vbo);
        glBindBuffer (GL_ARRAY_BUFFER, vbo);
        glBufferData (
            GL_ARRAY_BUFFER,
            3 * *point_count * sizeof (GLfloat),
            points,
            GL_STATIC_DRAW
        );
        glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
        glEnableVertexAttribArray (0);
        free (points);
    }

    // Copy all normal vectors in mesh data into VBOs
    {
        GLuint vbo;
        glGenBuffers (1, &vbo);
        glBindBuffer (GL_ARRAY_BUFFER, vbo);
        glBufferData (
            GL_ARRAY_BUFFER,
            3 * *point_count * sizeof (GLfloat),
            normals,
            GL_STATIC_DRAW
        );
        glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
        glEnableVertexAttribArray (1);
        free (normals);
    }

	aiReleaseImport (scene);
	printf ("mesh loaded\n");
	
	return true;
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

GLuint objVAO;
int ObjPointCount = 0;
	
float cam_speed = 1.0f; // 1 unit per second
float cam_yaw_speed = 10.0f; // 10 degrees per second
float cam_pos[] = {0.0f, 0.0f, 5.0f}; // don't start at zero, or we will be too close
float cam_xaw = 0.0f; // x-rotation in degrees
float cam_yaw = 0.0f; // y-rotation in degrees
float cam_zaw = 0.0f; // z-rotation in degrees

int view_mat_location, proj_mat_location;

double xLoc, yLoc;
double dx, dy;

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
    /*
    proj_mat_original = identity_mat4();
    proj_mat_original.m[0] = Sx;
    proj_mat_original.m[5] = Sy;
    proj_mat_original.m[10] = Sz;
    proj_mat_original.m[11] = -1.0f;
    proj_mat_original.m[14] = Pz;
    proj_mat_original.m[15] = 0.0f;
    dx = dy = 0.0f;
    */
	
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

inline void updateView() {
    mat4 T = translate (identity_mat4 (), vec3 (-cam_pos[0], -cam_pos[1], -cam_pos[2])); // cam translation
    mat4 R = rotate_z_deg(identity_mat4(), -cam_zaw)
            * rotate_y_deg(identity_mat4(), -cam_yaw)
            * rotate_x_deg(identity_mat4(), -cam_xaw);
    mat4 view_mat = T * R;
    glUniformMatrix4fv (view_mat_location, 1, GL_FALSE, view_mat.m);
}
void keyBoardEvent(double elapsed_seconds) {
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

void glfw_mouse_pos_callback(GLFWwindow* window, double xNow, double yNow) {
    float EPS = 1e-6;
    float distRatio = 1.0f;
    float angleRatio = 0.3f;
    dx = angleRatio*(xNow-xLoc), dy = angleRatio*(yNow-yLoc);
    //printf("[%f,%f] pos:(%f %f %f) angle:(%f %f %f)\n", dx, dy, cam_pos[0], cam_pos[1], cam_pos[2], cam_xaw, cam_yaw, cam_zaw);
    bool isMoved = (fabs(dx) > EPS) || (fabs(dy) > EPS);
    if(isMoved) {
        cam_xaw += dy;// * cos(dx);//distRatio * cos(dy) * cos(dx);// * cos(dy);
        //cam_yaw += dy * sin(dx);//distRatio * cos(dy) * sin(dx);
        cam_zaw += dx;//distRatio * sin(dy);// * cos(dy);
        updateView();
    }
    xLoc = xNow, yLoc = yNow;
    dx = dy = 0.0;
}

int main () {
	initializeOpenGL();
    glfwSetCursorPosCallback(g_window, glfw_mouse_pos_callback);
	// load the mesh using assimp

	assert (load_mesh (MESH_FILE, &objVAO, &ObjPointCount));
	GLuint shader_programme = create_programme_from_files (
		VERTEX_SHADER_FILE, FRAGMENT_SHADER_FILE
	);
	createShaders(shader_programme);	

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
        //mouseMoveEvent(proj_mat_location);
		keyBoardEvent(elapsed_seconds);
        //printf("[%f]\n", elapsed_seconds);
		// put the stuff we've been drawing onto the display
		glfwSwapBuffers (g_window);
	}
	
	// close GL context and any other GLFW resources
	glfwTerminate();
	return 0;
}
