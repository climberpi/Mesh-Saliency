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
#include <queue>
#include <iostream>
using namespace std;
#define GL_LOG_FILE "gl.log"
#define VERTEX_SHADER_FILE "shader/test_vs.glsl"
#define FRAGMENT_SHADER_FILE "shader/test_fs.glsl"
#define MESH_FILE "object/bunny.obj"
#define oo 88888888.0f

// keep track of window size for things like the viewport and the mouse cursor
int g_gl_width = 640;
int g_gl_height = 480;
GLFWwindow* g_window = NULL;

float fmin(float x, float y) {return x < y ? x : y;}
float fmax(float x, float y) {return x > y ? x : y;}

GLfloat* normals = NULL; // array of vertex normals
GLfloat* meanCurvature = NULL;
GLfloat* smoothSaliency = NULL;
int vertexCnt = 0.0f;

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
    vertexCnt = *point_count;
	
	/* generate a VAO, using the pass-by-reference parameter that we give to the
	function */
	glGenVertexArrays (1, vao);
	glBindVertexArray (*vao);
	
	/* we really need to copy out all the data from AssImp's funny little data
	structures into pure contiguous arrays before we copy it into data buffers
	because assimp's texture coordinates are not really contiguous in memory.
	i allocate some dynamic memory to do this. */
	GLfloat* points = NULL; // array of vertex points
    if (!mesh->HasPositions()) {
        fprintf(stderr, "ERROR: mesh %s don't have vertex data!\n", file_name);
        return false;
    }

    float xMin = oo, yMin = oo, zMin = oo;
    float xMax = -oo, yMax = -oo, zMax = -oo;
    // Get the mesh's verteices
    points = (GLfloat*)malloc (*point_count * 3 * sizeof (GLfloat));
    for (int i = 0; i < *point_count; i++) {
        const aiVector3D* vp = &(mesh->mVertices[i]);
        points[i * 3] = (GLfloat)vp->x;
        points[i * 3 + 1] = (GLfloat)vp->y;
        points[i * 3 + 2] = (GLfloat)vp->z;
        xMin = fmin(xMin, vp->x), xMax = fmax(xMax, vp->x);
        yMin = fmin(yMin, vp->y), yMax = fmax(yMax, vp->y);
        zMin = fmin(zMin, vp->z), zMax = fmax(zMax, vp->z);
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
        // vectors
        vec3 faceVec1 = vec3(v2->x - v1->x, v2->y - v1->y, v2->z - v1->z);
        vec3 faceVec2 = vec3(v3->x - v2->x, v3->y - v2->y, v3->z - v2->z);
        for(int k = 0; k < 3; k++) {
            vec3 crossProd = cross(faceVec1, faceVec2);
            normals[idx[k]*3+0] += (GLfloat)crossProd.v[0],
            normals[idx[k]*3+1] += (GLfloat)crossProd.v[1],
            normals[idx[k]*3+2] += (GLfloat)crossProd.v[2];
         }
    }

    for(int i = 0; i < *point_count; i++) {
        float norm = 0.0f;
        for(int k = 0; k < 3; k++)
            norm += normals[i*3+k]*normals[i*3+k];
        for(int k = 0; k < 3; k++)
            normals[i*3+k] /= sqrt(norm);
    }

    // Calculate each vertecies' shape operator
    mat3* shapeOperators = NULL;
    float* vertexArea = NULL;
    shapeOperators = (mat3*)malloc(*point_count * sizeof(mat3));
    vertexArea = (float*)malloc(*point_count * sizeof(float));
    for(int i = 0; i < *point_count; i++) {
        vertexArea[i] = 0.0f;
        for(int j = 0; j < 9; j++)
            shapeOperators[i].m[j] = 0.0f;
    }
    for(int k = 0; k < mesh->mNumFaces; k++) {
        // Calculate the face's area
        aiVector3D* aiVec[3];
        for(int idx = 0; idx < 3; idx++)
            aiVec[idx] = &(mesh->mVertices[mesh->mFaces[k].mIndices[idx]]);
        vec3 faceVec1 = vec3(aiVec[1]->x - aiVec[0]->x, 
                             aiVec[1]->y - aiVec[0]->y, 
                             aiVec[1]->z - aiVec[0]->z);
        vec3 faceVec2 = vec3(aiVec[2]->x - aiVec[1]->x,
                             aiVec[2]->y - aiVec[1]->y,
                             aiVec[2]->z - aiVec[1]->z);
        vec3 vecArea = cross(faceVec1, faceVec2);
        float faceArea = sqrt(vecArea.v[0]*vecArea.v[0] + vecArea.v[1]*vecArea.v[1] + vecArea.v[2]*vecArea.v[2]);

        for(int idx = 0; idx < 3; idx++) {
            int i = mesh->mFaces[k].mIndices[idx];
            int j = mesh->mFaces[k].mIndices[(idx+1)%3];
            // Get vertex i and j's normal vectors.
            vec3 Ni = vec3(normals[i*3], normals[i*3+1], normals[i*3+2]);
            vec3 Nj = vec3(normals[j*3], normals[j*3+1], normals[j*3+2]);
            // Get vertex i and j's location.
            const aiVector3D* aiVi = &(mesh->mVertices[i]);
            const aiVector3D* aiVj = &(mesh->mVertices[j]);
            vec3 Vi = vec3(aiVi->x, aiVi->y, aiVi->z);
            vec3 Vj = vec3(aiVj->x, aiVj->y, aiVj->z);
            
            // For vertex i, update the relative part of its shape operator
            vec3 Tij = (identity_mat3 - wedge(Ni, Ni)*(Vi-Vj);
            normalise(Tij);
            float kappa_ij = 2*dot(Ni, Vj-Vi);
            kappa_ij /= get_squared_dist(Vi, Vj);
            // Maintain vi's shape operator
            shapeOperators[i] = shapeOperators[i] + (wedge(Tij, Tij) * (kappa_ij * faceArea));
            vertexArea[i] += faceArea;

            // For vertex j, update the relative part of its shape operator
            vec3 Tji = (identity_mat3 - wedge(Nj, Nj)*(Vj-Vi);
            normalise(Tji);
            float kappa_ji = 2*dot(Nj, Vi-Vj);
            kappa_ji /= get_squared_dist(Vi, Vj);
            // Maintain vj's shape operator
            shapeOperators[j] = shapeOperators[j] + (wedge(Tji, Tji) * (kappa_ji * faceArea));
            vertexArea[j] += faceArea;
        }
    }

    for(int i = 0; i < *point_count; i++)
        shapeOperators[i] = shapeOperators[i] * (1.0f/vertexArea[i]);
    free(vertexArea);

    // Diagonalize the shape operator, and get the mean curvature
    meanCurvature = (float*)malloc(*point_count * sizeof(float));
    for(int k = 0; k < *point_count; k++) {
        vec3 E1 = vec3(1.0f, 0.0f, 0.0f);
        vec3 Nk = vec3(normals[k*3], normals[k*3+1], normals[k*3+2]);
        bool isMinus = get_squared_dist(E1, Nk) > get_squared_dist(E1, -1.0f*Nk);
        vec3 Wk;
        // Diagnoalization by the Householder transform
        if (!isMinus)
            Wk = E1 + Nk;
        else
            Wk = E1 - Nk;
        normalise(Wk);
        mat3 Qk = identity_mat3() - 2.0f * wedge(Wk, Wk);
        mat3 Mk = transpose(Qk) * shapeOperators[k] * Qk;
        // Calculate the mean curvature by M_k's trace;
        meanCurvature[k] = (GLfloat)(Mk.m[4] + Mk.m[8]);
    }
    free(shapeOperators);

    // Calculate the incident matrix ( as linked list )
    int* first = NULL;
    int* next = NULL;
    int* incidentVertex = NULL;
    first = (int*)malloc(*point_count * sizeof(int));
    memset(first, -1, sizeof(*point_count * sizeof(int)));
    next = (int*)malloc(mesh->mNumFaces * 6 * sizeof(int));
    incidentVertex = (int*)malloc(mesh->mNumFaces * 6 * sizeof(int));
    int edgeCnt = 0;
    for(int k = 0; k < mesh->mNumFaces; k++) {
        int idx[3];
        for(int i = 0; i < 3; i++)
            idx[i] = mesh->mFaces[k].mIndices[i];
        for(int i = 0; i < 3; i++) {
            int j1 = idx[(i+1)%3], j2 = idx[(i+2)%3];
            incidentVertex[++edgeCnt] = j1;
            next[edgeCnt] = first[idx[i]]; first[idx[i]] = edgeCnt;
            incidentVertex[++edgeCnt] = j2;
            next[edgeCnt] = first[idx[i]]; first[idx[i]] = edgeCnt;
        }
    }
    
    // Calculate the mesh saliency by BFS
    float diagonalLength = sqrt((xMax-xMin)*(xMax-xMin) + (yMax-yMin)*(yMax-yMin) + (zMax-zMin)*(zMax-zMin));
    float sigma = 0.003 * diagonalLength;
    float* saliency[7];
    float maxSaliency[7];
    for(int i = 2; i <= 6; i++) {
        saliency[i] = NULL;;
        saliency[i] = (float*)malloc(*point_count * sizeof(float));
        maxSaliency[i] = -oo;
    }
    // Labeled the vertecies whether covered or not.
    bool* used = NULL;
    used = (bool*)malloc(*point_count * sizeof(bool));
    for(int k = 0; k < *point_count; k++) {
        // Initialize the saliency and its local counter.
        for(int i = 2; i <= 6; i++)
            saliency[i][k] = 0.0f;
        // Initialize the saliency's Gaussian filter.
        float gaussianSigma1[7], gaussianSigma2[7], sumSigma1[7], sumSigma2[7];
        for(int i = 0; i < 7; i++)
            gaussianSigma1[i] = gaussianSigma2[i] = 0.0f,
            sumSigma1[i] = sumSigma2[i] = 0.0f;
        // Get the current vertex's information.
        aiVector3D* aiVec = &(mesh->mVertices[k]);
        vec3 vVec = vec3(aiVec->x, aiVec->y, aiVec->z);
        // Initialize the queue to find neighbourhood.
        queue<int> Q;
        memset(used, 0, *point_count * sizeof(bool));
        Q.push(k);
        used[k] = true;
        // Frsit BFS
        while(!Q.empty()) {
            // Get the front element in the queue.
            int idx = Q.front(); Q.pop();
            aiVec = &(mesh->mVertices[idx]);
            vec3 idxVec = vec3(aiVec->x, aiVec->y, aiVec->z);
            // Put the next level vertecies into the queue.
            for(int e = first[idx]; e != -1; e = next[e]) {
                int idxNext = incidentVertex[e];
                // Expand the next level vertecies.
                if(!used[idxNext]) {
                    aiVec = &(mesh->mVertices[idxNext]);
                    vec3 idxNextVec = vec3(aiVec->x, aiVec->y, aiVec->z);
                    if(get_squared_dist(vVec, idxNextVec) <= 36*sigma*sigma)
                        Q.push(incidentVertex[e]),
                        used[incidentVertex[e]] = 1;
                }
            }
            // Update Gaussian filter
            float dist = get_squared_dist(vVec, idxVec);
            for(int i = 2; i <= 6; i++) {
                float sigmaHere = i*i*sigma*sigma;
                if(dist <= sigmaHere) {
                    float factor = exp(-dist/(2*sigmaHere));
                    gaussianSigma1[i] += meanCurvature[idx] * factor;
                    sumSigma1[i] += factor;
                }
                if(dist <= 2*sigmaHere) {
                    float factor = exp(-dist/(8*sigma*sigma));
                    gaussianSigma2[i] += meanCurvature[idx] * factor;
                    sumSigma2[i] += factor;
                }
            }
        }
        for(int i = 2; i <= 6; i++) {
            saliency[i][k] = fabs(gaussianSigma1[i]/sumSigma1[i]
                                - gaussianSigma2[i]/sumSigma2[i]);
            maxSaliency[i] = fmax(maxSaliency[i], saliency[i][k]);
        }
    }

    // Second BFS and get the non-linear normailization of suppressian's saliency.
    smoothSaliency = (float*)malloc(*point_count * sizeof(float));
    for(int k = 0; k < *point_count; k++) {
        smoothSaliency[k] = 0.0f;
        float localMaxSaliency[7];//, localCntSaliency[7];
        for(int i = 2; i <= 6; i++)
            localMaxSaliency[i] = -oo;
        // Get the current vertex's information.
        aiVector3D* aiVec = &(mesh->mVertices[k]);
        vec3 vVec = vec3(aiVec->x, aiVec->y, aiVec->z);
        // Initialize the queue to find neighbourhood.
        queue<int> Q;
        memset(used, 0, *point_count * sizeof(bool));
        Q.push(k);
        used[k] = true;
        while(!Q.empty()) {
            // Get the front element in the queue.
            int idx = Q.front(); Q.pop();
            aiVec = &(mesh->mVertices[idx]);
            vec3 idxVec = vec3(aiVec->x, aiVec->y, aiVec->z);
            // Put the next level vertecies into the queue.
            for(int e = first[idx]; e != -1; e = next[e]) {
                int idxNext = incidentVertex[e];
                // Expand the next level vertecies.
                if(!used[idxNext]) {
                    aiVec = &(mesh->mVertices[idxNext]);
                    vec3 idxNextVec = vec3(aiVec->x, aiVec->y, aiVec->z);
                    if(get_squared_dist(vVec, idxNextVec) <= 36*sigma*sigma)
                        Q.push(incidentVertex[e]),
                        used[incidentVertex[e]] = 1;
                }
            }
            // Update Gaussian filter
            float dist = get_squared_dist(vVec, idxVec);
            for(int i = 2; i <= 6; i++) 
                localMaxSaliency[2] = fmax(localMaxSaliency[2], saliency[i][idx]);
        }
        // Calculate the weighted saliency
        float saliencySum = 0.0f;
        for(int i = 2; i <= 6; i++) {
            float factor = (maxSaliency[i]-localMaxSaliency[i])*(maxSaliency[i]-localMaxSaliency[i]);
            smoothSaliency[k] += (GLfloat)saliencySum[i][k] * factor;
            saliencySum += factor;
        }
        smoothSaliency[k] /= (GLfloat)saliencySum;
    }
    // Clean up resources
    free(first);
    free(next);
    free(incidentVertex);
    for(int i = 2; i <= 6; i++)
        free(saliency[i]);
	
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
        //free (normals);
    }

	aiReleaseImport (scene);
	printf ("mesh loaded\n");
	
	return true;
}

inline void checkLegal(float &x) {
    if(x > 0.9999f) x = 0.9999f;
    if(x < -0.9999f) x = -0.9999f;
}

// Change the YUV to RGB, and guarantee all R,G,B values in [-1, 1].
void getYUVtoRGB(float Y, float U, float V, float& R, float& G, float& B) {
    R = Y + 1.4075 * (V - 128.0);
    G = Y - 0.3455 * (U - 128.0) - (0.7169 * (V - 128.0f));
    B = Y + 1.7790 * (U - 128.0);
    checkLegal(R), checkLegal(G), checkLegal(B);   
}

void updateDisplayType(int type) {
    GLuint vbo;
    glGenBuffers (1, &vbo);
    glBindBuffer (GL_ARRAY_BUFFER, vbo);
    float* colors = NULL;
    colors = (float*)malloc(vertexCnt * 3 * sizeof(float));
    switch(type) {
        case 1: {
                    for(int i = 0; i < vertexCnt; i++)
                        colors[3*i+0] = normals[3*i+0],
                        colors[3*i+1] = normals[3*i+1],
                        colors[3*i+2] = normals[3*i+2];
                    break;
                }
        case 2: {
                    float xMin = oo, xMax = -oo;
                    for(int i = 0; i < vertexCnt; i++)
                        xMin = fmin(xMin, meanCurvature[i]),
                        xMax = fmax(xMax, meanCurvature[i]);
                    for(int i = 0; i < vertexCnt; i++) {
                        float Y = 2.0f*((meanCurvature[i]-xMin)/(xMax-xMin)-0.5f);
                        getYUVtoRGB(Y, 1.0f, 1.0f, colors[3*i], colors[3*i+1], colors[3*i+2]);
                    }
                }
        case 3: {
                    float xMin = oo, xMax = -oo;
                    for(int i = 0; i < vertexCnt; i++)
                        xMin = fmin(xMin, smoothSaliency[i]),
                        xMax = fmax(xMax, smoothSaliency[i]);
                    for(int i = 0; i < vertexCnt; i++) {
                        float Y = 2.0f*((smoothSaliency[i]-xMin)/(xMax-xMin)-0.5f);
                        getYUVtoRGB(Y, 1.0f, 1.0f, colors[3*i], colors[3*i+1], colors[3*i+2]);
                    }
                }
    }
    glBufferData (
        GL_ARRAY_BUFFER,
        3 * vertexCnt * sizeof (GLfloat),
        colors,
        GL_STATIC_DRAW
    );
    glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray (1);
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
    // Update the display mode.
    if (glfwGetKey(g_window, GLFW_KEY_1))
        updateDisplayType(1);
    else if (glfwGetKey(g_window, GLFW_KEY_2))
        updateDisplayType(2);
    else if (glfwGetKey(g_window, GLFW_KEY_3))
        updateDisplayType(3);
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
