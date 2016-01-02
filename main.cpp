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
#include "Mesh.h"
#include "QSlim.h"

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


static void qslim_init() {
    int i;
    cerr << "Reading input ..." << endl;
    cerr << "Cleaning up initial input ..." << endl;
	int initialVertCount = M0.vertCount();
	int initialEdgeCount = M0.edgeCount();
	int initialFaceCount = M0.faceCount();
    for(i=0; i<M0.faceCount(); i++)
	if( !M0.face(i)->plane().isValid() )
	    M0.killFace(M0.face(i));
    M0.removeDegeneracy(M0.allFaces());
    for(i=0; i<M0.vertCount(); i++)
    {
	if( M0.vertex(i)->edgeUses().length() == 0 )
	    M0.vertex(i)->kill();
    }
    cerr << "Input model summary:" << endl;
    cerr << "    Vertices    : " << initialVertCount << endl;
    cerr << "    Edges       : " << initialEdgeCount << endl;
    int man=0, non=0, bndry=0, bogus=0;
    for(i=0; i<M0.edgeCount(); i++)
        switch( M0.edge(i)->faceUses().length() )
        {
        case 0:
            bogus++;
            break;
        case 1:
            bndry++;
            break;
        case 2:
            man++;
            break;
        default:
            non++;
            break;
        }
    if( bogus )
        cerr << "        Bogus       : " << bogus << endl;
    cerr << "        Boundary    : " << bndry << endl;
    cerr << "        Manifold    : " << man << endl;
    cerr << "        Nonmanifold : " << non << endl;

    cerr << "    Faces       : " << initialFaceCount << endl;
}

static void qslim_run() {
    decimate_init(M0, pair_selection_tolerance);
    while( M0.validFaceCount > face_target&& decimate_min_error() < error_tolerance )
		decimate_contract(M0);
}

bool loadMeshQSlim(const char* fileName, Mesh& m) {
    // Load the mesh by assimp.
    const aiScene* scene = aiImportFile(fileName, aiProcess_Triangulate);
    if (!scene) {
        fprintf(stderr, "ERROR: reading mesh %s\n", fileName);
        return false;
    }
    const aiMesh* mesh = scene->mMeshes[0];
    const int vertexCntHere = mesh->mNumVertices;
    for(int i = 0; i < vertexCntHere; i++) {
        const aiVector3D* vp = &(mesh->mVertices[i]);
        Point3d p3d(vp->x, vp->y, vp->z);
        m.AddVertex(p3d);
        Vec3 v(vp->x, vp->y, vp->z);
        M0.in_Vertex(v);
    }
    const int faceCntHere = mesh->mNumFaces;
    for(int i = 0; i < faceCntHere; i++) {
        int idx[3];
        for(int k = 0; k < 3; k++)
            idx[k] = mesh->mFaces[i].mIndices[k];
        Triangle t(idx[0], idx[1], idx[2]);
        m.AddFace(t);
        M0.in_Face(idx[0], idx[1], idx[2]);
    }
    return true;
}

static void ReplaceM(Mesh& m) {
    vector<Point3d> newVertices;
	m.Vertices.swap(newVertices);
    vector<Triangle> newFaces;
	m.Faces.swap(newFaces);
	m.Vertices.reserve(M0.vertCount());
	m.Faces.reserve(M0.faceCount());
	int* map=new int[M0.vertCount()];
	for(int i=0;i<M0.vertCount();i++)
		map[i]=-1;
	for(int i=0;i<M0.vertCount();i++) {
		if(M0.vertex (i)->isValid()) {
			real* data=M0.vertex(i)->raw();
			Point3d p((float)data[0],(float)data[1],(float)data[2]);
			map[i]=m.AddVertex(p);
		}
	}
	for(int i=0;i<M0.faceCount();i++)
	{
		if(M0.face(i)->isValid())
		{
			Vertex* v0 = M0.face(i)->vertex(0);
			Vertex* v1 = M0.face(i)->vertex(1);
			Vertex* v2 = M0.face(i)->vertex(2);
			Triangle t(map[v0->uniqID], map[v1->uniqID], map[v2->uniqID]);
			m.AddFace(t);
		}
	}
	delete[] map;
}

void callQSlim(Mesh& m) {
    assert( loadMeshQSlim(MESH_FILE, m) );
	qslim_init();
    float ratio = 60.0f;
	face_target = (int)(1.0f*m.Faces.size()*ratio/100.0f);
	error_tolerance = oo;
	will_use_plane_constraint = true;
	will_use_vertex_constraint = false;
	will_preserve_boundaries = true;
	will_preserve_mesh_quality = true;
	will_constrain_boundaries = true;
	boundary_constraint_weight = 1.0;
	will_weight_by_area = false;
	placement_policy = 1;
	pair_selection_tolerance = 0.0;
	qslim_run();
	ReplaceM(m);
    printf("Simplification: %d / %lu\n", face_target, m.Faces.size());
    //updateDisplayMesh(m);
}

int main () {
	initializeOpenGL();
    glfwSetCursorPosCallback(g_window, mouseEventHandler);
	// load the mesh using assimp

	assert (load_mesh (MESH_FILE, &objVAO, &ObjPointCount));
	Mesh m;
    callQSlim(m);
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
		keyBoardEvent(elapsed_seconds);
		// put the stuff we've been drawing onto the display
		glfwSwapBuffers (g_window);
	}
	

	// close GL context and any other GLFW resources
	glfwTerminate();
	return 0;
}
