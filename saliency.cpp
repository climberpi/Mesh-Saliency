/*
 * author: Yupan Liu
 * date: Dec 27, 2015
 * brief: normal, mean curvature and mesh saliency
 */
#include "saliency.h"
// Load a mesh using the assimp library, and calculate its mesh saliency
bool load_mesh (const char* file_name, GLuint* vao, int* point_count) {
	const aiScene* scene = aiImportFile (file_name, aiProcess_Triangulate);
	if (!scene) {
		fprintf (stderr, "ERROR: reading mesh %s\n", file_name);
		return false;
	}
   
	// Get first mesh in file only
	const aiMesh* mesh = scene->mMeshes[0];
	printf ("    %i vertices in mesh[0]\n", mesh->mNumVertices);
	printf ("    %i faces in mesh[0]\n", mesh->mNumFaces);
	
	// Pass back number of vertex points in mesh 
	*point_count = mesh->mNumVertices;
    vertexCnt = *point_count;
	
	// Generate a VAO, using the pass-by-reference parameter that we give to the function 
    glGenVertexArrays (1, vao);
	glBindVertexArray (*vao);
	
	GLfloat* points = NULL; // array of vertex points
    if (!mesh->HasPositions()) {
        fprintf(stderr, "ERROR: mesh %s don't have vertex data!\n", file_name);
        return false;
    }

    float xMin = oo, yMin = oo, zMin = oo;
    float xMax = -oo, yMax = -oo, zMax = -oo;
    // Get the mesh's vertecies
    points = (GLfloat*)malloc (*point_count * 3 * sizeof (GLfloat));
    for (int i = 0; i < *point_count; i++) {
        const aiVector3D* vp = &(mesh->mVertices[i]);
        points[i * 3] = (GLfloat)vp->x;
        points[i * 3 + 1] = (GLfloat)vp->y;
        points[i * 3 + 2] = (GLfloat)vp->z;
        xMin = fMin(xMin, vp->x), xMax = fMax(xMax, vp->x);
        yMin = fMin(yMin, vp->y), yMax = fMax(yMax, vp->y);
        zMin = fMin(zMin, vp->z), zMax = fMax(zMax, vp->z);
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
            vec3 Tij = (identity_mat3() - wedge(Ni, Ni))*(Vi-Vj);
            Tij = normalise(Tij);
            float kappa_ij = 2*dot(Ni, Vj-Vi);
            kappa_ij /= get_squared_dist(Vi, Vj);
            // Maintain vi's shape operator
            shapeOperators[i] = shapeOperators[i] + (wedge(Tij, Tij) * (kappa_ij * faceArea));
            vertexArea[i] += faceArea;

            // For vertex j, update the relative part of its shape operator
            vec3 Tji = (identity_mat3() - wedge(Nj, Nj))*(Vj-Vi);
            Tji = normalise(Tji);
            float kappa_ji = 2*dot(Nj, Vi-Vj);
            kappa_ji /= get_squared_dist(Vi, Vj);
            // Maintain vj's shape operator
            shapeOperators[j] = shapeOperators[j] + (wedge(Tji, Tji) * (kappa_ji * faceArea));
            
            vertexArea[j] += faceArea;
        }
    }

    for(int i = 0; i < *point_count; i++) {
        shapeOperators[i] = shapeOperators[i] * (1.0f/vertexArea[i]);// * 10000000.0f;
        //print(shapeOperators[i]);
    }
    free(vertexArea);

    // Diagonalize the shape operator, and get the mean curvature
    meanCurvature = (float*)malloc(*point_count * sizeof(float));
    for(int k = 0; k < *point_count; k++) {
        vec3 E1 = vec3(1.0f, 0.0f, 0.0f);
        vec3 Nk = vec3(normals[k*3], normals[k*3+1], normals[k*3+2]);
        bool isMinus = get_squared_dist(E1, Nk) > get_squared_dist(E1 * (-1.0f), Nk);
        vec3 Wk;
        // Diagnoalization by the Householder transform
        if (!isMinus)
            Wk = E1 + Nk;
        else
            Wk = E1 - Nk;
        Wk = normalise(Wk);
        mat3 Qk = identity_mat3() - (wedge(Wk, Wk) * 2.0f);
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
    for(int i = 0; i < *point_count; i++)
        first[i] = -1;
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
    
    printf("BFS 1\n");
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
        if(k%1000 == 0)
            printf("#%d#\n", k);
        // Initialize the saliency and its local counter.
        for(int i = 2; i <= 6; i++)
            saliency[i][k] = 0.0f;
        // Initialize the saliency's Gaussian filter.
        float gaussianSigma1[7], gaussianSigma2[7], sumSigma1[7], sumSigma2[7];
        for(int i = 2; i <= 6; i++)
            gaussianSigma1[i] = gaussianSigma2[i] = 0.0f,
            sumSigma1[i] = sumSigma2[i] = 0.0f;
        // Get the current vertex's information.
        aiVector3D* aiVec = &(mesh->mVertices[k]);
        vec3 vVec = vec3(aiVec->x, aiVec->y, aiVec->z);
        // Initialize the queue to find neighbourhood.
        for(int i = 0; i < *point_count; i++)
            used[i] = false;
        queue<int> Q;
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
            maxSaliency[i] = fMax(maxSaliency[i], saliency[i][k]);
        }
    }

    printf("BFS 2\n");
    // Second BFS and get the non-linear normailization of suppressian's saliency.
    smoothSaliency = (float*)malloc(*point_count * sizeof(float));
    for(int k = 0; k < *point_count; k++) {
        if(k%1000 == 0)
            printf("[%d]\n", k);
        smoothSaliency[k] = 0.0f;
        float localMaxSaliency[7];//, localCntSaliency[7];
        for(int i = 2; i <= 6; i++)
            localMaxSaliency[i] = -oo;
        // Get the current vertex's information.
        aiVector3D* aiVec = &(mesh->mVertices[k]);
        vec3 vVec = vec3(aiVec->x, aiVec->y, aiVec->z);
        // Initialize the queue to find neighbourhood.
        for(int i = 0; i < *point_count; i++)
            used[i] = false;
        queue<int> Q;
        Q.push(k);
        used[k] = true;
        while(!Q.empty()) {
            // Get the front element in the queue.
            int idx = Q.front(); Q.pop();
            //aiVec = &(mesh->mVertices[idx]);
            //vec3 idxVec = vec3(aiVec->x, aiVec->y, aiVec->z);
            // Put the next level vertecies into the queue.
            for(int e = first[idx]; e != -1; e = next[e]) {
                int idxNext = incidentVertex[e]; 
                // Expand the next level vertecies.
                if(!used[idxNext]) {
                    aiVec = &(mesh->mVertices[idxNext]); 
                    vec3 idxNextVec = vec3(aiVec->x, aiVec->y, aiVec->z);
                    if(get_squared_dist(vVec, idxNextVec) <= 36*sigma*sigma)
                        Q.push(incidentVertex[ e]),
                        used[incidentVertex[e]] = 1;
                }
            }
            // Update Gaussian filter
            for(int i = 2; i <= 6; i++) 
                localMaxSaliency[i] = fMax(localMaxSaliency[i], saliency[i][idx]);
        }
        // Calculate the weighted saliency
        float saliencySum = 0.0f;
        for(int i = 2; i <= 6; i++) {
            float factor = (maxSaliency[i]-localMaxSaliency[i])*(maxSaliency[i]-localMaxSaliency[i]);
            smoothSaliency[k] += (GLfloat)saliency[i][k] * factor;
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
    updateDisplayType(2);
	aiReleaseImport (scene);
	printf ("mesh loaded\n");
	
	return true;
}
