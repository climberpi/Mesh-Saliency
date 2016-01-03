Mesh-Saliency
=======================================
The program based on siggraph 05 paper, mesh saliency(Lee et al, 2005), implemented by OpenGL 4 and Assimp. Using mean curvative(Taubin, 1995) to improve the performance of Surface simplification(Garland et al, 1997). Some codes modified by "Anton's OpenGL 4 Tutorials". And you can use the follow code to compile the program:
```
make -f Makefile
./mesh_saliency
```

### Reference
- Anton's OpenGL 4 Tutorials: <http://antongerdelan.net/opengl/>
- Lee C H, Varshney A, Jacobs D W. Mesh saliency[C]//ACM transactions on graphics (TOG). ACM, 2005, 24(3): 659-666.
- Taubin G. Estimating the tensor of curvature of a surface from a polyhedral approximation[C]//Computer Vision, 1995. Proceedings., Fifth International Conference on. IEEE, 1995: 902-907.
- Garland M, Heckbert P S. Surface simplification using quadric error metrics[C]//Proceedings of the 24th annual conference on Computer graphics and interactive techniques. ACM Press/Addison-Wesley Publishing Co., 1997: 209-216.

### QSlim
Using Michael Garland's code, and modified by Alec Jacobson:
<https://github.com/alecjacobson/qslim>
