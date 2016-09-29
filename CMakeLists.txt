cmake_minimum_required(VERSION 3.5)
project(Mesh_Saliency)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    lib/include/assimp/Compiler/poppack1.h
    lib/include/assimp/Compiler/pushpack1.h
    lib/include/assimp/ai_assert.h
    lib/include/assimp/anim.h
    lib/include/assimp/camera.h
    lib/include/assimp/cexport.h
    lib/include/assimp/cfileio.h
    lib/include/assimp/cimport.h
    lib/include/assimp/color4.h
    lib/include/assimp/config.h
    lib/include/assimp/DefaultLogger.hpp
    lib/include/assimp/defs.h
    lib/include/assimp/Exporter.hpp
    lib/include/assimp/Importer.hpp
    lib/include/assimp/importerdesc.h
    lib/include/assimp/IOStream.hpp
    lib/include/assimp/IOSystem.hpp
    lib/include/assimp/light.h
    lib/include/assimp/Logger.hpp
    lib/include/assimp/LogStream.hpp
    lib/include/assimp/material.h
    lib/include/assimp/matrix3x3.h
    lib/include/assimp/matrix4x4.h
    lib/include/assimp/mesh.h
    lib/include/assimp/NullLogger.hpp
    lib/include/assimp/postprocess.h
    lib/include/assimp/ProgressHandler.hpp
    lib/include/assimp/quaternion.h
    lib/include/assimp/scene.h
    lib/include/assimp/texture.h
    lib/include/assimp/types.h
    lib/include/assimp/vector2.h
    lib/include/assimp/vector3.h
    lib/include/assimp/version.h
    lib/include/GL/glew.h
    lib/include/GL/glxew.h
    lib/include/GL/wglew.h
    lib/include/GLFW/glfw3.h
    lib/include/GLFW/glfw3native.h
    lib/gl_utils.cpp
    lib/gl_utils.h
    lib/maths_funcs.cpp
    lib/maths_funcs.h
    display.cpp
    display.h
    main.cpp
    Mesh.h
    QSlim.cpp
    QSlim.h
    saliency.cpp
    saliency.h)

add_executable(Mesh_Saliency ${SOURCE_FILES})