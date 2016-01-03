BIN = mesh_saliency
CC = g++
FLAGS = -DAPPLE -Wall -pedantic -mmacosx-version-min=10.5 -arch x86_64 -fmessage-length=0 -UGLFW_CDECL -fprofile-arcs -ftest-coverage
INC = -I lib/include -I/sw/include -I/usr/local/include
LIB_PATH = lib/osx_64/
LOC_LIB = $(LIB_PATH)libGLEW.a $(LIB_PATH)libglfw3.a $(LIB_PATH)libassimp.a
SYS_LIB = -lz
FRAMEWORKS = -framework Cocoa -framework OpenGL -framework IOKit
SRC = main.cpp lib/maths_funcs.cpp lib/gl_utils.cpp saliency.cpp display.cpp

all:
	${CC} ${FLAGS} ${FRAMEWORKS} -o ${BIN} ${SRC} ${INC} ${LOC_LIB} ${SYS_LIB}

