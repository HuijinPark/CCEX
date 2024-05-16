#!/bin/sh

module purge 
module load icc18 icc18.impi fftw

#export LD_LIBRARY_PATH=/home/CQML/lib/gsl-2.7.1/lib:$LD_LIBRARY_PATH
#export INCLUDE=/home/CQML/lib/gsl-2.7.1/include:$INCLUDE

# $1 = idm
# $1 = mac

if [ $# -ne 1 ]; then
    echo "Error, no input value"
    echo "[1] Server : idm / mac"
    exit 1
fi


if [ ${1} == "idm" ]; then

cat << EOF > ./Makefile

CXX = mpiicpc

CXXFLAGS = -std=c++11 #-O2 -g #-Wall Higher level warning
#CXXFLAGS += -Wno-c++11-compat-deprecated-writable-strings 
CXXFLAGS += -Wno-deprecated-declarations
#CXXFLAGS += -Wno-writable-strings

SRC_DIR=./src
OBJ_DIR=./obj
BIN_DIR=./bin

TARGET = \$(BIN_DIR)/main.out

SRCS=\$(wildcard \$(SRC_DIR)/*.cpp)
SRCS += main.cpp
OBJS=\$(patsubst \$(SRC_DIR)/%.cpp,\$(OBJ_DIR)/%.o,\$(SRCS))

#INCLUDE_MPICH = -I/opt/homebrew/Cellar/mpich/4.2.1/include/
INCLUDE_EIGEN = -I./zlib/eigen-3.4.0/
INCLUDE_UTHASH = -I./zlib/uthash/include/
INCLUDE_MAIN = -I./include/

INCLUDE_MKL = -I /opt/intel/mkl/include
LIBRARY_MKL = -L /opt/intel/mkl/lib/intel64
LDFLAGS_MKL = -DMKL_ILP64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -m64 #-lgsl -lgslcblas 

INCLUDE = \$(INCLUDE_EIGEN) \$(INCLUDE_MKL) \$(INCLUDE_UTHASH) \$(INCLUDE_MAIN) #\$(INCLUDE_MPICH)
LIBRARY = \$(LIBRARY_MKL) \$(LDFLAGS_MKL)

all: \$(TARGET)

\$(TARGET): \$(OBJS) | \$(BIN_DIR)
	\$(CXX) \$(OBJS) -o \$(TARGET) \$(CXXFLAGS) \$(INCLUDE) \$(LIBRARY)

\$(OBJ_DIR)/%.o: \$(SRC_DIR)/%.cpp | \$(OBJ_DIR)
	\$(CXX) -c $< -o \$@ \$(CXXFLAGS) \$(INCLUDE) \$(LIBRARY)

\$(OBJ_DIR):
	mkdir -p \$(OBJ_DIR)

\$(BIN_DIR):
	mkdir -p \$(BIN_DIR)

clean:
	rm -rf \$(OBJ_DIR) \$(BIN_DIR)

EOF


elif [ ${1} == "mac" ]; then


cat << EOF > ./Makefile

CXX = mpicxx

CXXFLAGS = -std=c++11 #-O2 -g #-Wall Higher level warning
CXXFLAGS += -Wno-c++11-compat-deprecated-writable-strings 
CXXFLAGS += -Wno-deprecated-declarations
CXXFLAGS += -Wno-writable-strings

SRC_DIR=./src
OBJ_DIR=./obj
BIN_DIR=./bin

TARGET = \$(BIN_DIR)/main.out

SRCS=\$(wildcard \$(SRC_DIR)/*.cpp)
SRCS += main.cpp
OBJS=\$(patsubst \$(SRC_DIR)/%.cpp,\$(OBJ_DIR)/%.o,\$(SRCS))

INCLUDE_MPICH = -I/opt/homebrew/Cellar/mpich/4.2.1/include/
INCLUDE_EIGEN = -I./zlib/eigen-3.4.0/
INCLUDE_UTHASH = -I./zlib/uthash/include/
INCLUDE_MAIN = -I./include/

INCLUDE_MKL = -I/opt/intel/oneapi/mkl/latest/include/
LIBRARY_MKL = -L/opt/intel/oneapi/mkl/latest/lib/
LDFLAGS_MKL = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm #-ldl #-liomp5 

INCLUDE = \$(INCLUDE_EIGEN) \$(INCLUDE_MKL) \$(INCLUDE_UTHASH) \$(INCLUDE_MAIN) \$(INCLUDE_MPICH)
LIBRARY = \$(LIBRARY_MKL) \$(LDFLAGS_MKL)

all: \$(TARGET)

\$(TARGET): \$(OBJS) | \$(BIN_DIR)
	\$(CXX) \$(OBJS) -o \$(TARGET) \$(CXXFLAGS) \$(INCLUDE) \$(LIBRARY)

\$(OBJ_DIR)/%.o: \$(SRC_DIR)/%.cpp | \$(OBJ_DIR)
	\$(CXX) -c $< -o \$@ \$(CXXFLAGS) \$(INCLUDE) \$(LIBRARY)

\$(OBJ_DIR):
	mkdir -p \$(OBJ_DIR)

\$(BIN_DIR):
	mkdir -p \$(BIN_DIR)

clean:
	rm -rf \$(OBJ_DIR) \$(BIN_DIR)

EOF

else
 :
fi

make -j$(nproc)
#make 
