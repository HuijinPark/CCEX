#!/bin/sh

if [ $# -ne 1 ]; then
    echo "Error, no input value"
    echo "[1] Server : idm4 / mac / idm3"
    exit 1
fi


if [ ${1} == "idm4" ]; then

module purge 
#module load intel-23.2/fftw-3.3.10 intel-23.2/icc-23.2
module load intel-23.2/icc-23.2

cat << EOF > ./Makefile

CXX = mpiicpc

CXXFLAGS = -std=c++11 -O3 -g #-Wall Higher level warning
#CXXFLAGS += -Wno-c++11-compat-deprecated-writable-strings 
CXXFLAGS += -Wno-deprecated-declarations
CXXFLAGS += -diag-disable=2196
CXXFLAGS += -diag-disable=10441
#CXXFLAGS += -Wno-writable-strings

SRC_DIR=./src
OBJ_DIR=./obj
BIN_DIR=./bin

TARGET = \$(BIN_DIR)/main.out

SRCS=\$(wildcard \$(SRC_DIR)/*.cpp)
SRCS += main.cpp
OBJS=\$(patsubst \$(SRC_DIR)/%.cpp,\$(OBJ_DIR)/%.o,\$(SRCS))

#INCLUDE_MPICH = -I/opt/homebrew/Cellar/mpich/4.2.1/include/
INCLUDE_EIGEN = -I./zlib/eigen
INCLUDE_UTHASH = -I./zlib/uthash/include/
INCLUDE_MAIN = -I./include/

#INCLUDE_MKL = -I /opt/intel/mkl/include
#LIBRARY_MKL = -L /opt/intel/mkl/lib/intel64
LDFLAGS_MKL = -DMKL_ILP64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -liomp5 -m64 -xCORE-AVX512  #-lgsl -lgslcblas 

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

elif [ ${1} == "idm3" ]; then

module purge 
module load 22.2/gsl-2.7.1 22.2/icc-22.2 22.2/fftw-3.3.10

cat << EOF > ./Makefile

CXX = mpiicpc

CXXFLAGS = -std=gnu++17 -O3 -g #-Wall Higher level warning
#CXXFLAGS += -Wno-c++11-compat-deprecated-writable-strings 
CXXFLAGS += -Wno-deprecated-declarations
CXXFLAGS += -diag-disable=2196
#CXXFLAGS += -Wno-writable-strings

SRC_DIR=./src
OBJ_DIR=./obj
BIN_DIR=./bin

TARGET = \$(BIN_DIR)/main.out

SRCS=\$(wildcard \$(SRC_DIR)/*.cpp)
SRCS += main.cpp
OBJS=\$(patsubst \$(SRC_DIR)/%.cpp,\$(OBJ_DIR)/%.o,\$(SRCS))

#INCLUDE_MPICH = -I/opt/homebrew/Cellar/mpich/4.2.1/include/
INCLUDE_EIGEN = -I./zlib/eigen
INCLUDE_UTHASH = -I./zlib/uthash/include/
INCLUDE_MAIN = -I./include/

INCLUDE_MKL = -I /opt/intel/oneapi/mkl/latest/include/
LIBRARY_MKL = -L /opt/intel/oneapi/mkl/latest/lib/intel64/
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

CXXFLAGS = -std=c++11 -O2 -g #-Wall Higher level warning
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
INCLUDE_EIGEN = -I./zlib/eigen #-3.4.0/
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
