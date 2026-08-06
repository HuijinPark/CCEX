#!/bin/sh

if [ $# -ne 1 ]; then
    echo "Error, no input value"
    echo "[1] Server : idm4 / mac / idm3 / nurion"
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
CXXFLAGS += -Wno-writable-strings -fsanitize=address

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

elif [ ${1} == "nurion" ]; then

module purge
module load craype-network-opa craype-mic-knl intel/19.0.5 impi/19.0.5


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

INCLUDE_EIGEN = -I./zlib/eigen
INCLUDE_UTHASH = -I./zlib/uthash/include/
INCLUDE_MAIN = -I./include/

LIBRARY_MKL = -L /apps/compiler/intel/19.0.5/mkl/lib/intel64
LDFLAGS_MKL = -DMKL_ILP64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -liomp5 -m64 #-lgsl -lgslcblas 

INCLUDE = \$(INCLUDE_EIGEN) \$(INCLUDE_UTHASH) \$(INCLUDE_MAIN) #\$(INCLUDE_MPICH)
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
CXXFLAGS += -Wno-writable-strings -fsanitize=address

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

elif [ ${1} == "wsl" ]; then

# WSL/Ubuntu + Intel oneAPI 2026.
#
# Compiler: mpicxx (Intel MPI wrapper around g++; gcc backend, not icpx).
#   Avoids icpc/icpx deprecation noise; matches the reference WSL build that
#   is known to work with WSL networkingMode=mirrored. Critical: bash command
#   hash must be cleared (hash -r) after sourcing setvars.sh, otherwise a
#   previously-cached /usr/bin/mpicxx (broken system MPICH 4.2) shadows the
#   Intel MPI wrapper. As a defense-in-depth, the Makefile uses the absolute
#   path of mpicxx baked in at script time.
#
# MPI fabric: I_MPI_FABRICS=shm forces shared-memory transport, which
#   bypasses the WSL network namespace. Required for single-node MPI to work
#   correctly when networkingMode=mirrored is enabled.
#
# MKL: searches $(MKLROOT)/lib/intel64 first (oneAPI <=2025 layout), then
#   falls back to $(MKLROOT)/lib (oneAPI 2026+ layout). rpath embedded.

. /opt/intel/oneapi/setvars.sh > /dev/null 2>&1
hash -r
export I_MPI_FABRICS=shm

MPICXX_BIN="${I_MPI_ROOT}/bin/mpicxx"
if [ ! -x "$MPICXX_BIN" ]; then
    echo "[wsl] ERROR: Intel MPI's mpicxx not found at $MPICXX_BIN" >&2
    echo "[wsl] Check that Intel oneAPI is installed and setvars.sh succeeded." >&2
    exit 1
fi

echo "[wsl] I_MPI_ROOT      = $I_MPI_ROOT"
echo "[wsl] MKLROOT         = $MKLROOT"
echo "[wsl] MPICXX_BIN      = $MPICXX_BIN"
echo "[wsl] I_MPI_FABRICS   = $I_MPI_FABRICS  (shm = shared memory; mirrored-mode safe)"

cat << EOF > ./Makefile

# Auto-generated by do_compile.sh wsl on \$(date -u +"%Y-%m-%dT%H:%M:%SZ")
# WSL/Ubuntu + Intel oneAPI 2026 + Intel MPI (mpicxx, g++ backend)

CXX := ${MPICXX_BIN}

MKLROOT_FALLBACK ?= ${MKLROOT}
MKL_LIB_DIR := \$(firstword \$(wildcard \$(MKLROOT_FALLBACK)/lib/intel64 \$(MKLROOT_FALLBACK)/lib))
ifeq (\$(strip \$(MKL_LIB_DIR)),)
\$(error Could not find MKL library directory under \$(MKLROOT_FALLBACK)/lib)
endif

CXXFLAGS := -std=gnu++17 -O3 -g -march=native -mtune=native
CXXFLAGS += -D_FORTIFY_SOURCE=0
CXXFLAGS += -fpermissive
CXXFLAGS += -Wno-deprecated-declarations
CXXFLAGS += -Wno-writable-strings
CXXFLAGS += -Wno-unused-result
CXXFLAGS += -Wno-format-overflow
CXXFLAGS += -Wno-format-security

SRC_DIR := ./src
OBJ_DIR := ./obj
BIN_DIR := ./bin

TARGET := \$(BIN_DIR)/main.out

SRCS := \$(wildcard \$(SRC_DIR)/*.cpp)
SRCS += main.cpp
OBJS := \$(patsubst \$(SRC_DIR)/%.cpp,\$(OBJ_DIR)/%.o,\$(SRCS))

INCLUDE_EIGEN  := -I./zlib/eigen
INCLUDE_UTHASH := -I./zlib/uthash/include/
INCLUDE_MAIN   := -I./include/
INCLUDE_MKL    := -I\$(MKLROOT_FALLBACK)/include

LIBRARY_MKL := -L\$(MKL_LIB_DIR) -Wl,-rpath,\$(MKL_LIB_DIR)
LDFLAGS_MKL := -DMKL_ILP64 -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread -lm

INCLUDE := \$(INCLUDE_EIGEN) \$(INCLUDE_MKL) \$(INCLUDE_UTHASH) \$(INCLUDE_MAIN)
LIBRARY := \$(LIBRARY_MKL) \$(LDFLAGS_MKL)

.PHONY: all clean info
all: \$(TARGET)

\$(TARGET): \$(OBJS) | \$(BIN_DIR)
	\$(CXX) \$(OBJS) -o \$(TARGET) \$(CXXFLAGS) \$(INCLUDE) \$(LIBRARY)

\$(OBJ_DIR)/%.o: \$(SRC_DIR)/%.cpp | \$(OBJ_DIR)
	\$(CXX) -c \$< -o \$@ \$(CXXFLAGS) \$(INCLUDE) \$(LIBRARY)

\$(OBJ_DIR):
	mkdir -p \$(OBJ_DIR)

\$(BIN_DIR):
	mkdir -p \$(BIN_DIR)

info:
	@echo "CXX         = \$(CXX)"
	@echo "MKLROOT     = \$(MKLROOT_FALLBACK)"
	@echo "MKL_LIB_DIR = \$(MKL_LIB_DIR)"
	@\$(CXX) -show 2>/dev/null | head -1

clean:
	rm -rf \$(OBJ_DIR) \$(BIN_DIR)

EOF

else
 :
fi

make -j$(nproc)
#make 
