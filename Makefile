
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

TARGET = $(BIN_DIR)/main.out

SRCS=$(wildcard $(SRC_DIR)/*.cpp)
SRCS += main.cpp
OBJS=$(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

#INCLUDE_MPICH = -I/opt/homebrew/Cellar/mpich/4.2.1/include/
INCLUDE_EIGEN = -I./zlib/eigen
INCLUDE_UTHASH = -I./zlib/uthash/include/
INCLUDE_MAIN = -I./include/

#INCLUDE_MKL = -I /opt/intel/mkl/include
#LIBRARY_MKL = -L /opt/intel/mkl/lib/intel64
LDFLAGS_MKL = -DMKL_ILP64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -liomp5 -m64 -xCORE-AVX512  #-lgsl -lgslcblas 

INCLUDE = $(INCLUDE_EIGEN) $(INCLUDE_MKL) $(INCLUDE_UTHASH) $(INCLUDE_MAIN) #$(INCLUDE_MPICH)
LIBRARY = $(LIBRARY_MKL) $(LDFLAGS_MKL)

all: $(TARGET)

$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(OBJS) -o $(TARGET) $(CXXFLAGS) $(INCLUDE) $(LIBRARY)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDE) $(LIBRARY)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

