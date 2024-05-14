CXX=mpicxx

TARGET = test_defect

FLAGS = -std=c++11 -O2 -g -Wall -w -Wc++11-compat-deprecated-writable-strings


SRC_DIR=../src
OBJ_DIR=../obj

SRCS=$(wildcard $(SRC_DIR)/*.cpp) $(TARGET).cpp
OBJS=$(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))

INCLUDE_MPICH = -I/opt/homebrew/Cellar/mpich/4.2.1/include/
INCLUDE_EIGEN = -I../zlib/eigen-3.4.0/
INCLUDE_UTHASH = -I../zlib/uthash/include/
INCLUDE_MAIN = -I../include/

INCLUDE_MKL = -I/opt/intel/oneapi/mkl/latest/include/
LIBRARY_MKL = -L/opt/intel/oneapi/mkl/latest/lib/
#LDFLAGS_MKL = -Wl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm #-ldl #-liomp5 

INCLUDE = $(INCLUDE_EIGEN) $(INCLUDE_MKL) $(INCLUDE_UTHASH) $(INCLUDE_MAIN) $(INCLUDE_MPICH)
LIBRARY = $(LIBRARY_MKL) $(LDFLAGS_MKL)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(FLAG) $(INCLUDE) $(LIBRARY)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) -c $< -o $@ $(FLAG) $(INCLUDE) $(LIBRARY)

clean:
	rm -f $(OBJS) $(TARGET)

