# Compiler
CXX = g++
CXXFLAGS = -std=c++20 -O3 -march=native -I$(SRC_DIR)/../include -I$(SRC_DIR)/../utility -I$(CUR_DIR) -I$(HDF5_INCLUDE_DIR)
LDFLAGS= -pthread

# Define the path to your external parameters file.
PARAMS_FILE = ../params.txt

# If the file exists, include it.
ifneq ("$(wildcard $(PARAMS_FILE))", "")
  include $(PARAMS_FILE)
else
  # Default value if params.txt is not found.
  REAL_TYPE = double
endif

# Add the REAL_TYPE definition to the compiler flags.
CXXFLAGS += -DREAL_TYPE="$(REAL_TYPE)"

# Source and Include directories
SRC_DIR = $(CURDIR)/../src
CUR_DIR = $(CURDIR)
INCLUDE_DIRS = -I$(SRC_DIR)/../include -I$(SRC_DIR)/../utility -I$(CUR_DIR)

# Source files
SRCS = $(SRC_DIR)/main.cpp

# Object files (in the logs directory)
OBJDIR = $(CURDIR)/logs
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

# Executable
TARGET = Aphrodite

# Boost settings
BOOST_DIR = /sw-eb/software/Boost/1.81.0-GCC-12.2.0
BOOST_INCLUDE_DIR = $(BOOST_DIR)/include
BOOST_LIB_DIR = $(BOOST_DIR)/lib64

# Add Boost include directory to compiler flags
CXXFLAGS += -I$(BOOST_INCLUDE_DIR)

LIBS = -L$(HDF5_LIB_DIR) -lhdf5 -lhdf5_hl -pthread -L$(BOOST_LIB_DIR)

HDF5_INCLUDE_DIR = /sw-eb/software/HDF5/1.14.0-gompi-2022b/include
HDF5_LIB_DIR = /sw-eb/software/HDF5/1.14.0-gompi-2022b/lib

# Rules
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LIBS)

# Rule for building object files
$(OBJDIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJDIR)/*.o $(TARGET)

.PHONY: all clean