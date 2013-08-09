SHELL = /bin/sh

CXX     = mpic++

# Please ensure the following environmental variables are defined and point to
# the correct folders on your system before building.
# VTK_INCLUDE_DIR
# VTK_LIB_DIR
# BOOST_INCLUDE_DIR
# BOOST_LIB_DIR

CXXFLAGS = -I./ -I$(VTK_INCLUDE_DIR) -I$(PETSC_INCLUDE_DIR) $(OTHER_INCLUDE) -O3 -Wno-deprecated -std=c++11
CXXFLAGS_DEBUG = $(CXXFLAGS) -pg
LFLAGS = -L/groupvol/sjn/common/muparser/lib\
	 -L/groupvol/sjn/common/spud\
	 -L$(VTK_LIB_DIR)\
	 -L$(PETSC_LIB_DIR)

SPH_LIBS = -lboost_serialize -lboost_system -lboost_mpi
UTR_LIBS = -lvtkIO -lboost_serialize -lboost_system -lboost_mpi -lboost_program_options

MAKE    = make

# executable file names
SPH_TARGET = sph
UTR_TARGET = sphpp

SPH_OBJS = main.o \
	        core/Simulation.o\
	        utils/utils.o\
            kernels/WendlandQuintic.o

UTR_OBJS = main.o \
            Processor.o

SPH_2D_OBJS = $(SPH_OBJS:%=release/%_2D)
SPH_3D_OBJS = $(SPH_OBJS:%=release/%_3D)
UTR_2D_OBJS = $(UTR_OBJS:%=release/%_2D)
UTR_3D_OBJS = $(UTR_OBJS:%=release/%_3D)

SPH_SOURCE = $(SPH_OBJS:%.o=src/%.cpp)
UTR_SOURCE = $(UTR_OBJS:%.o=src/%.cpp)

#.SUFFIXES: .cpp .o
#
#.cpp.o:
#	$(CXX) $(CXXFLAGS) -c $< -o $@

dims2: $(SPH_2D_OBJS)
	@echo Linking 2D
	@$(CXX) $(LFLAGS) -o release/$(SPH_TARGET)_2D $(SPH_2D_OBJS) $(SPH_LIBS)

dims3: $(SPH_3D_OBJS)
	@echo Linking 2D
	@$(CXX) $(LFLAGS) -o release/$(SPH_TARGET)_2D $(SPH_3D_OBJS) $(SPH_LIBS)

release/%.o_2D: src/%.cpp
	@echo Compiling $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@ -DDIM=2

release/%.o_3D: src/%.cpp
	@echo Compiling debug $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@ -DDIM=3

clean:
	-find -name *.o | xargs rm
	-find -name *.o_2D | xargs rm
	-find -name *.o_3D | xargs rm

