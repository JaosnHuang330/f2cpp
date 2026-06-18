CXX = g++
CXXFLAGS = -O2 -std=c++17 -I. -Iruntime  -MMD -MP
LDFLAGS = 
TARGET = main
SRCS = mod_flux_kfvs.cpp mod_read_gmsh.cpp mod_vector_algebra.cpp mod_point_util.cpp mod_cell_2d.cpp mod_fvm_face_2d.cpp mod_write_vtk.cpp mod_objects.cpp mod_struct_to_array.cpp pre_processing.cpp mod_solver_kfvs.cpp solver.cpp euler2d.cpp runtime/fortran_io.cpp runtime/fortran_namelist.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
.PHONY: all clean
all: $(TARGET)
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
-include $(DEPS)
clean:
	rm -f $(TARGET) $(OBJS) $(DEPS)
