CXX = g++
CXXFLAGS = -O2 -std=c++17 -I. -Iruntime  -MMD -MP
LDFLAGS = 
TARGET = euler.x
SRCS = m_params.cpp m_mesh.cpp m_riemann.cpp m_io.cpp m_solver.cpp euler_shock_tube.cpp runtime/fortran_io.cpp
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
