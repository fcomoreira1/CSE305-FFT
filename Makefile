CXX=g++
CXXFLAGS=-std=c++17 -O3 -fopenmp

EXE=test
ODIR=build

SRCS = $(wildcard *.cpp)
DEPS = $(wildcard *.h)
_OBJS = $(subst .cpp,.o,$(SRCS))
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(EXE): $(SRCS) $(DEPS)
	$(CXX) -o $(EXE) $(SRCS) $(CXXFLAGS) $(LIBS)

.PHONY: incremental
incremental: $(OBJS)
	$(CXX) -o $(EXE) $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -r $(OBJS) $(EXE)
