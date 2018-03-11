
#CXX = clang++

CXXFLAGS += -O3 -fPIC -std=c++11
#CXXFLAGS += -O3 -g -fPIC -fsanitize=address -std=c++11

# LHAPDF
LHAPDFINCS = $(shell lhapdf-config --cppflags)
LHAPDFLIBS = $(shell lhapdf-config --ldflags)

# APFEL++
APFELPPINCS = $(shell apfelxx-config --cppflags)
APFELPPLIBS = $(shell apfelxx-config --ldflags)

# Now set up the compiler and link flags and libs
CXXFLAGS += $(LHAPDFINCS) $(APFELINCS) $(APFELPPINCS)
LDFLAGS  += $(CXXFLAGS)

CLIBS += $(LHAPDFLIBS) $(APFELLIBS) $(APFELPPLIBS)

install : all
all : interface_test

interface_test: interface_test.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cxx .o .f .c

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

.f.o:	 
	$(F77)  -c $< 

clean:
	rm -rf *.lo *.o *.la interface_test *~

