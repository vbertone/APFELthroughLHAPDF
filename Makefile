
CXX = clang++

CXXFLAGS += -O3 -fPIC -std=c++11

# LHAPDF
LHAPDFINCS = $(shell lhapdf-config --cppflags)
LHAPDFLIBS = $(shell lhapdf-config --ldflags)

# APFEL++
APFELPPINCS = $(shell apfelxx-config --cppflags)
APFELPPLIBS = $(shell apfelxx-config --ldflags)

# Now set up the compiler and link flags and libs
CXXFLAGS += $(LHAPDFINCS) $(APFELINCS) $(APFELPPINCS)
LDFLAGS  += $(LHAPDFINCS) $(APFELINCS) $(APFELPPINCS)

CLIBS += $(LHAPDFLIBS) $(APFELLIBS) $(APFELPPLIBS)

install : all
all : apfelxxTest interface_test

apfelxxTest: apfelxxTest.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

interface_test: interface_test.o 
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cxx .o .f .c

.cxx.o:	 
	$(CXX) $(CXXFLAGS) -c $< 

.f.o:	 
	$(F77)  -c $< 

clean:
	rm -rf *.lo *.o *.la apfelxxTest interface_test *~

