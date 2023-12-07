# %---------------------------------%
# |  SECTION 1: COMPILER C++        |
# %---------------------------------%

GPP = $(CXX)

UNAME := $(shell uname -s)
ifeq ($(UNAME),Linux)
        FLAGS= -O2 -Wall -pedantic -I /usr/include/eigen3/   
	FLAGS += -std=c++17 -fopenmp -fPIC
endif
ifeq ($(UNAME),Darwin)
        FLAGS= -O2 -Wall -pedantic -I /usr/local/include/eigen3/   
        FLAGS += -std=c++1z 
endif

#  %--------------------------------%
#  |  SECTION 2: THE ARCHIVER AR    |
#  %--------------------------------%

AR = ar
ARFLAGS = rcs

#  %--------------------------------%
#  |  SECTION 3: COMPILATION RULES  |
#  %--------------------------------%

libqoptkit.so: interface qoptkit.a
	$(GPP) $(FLAGS) -shared  pyqoptkit.o libqoptkit.a -o libqoptkit.so

interface:
	$(GPP)  $(FLAGS) -c pyqoptkit.cpp -o pyqoptkit.o	
	
	
qoptkit.a:   sim.o state.o qocircuit.o util.o
	     $(AR) $(ARFLAGS) libqoptkit.a sim.o state.o qocircuit.o util.o
	     
sim.o:       sim.cpp 
	     $(GPP)  $(FLAGS) -c sim.cpp  
	     	     
state.o:     state.cpp
	     $(GPP)  $(FLAGS) -c state.cpp

qocircuit.o: qocircuit.cpp
	     $(GPP)  $(FLAGS) -c qocircuit.cpp

util.o:      util.cpp
	     $(GPP)  $(FLAGS) -c util.cpp


clean:
	     rm *.so *.o *.a
