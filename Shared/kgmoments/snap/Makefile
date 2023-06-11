#
# Makefile for non-Microsoft compilers
#

## Linux  (uncomment the 2 lines below for compilation on Linux)
#CXXFLAGS += -std=c++98 -Wall
#LDFLAGS += -lrt

## CygWin (uncomment the 2 lines below for compilation on CygWin)
#CXXFLAGS += -Wall
#LDFLAGS += 

SNAP_DIR=$(HOME)/devextern/snap-2011-04/

## Main application file
MAIN = kronloglike

all: $(MAIN)

opt: CXXFLAGS += -O4
opt: LDFLAGS += -O4
opt: $(MAIN)

CXXFLAGS += -I$(SNAP_DIR)/glib -I$(SNAP_DIR)/snap

# COMPILE
$(MAIN): $(MAIN).cc Snap.o 
	g++ $(LDFLAGS) $(CXXFLAGS) -o $(MAIN) $(MAIN).cc $(SNAP_DIR)/snap/kronecker.cpp Snap.o -lrt
	
Snap.o: 
	g++ -c $(CXXFLAGS) $(SNAP_DIR)/snap/Snap.cpp 

clean:
	rm -f *.o  $(MAIN)  $(MAIN).exe
	rm -rf Debug Release

