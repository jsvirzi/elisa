ROOTINCS = -I$(shell root-config --incdir)
ROOTLIBS := $(shell root-config --libs) $(shell root-config --auxlibs)
UNFOLDLIBS = -L$(UNFOLDDIR)/lib -lResponseMatrix -lUnfold
UNFOLDINCS = -I$(UNFOLDDIR)/include
SYSLIBS = -lg2c -lm -ldl
INCLUDES = $(ROOTINCS) $(UNFOLDINCS)
LIBS = $(ROOTLIBS)

all : histogram_ensemble 

CFLAGS = -fPIC -Wall -c -g $(shell root-config --cflags) 
LFLAGS = 

install:

histogram_ensemble : histogram_ensemble.c 
	g++ histogram_ensemble.c -g -O2 $(INCLUDES) $(UNFOLDINCS) $(LIBS) $(UNFOLDLIBS) -o histogram_ensemble

clean:
	rm -f histogram_ensemble 

