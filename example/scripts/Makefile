ROOTINCS = -I$(shell root-config --incdir)
ROOTLIBS := $(shell root-config --libs) $(shell root-config --auxlibs)
UNFOLDLIBS = -L$(UNFOLDDIR)/lib -lResponseMatrix -lUnfold
UNFOLDINCS = -I$(UNFOLDDIR)/include
SYSLIBS = -lg2c -lm -ldl
INCLUDES = $(ROOTINCS) $(UNFOLDINCS)
LIBS = $(ROOTLIBS)

all : histogram_ensemble error_analysis_histograms

CFLAGS = -fPIC -Wall -c -g $(shell root-config --cflags) 
LFLAGS = 

install:

histogram_ensemble : histogram_ensemble.c 
	g++ histogram_ensemble.c -g -O2 $(INCLUDES) $(UNFOLDINCS) $(LIBS) $(UNFOLDLIBS) -o histogram_ensemble

error_analysis_histograms : error_analysis_histograms.c 
	g++ error_analysis_histograms.c -g -O2 $(INCLUDES) $(UNFOLDINCS) $(LIBS) $(UNFOLDLIBS) -o error_analysis_histograms

clean:
	rm -f histogram_ensemble 

