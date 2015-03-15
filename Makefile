ROOTINCS = -I$(shell root-config --incdir)
ROOTLIBS := $(shell root-config --libs) $(shell root-config --auxlibs)
UNFOLDLIBS = -L$(UNFOLDDIR)/lib -lResponseMatrix -lUnfold
UNFOLDINCS = -I$(UNFOLDDIR)/include
SYSLIBS = -lg2c -lm -ldl
INCLUDES = $(ROOTINCS) $(UNFOLDINCS)
LIBS = $(ROOTLIBS)

all : lib/libUnfold.so lib/libResponseMatrix.so example/toy_response_matrix example/unfold example/perform_measurement bin/poisson_extraction
# example/unfold example/ml example/bayes bin/analyze_bias bin/poisson_extractions

CFLAGS = -fPIC -Wall -c -g $(shell root-config --cflags) 
LFLAGS = 

install:
	@ mkdir -p lib
	@ mkdir -p bin 
	@ echo "export UNFOLDDIR=$(PWD)" >setup.sh
	@ echo "export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(PWD)/lib" >>setup.sh

bin/poisson_extraction : src/poisson_extraction.c 
	g++ example/poisson_extraction.c -g -O2 $(INCLUDES) $(UNFOLDINCS) $(LIBS) $(UNFOLDLIBS) -o bin/poisson_extraction

bin/analyze_bias : src/analyze_bias.c 
	g++ src/analyze_bias.c -g -O2 $(INCLUDES) $(UNFOLDINCS) $(LIBS) $(UNFOLDLIBS) -o bin/analyze_bias

example/unfold : example/unfold.c include/unfold.h include/response_matrix.h
	g++ example/unfold.c -g -O2 $(INCLUDES) $(UNFOLDINCS) $(LIBS) $(UNFOLDLIBS) -o example/unfold

example/ml : example/ml.c include/unfold.h include/response_matrix.h
	g++ example/ml.c -g -O2 $(INCLUDES) $(UNFOLDINCS) $(LIBS) $(UNFOLDLIBS) -o example/ml

example/perform_measurement : example/perform_measurement.c include/unfold.h include/response_matrix.h
	g++ example/perform_measurement.c -g -O2 $(ROOTINCS) $(ROOTLIBS) -o example/perform_measurement

example/unfold : example/unfold.c include/unfold.h include/response_matrix.h
	g++ example/unfold.c -g -O2 $(INCLUDES) $(UNFOLDINCS) $(LIBS) $(UNFOLDLIBS) -o example/unfold

example/bayes : example/bayes.c include/unfold.h include/response_matrix.h
	g++ example/bayes.c -g -O2 $(INCLUDES) $(UNFOLDINCS) $(LIBS) $(UNFOLDLIBS) -o example/bayes

example/toy_response_matrix : example/toy_response_matrix.c include/response_matrix.h
	g++ example/toy_response_matrix.c -g -O2 $(INCLUDES) $(UNFOLDINCS) $(LIBS) $(UNFOLDLIBS) -o example/toy_response_matrix

bin/bayes : src/bayes.c 
	g++ src/bayes.c -g -O2 $(INCLUDES) $(LIBS) -o bin/bayes

bin/unfold.o : src/unfold.c include/unfold.h
	g++ $(INCLUDES) $(CFLAGS) src/unfold.c -o bin/unfold.o 

bin/elisa.o : src/elisa.c include/unfold.h
	g++ $(INCLUDES) $(CFLAGS) src/elisa.c -o bin/elisa.o 

bin/response_matrix.o : src/response_matrix.c include/response_matrix.h
	g++ $(INCLUDES) $(CFLAGS) src/response_matrix.c -o bin/response_matrix.o 

lib/libUnfold.so : bin/unfold.o bin/elisa.o
	g++ -shared -Wall,-soname,libUnfold.so.1 -lc -o lib/libUnfold.so bin/unfold.o bin/elisa.o 

lib/libResponseMatrix.so : bin/response_matrix.o
	g++ -shared -Wall,-soname,libResponseMatrix.so.1 -lc -o lib/libResponseMatrix.so bin/response_matrix.o 

clean:
	rm -f lib/libUnfold.so
	rm -f lib/libResponseMatrix.so
	rm -f example/toy_response_matrix
	# rm -f example/unfold
	# rm -f example/ml
	# rm -f example/elisa
	# rm -f example/bayes
	rm -f bin/unfold.o
	rm -f bin/response_matrix.o
