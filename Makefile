# compiler to use
CC = c++
CFLAGS=-O2 -std=c++14 -g

LOCALPATH=/home/nirmalya/local/

INC = -I${LOCALPATH}/include  -I${LOCALPATH}/boost/include

PROG_OPT_LIB=${LOCALPATH}/boost/lib/libboost_program_options.a
LIBDIR=${LOCALPATH}/lib/
LIBS=$(PROG_OPT_LIB)

all: clean tools
    
tools:
	#$(CC) $(INC) $(STXXLINC) $(CFLAGS)  optimization_ex.cpp -o optimization_ex  $(LIBS) -lpthread -lcblas
	#$(CC) $(INC) $(STXXLINC) $(CFLAGS)  NBModelTest.cpp -o NBModelTest  $(LIBS) -lpthread -lcblas
	$(CC) $(INC) $(STXXLINC) $(CFLAGS)  MixtureModel.cpp -o MixtureModel  $(LIBS) -lpthread -lcblas
	#$(CC) $(INC) $(STXXLINC) $(CFLAGS)  NBinomEM.cpp -o NBinomEM  $(LIBS)
    
clean:
	rm -f NBinomEM

