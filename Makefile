# compiler to use
CC = c++
CFLAGS=-O2 -std=c++14 -g

LOCALPATH=/home/BROAD.MIT.EDU/nirmalya/local/

INC = -I${LOCALPATH}/include

PROG_OPT_LIB=${LOCALPATH}/lib/libboost_program_options.a
LIBDIR=${LOCALPATH}/lib/
LIBS=$(PROG_OPT_LIB)

all: clean tools
    
tools:
	$(CC) $(INC) $(STXXLINC) $(CFLAGS)  NBinomEM.cpp -o NBinomEM  $(LIBS)
    
clean:
	rm -f NBinomEM

