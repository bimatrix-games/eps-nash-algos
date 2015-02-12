SYSTEM= x86-64_Linux
LIBFORMAT=static_pic
CC=gcc
CPLEXDIR=/users/ecco/ptigwe/local/ibm/ILOG/CPLEX_Studio126/cplex
CPLEXINC=$(CPLEXDIR)/include
CPLEXLIB=$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)/
CFLAGS=-g -I$(CPLEXINC) -Wall
CPLEXLDFLAGS=-lcplex -lpthread
LDDIR=-L$(CPLEXLIB) -L$(HOME)/local/lib
LDFLAGS=$(CPLEXLDFLAGS) -lgmp -lm -lrt -lpthread

UTILS_DIR=../utils
UTILS_OBJ=$(UTILS_DIR)/cplp.o $(UTILS_DIR)/matrix.o $(UTILS_DIR)/io.o $(UTILS_DIR)/util.o 
OBJ=ts.o main.o

%.o : %.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
	
all: $(OBJ)
	$(CC) $(CFLAGS) $(LDDIR) $(OBJ) $(UTILS_OBJ) -o ts $(LDFLAGS)

clean:
	rm $(OBJ) ts
