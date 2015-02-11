CC=gcc
CPLEXDIR=/users/ecco/ptigwe/local/ibm/ILOG/CPLEX_Studio126/cplex
CPLEXINC=$(CPLEXDIR)/include
CPLEXLIB=$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)/
CFLAGS=-g -I$(CPLEXINC) -Wall
CPLEXLDFLAGS=-lcplex -lpthread
LDDIR=-L$(CPLEXLIB) -L$(HOME)/local/lib
LDFLAGS=$(CPLEXLDFLAGS) -lgmp -lm -lrt -lpthread

UTILS_DIR=../utils
UTILS_OBJ=$(UTILS_DIR)/matrix.o $(UTILS_DIR)/cplp.o
OBJ=lis.o checker.o main.o

%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@
	
all: $(OBJ)
	$(CC) $(CFLAGS) $(LDDIR) $(OBJ) $(UTILS_OBJ) -o supp_enum  -lm $(LDFLAGS) 

clean:
	rm $(OBJ) supp_enum
