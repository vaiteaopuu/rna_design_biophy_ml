CC=gcc
BIN=./bin

OBJS3 = \
	mcmc.o\
	seq_utils.o

OBJS4 = \
	mcmc_prof.o\
	seq_utils.o

# Adjust these paths to where ViennaRNA is installed
INCLUDE_DIR = ~/programs/vrna/include/
LIB_DIR = ~/programs/vrna/lib/

# Include directories for header files
INCLUDES = -I$(INCLUDE_DIR)

# Library directories for linking
LDFLAGS = -L$(LIB_DIR) -lRNA -lm -fopenmp -lstdc++ -g

# Compilation flags
CFLAGS = -Ofast -march=native -Wall -Wextra -DHAVE_INLINE $(INCLUDES)

all: mcmc

mcmc: $(OBJS3)
	$(CC) $(OBJS3) $(CFLAGS) $(LDFLAGS) -o $(BIN)/dca

mcmc_prof: $(OBJS4)
	$(CC) $(OBJS3) $(CFLAGS) $(LDFLAGS) -o $(BIN)/dca_prof

clean:
	rm -rf *.o
