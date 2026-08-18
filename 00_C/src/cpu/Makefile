HOME = ../..
include $(HOME)/Makefile_cpu.inc

EXE = poisson

OBJDIR = ./obj
SRC_DIR = ./src  

# List of source modules (without extension)
MODULES = timer para_range mpi_topology geometry matrix multigrid \
          poisson_matrix_operator rbgs_poisson_matrix poisson

# Generate object file paths
OBJS = $(OBJDIR)/timer.o \
	   $(OBJDIR)/para_range.o \
       $(OBJDIR)/mpi_topology.o \
       $(OBJDIR)/geometry.o \
       $(OBJDIR)/matrix.o \
       $(OBJDIR)/multigrid.o \
       $(OBJDIR)/poisson_matrix_operator.o \
       $(OBJDIR)/rbgs_poisson_matrix.o \
       $(OBJDIR)/poisson.o

all: directories $(EXE)

# Rule to create necessary directories
directories:
	mkdir -p $(OBJDIR)
	mkdir -p $(HOME)/run
	mkdir -p $(HOME)/run/result

# Rule for creating the executable
$(EXE): $(OBJS)
	$(CC) $(CFLAGS) -I$(HOME)/include -o $(HOME)/run/$@ $(OBJS) $(LDFLAGS)

# Compile each .c into obj/xxx.o，记得带上 CFLAGS
$(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -I$(HOME)/include -c $< -o $@ 

clean:
	rm -rf $(OBJDIR)/*.o
	rm -f $(HOME)/run/$(EXE)
	rm -rf $(HOME)/run/result

.PHONY: all clean directories

# module load gcc/15.2.0 
# module load mpi/openmpi-4.1.8 