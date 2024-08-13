HOME = ..
include $(HOME)/Makefile.inc
EXE = a.out  # Executable name without path

OBJDIR = ./obj
SRC_DIR = ./src  # Assuming your source files are in the 'src' directory

# Libraries
# Add CUDA, NVTX, and other libraries as required. Example:
# LIBS = -L/usr/local/cuda/lib64 -lnvToolsExt -lcufft -lcublas -lm

# List of source modules without file extensions
MODULES = para_range mpi_topology geometry matrix poisson_matrix_operator \
          cg_poisson_matrix rbgs_poisson_matrix multigrid_common multigrid_vcycle \
          multigrid poisson

# Generate object file paths
OBJS = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(MODULES)))

all: directories $(EXE)

# Rule to create necessary directories
directories:
	mkdir -p $(OBJDIR)
	mkdir -p $(HOME)/run/result

# Rule for creating the executable
$(EXE): $(OBJS)
	$(FC) -o $(HOME)/run/$@ $(FLAG) $^ $(LIBS)

# Generic rule for compiling Fortran source files, with module output directory
$(OBJDIR)/%.o: %.f90
	$(FC) -c $(FLAG) $(LIBS) $< -o $@ -module $(OBJDIR) -I$(OBJDIR)

# Additional dependencies
$(OBJDIR)/mpi_topology.o: $(OBJDIR)/para_range.o
$(OBJDIR)/geometry.o: $(OBJDIR)/mpi_topology.o
# ... other dependencies based on your module hierarchy ...

clean:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod
	rm -f $(HOME)/run/$(EXE)
	rm -rf $(HOME)/run/result

.PHONY: all clean directories
