ROOT = ../..
include $(ROOT)/Makefile_gpu.inc


# EXE = poisson
EXE = test_poisson

SRCS = poisson_gpu.cu \
	   rbgs_poisson_matrix_gpu.cu \
	   poisson_matrix_operator_gpu.cu \
	   multigrid_gpu.cu \
	   geometry_gpu.cu \
	   matrix_gpu.cu \
	   mpi_topology.c \
	   para_range.c 
# 	   timer.c


# $(EXE): 
# 	$(NVCC) $(NVCCFLAGS) $(SRCS) -I$(ROOT)/include -o $(ROOT)/run/$@ 
# # 	$< $(CUDA_LDFLAGS) 
# # 	-lmpi

# OBJDIR = ./obj
# SRC_DIR = ./src  

# # CUDA compiler
# NVCC = nvcc

# # List of source modules (without extension)
# MODULES = timer para_range mpi_topology geometry matrix multigrid \
#           poisson_matrix_operator rbgs_poisson_matrix poisson

# # Generate object file paths
# OBJS = $(OBJDIR)/timer.o \
# 	   $(OBJDIR)/para_range.o \
#        $(OBJDIR)/mpi_topology.o \
#        $(OBJDIR)/geometry.o \
#        $(OBJDIR)/matrix.o \
#        $(OBJDIR)/multigrid.o \
#        $(OBJDIR)/poisson_matrix_operator.o \
#        $(OBJDIR)/rbgs_poisson_matrix.o \
#        $(OBJDIR)/poisson.o

all: directories $(EXE)

directories:
	mkdir -p $(ROOT)/run

$(EXE): $(SRCS)
	$(NVCC) $(NVCCFLAGS) $(SRCS) \
	-I$(ROOT)/include \
	-o $(ROOT)/run/$@

# # Rule to create necessary directories
# directories:
# 	mkdir -p $(OBJDIR)
# 	mkdir -p $(HOME)/run
# 	mkdir -p $(HOME)/run/result

# # Rule for creating the executable
# $(EXE): $(OBJS)
# 	$(CC) $(CFLAGS) -I$(HOME)/include -o $(HOME)/run/$@ $(OBJS) $(LDFLAGS)

# # Compile each .c into obj/xxx.o，
# $(OBJDIR)/%.o: %.c
# 	$(CC) $(CFLAGS) -I$(HOME)/include -c $< -o $@ 

clean:
# 	rm -rf $(OBJDIR)/*.o
	rm -f $(ROOT)/run/$(EXE)
	rm -rf $(ROOT)/run/result

.PHONY: all clean directories
