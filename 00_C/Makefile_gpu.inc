# ===== CUDA compiler =====
NVCC = nvcc

# ===== Architecture =====
ARCH = -arch=sm_80

# ===== Host compiler =====
# HOST = -ccbin nvc++
HOST = -ccbin mpicxx

# ===== Optimization =====
OPT = -O0

# ===== Warning suppression =====
SUPPRESS = -diag-suppress=177 \
           -diag-suppress=550 \
           -diag-suppress=780 \
           -Xcompiler=--diag_suppress=declared_but_not_referenced \
           -Xcompiler=--diag_suppress=set_but_not_used \
           -Xcompiler=--diag_suppress=hidden_by_old_for_init

# ===== Combine flags =====
NVCCFLAGS = $(OPT) $(ARCH) $(HOST) $(SUPPRESS)
