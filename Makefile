# settings
DIALECT      = -std=c++17
#OPTIMIZATION = -O0 -g
OPTIMIZATION = -O0 -g
WARNINGS     = -Xcompiler="-Wall -Wextra"
# NVCC_FLAGS   = -arch=sm_61 -lineinfo --expt-relaxed-constexpr -rdc=true
NVCC_FLAGS   = -arch=native -lineinfo --expt-relaxed-constexpr -rdc=true --extended-lambda -lnvToolsExt -Xcompiler="-fopenmp" #-res-usage #-Xptxas "-v"
#LDFLAGS      = -Xcompiler="-pthread -s"  $(NVCC_FLAGS) -lz
LDFLAGS      = -Xcompiler="-pthread"  $(NVCC_FLAGS) -lz
COMPILER     = nvcc
ARTIFACT     = align

MAKEDB = makedb
MODIFYDB = modifydb
GRIDSEARCH = gridsearch

# make targets
.PHONY: clean

release: $(ARTIFACT) $(MAKEDB) #$(MODIFYDB) $(GRIDSEARCH)

clean :
	rm -f *.o
	rm -f $(ARTIFACT)
	rm -f $(MAKEDB)


# compiler call
COMPILE = $(COMPILER) $(NVCC_FLAGS) $(DIALECT) $(OPTIMIZATION) $(WARNINGS) -c $< -o $@


# link object files into executable
$(ARTIFACT): main.o sequence_io.o dbdata.o options.o blosum.o half2_kernel_instantiations.o float_kernel_instantiations.o dpx_s32_kernel_instantiations.o dpx_s16_kernel_instantiations.o
	$(COMPILER) $^ -o $(ARTIFACT) $(LDFLAGS)

$(MAKEDB): makedb.o sequence_io.o dbdata.o
	$(COMPILER) $^ -o $(MAKEDB) $(LDFLAGS)

$(MODIFYDB): modifydb.o sequence_io.o dbdata.o
	$(COMPILER) $^ -o $(MODIFYDB) $(LDFLAGS)

$(GRIDSEARCH): gridsearch.o sequence_io.o dbdata.o blosum.o
	$(COMPILER) $^ -o $(GRIDSEARCH) $(LDFLAGS)

# compile CUDA files
main.o : main.cu sequence_io.h length_partitions.hpp dbdata.hpp cudasw4.cuh new_kernels.cuh convert.cuh float_kernels.cuh half2_kernels.cuh dpx_s16_kernels.cuh blosum.hpp types.hpp
	$(COMPILE)

# compile pure C++ files
sequence_io.o : sequence_io.cpp sequence_io.h
	$(COMPILE)

# compile pure C++ files
dbdata.o : dbdata.cpp dbdata.hpp mapped_file.hpp sequence_io.h length_partitions.hpp
	$(COMPILE)

# compile pure C++ files
options.o : options.cpp options.hpp types.hpp
	$(COMPILE)

blosum.o : blosum.cu blosum.hpp
	$(COMPILE)

half2_kernel_instantiations.o: half2_kernel_instantiations.cu half2_kernels.cuh
	$(COMPILE)

float_kernel_instantiations.o: float_kernel_instantiations.cu float_kernels.cuh
	$(COMPILE)

dpx_s16_kernel_instantiations.o: dpx_s16_kernel_instantiations.cu dpx_s16_kernels.cuh
	$(COMPILE)

dpx_s32_kernel_instantiations.o: dpx_s32_kernel_instantiations.cu dpx_s32_kernels.cuh
	$(COMPILE)

# compile pure C++ files
makedb.o : makedb.cpp dbdata.hpp sequence_io.h
	$(COMPILE)

modifydb.o : modifydb.cpp dbdata.hpp sequence_io.h
	$(COMPILE)

gridsearch.o : gridsearch.cu length_partitions.hpp dbdata.hpp
	$(COMPILE)

