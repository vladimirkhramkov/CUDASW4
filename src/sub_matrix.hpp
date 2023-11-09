#ifndef SUB_MATRIX_HPP
#define SUB_MATRIX_HPP

#include "types.hpp"
#include "util.cuh"
#include <array>
#include <string>
#include <vector>

namespace cudasw4{

#ifdef __CUDACC__

extern __constant__ char deviceSubMatrix[25*25];
extern __constant__ int deviceSubMatrixDim;
extern __constant__ int deviceSubMatrixDimSquared;

#endif

extern char hostSubMatrix[25*25];
extern int hostSubMatrixDim;
extern int hostSubMatrixDimSquared;

//set host and device global substitution matrix variables
void setProgramWideSubMatrix(SubMatrixType subMatrixType, const std::vector<int>& deviceIds);

} //namespace cudasw4

#endif