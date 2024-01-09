#include "sub_matrix.hpp"
#include "util.cuh"

#include <cassert>

namespace cudasw4{

    #ifdef __CUDACC__
    __constant__ char deviceSubMatrix[25*25];
    __constant__ int deviceSubMatrixDim;
    __constant__ int deviceSubMatrixDimSquared;
    #endif
    
    char hostSubMatrix[25*25];
    int hostSubMatrixDim;
    int hostSubMatrixDimSquared;
    
    //set host and device global variables
    
    
    void setProgramWideSubMatrix(SubMatrixType subMatrixType, const std::vector<int>& deviceIds){
        switch(subMatrixType){
            case SubMatrixType::DNA:
                {
                    const auto subMatrix = DNA::get1D();
                    const int dim = DNA::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);
                }
                break;               
            case SubMatrixType::NUC44:
                {
                    const auto subMatrix = NUC44::get1D();
                    const int dim = NUC44::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;            
            case SubMatrixType::PAM30:
                {
                    const auto subMatrix = PAM30::get1D();
                    const int dim = PAM30::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;            
            case SubMatrixType::PAM70:
                {
                    const auto subMatrix = PAM70::get1D();
                    const int dim = PAM70::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;            
            case SubMatrixType::BLOSUM45:
                {
                    const auto subMatrix = BLOSUM45::get1D();
                    const int dim = BLOSUM45::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;
            case SubMatrixType::BLOSUM50:
                {
                    const auto subMatrix = BLOSUM50::get1D();
                    const int dim = BLOSUM50::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;
            case SubMatrixType::BLOSUM62:
                {
                    const auto subMatrix = BLOSUM62::get1D();
                    const int dim = BLOSUM62::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;
            case SubMatrixType::BLOSUM80:
                {
                    const auto subMatrix = BLOSUM80::get1D();
                    const int dim = BLOSUM80::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;
            case SubMatrixType::PAM30_20:
                {
                    const auto subMatrix = PAM30_20::get1D();
                    const int dim = PAM30_20::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;                
            case SubMatrixType::PAM70_20:
                {
                    const auto subMatrix = PAM70_20::get1D();
                    const int dim = PAM70_20::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;                
            case SubMatrixType::BLOSUM45_20:
                {
                    const auto subMatrix = BLOSUM45_20::get1D();
                    const int dim = BLOSUM45_20::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;
            case SubMatrixType::BLOSUM50_20:
                {
                    const auto subMatrix = BLOSUM50_20::get1D();
                    const int dim = BLOSUM50_20::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;
            case SubMatrixType::BLOSUM62_20:
                {
                    const auto subMatrix = BLOSUM62_20::get1D();
                    const int dim = BLOSUM62_20::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;
            case SubMatrixType::BLOSUM80_20:
                {
                    const auto subMatrix = BLOSUM80_20::get1D();
                    const int dim = BLOSUM80_20::dim;
                    hostSubMatrixDim = dim;
                    hostSubMatrixDimSquared = dim * dim;
                    auto it = std::copy(subMatrix.begin(), subMatrix.end(), hostSubMatrix);
                    assert(std::distance(hostSubMatrix, it) <= 25 * 25);                
                }
                break;
            default:
                assert(false && "unimplemented substitution matrix copy");
                break;
        }
    #ifdef __CUDACC__
        RevertDeviceId rdi{};
    
        int numGpus = deviceIds.size();
    
        for(int gpu = 0; gpu < numGpus; gpu++){
            cudaSetDevice(deviceIds[gpu]); CUERR;
            cudaMemcpyToSymbol(deviceSubMatrix, &(hostSubMatrix[0]), sizeof(char) * hostSubMatrixDim * hostSubMatrixDim); CUERR;
            cudaMemcpyToSymbol(deviceSubMatrixDim, &hostSubMatrixDim, sizeof(int)); CUERR;
            cudaMemcpyToSymbol(deviceSubMatrixDimSquared, &hostSubMatrixDimSquared, sizeof(int)); CUERR;
        }
    #endif    
    }

} //namespace cudasw4