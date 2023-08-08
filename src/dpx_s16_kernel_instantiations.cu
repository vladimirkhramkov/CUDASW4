//#include "util.cuh"
#include "new_kernels.cuh"

//#include "thrust/iterator/counting_iterator.h"

namespace cudasw4{

// using ScoreIterInst = TopNMaximaArray;
// using PosIterInst = decltype(thrust::make_counting_iterator<size_t>(0));
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 2, 24, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 4, 16, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 10, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 12, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 14, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 16, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 18, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 20, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 22, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 24, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 26, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 28, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 30, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 8, 32, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 16, 18, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 16, 20, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 16, 22, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 16, 24, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 16, 26, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 16, 28, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 16, 30, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 16, 32, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 18, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 20, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 22, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 24, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 26, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 28, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 30, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 32, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 34, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 36, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 38, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);
template void call_NW_local_affine_single_pass_s16_DPX_new<256, 32, 40, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);

template void call_NW_local_affine_many_pass_s16_DPX_new<256, 32, 22, ScoreIterInst, PosIterInst>(BlosumType, const char * const, ScoreIterInst const, short2*, short2*, const size_t* const, const SequenceLengthT* const, PosIterInst const, const int, ReferenceIdT* const, int* const, const bool, const char4*, const SequenceLengthT, const int, const int, cudaStream_t);

} // namespace cudasw4