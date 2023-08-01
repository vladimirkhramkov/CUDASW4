#ifndef CUDASW4_CUH
#define CUDASW4_CUH

#include "hpc_helpers/cuda_raiiwrappers.cuh"
#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/nvtx_markers.cuh"
#include "hpc_helpers/simple_allocation.cuh"

#include "sequence_io.h"
#include "kseqpp/kseqpp.hpp"
#include "dbdata.hpp"
#include "length_partitions.hpp"
#include "util.cuh"
#include "new_kernels.cuh"
#include "blosum.hpp"
#include "types.hpp"
#include "dbbatching.cuh"
#include "convert.cuh"

#include <thrust/binary_search.h>
#include <thrust/sort.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <string_view>

namespace cudasw4{

    template<class T, T factor>
    struct RoundToNextMultiple{
        __host__ __device__ 
        T operator()(const T& value){
            return SDIV(value, factor) * factor;
        }        
    };

    __global__
    void addKernel(int* output, const int* input1, const int* input2){
        *output = *input1 + *input2;
    }


    __global__
    void sumNumOverflowsKernel(int* output, const int* input, int numGpus){
        int sum = input[0];
        for(int gpu = 1; gpu < numGpus; gpu++){
            sum += input[gpu];
        }
        output[0] = sum;
    }

    template<class PartitionOffsets, class Indices>
    __global__
    void transformLocalSequenceIndicesToGlobalIndices(
        int gpu,
        int N,
        PartitionOffsets partitionOffsets,
        Indices maxReduceArrayIndices
    ){
        const int tid = threadIdx.x + blockIdx.x * blockDim.x;

        if(tid < N){
            maxReduceArrayIndices[tid] = partitionOffsets.getGlobalIndex(gpu, maxReduceArrayIndices[tid]);
        }
    }

    struct BenchmarkStats{
        int numOverflows{};
        double seconds{};
        double gcups{};
    };

    struct ScanResult{
        std::vector<int> scores{};
        std::vector<size_t> referenceIds{};
        BenchmarkStats stats{};
    };

    struct KernelTypeConfig{
        KernelType singlePassType = KernelType::Half2;
        KernelType manyPassType_small = KernelType::Half2;
        KernelType manyPassType_large = KernelType::Float;
        KernelType overflowType = KernelType::Float;
    };

    struct MemoryConfig{
        size_t maxBatchBytes = 128ull * 1024ull * 1024ull;
        size_t maxBatchSequences = 10'000'000;
        size_t maxTempBytes = 4ull * 1024ull * 1024ull * 1024ull;
        size_t maxGpuMem = std::numeric_limits<size_t>::max();
    };

    template<size_t size>
    struct TopNMaximaArray{
        struct Ref{
            size_t index;
            size_t indexOffset;
            int* d_locks;
            volatile float* d_scores;
            size_t* d_indices;

            __device__
            Ref& operator=(float newscore){            

                const size_t slot = (indexOffset + index) % size;

                // if(index + indexOffset == 51766021){
                //     printf("Ref operator=(%f), index %lu indexOffset %lu, slot %lu, griddimx %d, blockdimx %d, blockIdxx %d, threadidxx %d\n", 
                //         newscore, index, indexOffset, slot, gridDim.x, blockDim.x, blockIdx.x, threadIdx.x);
                // }

                int* const lock = &d_locks[slot];

                while (0 != (atomicCAS(lock, 0, 1))) {}

                const float currentScore = d_scores[slot];
                if(currentScore < newscore){
                    d_scores[slot] = newscore;
                    d_indices[slot] = indexOffset + index;
                }

                atomicExch(lock, 0);
            }
        };

        TopNMaximaArray(float* d_scores_, size_t* d_indices_, int* d_locks_, size_t offset)
            : indexOffset(offset), d_locks(d_locks_), d_scores(d_scores_), d_indices(d_indices_){}

        template<class Index>
        __device__
        Ref operator[](Index index) const{
            Ref r;
            r.index = index;
            r.indexOffset = indexOffset;
            r.d_locks = d_locks;
            r.d_scores = d_scores;
            r.d_indices = d_indices;
            return r;
        }

        size_t indexOffset = 0;
        int* d_locks;
        volatile float* d_scores;
        size_t* d_indices;
    };

    struct HostGpuPartitionOffsets{
        int numGpus;
        int numLengthPartitions;
        std::vector<size_t> partitionSizes;
        std::vector<size_t> horizontalPS;
        std::vector<size_t> verticalPS;
        std::vector<size_t> totalPerLengthPartitionPS;

        HostGpuPartitionOffsets(int numGpus_, int numLengthpartitions_, std::vector<size_t> partitionSizes_)
            : numGpus(numGpus_), 
            numLengthPartitions(numLengthpartitions_), 
            partitionSizes(std::move(partitionSizes_)),
            horizontalPS(numGpus * numLengthPartitions, 0),
            verticalPS(numGpus * numLengthPartitions, 0),
            totalPerLengthPartitionPS(numLengthPartitions, 0)
        {
            assert(partitionSizes.size() == numGpus * numLengthPartitions);

            for(int gpu = 0; gpu < numGpus; gpu++){
                for(int l = 1; l < numLengthPartitions; l++){
                    horizontalPS[gpu * numLengthPartitions + l] = horizontalPS[gpu * numLengthPartitions + l-1] + partitionSizes[gpu * numLengthPartitions + l-1];
                }
            }
            for(int l = 0; l < numLengthPartitions; l++){
                for(int gpu = 1; gpu < numGpus; gpu++){
                    verticalPS[gpu * numLengthPartitions + l] = verticalPS[(gpu-1) * numLengthPartitions + l] + partitionSizes[(gpu-1) * numLengthPartitions + l];
                }
            }
            for(int l = 1; l < numLengthPartitions; l++){
                totalPerLengthPartitionPS[l] = totalPerLengthPartitionPS[l-1] 
                    + (verticalPS[(numGpus - 1) * numLengthPartitions + (l-1)] + partitionSizes[(numGpus-1) * numLengthPartitions + (l-1)]);
            }
        }

        size_t getGlobalIndex(int gpu, size_t localIndex) const {
            const size_t* const myHorizontalPS = horizontalPS.data() + gpu * numLengthPartitions;
            const auto it = std::lower_bound(myHorizontalPS, myHorizontalPS + numLengthPartitions, localIndex+1);
            const int whichPartition = std::distance(myHorizontalPS, it) - 1;
            const size_t occurenceInPartition = localIndex - myHorizontalPS[whichPartition];
            const size_t globalPartitionBegin = totalPerLengthPartitionPS[whichPartition];
            const size_t elementsOfOtherPreviousGpusInPartition = verticalPS[gpu * numLengthPartitions + whichPartition];
            //std::cout << "whichPartition " << whichPartition << ", occurenceInPartition " << occurenceInPartition 
            //    << ", globalPartitionBegin " << globalPartitionBegin << ", elementsOfOtherPreviousGpusInPartition " << elementsOfOtherPreviousGpusInPartition << "\n";
            return globalPartitionBegin + elementsOfOtherPreviousGpusInPartition + occurenceInPartition;
        };

        void print(std::ostream& os){
            os << "numGpus " << numGpus << "\n";
            os << "numLengthPartitions " << numLengthPartitions << "\n";
            os << "partitionSizes\n";
            for(int gpu = 0; gpu < numGpus; gpu++){
                for(int l = 0; l < numLengthPartitions; l++){
                    os << partitionSizes[gpu * numLengthPartitions + l] << " ";
                }
                os << "\n";
            }
            os << "\n";
            os << "horizontalPS\n";
            for(int gpu = 0; gpu < numGpus; gpu++){
                for(int l = 0; l < numLengthPartitions; l++){
                    os << horizontalPS[gpu * numLengthPartitions + l] << " ";
                }
                os << "\n";
            }
            os << "\n";
            os << "verticalPS\n";
            for(int gpu = 0; gpu < numGpus; gpu++){
                for(int l = 0; l < numLengthPartitions; l++){
                    os << verticalPS[gpu * numLengthPartitions + l] << " ";
                }
                os << "\n";
            }
            os << "\n";
            os << "totalPerLengthPartitionPS\n";
            for(int l = 0; l < numLengthPartitions; l++){
                os << totalPerLengthPartitionPS[l] << " ";
            }
            os << "\n";
        }
    };

    struct DeviceGpuPartitionOffsets{
        template<class T>
        using MyDeviceBuffer = helpers::SimpleAllocationDevice<T, 0>;

        int numGpus;
        int numLengthPartitions;
        MyDeviceBuffer<size_t> partitionSizes;
        MyDeviceBuffer<size_t> horizontalPS;
        MyDeviceBuffer<size_t> verticalPS;
        MyDeviceBuffer<size_t> totalPerLengthPartitionPS;

        struct View{
            int numGpus;
            int numLengthPartitions;
            const size_t* partitionSizes;
            const size_t* horizontalPS;
            const size_t* verticalPS;
            const size_t* totalPerLengthPartitionPS;

            __device__
            size_t getGlobalIndex(int gpu, size_t localIndex) const {
                const size_t* const myHorizontalPS = horizontalPS + gpu * numLengthPartitions;
                const auto it = thrust::lower_bound(thrust::seq, myHorizontalPS, myHorizontalPS + numLengthPartitions, localIndex+1);
                const int whichPartition = thrust::distance(myHorizontalPS, it) - 1;
                const size_t occurenceInPartition = localIndex - myHorizontalPS[whichPartition];
                const size_t globalPartitionBegin = totalPerLengthPartitionPS[whichPartition];
                const size_t elementsOfOtherPreviousGpusInPartition = verticalPS[gpu * numLengthPartitions + whichPartition];
                return globalPartitionBegin + elementsOfOtherPreviousGpusInPartition + occurenceInPartition;
            };
        };

        DeviceGpuPartitionOffsets() = default;
        DeviceGpuPartitionOffsets(const HostGpuPartitionOffsets& hostData)
            : numGpus(hostData.numGpus),
            numLengthPartitions(hostData.numLengthPartitions),
            partitionSizes(numGpus * numLengthPartitions),
            horizontalPS(numGpus * numLengthPartitions),
            verticalPS(numGpus * numLengthPartitions),
            totalPerLengthPartitionPS(numLengthPartitions)
        {
            cudaMemcpyAsync(partitionSizes.data(), hostData.partitionSizes.data(), sizeof(size_t) * numGpus * numLengthPartitions, H2D, cudaStreamLegacy); CUERR;
            cudaMemcpyAsync(horizontalPS.data(), hostData.horizontalPS.data(), sizeof(size_t) * numGpus * numLengthPartitions, H2D, cudaStreamLegacy); CUERR;
            cudaMemcpyAsync(verticalPS.data(), hostData.verticalPS.data(), sizeof(size_t) * numGpus * numLengthPartitions, H2D, cudaStreamLegacy); CUERR;
            cudaMemcpyAsync(totalPerLengthPartitionPS.data(), hostData.totalPerLengthPartitionPS.data(), sizeof(size_t) * numLengthPartitions, H2D, cudaStreamLegacy); CUERR;
        }

        View getDeviceView() const{
            View view;
            view.numGpus = numGpus;
            view.numLengthPartitions = numLengthPartitions;
            view.partitionSizes = partitionSizes.data();
            view.horizontalPS = horizontalPS.data();
            view.verticalPS = verticalPS.data();
            view.totalPerLengthPartitionPS = totalPerLengthPartitionPS.data();
            return view;
        }
    };


    class CudaSW4{
    public:
        template<class T>
        using MyPinnedBuffer = helpers::SimpleAllocationPinnedHost<T, 0>;
        template<class T>
        using MyDeviceBuffer = helpers::SimpleAllocationDevice<T, 0>;

        static constexpr int maxReduceArraySize = 512 * 1024;

        struct GpuWorkingSet{

            using MaxReduceArray = TopNMaximaArray<maxReduceArraySize>;
        
            GpuWorkingSet(
                size_t gpumemlimit,
                size_t maxBatchBytes,
                size_t maxBatchSequences,
                size_t maxTempBytes,
                const std::vector<DBdataView>& dbPartitions
            ){
                cudaGetDevice(&deviceId);
        
                size_t numSubjects = 0;
                size_t numSubjectBytes = 0;
                for(const auto& p : dbPartitions){
                    numSubjects += p.numSequences();
                    numSubjectBytes += p.numChars();
                }
        
                d_query.resize(1024*1024); CUERR
        
                size_t usedGpuMem = 0;
                usedGpuMem += sizeof(int) * maxReduceArraySize; // d_maxReduceArrayLocks
                usedGpuMem += sizeof(float) * maxReduceArraySize; // d_maxReduceArrayScores
                usedGpuMem += sizeof(size_t) * maxReduceArraySize; // d_maxReduceArrayIndices
        
                d_maxReduceArrayLocks.resize(maxReduceArraySize);
                d_maxReduceArrayScores.resize(maxReduceArraySize);
                d_maxReduceArrayIndices.resize(maxReduceArraySize);
        
                cudaMemset(d_maxReduceArrayLocks.data(), 0, sizeof(int) * maxReduceArraySize);
        
                //devAlignmentScoresFloat.resize(numSubjects);
                Fillchar.resize(16*512);
                cudaMemset(Fillchar.data(), 20, 16*512);
        
                forkStreamEvent = CudaEvent{cudaEventDisableTiming}; CUERR;
                numWorkStreamsWithoutTemp = 10;
                workstreamIndex = 0;
                workStreamsWithoutTemp.resize(numWorkStreamsWithoutTemp);
        
                size_t memoryRequiredForFullDB = 0;
                memoryRequiredForFullDB += numSubjectBytes; // d_fulldb_chardata
                memoryRequiredForFullDB += sizeof(size_t) * numSubjects; //d_fulldb_lengthdata
                memoryRequiredForFullDB += sizeof(size_t) * (numSubjects+1); //d_fulldb_offsetdata
                memoryRequiredForFullDB += sizeof(size_t) * numSubjects * 2; //d_overflow_positions_vec
        
                size_t memoryRequiredForBatchedProcessing = 0;
                memoryRequiredForBatchedProcessing += maxBatchBytes * 2; // d_chardata_vec
                memoryRequiredForBatchedProcessing += sizeof(size_t) * maxBatchSequences * 2; //d_lengthdata_vec
                memoryRequiredForBatchedProcessing += sizeof(size_t) * (maxBatchSequences+1) * 2; //d_offsetdata_vec
                memoryRequiredForBatchedProcessing += sizeof(size_t) * maxBatchSequences * 2; //d_overflow_positions_vec
        
                const bool hasEnoughMemoryForFullDB = usedGpuMem + memoryRequiredForFullDB + maxTempBytes <= gpumemlimit;
                const bool hasEnoughMemoryForBatchedDB = usedGpuMem + memoryRequiredForBatchedProcessing + maxTempBytes <= gpumemlimit;
        
                numCopyBuffers = 2;
        
                h_chardata_vec.resize(numCopyBuffers);
                h_lengthdata_vec.resize(numCopyBuffers);
                h_offsetdata_vec.resize(numCopyBuffers);
                d_chardata_vec.resize(numCopyBuffers);
                d_lengthdata_vec.resize(numCopyBuffers);
                d_offsetdata_vec.resize(numCopyBuffers);
                copyStreams.resize(numCopyBuffers);
                pinnedBufferEvents.resize(numCopyBuffers);
                deviceBufferEvents.resize(numCopyBuffers);
                d_total_overflow_number.resize(1);
                d_overflow_number.resize(numCopyBuffers);
                h_overflow_number.resize(numCopyBuffers);
                d_overflow_positions_vec.resize(numCopyBuffers);
        
                auto allocateForFullDBProcessing = [&](){
                    d_fulldb_chardata.resize(numSubjectBytes);
                    d_fulldb_lengthdata.resize(numSubjects);
                    d_fulldb_offsetdata.resize(numSubjects+1);
                    
                    //d_selectedPositions.resize(numSubjects);
                    
                    for(int i = 0; i < numCopyBuffers; i++){
                        h_chardata_vec[i].resize(maxBatchBytes);
                        h_lengthdata_vec[i].resize(maxBatchSequences);
                        h_offsetdata_vec[i].resize(maxBatchSequences+1);               
                        pinnedBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        deviceBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        d_overflow_positions_vec[i].resize(numSubjects);
                    }
                    canStoreFullDB = true;
                    usedGpuMem += memoryRequiredForFullDB;
                };
        
                auto allocateForBatchedProcessing = [&](){
                    //d_selectedPositions.resize(maxBatchSequences);
        
                    for(int i = 0; i < numCopyBuffers; i++){
                        h_chardata_vec[i].resize(maxBatchBytes);
                        h_lengthdata_vec[i].resize(maxBatchSequences);
                        h_offsetdata_vec[i].resize(maxBatchSequences+1);
                        d_chardata_vec[i].resize(maxBatchBytes);
                        d_lengthdata_vec[i].resize(maxBatchSequences);
                        d_offsetdata_vec[i].resize(maxBatchSequences+1);
                        pinnedBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        deviceBufferEvents[i] = CudaEvent{cudaEventDisableTiming}; CUERR;
                        d_overflow_positions_vec[i].resize(maxBatchSequences);
                    }
                    canStoreFullDB = false;
                    usedGpuMem += memoryRequiredForBatchedProcessing;
                };
                
        
                if(hasEnoughMemoryForFullDB){
                    try{
                        allocateForFullDBProcessing();
                    }catch(...){
                        if(hasEnoughMemoryForBatchedDB){
                            allocateForBatchedProcessing();
                        }else{
                            throw std::runtime_error("Not enough GPU memory available on device id " + std::to_string(deviceId));
                        }   
                    }
                }else{
                    if(hasEnoughMemoryForBatchedDB){
                        allocateForBatchedProcessing();
                    }else{
                        throw std::runtime_error("Not enough GPU memory available on device id " + std::to_string(deviceId));
                    }
                }        
                // thrust::sequence(
                //     thrust::device,
                //     d_selectedPositions.begin(),
                //     d_selectedPositions.end(),
                //     size_t(0)
                // );
                //allocations for db processing did not actually
        
                assert(usedGpuMem <= gpumemlimit);
        
                numTempBytes = std::min(maxTempBytes, gpumemlimit - usedGpuMem);
        
                d_tempStorageHE.resize(numTempBytes);
            }
        
            MaxReduceArray getMaxReduceArray(size_t offset){
                return MaxReduceArray(
                    d_maxReduceArrayScores.data(), 
                    d_maxReduceArrayIndices.data(), 
                    d_maxReduceArrayLocks.data(), 
                    offset
                );
            }
        
            void resetMaxReduceArray(cudaStream_t stream){
                thrust::fill(thrust::cuda::par_nosync.on(stream),
                    d_maxReduceArrayScores.data(),
                    d_maxReduceArrayScores.data() + maxReduceArraySize,
                    -1.f
                );
                cudaMemsetAsync(d_maxReduceArrayIndices.data(), 0, sizeof(int) * maxReduceArraySize, stream);
            }
        
            void setPartitionOffsets(const HostGpuPartitionOffsets& offsets){
                deviceGpuPartitionOffsets = DeviceGpuPartitionOffsets(offsets);
            }        
        
            bool singleBatchDBisOnGpu = false;
            int deviceId;
            int numCopyBuffers;
            int numWorkStreamsWithoutTemp = 1;
            int workstreamIndex;
            int copyBufferIndex = 0;
            size_t numTempBytes;
        
            MyDeviceBuffer<int> d_maxReduceArrayLocks;
            MyDeviceBuffer<float> d_maxReduceArrayScores;
            MyDeviceBuffer<size_t> d_maxReduceArrayIndices;
        
            MyDeviceBuffer<char> d_query;
            MyDeviceBuffer<char> d_tempStorageHE;
            //MyDeviceBuffer<float> devAlignmentScoresFloat;
            MyDeviceBuffer<char> Fillchar;
            MyDeviceBuffer<size_t> d_selectedPositions;
            MyDeviceBuffer<int> d_total_overflow_number;
            MyDeviceBuffer<int> d_overflow_number;
            MyPinnedBuffer<int> h_overflow_number;
            CudaStream hostFuncStream;
            CudaStream workStreamForTempUsage;
            CudaEvent forkStreamEvent;
        
            bool canStoreFullDB = false;
            bool fullDBisUploaded = false;
            MyDeviceBuffer<char> d_fulldb_chardata;
            MyDeviceBuffer<size_t> d_fulldb_lengthdata;
            MyDeviceBuffer<size_t> d_fulldb_offsetdata;
            
            std::vector<MyPinnedBuffer<char>> h_chardata_vec;
            std::vector<MyPinnedBuffer<size_t>> h_lengthdata_vec;
            std::vector<MyPinnedBuffer<size_t>> h_offsetdata_vec;
            std::vector<MyDeviceBuffer<char>> d_chardata_vec;
            std::vector<MyDeviceBuffer<size_t>> d_lengthdata_vec;
            std::vector<MyDeviceBuffer<size_t>> d_offsetdata_vec;
            std::vector<CudaStream> copyStreams;
            std::vector<CudaEvent> pinnedBufferEvents;
            std::vector<CudaEvent> deviceBufferEvents;
            std::vector<CudaStream> workStreamsWithoutTemp;
            std::vector<MyDeviceBuffer<size_t>> d_overflow_positions_vec;
        
            DeviceGpuPartitionOffsets deviceGpuPartitionOffsets;
        
        };

        struct SequenceLengthStatistics{
            size_t max_length = 0;
            size_t min_length = std::numeric_limits<size_t>::max();
            size_t sumOfLengths = 0;
        };

    public:
        CudaSW4(
            std::vector<int> deviceIds_, 
            int numTop,
            BlosumType blosumType,
            const KernelTypeConfig& kernelTypeConfig,
            const MemoryConfig& memoryConfig,
            bool verbose_
        ) : deviceIds(std::move(deviceIds_)), verbose(verbose_)
        {
            if(deviceIds.size() == 0){ 
                throw std::runtime_error("No device selected");
            
            }
            RevertDeviceId rdi{};

            initializeGpus();

            const int numDBChunks = 1;
            resultNumOverflows.resize(numDBChunks);

            const int numGpus = deviceIds.size();
            cudaSetDevice(deviceIds[0]);
            devAllAlignmentScoresFloat.resize(maxReduceArraySize * numGpus);
            dev_sorted_indices.resize(maxReduceArraySize * numGpus);
            d_resultNumOverflows.resize(numGpus);
            scanTimer = std::make_unique<helpers::GpuTimer>("Scan");
            totalTimer = std::make_unique<helpers::GpuTimer>("Total");

            setBlosum(blosumType);
            setNumTop(numTop);
            setKernelTypeConfig(kernelTypeConfig);
            setMemoryConfig(memoryConfig);

            dbIsReady = false;
        }

        CudaSW4() = delete;
        CudaSW4(const CudaSW4&) = delete;
        CudaSW4(CudaSW4&&) = default;
        CudaSW4& operator=(const CudaSW4&) = delete;
        CudaSW4& operator=(CudaSW4&&) = default;

        void setGapOpenScore(int score){
            if(score >= 0){
                std::cout << "Warning, gap open score set to non-negative value. Is this intended?\n";
            }
            gop = score;
        }
        void setGapExtendScore(int score){
            if(score >= 0){
                std::cout << "Warning, gap extend score set to non-negative value. Is this intended?\n";
            }
            gex = score;
        }

        void setDatabase(std::shared_ptr<DB> dbPtr){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReady();
        }

        void setDatabase(std::shared_ptr<PseudoDB> dbPtr){
            RevertDeviceId rdi{};
            fullDB = AnyDBWrapper(dbPtr);
            makeReady();
        }

        void setBlosum(BlosumType blosumType){
            setProgramWideBlosum(blosumType, deviceIds);
        }

        void setNumTop(int value){
            if(value >= 0){
                numTop = value;
                updateNumResultsPerQuery();

                const int numDBChunks = 1;
                alignment_scores_float.resize(numDBChunks * results_per_query);
                sorted_indices.resize(numDBChunks * results_per_query);
                resultDbChunkIndices.resize(numDBChunks * results_per_query);
            }
        }

        void setKernelTypeConfig(const KernelTypeConfig& val){
            if(!isValidSinglePassType(val.singlePassType)){
                throw std::runtime_error("Invalid singlepass kernel type");
            }

            if(!isValidMultiPassType_small(val.manyPassType_small)){
                throw std::runtime_error("Invalid manyPassType_small kernel type");
            }

            if(!isValidMultiPassType_large(val.manyPassType_large)){
                throw std::runtime_error("Invalid manyPassType_large kernel type");
            }

            if(!isValidOverflowType(val.overflowType)){
                throw std::runtime_error("Invalid overflow kernel type");
            }

            kernelTypeConfig = val;
        }

        void setMemoryConfig(const MemoryConfig& val){
            memoryConfig = val;
        }

        std::string_view getReferenceHeader(size_t referenceId) const{
            const int dbChunkIndex = 0;
            const auto& chunkData = fullDB.getChunks()[dbChunkIndex];
            const char* const headerBegin = chunkData.headers() + chunkData.headerOffsets()[referenceId];
            const char* const headerEnd = chunkData.headers() + chunkData.headerOffsets()[referenceId+1];
            return std::string_view(headerBegin, std::distance(headerBegin, headerEnd));
        }

        int getReferenceLength(size_t referenceId) const{
            const int dbChunkIndex = 0;
            const auto& chunkData = fullDB.getChunks()[dbChunkIndex];
            return chunkData.lengths()[referenceId];
        }

        std::string getReferenceSequence(size_t referenceId) const{
            const int dbChunkIndex = 0;
            const auto& chunkData = fullDB.getChunks()[dbChunkIndex];
            const char* const begin = chunkData.chars() + chunkData.offsets()[referenceId];
            const char* const end = chunkData.chars() + chunkData.offsets()[referenceId+1];

            std::string sequence(end - begin, '\0');
            std::transform(
                begin, 
                end,
                sequence.begin(),
                inverse_convert_AA
            );

            return sequence;
        }

        void prefetchFullDBToGpus(){
            RevertDeviceId rdi{};

            const int numGpus = deviceIds.size();
            std::vector<int> copyIds;

            helpers::CpuTimer copyTimer("transfer DB to GPUs");
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]);
                auto& ws = *workingSets[gpu];
                const int chunkId = 0;
                if(ws.canStoreFullDB){
                    const auto& plan = batchPlans_fulldb_perChunk[chunkId][gpu][0];

                    const int currentBuffer = 0;

                    cudaStream_t H2DcopyStream = ws.copyStreams[currentBuffer];

                    executeCopyPlanH2DDirect(
                        plan,
                        ws.d_fulldb_chardata.data(),
                        ws.d_fulldb_lengthdata.data(),
                        ws.d_fulldb_offsetdata.data(),
                        subPartitionsForGpus_perDBchunk[chunkId][gpu],
                        H2DcopyStream
                    );
                    
                    ws.fullDBisUploaded = true;
                    copyIds.push_back(gpu);
                }
            }
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]);
                cudaDeviceSynchronize(); CUERR;
            }
            copyTimer.stop();
            if(copyIds.size() > 0){
                if(verbose){
                    std::cout << "Transferred DB data in advance to GPU(s) ";
                    for(auto x : copyIds){
                        std::cout << x << " ";
                    }
                    std::cout << "\n";
                    copyTimer.print();
                }
            }
        }

        ScanResult scan(const char* query, int queryLength){
            if(!dbIsReady){
                throw std::runtime_error("DB not set correctly");
            }
            RevertDeviceId rdi{};

            const int masterDeviceId = deviceIds[0];
            cudaSetDevice(masterDeviceId);

            scanTimer->reset();
            scanTimer->start();

            setQuery(query, queryLength);

            scanDatabaseForQuery();

            totalProcessedQueryLengths += queryLength;
            totalNumOverflows += resultNumOverflows[0];

            const auto& sequenceLengthStatistics = getSequenceLengthStatistics();

            ScanResult result;
            //std::cout << resultNumOverflows[0] << " overflows\n";
            //std::vector<float> final_alignment_scores_float(alignment_scores_float.begin(), alignment_scores_float.begin() + results_per_query);
            //std::vector<size_t> final_sorted_indices(sorted_indices.begin(), sorted_indices.begin() + results_per_query);
            //std::vector<int> final_resultDbChunkIndices(resultDbChunkIndices.begin(), resultDbChunkIndices.begin() + results_per_query);
            
            //sort the chunk results per query to find overall top results
            // for (int query_num=0; query_num < numQueries; query_num++) {
            //     float* scores =  &alignment_scores_float[query_num * numDBChunks * results_per_query];
            //     size_t* indices =  &sorted_indices[query_num * numDBChunks * results_per_query];
            //     int* chunkIds =  &resultDbChunkIndices[query_num * numDBChunks * results_per_query];

            //     std::vector<int> permutation(results_per_query * numDBChunks);
            //     std::iota(permutation.begin(), permutation.end(), 0);
            //     std::sort(permutation.begin(), permutation.end(),
            //         [&](const auto& l, const auto& r){
            //             return scores[l] > scores[r];
            //         }
            //     );

            //     for(int i = 0; i < results_per_query; i++){
            //         final_alignment_scores_float[query_num * results_per_query + i] = scores[permutation[i]];
            //         final_sorted_indices[query_num * results_per_query + i] = indices[permutation[i]];
            //         final_resultDbChunkIndices[query_num * results_per_query + i] = chunkIds[permutation[i]];
            //     }        
            // }

            result.stats = makeBenchmarkStats(
                scanTimer->elapsed() / 1000, 
                sequenceLengthStatistics.sumOfLengths * queryLength, 
                resultNumOverflows[0]
            );
            result.scores.insert(result.scores.end(), alignment_scores_float.begin(), alignment_scores_float.begin() + results_per_query);
            result.referenceIds.insert(result.referenceIds.end(), sorted_indices.begin(), sorted_indices.begin() + results_per_query);

            // for(int i = 0; i < results_per_query; i++){
            //     const int sortedIndex = final_sorted_indices[i];
            //     const int dbChunkIndex = 0;
            //     const auto& chunkData = fullDB.getChunks()[dbChunkIndex];

            //     const char* headerBegin = chunkData.headers() + chunkData.headerOffsets()[sortedIndex];
            //     const char* headerEnd = chunkData.headers() + chunkData.headerOffsets()[sortedIndex+1];
            //     result.headers.emplace_back(headerBegin, std::distance(headerBegin, headerEnd));
            // }

            return result;
        }


        void printDBInfo() const{
            const size_t numSequences = fullDB.getChunks()[0].numSequences();
            std::cout << numSequences << " sequences, " << fullDB.getChunks()[0].numChars() << " characters\n";

            SequenceLengthStatistics stats = getSequenceLengthStatistics();

            std::cout << "Min length " << stats.min_length << ", max length " << stats.max_length 
                << ", avg length " << stats.sumOfLengths / numSequences << "\n";
        }

        void printDBLengthPartitions() const{
            auto lengthBoundaries = getLengthPartitionBoundaries();
            const int numLengthPartitions = getLengthPartitionBoundaries().size();

            for(int i = 0; i < numLengthPartitions; i++){
                std::cout << "<= " << lengthBoundaries[i] << ": " << fullDB_numSequencesPerLengthPartition[i] << "\n";
            }
        }

        void totalTimerStart(){
            totalProcessedQueryLengths = 0;
            totalNumOverflows = 0;
            totalTimer->start();
        }

        BenchmarkStats totalTimerStop(){
            totalTimer->stop();

            const auto& sequenceLengthStatistics = getSequenceLengthStatistics();
            BenchmarkStats stats = makeBenchmarkStats(
                totalTimer->elapsed() / 1000,
                totalProcessedQueryLengths * sequenceLengthStatistics.sumOfLengths,
                totalNumOverflows
            );

            return stats;
        }

        bool isValidSinglePassType(KernelType type) const{
            return (type == KernelType::Half2 || type == KernelType::DPXs16);
        }

        bool isValidMultiPassType_small(KernelType type) const{
            return (type == KernelType::Half2 || type == KernelType::DPXs16);
        }

        bool isValidMultiPassType_large(KernelType type) const{
            return (type == KernelType::Float);
        }

        bool isValidOverflowType(KernelType type) const{
            return (type == KernelType::Float);
        }

    private:
        void initializeGpus(){
            const int numGpus = deviceIds.size();

            for(int i = 0; i < numGpus; i++){
                cudaSetDevice(deviceIds[i]); CUERR
                helpers::init_cuda_context(); CUERR
                cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);CUERR
                cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte); CUERR
        
                cudaMemPool_t mempool;
                cudaDeviceGetDefaultMemPool(&mempool, deviceIds[i]); CUERR
                uint64_t threshold = UINT64_MAX;
                cudaMemPoolSetAttribute(mempool, cudaMemPoolAttrReleaseThreshold, &threshold);CUERR
            
                gpuStreams.emplace_back();
                gpuEvents.emplace_back(cudaEventDisableTiming);
            }
        }

        void makeReady(){
            const int numDBChunks = fullDB.getInfo().numChunks;
            assert(numDBChunks == 1);

            dbSequenceLengthStatistics = nullptr;

            computeTotalNumSequencePerLengthPartition();
            partitionDBAmongstGpus();

            allocateGpuWorkingSets();

            createDBBatchesForGpus();
            
            dbIsReady = true;
            updateNumResultsPerQuery();
        }

        void computeTotalNumSequencePerLengthPartition(){
            auto lengthBoundaries = getLengthPartitionBoundaries();
            const int numLengthPartitions = getLengthPartitionBoundaries().size();

            fullDB_numSequencesPerLengthPartition.resize(numLengthPartitions);

            const auto& dbData = fullDB.getChunks()[0];
            auto partitionBegin = dbData.lengths();
            for(int i = 0; i < numLengthPartitions; i++){
                //length k is in partition i if boundaries[i-1] < k <= boundaries[i]
                int searchFor = lengthBoundaries[i];
                if(searchFor < std::numeric_limits<int>::max()){
                    searchFor += 1;
                }
                auto partitionEnd = std::lower_bound(
                    partitionBegin, 
                    dbData.lengths() + dbData.numSequences(), 
                    searchFor
                );
                fullDB_numSequencesPerLengthPartition[i] = std::distance(partitionBegin, partitionEnd);
                partitionBegin = partitionEnd;
            }
        }

        void partitionDBAmongstGpus(){
            const int numDBChunks = fullDB.getInfo().numChunks;
            const int numGpus = deviceIds.size();
            const int numLengthPartitions = getLengthPartitionBoundaries().size();

            numSequencesPerLengthPartitionPrefixSum_perDBchunk.clear();
            dbPartitionsByLengthPartitioning_perDBchunk.clear();
            subPartitionsForGpus_perDBchunk.clear();
            lengthPartitionIdsForGpus_perDBchunk.clear();
            numSequencesPerGpu_perDBchunk.clear();
            numSequencesPerGpuPrefixSum_perDBchunk.clear();

            numSequencesPerLengthPartitionPrefixSum_perDBchunk.resize(numDBChunks);
            dbPartitionsByLengthPartitioning_perDBchunk.resize(numDBChunks);
            subPartitionsForGpus_perDBchunk.resize(numDBChunks);
            lengthPartitionIdsForGpus_perDBchunk.resize(numDBChunks);
            numSequencesPerGpu_perDBchunk.resize(numDBChunks);
            numSequencesPerGpuPrefixSum_perDBchunk.resize(numDBChunks);

            for(int chunkId = 0; chunkId < numDBChunks; chunkId++){
                const auto& dbChunk = fullDB.getChunks()[chunkId];
        
                auto& numSequencesPerLengthPartitionPrefixSum = numSequencesPerLengthPartitionPrefixSum_perDBchunk[chunkId];
                auto& dbPartitionsByLengthPartitioning = dbPartitionsByLengthPartitioning_perDBchunk[chunkId];
                auto& subPartitionsForGpus = subPartitionsForGpus_perDBchunk[chunkId];
                auto& lengthPartitionIdsForGpus = lengthPartitionIdsForGpus_perDBchunk[chunkId];
                auto& numSequencesPerGpu = numSequencesPerGpu_perDBchunk[chunkId];
                auto& numSequencesPerGpuPrefixSum = numSequencesPerGpuPrefixSum_perDBchunk[chunkId];
        
                subPartitionsForGpus.resize(numGpus);
                lengthPartitionIdsForGpus.resize(numGpus);
                numSequencesPerGpu.resize(numGpus, 0);
                numSequencesPerGpuPrefixSum.resize(numGpus, 0);
        
                numSequencesPerLengthPartitionPrefixSum.resize(numLengthPartitions, 0);
                for(int i = 0; i < numLengthPartitions-1; i++){
                    numSequencesPerLengthPartitionPrefixSum[i+1] = numSequencesPerLengthPartitionPrefixSum[i] + fullDB_numSequencesPerLengthPartition[i];
                }
        
                for(int i = 0; i < numLengthPartitions; i++){
                    size_t begin = numSequencesPerLengthPartitionPrefixSum[i];
                    size_t end = begin + fullDB_numSequencesPerLengthPartition[i];
                    dbPartitionsByLengthPartitioning.emplace_back(dbChunk, begin, end);        
                }
        
                std::vector<std::vector<int>> numSubPartitionsPerLengthPerGpu(numGpus, std::vector<int>(numLengthPartitions, 0));
            
                for(int lengthPartitionId = 0; lengthPartitionId < numLengthPartitions; lengthPartitionId++){
                    const auto& lengthPartition = dbPartitionsByLengthPartitioning[lengthPartitionId];        
                    const auto partitionedByGpu = partitionDBdata_by_numberOfChars(lengthPartition, lengthPartition.numChars() / numGpus);
            
                    assert(int(partitionedByGpu.size()) <= numGpus);
                    for(int gpu = 0; gpu < numGpus; gpu++){
                        if(gpu < int(partitionedByGpu.size())){
                            subPartitionsForGpus[gpu].push_back(partitionedByGpu[gpu]);
                            lengthPartitionIdsForGpus[gpu].push_back(lengthPartitionId);
                        }else{
                            //add empty partition
                            subPartitionsForGpus[gpu].push_back(DBdataView(dbChunk, 0, 0));
                            lengthPartitionIdsForGpus[gpu].push_back(0);
                        }
                    }
                }
            
                for(int i = 0; i < numGpus; i++){
                    for(const auto& p : subPartitionsForGpus[i]){
                        numSequencesPerGpu[i] += p.numSequences();
                    }
                }
                for(int i = 0; i < numGpus-1; i++){
                    numSequencesPerGpuPrefixSum[i+1] = numSequencesPerGpuPrefixSum[i] + numSequencesPerGpu[i];
                }
            }
        
            numSequencesPerGpu_total.resize(numGpus);
            numSequencesPerGpuPrefixSum_total.resize(numGpus);
            numSequencesPerGpuPrefixSum_total[0] = 0;

        
            for(int i = 0; i < numGpus; i++){
                size_t num = 0;
                for(int chunkId = 0; chunkId < numDBChunks; chunkId++){
                    num += numSequencesPerGpu_perDBchunk[chunkId][i];
                }
                numSequencesPerGpu_total[i] = num;
                if(i < numGpus - 1){
                    numSequencesPerGpuPrefixSum_total[i+1] = numSequencesPerGpuPrefixSum_total[i] + num;
                }
            }
        }

        void allocateGpuWorkingSets(){
            const int numGpus = deviceIds.size();

            workingSets.clear();
            workingSets.resize(numGpus);

            if(verbose){
                std::cout << "Allocate Memory: \n";
            }
            //nvtx::push_range("ALLOC_MEM", 0);
            helpers::CpuTimer allocTimer("ALLOC_MEM");

            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]);

                size_t freeMem, totalMem;
                cudaMemGetInfo(&freeMem, &totalMem);
                constexpr size_t safety = 256*1024*1024;
                size_t memlimit = std::min(freeMem, memoryConfig.maxGpuMem);
                if(memlimit > safety){
                    memlimit -= safety;
                }

                if(verbose){
                    std::cout << "gpu " << gpu << " may use " << memlimit << " bytes. ";
                }

                workingSets[gpu] = std::make_unique<GpuWorkingSet>(
                    memlimit,
                    memoryConfig.maxBatchBytes,
                    memoryConfig.maxBatchSequences,
                    memoryConfig.maxTempBytes,
                    subPartitionsForGpus_perDBchunk[0][gpu]
                );

                if(verbose){
                    std::cout << "Using " << workingSets[gpu]->numTempBytes << " temp bytes. ";

                    if(workingSets[gpu]->canStoreFullDB){
                        std::cout << "It can store its DB in memory\n";
                    }else{
                        std::cout << "It will process its DB in batches\n";
                    }
                }

                //spin up the host callback thread
                auto noop = [](void*){};
                cudaLaunchHostFunc(
                    gpuStreams[gpu], 
                    noop, 
                    nullptr
                ); CUERR
            }    

            if(verbose){
                allocTimer.print();
            }
        }

        void createDBBatchesForGpus(){
            const int numDBChunks = fullDB.getInfo().numChunks;
            const int numGpus = deviceIds.size();

            batchPlans_perChunk.clear();
            batchPlans_fulldb_perChunk.clear();
            batchPlans_perChunk.resize(numDBChunks);
            batchPlans_fulldb_perChunk.resize(numDBChunks);

            for(int chunkId = 0; chunkId < numDBChunks; chunkId++){
                batchPlans_perChunk[chunkId].resize(numGpus);
                batchPlans_fulldb_perChunk[chunkId].resize(numGpus);
        
                for(int gpu = 0; gpu < numGpus; gpu++){
                    const auto& ws = *workingSets[gpu];
                    
                    batchPlans_perChunk[chunkId][gpu] = computeDbCopyPlan(
                        subPartitionsForGpus_perDBchunk[chunkId][gpu],
                        lengthPartitionIdsForGpus_perDBchunk[chunkId][gpu],
                        sizeof(char) * ws.h_chardata_vec[0].size(),
                        ws.h_lengthdata_vec[0].size()
                    );
                    if(verbose){
                        std::cout << "Batch plan chunk " << chunkId << ", gpu " << gpu << ": " << batchPlans_perChunk[chunkId][gpu].size() << " batches\n";
                    }
        
                    if(ws.canStoreFullDB){
                        batchPlans_fulldb_perChunk[chunkId][gpu] = computeDbCopyPlan(
                            subPartitionsForGpus_perDBchunk[chunkId][gpu],
                            lengthPartitionIdsForGpus_perDBchunk[chunkId][gpu],
                            sizeof(char) * ws.d_fulldb_chardata.size(),
                            ws.d_fulldb_lengthdata.size()
                        );
                        assert(batchPlans_fulldb_perChunk[chunkId][gpu].size() == 1);
                        //std::cout << "Batch plan fulldb chunk " << chunkId << ", gpu " << gpu << ": " << batchPlans_fulldb_perChunk[chunkId][gpu].size() << " batches\n";
                    }else{
                        batchPlans_fulldb_perChunk[chunkId][gpu] = batchPlans_perChunk[chunkId][gpu]; //won't be used in this case, simply set it to batched plan
                    }
        
                    // for(int i = 0; i < std::min(5, int(batchPlans_perChunk[chunkId][gpu].size())); i++){
                    //     std::cout << batchPlans_perChunk[chunkId][gpu][i] << "\n";
                    // }
                }
            }
        }

        void printDBDataView(const DBdataView& view) const{
            std::cout << "Sequences: " << view.numSequences() << "\n";
            std::cout << "Chars: " << view.offsets()[0] << " - " << view.offsets()[view.numSequences()] << " (" << (view.offsets()[view.numSequences()] - view.offsets()[0]) << ")"
                << " " << view.numChars() << "\n";
        }

        void printDBDataViews(const std::vector<DBdataView>& views) const {
            size_t numViews = views.size();
            for(size_t p = 0; p < numViews; p++){
                const DBdataView& view = views[p];
        
                std::cout << "View " << p << "\n";
                printDBDataView(view);
            }
        }

        SequenceLengthStatistics getSequenceLengthStatistics() const{
            if(dbSequenceLengthStatistics == nullptr){
                dbSequenceLengthStatistics = std::make_unique<SequenceLengthStatistics>();

                const int numDBChunks = fullDB.getInfo().numChunks;
                for(int i = 0; i < numDBChunks; i++){
                    const auto& chunkData = fullDB.getChunks()[i];
                    size_t numSeq = chunkData.numSequences();

                    for (size_t i=0; i < numSeq; i++) {
                        if (chunkData.lengths()[i] > dbSequenceLengthStatistics->max_length) dbSequenceLengthStatistics->max_length = chunkData.lengths()[i];
                        if (chunkData.lengths()[i] < dbSequenceLengthStatistics->min_length) dbSequenceLengthStatistics->min_length = chunkData.lengths()[i];
                        dbSequenceLengthStatistics->sumOfLengths += chunkData.lengths()[i];
                    }
                }
            }
            return *dbSequenceLengthStatistics;
        }

        std::vector<DeviceBatchCopyToPinnedPlan> computeDbCopyPlan(
            const std::vector<DBdataView>& dbPartitions,
            const std::vector<int>& lengthPartitionIds,
            size_t MAX_CHARDATA_BYTES,
            size_t MAX_SEQ
        ) const {
            std::vector<DeviceBatchCopyToPinnedPlan> result;
        
            size_t currentCopyPartition = 0;
            size_t currentCopySeqInPartition = 0;
        
            //size_t processedSequences = 0;
            while(currentCopyPartition < dbPartitions.size()){
                
                size_t usedBytes = 0;
                size_t usedSeq = 0;
        
                DeviceBatchCopyToPinnedPlan plan;
        
                while(currentCopyPartition < dbPartitions.size()){
                    if(dbPartitions[currentCopyPartition].numSequences() == 0){
                        currentCopyPartition++;
                        continue;
                    }
        
                    //figure out how many sequences to copy to pinned
                    size_t remainingBytes = MAX_CHARDATA_BYTES - usedBytes;
                    
                    auto dboffsetsBegin = dbPartitions[currentCopyPartition].offsets() + currentCopySeqInPartition;
                    auto dboffsetsEnd = dbPartitions[currentCopyPartition].offsets() + dbPartitions[currentCopyPartition].numSequences() + 1;
                    
                    auto searchFor = dbPartitions[currentCopyPartition].offsets()[currentCopySeqInPartition] + remainingBytes + 1; // +1 because remainingBytes is inclusive
                    auto it = std::lower_bound(
                        dboffsetsBegin,
                        dboffsetsEnd,
                        searchFor
                    );
        
                    size_t numToCopyByBytes = 0;
                    if(it != dboffsetsBegin){
                        numToCopyByBytes = std::distance(dboffsetsBegin, it) - 1;
                    }
                    if(numToCopyByBytes == 0 && currentCopySeqInPartition == 0){
                        std::cout << "Warning. copy buffer size too small. skipped a db portion\n";
                        break;
                    }
                    
                    size_t remainingSeq = MAX_SEQ - usedSeq;            
                    size_t numToCopyBySeq = std::min(dbPartitions[currentCopyPartition].numSequences() - currentCopySeqInPartition, remainingSeq);
                    size_t numToCopy = std::min(numToCopyByBytes,numToCopyBySeq);
        
                    if(numToCopy > 0){
                        DeviceBatchCopyToPinnedPlan::CopyRange copyRange;
                        copyRange.lengthPartitionId = lengthPartitionIds[currentCopyPartition];
                        copyRange.currentCopyPartition = currentCopyPartition;
                        copyRange.currentCopySeqInPartition = currentCopySeqInPartition;
                        copyRange.numToCopy = numToCopy;
                        plan.copyRanges.push_back(copyRange);
        
                        if(usedSeq == 0){
                            plan.h_partitionIds.push_back(lengthPartitionIds[currentCopyPartition]);
                            plan.h_numPerPartition.push_back(numToCopy);
                        }else{
                            //if is same length partition as previous copy 
                            if(plan.h_partitionIds.back() == lengthPartitionIds[currentCopyPartition]){
                                plan.h_numPerPartition.back() += numToCopy;
                            }else{
                                //new length partition
                                plan.h_partitionIds.push_back(lengthPartitionIds[currentCopyPartition]);
                                plan.h_numPerPartition.push_back(numToCopy);
                            }
                        }
                        usedBytes += (dbPartitions[currentCopyPartition].offsets()[currentCopySeqInPartition+numToCopy] 
                            - dbPartitions[currentCopyPartition].offsets()[currentCopySeqInPartition]);
                        usedSeq += numToCopy;
        
                        currentCopySeqInPartition += numToCopy;
                        if(currentCopySeqInPartition == dbPartitions[currentCopyPartition].numSequences()){
                            currentCopySeqInPartition = 0;
                            currentCopyPartition++;
                        }
                    }else{
                        break;
                    }
                }
        
                plan.usedBytes = usedBytes;
                plan.usedSeq = usedSeq;    
                
                if(usedSeq == 0 && currentCopyPartition < dbPartitions.size() && dbPartitions[currentCopyPartition].numSequences() > 0){
                    std::cout << "Warning. copy buffer size too small. skipped a db portion. stop\n";
                    break;
                }
        
                if(plan.usedSeq > 0){
                    result.push_back(plan);
                }
            }
        
            return result;
        }


        void setQuery(const char* query, int queryLength){
            currentQueryLength = queryLength;

            const int numGpus = deviceIds.size();
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                auto& ws = *workingSets[gpu];
                const int roundedLength = SDIV(queryLength, 128) * 128;
                ws.d_query.resize(roundedLength);
                cudaMemsetAsync(ws.d_query.data(), 20, roundedLength, gpuStreams[gpu]);
                cudaMemcpyAsync(ws.d_query.data(), query, queryLength, cudaMemcpyDefault, gpuStreams[gpu]); CUERR
                NW_convert_protein_single<<<SDIV(queryLength, 128), 128, 0, gpuStreams[gpu]>>>(ws.d_query.data(), queryLength); CUERR
            }
        }

        void scanDatabaseForQuery(){
            const int numGpus = deviceIds.size();
            const int masterDeviceId = deviceIds[0];
            const auto& masterStream1 = gpuStreams[0];
            auto& masterevent1 = gpuEvents[0];

            cudaSetDevice(masterDeviceId);
            // scanTimer->reset();
            // scanTimer->start();

            thrust::fill(
                thrust::cuda::par_nosync.on(masterStream1),
                devAllAlignmentScoresFloat.begin(),
                devAllAlignmentScoresFloat.end(),
                0
            );

            //std::cout << "Starting NW_local_affine_half2 for Query " << query_num << "\n";

            const int numDBChunks = 1;

            for(int chunkId = 0; chunkId < numDBChunks; chunkId++){
                cudaSetDevice(masterDeviceId);           

                cudaEventRecord(masterevent1, masterStream1); CUERR;

                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    cudaStreamWaitEvent(gpuStreams[gpu], masterevent1, 0); CUERR;
                }

                processQueryOnGpus();

                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];

                    //we could sort the maxReduceArray here as well and only send the best results_per_query entries

                    if(numGpus > 1){
                        //transform per gpu local sequence indices into global sequence indices
                        transformLocalSequenceIndicesToGlobalIndices<<<SDIV(maxReduceArraySize, 128), 128, 0, gpuStreams[gpu]>>>(
                            gpu,
                            maxReduceArraySize,
                            ws.deviceGpuPartitionOffsets.getDeviceView(),
                            ws.d_maxReduceArrayIndices.data()
                        ); CUERR;
                    }

                    cudaMemcpyAsync(
                        devAllAlignmentScoresFloat.data() + maxReduceArraySize*gpu,
                        ws.d_maxReduceArrayScores.data(),
                        sizeof(float) * maxReduceArraySize,
                        cudaMemcpyDeviceToDevice,
                        gpuStreams[gpu]
                    ); CUERR;
                    cudaMemcpyAsync(
                        dev_sorted_indices.data() + maxReduceArraySize*gpu,
                        ws.d_maxReduceArrayIndices.data(),
                        sizeof(size_t) * maxReduceArraySize,
                        cudaMemcpyDeviceToDevice,
                        gpuStreams[gpu]
                    ); CUERR;                
                    cudaMemcpyAsync(
                        d_resultNumOverflows.data() + gpu,
                        ws.d_total_overflow_number.data(),
                        sizeof(int),
                        cudaMemcpyDeviceToDevice,
                        gpuStreams[gpu]
                    ); CUERR;                

                    cudaEventRecord(ws.forkStreamEvent, gpuStreams[gpu]); CUERR;

                    cudaSetDevice(masterDeviceId);
                    cudaStreamWaitEvent(masterStream1, ws.forkStreamEvent, 0); CUERR;
                }

                cudaSetDevice(masterDeviceId);

                // thrust::sequence(
                //     thrust::cuda::par_nosync.on(masterStream1),
                //     dev_sorted_indices.begin(), 
                //     dev_sorted_indices.end(),
                //     0
                // );
                // thrust::sort_by_key(
                //     thrust::cuda::par_nosync(thrust_async_allocator<char>(masterStream1)).on(masterStream1),
                //     // thrust::cuda::par_nosync(thrust_preallocated_single_allocator<char>((void*)workingSets[0]->d_tempStorageHE, 
                //     //     workingSets[0]->numTempBytes)).on(masterStream1),
                //     devAllAlignmentScoresFloat.begin(),
                //     devAllAlignmentScoresFloat.end(),
                //     dev_sorted_indices.begin(),
                //     thrust::greater<float>()
                // );

                thrust::sort_by_key(
                    thrust::cuda::par_nosync(thrust_async_allocator<char>(masterStream1)).on(masterStream1),
                    // thrust::cuda::par_nosync(thrust_preallocated_single_allocator<char>((void*)workingSets[0]->d_tempStorageHE, 
                    //     workingSets[0]->numTempBytes)).on(masterStream1),
                    devAllAlignmentScoresFloat.begin(),
                    devAllAlignmentScoresFloat.begin() + maxReduceArraySize * numGpus,
                    dev_sorted_indices.begin(),
                    thrust::greater<float>()
                );

                if(numGpus > 1){
                    //sum the overflows per gpu
                    sumNumOverflowsKernel<<<1,1,0,masterStream1>>>(d_resultNumOverflows.data(), d_resultNumOverflows.data(), numGpus); CUERR;                
                }

                std::fill(
                    &resultDbChunkIndices[numDBChunks*results_per_query + chunkId * results_per_query],
                    &resultDbChunkIndices[numDBChunks*results_per_query + chunkId * results_per_query] + results_per_query,
                    chunkId
                );

                cudaMemcpyAsync(
                    &(alignment_scores_float[chunkId * results_per_query]), 
                    devAllAlignmentScoresFloat.data(), 
                    sizeof(float) * results_per_query, 
                    cudaMemcpyDeviceToHost, 
                    masterStream1
                );  CUERR
                cudaMemcpyAsync(
                    &(sorted_indices[chunkId * results_per_query]), 
                    dev_sorted_indices.data(), 
                    sizeof(size_t) * results_per_query, 
                    cudaMemcpyDeviceToHost, 
                    masterStream1
                );  CUERR
                cudaMemcpyAsync(
                    &(resultNumOverflows[chunkId]), 
                    d_resultNumOverflows.data(), 
                    sizeof(int), 
                    cudaMemcpyDeviceToHost, 
                    masterStream1
                );  CUERR

                

            }

            scanTimer->stop();
            cudaStreamSynchronize(masterStream1); CUERR;
        }

        void processQueryOnGpus(){

            const int chunkId = 0;
            const std::vector<std::vector<DBdataView>>& dbPartitionsPerGpu = subPartitionsForGpus_perDBchunk[chunkId];
            const std::vector<std::vector<DeviceBatchCopyToPinnedPlan>>& batchPlansPerGpu_batched = batchPlans_perChunk[chunkId];
            const std::vector<std::vector<DeviceBatchCopyToPinnedPlan>>& batchPlansPerGpu_full = batchPlans_fulldb_perChunk[chunkId];
            const std::vector<size_t>& numberOfSequencesPerGpu = numSequencesPerGpu_perDBchunk[chunkId];

            
            constexpr auto boundaries = getLengthPartitionBoundaries();
            constexpr int numLengthPartitions = boundaries.size();
            const int numGpus = deviceIds.size();
            const bool useExtraThreadForBatchTransfer = numGpus > 1;
        
            size_t totalNumberOfSequencesToProcess = std::reduce(numberOfSequencesPerGpu.begin(), numberOfSequencesPerGpu.end());
            
            size_t totalNumberOfProcessedSequences = 0;
            // std::vector<size_t> processedSequencesPerGpu(numGpus, 0);
        
            // std::vector<size_t> processedBatchesPerGpu(numGpus, 0);
        
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                auto& ws = *workingSets[gpu];
        
                cudaMemsetAsync(ws.d_total_overflow_number.data(), 0, sizeof(int), gpuStreams[gpu]);
                
                ws.resetMaxReduceArray(gpuStreams[gpu]);
        
                //create dependency on mainStream
                cudaEventRecord(ws.forkStreamEvent, gpuStreams[gpu]); CUERR;
                cudaStreamWaitEvent(ws.workStreamForTempUsage, ws.forkStreamEvent, 0); CUERR;
                for(auto& stream : ws.workStreamsWithoutTemp){
                    cudaStreamWaitEvent(stream, ws.forkStreamEvent, 0); CUERR;
                }
                cudaStreamWaitEvent(ws.hostFuncStream, ws.forkStreamEvent, 0); CUERR;
            }       
        
            //variables per gpu to keep between loops
            struct Variables{
                int currentBuffer = 0;
                int previousBuffer = 0;
                cudaStream_t H2DcopyStream = cudaStreamLegacy;
                char* h_inputChars = nullptr;
                size_t* h_inputLengths = nullptr;
                size_t* h_inputOffsets = nullptr;
                char* d_inputChars = nullptr;
                size_t* d_inputLengths = nullptr;
                size_t* d_inputOffsets = nullptr;
                int* d_overflow_number = nullptr;
                size_t* d_overflow_positions = nullptr;
                size_t pointerSequencesOffset = 0;
                size_t pointerBytesOffset = 0;
                const std::vector<DeviceBatchCopyToPinnedPlan>* batchPlansPtr;
                size_t processedSequences = 0;
                size_t processedBatches = 0;
            };
        
            std::vector<Variables> variables_vec(numGpus);
            //init variables
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                const auto& ws = *workingSets[gpu];
                auto& variables = variables_vec[gpu];
                variables.processedSequences = 0;
                variables.processedBatches = 0;
                if(ws.canStoreFullDB){
                    if(ws.fullDBisUploaded){
                        assert(batchPlansPerGpu_full[gpu].size() == 1);
                        variables.batchPlansPtr = &batchPlansPerGpu_full[gpu];
                    }else{
                        variables.batchPlansPtr = &batchPlansPerGpu_batched[gpu];
                    }
                }else{
                    variables.batchPlansPtr = &batchPlansPerGpu_batched[gpu];
                }
            }
            
            while(totalNumberOfProcessedSequences < totalNumberOfSequencesToProcess){
                //set up gpu variables for current iteration
                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];
                    auto& variables = variables_vec[gpu];
                    if(variables.processedBatches < variables.batchPlansPtr->size()){
                        if(ws.canStoreFullDB){
                            if(ws.fullDBisUploaded){
                                //single buffer, no transfer, 1 batch
                                variables.currentBuffer = 0;
                                variables.previousBuffer = 0;
                                variables.H2DcopyStream = ws.copyStreams[0];
                                variables.h_inputChars = nullptr;
                                variables.h_inputLengths = nullptr;
                                variables.h_inputOffsets = nullptr;
                                variables.d_inputChars = ws.d_fulldb_chardata.data();
                                variables.d_inputLengths = ws.d_fulldb_lengthdata.data();
                                variables.d_inputOffsets = ws.d_fulldb_offsetdata.data();
                                variables.d_overflow_number = ws.d_overflow_number.data() + 0;
                                variables.d_overflow_positions = ws.d_overflow_positions_vec[0].data();
                            }else{
                                //first query, double buffer batched transfer into full db buffer
                                variables.currentBuffer = ws.copyBufferIndex;
                                if(variables.currentBuffer == 0){
                                    variables.previousBuffer = ws.numCopyBuffers - 1;
                                }else{
                                    variables.previousBuffer = (variables.currentBuffer - 1);
                                } 
                                variables.H2DcopyStream = ws.copyStreams[variables.currentBuffer];
                                variables.h_inputChars = ws.h_chardata_vec[variables.currentBuffer].data();
                                variables.h_inputLengths = ws.h_lengthdata_vec[variables.currentBuffer].data();
                                variables.h_inputOffsets = ws.h_offsetdata_vec[variables.currentBuffer].data();
                                variables.d_inputChars = ws.d_fulldb_chardata.data() + variables.pointerBytesOffset;
                                variables.d_inputLengths = ws.d_fulldb_lengthdata.data() + variables.pointerSequencesOffset;
                                variables.d_inputOffsets = ws.d_fulldb_offsetdata.data() + variables.pointerSequencesOffset;
                                variables.d_overflow_number = ws.d_overflow_number.data() + variables.currentBuffer;
                                variables.d_overflow_positions = ws.d_overflow_positions_vec[variables.currentBuffer].data();
                            }
                        }else{
                            //full db not possible, any query, double buffer batched transfer
                            variables.currentBuffer = ws.copyBufferIndex;
                            if(variables.currentBuffer == 0){
                                variables.previousBuffer = ws.numCopyBuffers - 1;
                            }else{
                                variables.previousBuffer = (variables.currentBuffer - 1);
                            } 
                            variables.H2DcopyStream = ws.copyStreams[variables.currentBuffer];
                            variables.h_inputChars = ws.h_chardata_vec[variables.currentBuffer].data();
                            variables.h_inputLengths = ws.h_lengthdata_vec[variables.currentBuffer].data();
                            variables.h_inputOffsets = ws.h_offsetdata_vec[variables.currentBuffer].data();
                            variables.d_inputChars = ws.d_chardata_vec[variables.currentBuffer].data();
                            variables.d_inputLengths = ws.d_lengthdata_vec[variables.currentBuffer].data();
                            variables.d_inputOffsets = ws.d_offsetdata_vec[variables.currentBuffer].data();
                            variables.d_overflow_number = ws.d_overflow_number.data() + variables.currentBuffer;
                            variables.d_overflow_positions = ws.d_overflow_positions_vec[variables.currentBuffer].data();
                        }
                    }
                }
                //upload batch
                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];
                    auto& variables = variables_vec[gpu];
                    if(variables.processedBatches < variables.batchPlansPtr->size()){
                        const auto& plan = (*variables.batchPlansPtr)[variables.processedBatches];
        
                        if((ws.canStoreFullDB && !ws.fullDBisUploaded) || !ws.canStoreFullDB){
                            //transfer data
                            //can only overwrite device buffer if it is no longer in use on workstream
                            cudaStreamWaitEvent(variables.H2DcopyStream, ws.deviceBufferEvents[variables.currentBuffer], 0); CUERR;
        
                            if(useExtraThreadForBatchTransfer){
                                cudaStreamWaitEvent(ws.hostFuncStream, ws.pinnedBufferEvents[variables.currentBuffer]); CUERR;
                                executePinnedCopyPlanWithHostCallback(
                                    plan, 
                                    variables.h_inputChars,
                                    variables.h_inputLengths,
                                    variables.h_inputOffsets,
                                    dbPartitionsPerGpu[gpu], 
                                    ws.hostFuncStream
                                );
                                cudaEventRecord(ws.forkStreamEvent, ws.hostFuncStream); CUERR;
                                cudaStreamWaitEvent(variables.H2DcopyStream, ws.forkStreamEvent, 0);
        
                                cudaMemcpyAsync(
                                    variables.d_inputChars,
                                    variables.h_inputChars,
                                    plan.usedBytes,
                                    H2D,
                                    variables.H2DcopyStream
                                ); CUERR;
                                cudaMemcpyAsync(
                                    variables.d_inputLengths,
                                    variables.h_inputLengths,
                                    sizeof(size_t) * plan.usedSeq,
                                    H2D,
                                    variables.H2DcopyStream
                                ); CUERR;
                                cudaMemcpyAsync(
                                    variables.d_inputOffsets,
                                    variables.h_inputOffsets,
                                    sizeof(size_t) * (plan.usedSeq+1),
                                    H2D,
                                    variables.H2DcopyStream
                                ); CUERR;
                            }else{
                                //synchronize to avoid overwriting pinned buffer of target before it has been fully transferred
                                cudaEventSynchronize(ws.pinnedBufferEvents[variables.currentBuffer]); CUERR;
                                //executePinnedCopyPlanSerial(plan, ws, currentBuffer, dbPartitionsPerGpu[gpu]);
        
                                executePinnedCopyPlanSerialAndTransferToGpu(
                                    plan, 
                                    variables.h_inputChars,
                                    variables.h_inputLengths,
                                    variables.h_inputOffsets,
                                    variables.d_inputChars,
                                    variables.d_inputLengths,
                                    variables.d_inputOffsets,
                                    dbPartitionsPerGpu[gpu], 
                                    variables.H2DcopyStream
                                );
                            }
                            
                            cudaEventRecord(ws.pinnedBufferEvents[variables.currentBuffer], variables.H2DcopyStream); CUERR;
                        }
        
                        cudaMemsetAsync(variables.d_overflow_number, 0, sizeof(int), variables.H2DcopyStream);
                        // cudaMemsetAsync(d_overflow_positions, 0, sizeof(size_t) * 10000, H2DcopyStream);
                        
                        //all data is ready for alignments. create dependencies for work streams
                        cudaEventRecord(ws.forkStreamEvent, variables.H2DcopyStream); CUERR;
                        cudaStreamWaitEvent(ws.workStreamForTempUsage, ws.forkStreamEvent, 0); CUERR;
                        for(auto& stream : ws.workStreamsWithoutTemp){
                            cudaStreamWaitEvent(stream, ws.forkStreamEvent, 0); CUERR;
                        }
                        //wait for previous batch to finish
                        cudaStreamWaitEvent(ws.workStreamForTempUsage, ws.deviceBufferEvents[variables.previousBuffer], 0); CUERR;
                        for(auto& stream : ws.workStreamsWithoutTemp){
                            cudaStreamWaitEvent(stream, ws.deviceBufferEvents[variables.previousBuffer], 0); CUERR;
                        }
                //     }
                // }
        
                //launch alignment kernels
        
                // for(int gpu = 0; gpu < numGpus; gpu++){
                //     if(processedBatchesPerGpu[gpu] < batchPlanPerGpu[gpu].size()){
                //         cudaSetDevice(deviceIds[gpu]); CUERR;
                //         auto& ws = *workingSets[gpu];
                //         const auto& plan = batchPlanPerGpu[gpu][processedBatchesPerGpu[gpu]];
                //         int currentBuffer = 0;
                //         //int previousBuffer = 0;
                //         //cudaStream_t H2DcopyStream = ws.copyStreams[currentBuffer];
                //         if(!ws.singleBatchDBisOnGpu){
                //             currentBuffer = ws.copyBufferIndex;
                //             //H2DcopyStream = ws.copyStreams[currentBuffer];
                //         }
        
                        // const char* const inputChars = ws.d_chardata_vec[currentBuffer].data();
                        // const size_t* const inputOffsets = ws.d_offsetdata_vec[currentBuffer].data();
                        // const size_t* const inputLengths = ws.d_lengthdata_vec[currentBuffer].data();
                        // int* const d_overflow_number = ws.d_overflow_number.data() + currentBuffer;
                        // size_t* const d_overflow_positions = ws.d_overflow_positions_vec[currentBuffer].data();
        
                        const char* const inputChars = variables.d_inputChars;
                        const size_t* const inputLengths = variables.d_inputLengths;
                        const size_t* const inputOffsets = variables.d_inputOffsets;
                        int* const d_overflow_number = variables.d_overflow_number;
                        size_t* const d_overflow_positions = variables.d_overflow_positions;
        
                        auto runAlignmentKernels = [&](auto& d_scores, size_t* d_overflow_positions, int* d_overflow_number){
                            const char4* const d_query = reinterpret_cast<char4*>(ws.d_query.data());
        
                            auto nextWorkStreamNoTemp = [&](){
                                ws.workstreamIndex = (ws.workstreamIndex + 1) % ws.numWorkStreamsWithoutTemp;
                                return (cudaStream_t)ws.workStreamsWithoutTemp[ws.workstreamIndex];
                            };
                            std::vector<int> numPerPartitionPrefixSum(plan.h_numPerPartition.size());
                            for(int i = 0; i < int(plan.h_numPerPartition.size())-1; i++){
                                numPerPartitionPrefixSum[i+1] = numPerPartitionPrefixSum[i] + plan.h_numPerPartition[i];
                            }
                            //size_t exclPs = 0;
                            //for(int lp = 0; lp < int(plan.h_partitionIds.size()); lp++){
                            for(int lp = plan.h_partitionIds.size() - 1; lp >= 0; lp--){
                                const int partId = plan.h_partitionIds[lp];
                                const int numSeq = plan.h_numPerPartition[lp];
                                const int start = numPerPartitionPrefixSum[lp];
                                //std::cout << "partId " << partId << " numSeq " << numSeq << "\n";
                                
                                //const size_t* const d_selectedPositions = ws.d_selectedPositions.data() + start;
                                auto d_selectedPositions = thrust::make_counting_iterator<size_t>(start);
        
                                #if 1
                                if(kernelTypeConfig.singlePassType == KernelType::Half2){
                                    
                                    constexpr int sh2bs = 256; // single half2 blocksize 
                                    if (partId == 0){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 2, 24>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 1){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 4, 16>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 2){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 10>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 3){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 12>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 4){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 14>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 5){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 16>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 6){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 18>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 7){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 20>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 8){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 22>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 9){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 24>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 10){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 26>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 11){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 28>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 12){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 30>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 13){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 8, 32>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 14){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 16, 18>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 15){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 16, 20>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 16){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 16, 22>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 17){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 16, 24>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 18){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 16, 26>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 19){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 16, 28>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 20){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 16, 30>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 21){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 16, 32>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 22){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 18>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 23){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 20>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 24){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 22>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 25){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 24>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 26){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 26>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 27){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 28>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 28){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 30>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 29){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 32>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 30){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 34>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 31){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 36>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 32){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 38>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 33){call_NW_local_affine_Protein_single_pass_half2_new<sh2bs, 32, 40>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                
                                }else if(kernelTypeConfig.singlePassType == KernelType::DPXs16){
        
                                    constexpr int sh2bs = 256; // dpx s16 blocksize 
                                    if (partId == 0){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 2, 24>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 1){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 4, 16>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 2){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 10>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 3){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 12>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 4){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 14>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 5){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 16>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 6){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 18>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 7){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 20>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 8){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 22>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 9){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 24>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 10){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 26>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 11){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 28>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 12){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 30>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 0, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 13){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 8, 32>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 14){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 16, 18>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 15){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 16, 20>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 16){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 16, 22>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 17){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 16, 24>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 18){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 16, 26>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 19){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 16, 28>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 20){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 16, 30>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 21){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 16, 32>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 22){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 18>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 23){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 20>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 24){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 22>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 25){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 24>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 26){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 26>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 27){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 28>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 28){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 30>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 29){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 32>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 30){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 34>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 31){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 36>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 32){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 38>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }
                                    if (partId == 33){call_NW_local_affine_single_pass_s16_DPX_new<sh2bs, 32, 40>(blosumType, inputChars, d_scores, inputOffsets , inputLengths, d_selectedPositions, numSeq, d_overflow_positions, d_overflow_number, 1, d_query, currentQueryLength, gop, gex, nextWorkStreamNoTemp()); CUERR }

                                }
                                #endif
        
        
                                if(partId == numLengthPartitions - 2){
                                    if(kernelTypeConfig.manyPassType_small == KernelType::Half2){
                                        constexpr int blocksize = 32 * 8;
                                        constexpr int groupsize = 32;
                                        constexpr int groupsPerBlock = blocksize / groupsize;
                                        constexpr int alignmentsPerGroup = 2;
                                        constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
                                        
                                        const size_t tempBytesPerBlockPerBuffer = sizeof(__half2) * alignmentsPerBlock * currentQueryLength;
        
                                        const size_t maxNumBlocks = ws.numTempBytes / (tempBytesPerBlockPerBuffer * 2);
                                        const size_t maxSubjectsPerIteration = std::min(maxNumBlocks * alignmentsPerBlock, size_t(numSeq));
        
                                        const size_t numBlocksPerIteration = SDIV(maxSubjectsPerIteration, alignmentsPerBlock);
                                        const size_t requiredTempBytes = tempBytesPerBlockPerBuffer * 2 * numBlocksPerIteration;
        
                                        __half2* d_temp = (__half2*)ws.d_tempStorageHE.data();
                                        __half2* d_tempHcol2 = d_temp;
                                        __half2* d_tempEcol2 = (__half2*)(((char*)d_tempHcol2) + requiredTempBytes / 2);
        
                                        const int numIters =  SDIV(numSeq, maxSubjectsPerIteration);
        
                                        for(int iter = 0; iter < numIters; iter++){
                                            const size_t begin = iter * maxSubjectsPerIteration;
                                            const size_t end = iter < numIters-1 ? (iter+1) * maxSubjectsPerIteration : numSeq;
                                            const size_t num = end - begin;                      
        
                                            cudaMemsetAsync(d_temp, 0, requiredTempBytes, ws.workStreamForTempUsage); CUERR;
                                            //std::cout << "iter " << iter << " / " << numIters << " gridsize " << SDIV(num, alignmentsPerBlock) << "\n";
        
                                            //cudaDeviceSynchronize(); CUERR;
                                            
                                            //NW_local_affine_Protein_many_pass_half2<groupsize, 12><<<SDIV(num, alignmentsPerBlock), blocksize, 0, ws.workStreamForTempUsage>>>(
                                            call_NW_local_affine_Protein_many_pass_half2_new<blocksize, groupsize, 22>(
                                                blosumType,
                                                inputChars, 
                                                d_scores, 
                                                d_tempHcol2, 
                                                d_tempEcol2, 
                                                inputOffsets, 
                                                inputLengths, 
                                                d_selectedPositions + begin, 
                                                num, 
                                                d_overflow_positions, 
                                                d_overflow_number, 
                                                1, 
                                                d_query,
                                                currentQueryLength, 
                                                gop, 
                                                gex,
                                                ws.workStreamForTempUsage
                                            ); CUERR
                                        }
                                    }else if(kernelTypeConfig.manyPassType_small == KernelType::DPXs16){
                                        constexpr int blocksize = 32 * 8;
                                        constexpr int groupsize = 32;
                                        constexpr int groupsPerBlock = blocksize / groupsize;
                                        constexpr int alignmentsPerGroup = 2;
                                        constexpr int alignmentsPerBlock = groupsPerBlock * alignmentsPerGroup;
                                        
                                        const size_t tempBytesPerBlockPerBuffer = sizeof(short2) * alignmentsPerBlock * currentQueryLength;
        
                                        const size_t maxNumBlocks = ws.numTempBytes / (tempBytesPerBlockPerBuffer * 2);
                                        const size_t maxSubjectsPerIteration = std::min(maxNumBlocks * alignmentsPerBlock, size_t(numSeq));
        
                                        const size_t numBlocksPerIteration = SDIV(maxSubjectsPerIteration, alignmentsPerBlock);
                                        const size_t requiredTempBytes = tempBytesPerBlockPerBuffer * 2 * numBlocksPerIteration;
        
                                        short2* d_temp = (short2*)ws.d_tempStorageHE.data();
                                        short2* d_tempHcol2 = d_temp;
                                        short2* d_tempEcol2 = (short2*)(((char*)d_tempHcol2) + requiredTempBytes / 2);
        
                                        const int numIters =  SDIV(numSeq, maxSubjectsPerIteration);
        
                                        for(int iter = 0; iter < numIters; iter++){
                                            const size_t begin = iter * maxSubjectsPerIteration;
                                            const size_t end = iter < numIters-1 ? (iter+1) * maxSubjectsPerIteration : numSeq;
                                            const size_t num = end - begin;                      
        
                                            cudaMemsetAsync(d_temp, 0, requiredTempBytes, ws.workStreamForTempUsage); CUERR;
                                            //std::cout << "iter " << iter << " / " << numIters << " gridsize " << SDIV(num, alignmentsPerBlock) << "\n";
        
                                            //cudaDeviceSynchronize(); CUERR;
                                            
                                            call_NW_local_affine_many_pass_s16_DPX_new<blocksize, groupsize, 22>(
                                                blosumType,
                                                inputChars, 
                                                d_scores, 
                                                d_tempHcol2, 
                                                d_tempEcol2, 
                                                inputOffsets, 
                                                inputLengths, 
                                                d_selectedPositions + begin, 
                                                num, 
                                                d_overflow_positions, 
                                                d_overflow_number, 
                                                1, 
                                                d_query,
                                                currentQueryLength, 
                                                gop, 
                                                gex,
                                                ws.workStreamForTempUsage
                                            ); CUERR
                                        }
                                    }else{
                                        assert(false);
                                    }
                                }
        
        
                                if(partId == numLengthPartitions - 1){
                                    if(kernelTypeConfig.manyPassType_large == KernelType::Float){
                                        const size_t tempBytesPerSubjectPerBuffer = sizeof(float2) * SDIV(currentQueryLength,32) * 32;
                                        const size_t maxSubjectsPerIteration = std::min(size_t(numSeq), ws.numTempBytes / (tempBytesPerSubjectPerBuffer * 2));
        
                                        float2* d_temp = (float2*)ws.d_tempStorageHE.data();
                                        float2* d_tempHcol2 = d_temp;
                                        float2* d_tempEcol2 = (float2*)(((char*)d_tempHcol2) + maxSubjectsPerIteration * tempBytesPerSubjectPerBuffer);
        
                                        const int numIters =  SDIV(numSeq, maxSubjectsPerIteration);
                                        for(int iter = 0; iter < numIters; iter++){
                                            const size_t begin = iter * maxSubjectsPerIteration;
                                            const size_t end = iter < numIters-1 ? (iter+1) * maxSubjectsPerIteration : numSeq;
                                            const size_t num = end - begin;
        
                                            cudaMemsetAsync(d_temp, 0, tempBytesPerSubjectPerBuffer * 2 * num, ws.workStreamForTempUsage); CUERR;
        
                                            //cudaDeviceSynchronize(); CUERR;
        
                                            //NW_local_affine_read4_float_query_Protein<32, 32><<<num, 32, 0, ws.workStreamForTempUsage>>>(
                                            call_NW_local_affine_read4_float_query_Protein_new<20>(
                                                blosumType,
                                                inputChars, 
                                                d_scores, 
                                                d_tempHcol2, 
                                                d_tempEcol2, 
                                                inputOffsets, 
                                                inputLengths, 
                                                d_selectedPositions + begin, 
                                                num,
                                                d_query,
                                                currentQueryLength, 
                                                gop, 
                                                gex,
                                                ws.workStreamForTempUsage
                                            ); CUERR 
                                        }
                                    }else{
                                        assert(false);
                                    }
                                }
        
                                // std::cout << "partId = " << partId << ",numSeq " << numSeq << "\n";
                                //cudaDeviceSynchronize(); CUERR;
        
                                //exclPs += numSeq;
                            }
                        };
        
                        auto maxReduceArray = ws.getMaxReduceArray(variables.processedSequences);
                        //runAlignmentKernels(ws.devAlignmentScoresFloat.data() + variables.processedSequences, d_overflow_positions, d_overflow_number);
                        runAlignmentKernels(maxReduceArray, d_overflow_positions, d_overflow_number);
        
        
                        //alignments are done in workstreams. now, join all workstreams into workStreamForTempUsage to process overflow alignments
                        for(auto& stream : ws.workStreamsWithoutTemp){
                            cudaEventRecord(ws.forkStreamEvent, stream); CUERR;
                            cudaStreamWaitEvent(ws.workStreamForTempUsage, ws.forkStreamEvent, 0); CUERR;    
                        }
                    }
                }
        
                //process overflow alignments and finish processing of batch
                for(int gpu = 0; gpu < numGpus; gpu++){
                    cudaSetDevice(deviceIds[gpu]); CUERR;
                    auto& ws = *workingSets[gpu];
                    const auto& variables = variables_vec[gpu];
                    if(variables.processedBatches < variables.batchPlansPtr->size()){
        
                        const char* const inputChars = variables.d_inputChars;
                        const size_t* const inputLengths = variables.d_inputLengths;
                        const size_t* const inputOffsets = variables.d_inputOffsets;
                        int* const d_overflow_number = variables.d_overflow_number;
                        size_t* const d_overflow_positions = variables.d_overflow_positions;
        
                        const char4* const d_query = reinterpret_cast<char4*>(ws.d_query.data());
        
                        auto maxReduceArray = ws.getMaxReduceArray(variables.processedSequences);
        
                        if(kernelTypeConfig.overflowType == KernelType::Float){
                            //std::cerr << "overflow processing\n";
                            float2* d_temp = (float2*)ws.d_tempStorageHE.data();
        
                            //launch_process_overflow_alignments_kernel_NW_local_affine_read4_float_query_Protein<32,32><<<1,1,0, ws.workStreamForTempUsage>>>(
                            if(hostBlosumDim == 21){
                                launch_process_overflow_alignments_kernel_NW_local_affine_read4_float_query_Protein_new<20,21><<<1,1,0, ws.workStreamForTempUsage>>>(
                                    d_overflow_number,
                                    d_temp, 
                                    ws.numTempBytes,
                                    inputChars, 
                                    //ws.devAlignmentScoresFloat.data() + variables.processedSequences, 
                                    maxReduceArray,
                                    inputOffsets, 
                                    inputLengths, 
                                    d_overflow_positions, 
                                    d_query,
                                    currentQueryLength, 
                                    gop, 
                                    gex
                                ); CUERR
                            #ifdef CAN_USE_FULL_BLOSUM
                            }else if(hostBlosumDim == 25){
                                launch_process_overflow_alignments_kernel_NW_local_affine_read4_float_query_Protein_new<20,25><<<1,1,0, ws.workStreamForTempUsage>>>(
                                    d_overflow_number,
                                    d_temp, 
                                    ws.numTempBytes,
                                    inputChars, 
                                    //ws.devAlignmentScoresFloat.data() + variables.processedSequences, 
                                    maxReduceArray,
                                    inputOffsets, 
                                    inputLengths, 
                                    d_overflow_positions, 
                                    d_query,
                                    currentQueryLength, 
                                    gop, 
                                    gex
                                ); CUERR
                            #endif
                            }else{
                                assert(false);
                            }
                        }else{
                            assert(false);
                        }
        
                        //update total num overflows for query
                        addKernel<<<1,1,0, ws.workStreamForTempUsage>>>(ws.d_total_overflow_number.data(), ws.d_total_overflow_number.data(), d_overflow_number); CUERR;
        
                        //after processing overflow alignments, the batch is done and its data can be resused
                        cudaEventRecord(ws.deviceBufferEvents[variables.currentBuffer], ws.workStreamForTempUsage); CUERR;
        
                        //let other workstreams depend on temp usage stream
                        for(auto& stream : ws.workStreamsWithoutTemp){
                            cudaStreamWaitEvent(ws.workStreamForTempUsage, ws.deviceBufferEvents[variables.currentBuffer], 0); CUERR;    
                        }
        
                        ws.copyBufferIndex = (ws.copyBufferIndex+1) % ws.numCopyBuffers;
                    }
                }
        
                //update running numbers
                for(int gpu = 0; gpu < numGpus; gpu++){
                    auto& variables = variables_vec[gpu];
                    if(variables.processedBatches < variables.batchPlansPtr->size()){
                        const auto& plan = (*variables.batchPlansPtr)[variables.processedBatches];
                        variables.pointerSequencesOffset += plan.usedSeq;
                        variables.pointerBytesOffset += plan.usedBytes;
                        variables.processedSequences += plan.usedSeq;
                        variables.processedBatches++;
        
                        totalNumberOfProcessedSequences += plan.usedSeq;
                    } 
                }
        
            } //while not done
        
        
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                auto& ws = *workingSets[gpu];
                if(ws.canStoreFullDB && !ws.fullDBisUploaded){
                    ws.fullDBisUploaded = true;
        
                    // current offsets in d_fulldb_offsetdata store the offsets for each batch, i.e. for each batch the offsets will start again at 0
                    // compute prefix sum of d_fulldb_lengthdata to obtain the single-batch offsets
        
                    cudaMemsetAsync(ws.d_fulldb_offsetdata.data(), 0, sizeof(size_t), ws.workStreamForTempUsage); CUERR;
        
                    auto d_paddedLengths = thrust::make_transform_iterator(
                        ws.d_fulldb_lengthdata.data(),
                        RoundToNextMultiple<size_t, 4>{}
                    );
        
                    thrust::inclusive_scan(
                        thrust::cuda::par_nosync(thrust_async_allocator<char>(ws.workStreamForTempUsage)).on(ws.workStreamForTempUsage),
                        d_paddedLengths,
                        d_paddedLengths + ws.d_fulldb_lengthdata.size(),
                        ws.d_fulldb_offsetdata.data() + 1
                    );
                }
            }
        
        
        
            for(int gpu = 0; gpu < numGpus; gpu++){
                cudaSetDevice(deviceIds[gpu]); CUERR;
                auto& ws = *workingSets[gpu];
                //create dependency for gpuStreams[gpu]
                cudaEventRecord(ws.forkStreamEvent, ws.workStreamForTempUsage); CUERR;
                cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
        
                for(auto& stream : ws.workStreamsWithoutTemp){
                    cudaEventRecord(ws.forkStreamEvent, stream); CUERR;
                    cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
                }
        
                // for(auto& stream : ws.copyStreams){
                //     cudaEventRecord(ws.forkStreamEvent, stream); CUERR;
                //     cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
                // }
        
                cudaEventRecord(ws.forkStreamEvent, ws.hostFuncStream); CUERR;
                cudaStreamWaitEvent(gpuStreams[gpu], ws.forkStreamEvent, 0); CUERR;
            }
            
        
        }

        BenchmarkStats makeBenchmarkStats(double seconds, double cells, int overflows) const{
            BenchmarkStats stats;
            stats.seconds = seconds;
            stats.gcups = cells / 1000. / 1000. / 1000.;
            stats.gcups = stats.gcups / stats.seconds;
            stats.numOverflows = overflows;
            return stats;
        }

        void updateNumResultsPerQuery(){

            results_per_query = std::min(size_t(numTop), size_t(maxReduceArraySize));
            if(dbIsReady){
                results_per_query = std::min(size_t(results_per_query), fullDB.getChunks()[0].numSequences());
            }
        }

        std::vector<size_t> fullDB_numSequencesPerLengthPartition;
        std::vector<size_t> numSequencesPerGpu_total;
        std::vector<size_t> numSequencesPerGpuPrefixSum_total;

        //partition chars of whole DB amongst the gpus
        std::vector<std::vector<size_t>> numSequencesPerLengthPartitionPrefixSum_perDBchunk;
        std::vector<std::vector<DBdataView>> dbPartitionsByLengthPartitioning_perDBchunk;
        std::vector<std::vector<std::vector<DBdataView>>> subPartitionsForGpus_perDBchunk;
        std::vector<std::vector<std::vector<int>>> lengthPartitionIdsForGpus_perDBchunk;
        std::vector<std::vector<size_t>> numSequencesPerGpu_perDBchunk;
        std::vector<std::vector<size_t>> numSequencesPerGpuPrefixSum_perDBchunk;
        std::vector<CudaStream> gpuStreams;
        std::vector<CudaEvent> gpuEvents;
        std::vector<std::unique_ptr<GpuWorkingSet>> workingSets;  

        std::vector<std::vector<std::vector<DeviceBatchCopyToPinnedPlan>>> batchPlans_perChunk;
        std::vector<std::vector<std::vector<DeviceBatchCopyToPinnedPlan>>> batchPlans_fulldb_perChunk;
        int results_per_query;
        int currentQueryLength;

        bool dbIsReady{};
        AnyDBWrapper fullDB;

        mutable std::unique_ptr<SequenceLengthStatistics> dbSequenceLengthStatistics;

        //final scan results. device data resides on gpu deviceIds[0]
        MyPinnedBuffer<float> alignment_scores_float;
        MyPinnedBuffer<size_t> sorted_indices;
        MyPinnedBuffer<int> resultDbChunkIndices;
        MyPinnedBuffer<int> resultNumOverflows;
        MyDeviceBuffer<float> devAllAlignmentScoresFloat;
        MyDeviceBuffer<size_t> dev_sorted_indices;
        MyDeviceBuffer<int> d_resultNumOverflows;
        std::unique_ptr<helpers::GpuTimer> scanTimer;

        size_t totalProcessedQueryLengths{};
        size_t totalNumOverflows{};
        std::unique_ptr<helpers::GpuTimer> totalTimer;



        

        //--------------------------------------
        bool verbose = false;
        int gop = -11;
        int gex = -1;
        int numTop = 10;
        BlosumType blosumType = BlosumType::BLOSUM62_20;
        // KernelType singlePassType = KernelType::Half2;
        // KernelType manyPassType_small = KernelType::Half2;
        // KernelType manyPassType_large = KernelType::Float;
        // KernelType overflowType = KernelType::Float;

        KernelTypeConfig kernelTypeConfig;
        MemoryConfig memoryConfig;

        // size_t maxBatchBytes = 128ull * 1024ull * 1024ull;
        // size_t maxBatchSequences = 10'000'000;
        // size_t maxTempBytes = 4ull * 1024ull * 1024ull * 1024ull;
        // size_t maxGpuMem = std::numeric_limits<size_t>::max();

        
        std::vector<int> deviceIds;

    };


} //namespace cudasw4

#endif


