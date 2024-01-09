#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>

#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/peer_access.cuh"

// #include "kseqpp/kseqpp.hpp"
// #include "sequence_io.h"
#include "types.hpp"
#include "options.hpp"
#include "dbdata.hpp"
#include "cudasw4.cuh"
#include "config.hpp"
#include "reverse.hpp"
#include "alignments.hpp"

int main(int argc, char* argv[])
{
    ProgramOptions options;
    bool parseSuccess = parseArgs(argc, argv, options);

    if (options.version) {
        printVersion();
        return 0;
    }

    if(!parseSuccess || options.help) {
        printHelp(argc, argv);
        return 0;
    }

    if (options.reverseComplement) {
        if(options.sequenceType == cudasw4::SequenceType::Nucleotide) {
            const size_t numQueries = options.queries.size();
            for (size_t i = 0; i < numQueries; i++) {
                std::string reverseComplement = getReverseComplement(options.queries[i].sequence);
                options.queries.push_back(cudasw4::QuerySequence(options.queries[i].header + "_reversed", reverseComplement));
            }
        } else {
            options.reverseComplement = false;
        }
    }

    printOptions(options);

    std::vector<int> deviceIds;
    {
        int num = 0;
        cudaGetDeviceCount(&num); CUERR
        for(int i = 0; i < num; i++){
            deviceIds.push_back(i);
        }
        if(deviceIds.size() > 0){
            if(options.verbose){
                std::cout << "Will use GPU";
                for(auto x : deviceIds){
                    std::cout << " " << x;
                }
                std::cout << "\n";
            }
        }else{
            throw std::runtime_error("No GPU found");
        }
    }

    helpers::PeerAccess peerAccess(deviceIds, false);
 
    using KernelTypeConfig = cudasw4::KernelTypeConfig;
    using MemoryConfig = cudasw4::MemoryConfig;
    using ScanResult = cudasw4::ScanResult;

    KernelTypeConfig kernelTypeConfig;
    kernelTypeConfig.singlePassType = options.singlePassType;
    kernelTypeConfig.manyPassType_small = options.manyPassType_small;
    kernelTypeConfig.manyPassType_large = options.manyPassType_large;
    kernelTypeConfig.overflowType = options.overflowType;

    MemoryConfig memoryConfig;
    memoryConfig.maxBatchBytes = options.maxBatchBytes;
    memoryConfig.maxBatchSequences = options.maxBatchSequences;
    memoryConfig.maxTempBytes = options.maxTempBytes;
    memoryConfig.maxGpuMem = options.maxGpuMem;

    std::ofstream outputfile(options.outputfile);
    if(!bool(outputfile)){
        throw std::runtime_error("Cannot open file " + options.outputfile);
    }

    int progressFileDescriptor = options.progressPipePath.length() == 0 ? -1 : open(options.progressPipePath.c_str(), O_WRONLY | O_NONBLOCK);

    cudasw4::CudaSW4 cudaSW4(
        deviceIds, 
        options.numTopOutputs,
        options.sequenceType,
        options.subMatrixType, 
        kernelTypeConfig, 
        memoryConfig, 
        options.verbose,
        options.progressKey,
        progressFileDescriptor
    );

    size_t db_count = options.databases.size();

    std::unordered_set<std::string> uniqueHeaders;
    std::vector<cudasw4::HitResult> totalResults;

    for (size_t db_index = 0; db_index < db_count; db_index ++) {

        if(options.verbose) {
            std::cout << "Reading Database: \n";
        }
        try {
            helpers::CpuTimer timer_read_db("Read DB");
            constexpr bool writeAccess = false;
            const bool prefetchSeq = options.prefetchDBFile;

            auto fullDB_tmp = std::make_shared<cudasw4::DB>(cudasw4::loadDB(options.databases[db_index], writeAccess, prefetchSeq));
            if(options.verbose){
                timer_read_db.print();
            }

            cudaSW4.setDatabase(fullDB_tmp);
        } catch(cudasw4::LoadDBException& ex) {
            if(options.verbose){
                std::cout << "Failed to map db files. Using fallback db. Error message: " << ex.what() << "\n";
            }
            helpers::CpuTimer timer_read_db("Read DB");
            auto fullDB_tmp = std::make_shared<cudasw4::DBWithVectors>(cudasw4::loadDBWithVectors(options.databases[db_index]));
            if(options.verbose){
                timer_read_db.print();
            }

            cudaSW4.setDatabase(fullDB_tmp);
        }

        if (options.verbose) {
            cudaSW4.printDBInfo();
            if(options.printLengthPartitions){
                cudaSW4.printDBLengthPartitions();
            }
        }

        if (options.loadFullDBToGpu) {
            cudaSW4.prefetchFullDBToGpus();
        }

        size_t queries_count = options.queries.size();
        size_t query_num = 0;
        for (const auto& query : options.queries) {

            cudaSW4.totalTimerStart();

            cudasw4::DbQueryIndex scan_index(db_index, db_count, query_num, queries_count);

            std::cout << "Processing query " << query_num << " ... ";
            std::cout.flush();
            const std::string& header   = query.header;
            const std::string& sequence = query.sequence;

            ScanResult scanResult = cudaSW4.scan(scan_index, sequence.data(), sequence.size());

            if(options.verbose){
                std::cout << "Done. Scan time: " << scanResult.stats.seconds << " s, " << scanResult.stats.gcups << " GCUPS\n";
            }else{
                std::cout << "Done.\n";
            }

            for (size_t i = 0; i < scanResult.scores.size(); ++i) {
                if (scanResult.scores[i] < options.minScore) break;

                std::string refHeader = (std::string)cudaSW4.getReferenceHeader(scanResult.referenceIds[i]);

                if (uniqueHeaders.find(refHeader) == uniqueHeaders.end()) {
                    // The header is not in the set, so it's a new unique element
                    // std::string refSequence = (std::string)cudaSW4.getReferenceSequence(scanResult.referenceIds[i]);
                    totalResults.push_back(
                        cudasw4::HitResult(
                            query_num, 
                            header,
                            sequence.size(),
                            sequence,
                            refHeader,
                            cudaSW4.getReferenceLength(scanResult.referenceIds[i]),
                            //TODO: think about how to store subject sequence for final alignment calculation
                            cudaSW4.getReferenceSequence(scanResult.referenceIds[i]), 
                            scanResult.scores[i]
                        )
                    );
                    uniqueHeaders.insert(refHeader); // Add the header to the set
                }
            }

            query_num++;

            auto totalBenchmarkStats = cudaSW4.totalTimerStop();
            if(options.verbose){
                std::cout << "Total time: " << totalBenchmarkStats.seconds << " s, " << totalBenchmarkStats.gcups << " GCUPS\n";
            }
        }
    }

    if(options.numTopOutputs > 0){
        std::sort(totalResults.begin(), totalResults.end(), [](const auto& a, const auto& b) {
            return a.score > b.score;
        });

        cudasw2::calculateAlignments(
            outputfile, 
            totalResults, 
            options,
            cudasw4::hostSubMatrix, 
            cudasw4::hostSubMatrixDim
        );

        // cudasw2::printScanResultCSV(outputfile, totalResults, options);
        
        outputfile.flush();
    }

    if (progressFileDescriptor != -1) close(progressFileDescriptor);
}
