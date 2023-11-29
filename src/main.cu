#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "hpc_helpers/all_helpers.cuh"
#include "hpc_helpers/peer_access.cuh"

// #include "kseqpp/kseqpp.hpp"
// #include "sequence_io.h"
#include "options.hpp"
#include "dbdata.hpp"
#include "cudasw4.cuh"
#include "config.hpp"
#include "reverse.hpp"

std::vector<std::string> split(const std::string& str, char c){
	std::vector<std::string> result;

	std::stringstream ss(str);
	std::string s;

	while (std::getline(ss, s, c)) {
		result.emplace_back(s);
	}

	return result;
}

void printScanResultPlain(
    std::ostream& os, 
    const cudasw4::ScanResult& scanResult, 
    const cudasw4::CudaSW4& cudaSW4, 
    const ProgramOptions& options
){
    const int n = scanResult.scores.size();
    for(int i = 0; i < n; i++){
        if (scanResult.scores[i] < options.minScore) continue;

        const auto referenceId = scanResult.referenceIds[i];
        os << "Result " << i << ".";
        os << " Score: " << scanResult.scores[i] << ".";
        os << " Length: " << cudaSW4.getReferenceLength(referenceId) << ".";
        os << " Header " << cudaSW4.getReferenceHeader(referenceId) << ".";
        os << " referenceId " << referenceId;
        os << "\n";
        //std::cout << " Sequence " << cudaSW4.getReferenceSequence(referenceId) << "\n";

    }
}

void printScanResultCSV(
    std::ostream& os, 
    const cudasw4::ScanResult& scanResult, 
    const cudasw4::CudaSW4& cudaSW4, 
    const ProgramOptions& options,
    int64_t queryId,
    cudasw4::SequenceLengthT queryLength,
    std::string_view queryHeader
){
    constexpr char sep = ',';

    const int n = scanResult.scores.size();
    for(int i = 0; i < n; i++){
        if (scanResult.scores[i] < options.minScore) continue;

        const auto referenceId = scanResult.referenceIds[i];
        
        // os << queryId << sep 
        //     << queryLength << sep
        //     << queryHeader << sep
        //     << i << sep
        //     << scanResult.scores[i] << sep
        //     << cudaSW4.getReferenceLength(referenceId) << sep
        //     << cudaSW4.getReferenceHeader(referenceId) << sep
        //     << referenceId << "\n";

        //std::cout << " Sequence " << cudaSW4.getReferenceSequence(referenceId) << "\n";

        for (unsigned int paramIndex = 0; paramIndex < options.csvColumns.size(); ++paramIndex) {
            std::string parameter = options.csvColumns[paramIndex];
            if (paramIndex)
                os << sep;

            if (parameter == "qacc") {
                os << queryHeader || queryId; //initialQueries[answerEntry.queryIndex].first;
            } else if (parameter == "qlen") {
                os << queryLength; //initialQueries[answerEntry.queryIndex].second.length();
            } else if (parameter == "sacc") {
                os << cudaSW4.getReferenceHeader(referenceId);//answerEntry.getSacc();
            } else if (parameter == "slen") {
                os << cudaSW4.getReferenceLength(referenceId); //answerEntry.getSlen();
            } else if (parameter == "score") {
                os << scanResult.scores[i]; //answerEntry.getScore();
            } else if (parameter == "length") {
                // os << stringAlignment[0].length();
            } else if (parameter == "nident") {
                // int matches = 0;
                // for (int i = 0; i < stringAlignment[0].length(); ++i)
                //     if (stringAlignment[0][i] != '-' && (stringAlignment[0][i] == stringAlignment[2][i] || stringAlignment[1][i] == '|' && nucleotideSearch))
                //         ++matches;
                // os << matches;
            } else if (parameter == "gaps") {
                // int gaps = 0;
                // for (int i = 0; i < stringAlignment[0].length(); ++i)
                //     gaps += (stringAlignment[0][i] == '-')
                //             + (stringAlignment[2][i] == '-');
                // os << gaps;
            } else if (parameter == "qstart") {
                // os << answerEntry.getQstart() + 1;
            } else if (parameter == "qend") {
                // os << answerEntry.getQend() + 1;
            } else if (parameter == "sstart") {
                // os << answerEntry.getSstart() + 1;
            } else if (parameter == "send") {
                // os << answerEntry.getSend() + 1;
            } else if (parameter == "positive") {
                // int positiveMatches = 0;
                // for (int i = 0; i < stringAlignment[0].length(); ++i) {
                //     if (stringAlignment[0][i] != '-'
                //             && stringAlignment[2][i] != '-') {
                //         positiveMatches +=
                //                 scoringMatrix[CFastaFile::amino_acids_trans[stringAlignment[0][i]]][CFastaFile::amino_acids_trans[stringAlignment[2][i]]]
                //                         > 0;
                //     }
                // }
                // os << positiveMatches;
            } else if (parameter == "btop") {
                // os << BTOP;
            } else if (parameter == "topline") {
                // os << stringAlignment[0];
            } else if (parameter == "middleline") {
                // os << stringAlignment[1];
            } else if (parameter == "bottomline") {
                // os << stringAlignment[2];
            } else if (parameter == "reversed") {
                // os << answerEntry.reversed;
            } else if (parameter == "qcovs") {
                // int coverageLength = queryCoverageLength[answerEntry.sequenceName];
                // os << fixed << setprecision(2) << 100.0 * coverageLength / query.length();
                // rounded percent
                //os << (200 * coverageLength + query.length()) / (2 * query.length());
            }
        }
        os << "\n";

    }
}

struct BatchOfQueries{
    std::vector<char> chars;               
    std::vector<std::size_t> offsets;  
    std::vector<cudasw4::SequenceLengthT> lengths;  
    std::vector<std::string> headers;  
};

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

    if(!options.usePseudoDB){
        if(options.verbose){
            std::cout << "Reading Database: \n";
        }
        try{
            helpers::CpuTimer timer_read_db("Read DB");
            constexpr bool writeAccess = false;
            const bool prefetchSeq = options.prefetchDBFile;

            auto fullDB_tmp = std::make_shared<cudasw4::DB>(cudasw4::loadDB(options.dbPrefix, writeAccess, prefetchSeq));
            if(options.verbose){
                timer_read_db.print();
            }

            cudaSW4.setDatabase(fullDB_tmp);
        }catch(cudasw4::LoadDBException& ex){
            if(options.verbose){
                std::cout << "Failed to map db files. Using fallback db. Error message: " << ex.what() << "\n";
            }
            helpers::CpuTimer timer_read_db("Read DB");
            auto fullDB_tmp = std::make_shared<cudasw4::DBWithVectors>(cudasw4::loadDBWithVectors(options.dbPrefix));
            if(options.verbose){
                timer_read_db.print();
            }

            cudaSW4.setDatabase(fullDB_tmp);
        }
    }

    if(options.verbose){
        cudaSW4.printDBInfo();
        if(options.printLengthPartitions){
            cudaSW4.printDBLengthPartitions();
        }
    }

    if(options.loadFullDBToGpu){
        cudaSW4.prefetchFullDBToGpus();
    }

    int64_t query_num = 0;
    for(const auto& query : options.queries){

        cudaSW4.totalTimerStart();

        std::cout << "Processing query " << query_num << " ... ";
        std::cout.flush();
        const std::string& header   = query.header;
        const std::string& sequence = query.sequence;

        ScanResult scanResult = cudaSW4.scan(sequence.data(), sequence.size());
        if(options.verbose){
            std::cout << "Done. Scan time: " << scanResult.stats.seconds << " s, " << scanResult.stats.gcups << " GCUPS\n";
        }else{
            std::cout << "Done.\n";
        }

        if(options.numTopOutputs > 0){
            if(options.outputMode == ProgramOptions::OutputMode::Plain){
                printScanResultPlain(outputfile, scanResult, cudaSW4, options);
            }else{
                printScanResultCSV(outputfile, scanResult, cudaSW4, options, query_num, sequence.size(), header);
            }
            outputfile.flush();
        }

        query_num++;

        auto totalBenchmarkStats = cudaSW4.totalTimerStop();
        if(options.verbose){
            std::cout << "Total time: " << totalBenchmarkStats.seconds << " s, " << totalBenchmarkStats.gcups << " GCUPS\n";
        }
    }

    if (progressFileDescriptor != -1) close(progressFileDescriptor);
}
