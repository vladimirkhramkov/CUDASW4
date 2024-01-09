#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include "types.hpp"
#include <string>
#include <iostream>


namespace cudasw4 {
    
    struct QuerySequence {
        std::string header;
        std::string sequence;

        // Default constructor
        QuerySequence() : header(""), sequence("") {}

        // Constructor with parameters
        QuerySequence(const std::string& h, const std::string& s) : header(h), sequence(s) {}
    };  

}

struct ProgramOptions{
    enum class OutputMode{
        Plain,
        CSV
    };

    bool version = false;
    bool help = false;
    bool loadFullDBToGpu = false;
    bool usePseudoDB = false;
    bool printLengthPartitions = false;
    bool verbose = false;
    bool prefetchDBFile = false;
    bool reverseComplement = false;

    int numTopOutputs = 10;
    int gop = -11;
    int gex = -1;
    int pseudoDBLength = 0;
    int pseudoDBSize = 0;

    int minScore = 1;

    cudasw4::SequenceType sequenceType = cudasw4::SequenceType::Protein;
    cudasw4::SubMatrixType subMatrixType = cudasw4::SubMatrixType::BLOSUM62_20;
    cudasw4::KernelType singlePassType = cudasw4::KernelType::Half2;
    cudasw4::KernelType manyPassType_small = cudasw4::KernelType::Half2;
    cudasw4::KernelType manyPassType_large = cudasw4::KernelType::Float;
    cudasw4::KernelType overflowType = cudasw4::KernelType::Float;
    OutputMode outputMode = OutputMode::CSV;

    size_t maxBatchBytes = 128ull * 1024ull * 1024ull;
    size_t maxBatchSequences = 10'000'000;
    size_t maxTempBytes = 4ull * 1024ull * 1024ull * 1024ull;

    size_t maxGpuMem = std::numeric_limits<size_t>::max();

    std::string outputfile = "/dev/stdout";
    std::vector<std::string> databases;

    std::string progressKey = "";
    std::string progressPipePath = "";

    std::vector<std::string> csvColumns;

    std::vector<cudasw4::QuerySequence> queries;

    std::string outputModeString() const{
        switch(outputMode){
            case OutputMode::Plain: return "Plain";
            case OutputMode::CSV: return "CSV";
            default: return "Unnamed output mode";
        }
    }

    bool isDnaDatabase() const {
        return (sequenceType == cudasw4::SequenceType::Nucleotide);
    }
};

void printOptions(const ProgramOptions& options);

bool parseArgs(int argc, char** argv, ProgramOptions& options);

void printVersion();

void printHelp(int argc, char** argv);

#endif