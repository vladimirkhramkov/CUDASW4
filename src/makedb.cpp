
#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP

#include <algorithm>
#include <iterator>
#include <iostream>
#include <chrono>

#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "config.hpp"

#include "hpc_helpers/all_helpers.cuh"
#include "sequence_io.h"
#include "dbdata.hpp"
#include "types.hpp"
#include "convert.cuh"
#include "kseqpp/kseqpp.hpp"
#include "length_partitions.hpp"
#include "mmapbuffer.hpp"

#include <thrust/iterator/transform_iterator.h>
#include <thrust/execution_policy.h>
#include <thrust/transform.h>

std::size_t getAvailableMemoryInKB_linux(){
    //https://stackoverflow.com/questions/349889/how-do-you-determine-the-amount-of-linux-system-ram-in-c
    std::string token;
    std::ifstream file("/proc/meminfo");
    assert(bool(file));
    while(file >> token) {
        if(token == "MemAvailable:") {
            std::size_t mem;
            if(file >> mem) {
                return mem;
            } else {
                return 0;       
            }
        }
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return 0;
}


std::size_t getCurrentRSS_linux(){
    std::ifstream in("/proc/self/statm");
    std::size_t tmp, rss;
    in >> tmp >> rss;
    
    return rss * sysconf(_SC_PAGESIZE);
}

std::size_t getRSSLimit_linux(){
    rlimit rlim;
    int ret = getrlimit(RLIMIT_RSS, &rlim);
    if(ret != 0){
        std::perror("Could not get RSS limit!");
        return 0;
    }
    return rlim.rlim_cur;    
}


std::size_t getAvailableMemoryInKB(){
    return std::min(getAvailableMemoryInKB_linux(), (getRSSLimit_linux() - getCurrentRSS_linux()) / 1024);
}


struct InMemoryBatch{
    std::vector<char> chars;               
    std::vector<std::size_t> offsets;  
    std::vector<cudasw4::SequenceLengthT> lengths;  
    std::vector<char> headers;  
    std::vector<std::size_t> headerOffsets;  
};

struct HybridBatch{
    static constexpr size_t charweight = 50;
    static constexpr size_t offsetweight = 7;
    static constexpr size_t lengthweight = 7;
    static constexpr size_t headersweight = 29;
    static constexpr size_t headeroffsetweight = 7;

    static_assert(charweight + offsetweight + lengthweight + headersweight + headeroffsetweight == 100);

    HybridBatch(const std::string& temppath, size_t memoryLimit) : 
        chars(0, (memoryLimit / 100) * charweight, temppath + "_cudasw4tmpchars"),
        offsets(0, (memoryLimit / 100) * offsetweight, temppath + "_cudasw4tmpoffsets"),
        lengths(0, (memoryLimit / 100) * lengthweight, temppath + "_cudasw4tmplengths"),
        headers(0, (memoryLimit / 100) * headersweight, temppath + "_cudasw4tmpheaders"),
        headerOffsets(0, (memoryLimit / 100) * headeroffsetweight, temppath + "_cudasw4tmpheaderOffsets")
    {

    }
    cudasw4::FileBackedUVector<char> chars;               
    cudasw4::FileBackedUVector<std::size_t> offsets;  
    cudasw4::FileBackedUVector<cudasw4::SequenceLengthT> lengths;  
    cudasw4::FileBackedUVector<char> headers;  
    cudasw4::FileBackedUVector<std::size_t> headerOffsets;  
};


template<class Batch>
void loadWholeFileIntoBatch_withPaddedSequences(const std::string& inputfilename, Batch& batch){
    constexpr int ALIGN = 4;

    batch.chars.clear();
    batch.offsets.clear();
    batch.lengths.clear();
    batch.headers.clear();
    batch.headerOffsets.clear();
    batch.offsets.push_back(0);
    batch.headerOffsets.push_back(0);

    kseqpp::KseqPP reader(inputfilename);
    while(reader.next() >= 0){
        const std::string& header = reader.getCurrentHeader();
        const std::string& sequence = reader.getCurrentSequence();

        const size_t sequencepadding = (sequence.size() % ALIGN == 0) ? 0 : ALIGN - sequence.size() % ALIGN;

        batch.chars.insert(batch.chars.end(), sequence.begin(), sequence.end());
        batch.chars.insert(batch.chars.end(), sequencepadding, ' ');
        batch.offsets.push_back(batch.chars.size());
        batch.lengths.push_back(sequence.size());

        batch.headers.insert(batch.headers.end(), header.begin(), header.end());
        batch.headerOffsets.push_back(batch.headers.size());
    }
}

template<class Batch>
void loadWholeFileIntoBatch_withPaddedConvertedSequences(const std::string& inputfilename, Batch& batch, const cudasw4::SequenceType& sequenceType){
    constexpr int ALIGN = 4;

    batch.chars.clear();
    batch.offsets.clear();
    batch.lengths.clear();
    batch.headers.clear();
    batch.headerOffsets.clear();
    batch.offsets.push_back(0);
    batch.headerOffsets.push_back(0);

    kseqpp::KseqPP reader(inputfilename);
    while(reader.next() >= 0){
        if(batch.lengths.size() > cudasw4::MaxSequencesInDB::value()){
            std::string msg = "File contains at least " + std::to_string(batch.lengths.size()+1) 
                    + " sequences, but config allows at most " + std::to_string(cudasw4::MaxSequencesInDB::value());
            throw std::runtime_error(msg);
        }

        const std::string& header = reader.getCurrentHeader();
        const std::string& sequence = reader.getCurrentSequence();

        size_t sequenceLength = sequence.size();
        if(sequenceLength > cudasw4::MaxSequenceLength::value()){
            std::string msg = "Got sequence of length " + std::to_string(sequenceLength) 
                + ", but config allows only lengths <= " + std::to_string(cudasw4::MaxSequenceLength::value());
            throw std::runtime_error(msg);
        }


        const size_t sequencepadding = (sequence.size() % ALIGN == 0) ? 0 : ALIGN - sequence.size() % ALIGN;

        const size_t oldCharsSize = batch.chars.size();
        const size_t newCharsSize = oldCharsSize + sequence.size() + sequencepadding;
        batch.chars.resize(newCharsSize);

        if (sequenceType == cudasw4::SequenceType::Protein) {
            auto convert = cudasw4::ConvertAA_20{};
            auto it = std::transform(sequence.begin(), sequence.end(), batch.chars.begin() + oldCharsSize, convert);
            std::fill(it, batch.chars.end(), convert(' ')); // add converted padding
        } else {
            auto convert = cudasw4::ConvertNA{};
            auto it = std::transform(sequence.begin(), sequence.end(), batch.chars.begin() + oldCharsSize, convert);
            std::fill(it, batch.chars.end(), convert(' ')); // add converted padding
        }

        batch.offsets.push_back(newCharsSize);
        batch.lengths.push_back(sequence.size());

        batch.headers.insert(batch.headers.end(), header.begin(), header.end());
        batch.headerOffsets.push_back(batch.headers.size());
    }
}

template<class Batch>
void createDBfilesFromSequenceBatch(const std::string& outputPrefix, const Batch& batch){
    using cudasw4::DBdataIoConfig;

    const size_t numSequences = batch.lengths.size();

    std::vector<cudasw4::ReferenceIdT> indices(numSequences);
    std::iota(indices.begin(), indices.end(), cudasw4::ReferenceIdT(0));

    auto compareIndicesByLength = [&](const auto& l, const auto& r){
        return batch.lengths[l] < batch.lengths[r];
    };

    std::sort(indices.begin(), indices.end(), compareIndicesByLength);

    const auto sortedLengthsIterator = thrust::make_transform_iterator(
        indices.begin(),
        [&](auto i){ return batch.lengths[i]; }
    );
    const auto sortedLengths_endIterator = sortedLengthsIterator + numSequences;

    auto lengthBoundaries = cudasw4::getLengthPartitionBoundaries();
    const int numPartitions = lengthBoundaries.size();

    std::vector<size_t> numSequencesPerPartition(numPartitions);

    auto partitionBegin = sortedLengthsIterator;
    for(int i = 0; i < numPartitions; i++){
        //length k is in partition i if boundaries[i-1] < k <= boundaries[i]
        int searchFor = lengthBoundaries[i];
        if(searchFor < std::numeric_limits<int>::max()){
            searchFor += 1;
        }
        auto partitionEnd = std::lower_bound(
            partitionBegin, 
            sortedLengths_endIterator, 
            searchFor
        );
        numSequencesPerPartition[i] = std::distance(partitionBegin, partitionEnd);
        partitionBegin = partitionEnd;
    }
    // for(int i = 0; i < numPartitions; i++){
    //     std::cout << "numInPartition " << i << " (<= " << lengthBoundaries[i] << " ) : " << numSequencesPerPartition[i] << "\n";
    // }

    //write partition data to metadata file
    std::ofstream metadataout(outputPrefix + DBdataIoConfig::metadatafilename(), std::ios::binary);
    if(!metadataout) throw std::runtime_error("Cannot open output file " + outputPrefix + DBdataIoConfig::metadatafilename());

    metadataout.write((const char*)&numPartitions, sizeof(int));
    for(int i = 0; i < numPartitions; i++){
        const int limit = lengthBoundaries[i];
        metadataout.write((const char*)&limit, sizeof(int));
    }
    metadataout.write((const char*)numSequencesPerPartition.data(), sizeof(size_t) * numPartitions);


    //write db files with sequences sorted by length

    std::ofstream headersout(outputPrefix + DBdataIoConfig::headerfilename(), std::ios::binary);
    if(!headersout) throw std::runtime_error("Cannot open output file " + outputPrefix + DBdataIoConfig::headerfilename());
    std::ofstream headersoffsetsout(outputPrefix + DBdataIoConfig::headeroffsetsfilename(), std::ios::binary);
    if(!headersoffsetsout) throw std::runtime_error("Cannot open output file " + outputPrefix + DBdataIoConfig::headeroffsetsfilename());
    std::ofstream charsout(outputPrefix + DBdataIoConfig::sequencesfilename(), std::ios::binary);
    if(!charsout) throw std::runtime_error("Cannot open output file " + outputPrefix + DBdataIoConfig::sequencesfilename());
    std::ofstream offsetsout(outputPrefix + DBdataIoConfig::sequenceoffsetsfilename(), std::ios::binary);
    if(!offsetsout) throw std::runtime_error("Cannot open output file " + outputPrefix + DBdataIoConfig::sequenceoffsetsfilename());
    std::ofstream lengthsout(outputPrefix + DBdataIoConfig::sequencelengthsfilename(), std::ios::binary);
    if(!lengthsout) throw std::runtime_error("Cannot open output file " + outputPrefix + DBdataIoConfig::sequencelengthsfilename());

    size_t currentHeaderOffset = 0;
    size_t currentCharOffset = 0;
    headersoffsetsout.write((const char*)&currentHeaderOffset, sizeof(size_t));
    offsetsout.write((const char*)&currentCharOffset, sizeof(size_t));
    for(size_t i = 0; i < numSequences; i++){
        const size_t sortedIndex = indices[i];

        const char* const header = batch.headers.data() + batch.headerOffsets[sortedIndex];
        const int headerLength = batch.headerOffsets[sortedIndex+1] - batch.headerOffsets[sortedIndex];

        headersout.write(header, headerLength);
        currentHeaderOffset += headerLength;
        headersoffsetsout.write((const char*)&currentHeaderOffset, sizeof(size_t));

        const size_t numChars = batch.offsets[sortedIndex+1] - batch.offsets[sortedIndex];
        const cudasw4::SequenceLengthT length = batch.lengths[sortedIndex];
        const char* const sequence = batch.chars.data() + batch.offsets[sortedIndex];


        charsout.write(sequence, numChars);
        lengthsout.write((const char*)&length, sizeof(cudasw4::SequenceLengthT));
        currentCharOffset += numChars;
        offsetsout.write((const char*)&currentCharOffset, sizeof(size_t));
    }
}


int main(int argc, char* argv[])
{
    if(argc < 3) {
        std::cout << "Usage:\n  " << argv[0] << " <FASTA/FASTQ filename> pathtodb/dbname [options]\n";
        std::cout << "Input file may be gzip'ed. pathtodb must exist.\n";
        std::cout << "Options:\n";
        std::cout << "    --type val : Sequence type: (n|p). Nucleotides (n) or proteins (p) . Default value: p.\n";
        std::cout << "    --mem val : Memory limit. Can use suffix K,M,G. If makedb requires more memory, temp files in temp directory will be used. Default all available memory.\n";
        std::cout << "    --tempdir val : Temp directory for temporary files. Must exist. Default is db output directory.\n";
        return 0;
    }

    auto parseMemoryString = [](const std::string& string){
        std::size_t result = 0;
        if(string.length() > 0){
            std::size_t factor = 1;
            bool foundSuffix = false;
            switch(string.back()){
                case 'K':{
                    factor = std::size_t(1) << 10; 
                    foundSuffix = true;
                }break;
                case 'M':{
                    factor = std::size_t(1) << 20;
                    foundSuffix = true;
                }break;
                case 'G':{
                    factor = std::size_t(1) << 30;
                    foundSuffix = true;
                }break;
            }
            if(foundSuffix){
                const auto numberString = string.substr(0, string.size()-1);
                result = factor * std::stoull(numberString);
            }else{
                result = std::stoull(string);
            }
        }else{
            result = 0;
        }
        return result;
    };

    const std::string fastafilename = argv[1];
    const std::string outputPrefix = argv[2];
    std::string temppath = outputPrefix;
    cudasw4::SequenceType sequenceType = cudasw4::SequenceType::Protein;
    size_t availableMem = getAvailableMemoryInKB() * 1024;
    constexpr size_t GB = 1024*1024*1024;
    if(availableMem > 1*GB){
        availableMem -= 1*GB;
    }

    for(int i = 3; i < argc; i++){
        const std::string arg = argv[i];
        if(arg == "--type"){
            if(i + 1 < argc) { // Make sure there is a value after --type
                std::string typeArg = argv[++i];
                if(typeArg == "n") {
                    sequenceType = cudasw4::SequenceType::Nucleotide;
                }
            } else {
                std::cerr << "No sequence type specified after --type." << std::endl;
                return 1;
            }
        }else if(arg == "--mem"){
            availableMem = parseMemoryString(argv[++i]);
        }else if(arg == "--tempdir"){
           temppath = argv[++i];
           if(temppath.back() != '/'){
                temppath += '/';
           }
        }else{
            std::cout << "Unexpected arg " << arg << "\n";
        }
    }
    std::cout << "availableMem: " << availableMem << "\n";

    //InMemoryBatch batch;
    HybridBatch batch(temppath, availableMem);

    std::cout << "Parsing file\n";
    helpers::CpuTimer timer1("file parsing");
    //loadWholeFileIntoBatch_withPaddedConvertedSequences(fastafilename, batch, sequenceType);
    loadWholeFileIntoBatch_withPaddedSequences(fastafilename, batch);
    timer1.print();

    std::cout << "Number of input sequences:  " << batch.offsets.size() - 1 << '\n';
    std::cout << "Number of input characters: " << batch.chars.size() << '\n';

    std::cout << "Converting sequences\n";
    helpers::CpuTimer timer2("sequences conversion");

    if (sequenceType == cudasw4::SequenceType::Nucleotide) {
        thrust::transform(thrust::omp::par, batch.chars.begin(), batch.chars.end(), batch.chars.begin(), cudasw4::ConvertNA{});
    } else {
        thrust::transform(thrust::omp::par, batch.chars.begin(), batch.chars.end(), batch.chars.begin(), cudasw4::ConvertAA_20{});
    }
    timer2.print();

    std::cout << "Creating DB files\n";
    const std::string batchOutputPrefix = outputPrefix + std::to_string(0);
    helpers::CpuTimer timer3("db creation");
    createDBfilesFromSequenceBatch(batchOutputPrefix, batch);
    timer3.print();

    cudasw4::DBGlobalInfo info;

    cudasw4::writeGlobalDbInfo(outputPrefix, info);

}
