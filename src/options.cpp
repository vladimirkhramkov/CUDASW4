#include "options.hpp"
#include "types.hpp"
#include "hpc_helpers/all_helpers.cuh"
#include "version.hpp"
#include "kseqpp/kseqpp.hpp"

#include <string>
#include <iostream>
#include <sstream>
#include <numeric> // For std::accumulate

void printOptions(const ProgramOptions& options) {
    std::string delimiter = " ";

    std::cout << "Selected options:\n";
    std::cout << "verbose: " << options.verbose << "\n";
    std::cout << "numTopOutputs: " << options.numTopOutputs << "\n";
    std::cout << "gop: " << options.gop << "\n";
    std::cout << "gex: " << options.gex << "\n";

    for(size_t i = 0; i < options.queries.size(); i++){
        std::cout << "query " << i  << " : " << options.queries[i].header << "\n";
    }
    std::cout << "sequence type: " << to_string(options.sequenceType) << "\n";
    #ifdef CAN_USE_FULL_BLOSUM
    std::cout << "matrix: " << to_string(options.subMatrixType) << "\n";
    #else
    std::cout << "matrix: " << to_string_nodim(options.subMatrixType) << "\n";
    #endif
    std::cout << "reverse complement: " << options.reverseComplement << "\n";

    std::cout << "using db file: " << std::accumulate(
        std::next(options.databases.begin()), options.databases.end(),
        options.databases[0], // Start with the first element
        [&delimiter](const std::string &a, const std::string &b) {
            return a + delimiter + b;
        }
    ) << "\n";

    std::cout << "output mode: " << options.outputModeString() << "\n";
    std::cout << "output file: " << options.outputfile << "\n";
}

bool parseArgs(int argc, char** argv, ProgramOptions& options) {

    bool gotQuery = false;
    bool gotDB = false;
    bool gotGex = false;
    bool gotGop = false;
    bool gotDPX = false;

    options.queries.clear();

    for(int i = 1; i < argc; i++) {
        const std::string arg = argv[i];
        if(arg == "-help") {
            options.help = true;
        } else if(arg == "-version") {
            options.version = true;
        } else if(arg == "-verbose") {
            options.verbose = true;            
        } else if(arg == "-min_score") {
            options.minScore = std::atoi(argv[++i]);
        } else if(arg == "-top" || arg == "-topscore_num") {
            options.numTopOutputs = std::atoi(argv[++i]);
        } else if(arg == "-gop" || arg == "-gapo") {
            options.gop = - std::abs(std::atoi(argv[++i]));
            gotGop = true;
        } else if(arg == "-gex" || arg == "-gape") {
            options.gex = - std::abs(std::atoi(argv[++i]));
            gotGex = true;
        } else if(arg == "-query") {
            kseqpp::KseqPP reader(argv[++i]);

            while(reader.next() >= 0) {
                cudasw4::QuerySequence query(reader.getCurrentHeader(), reader.getCurrentSequence());

                options.queries.push_back(query);
            }

            if (options.queries.size() > 0) gotQuery = true;
        } else if(arg == "-db") {
            options.databases.clear();
            std::string names = "";

            for (size_t j = 0; j < strlen(argv[i + 1]); j++) {
                if (argv[i+1][j] == ' ') {
                    if (names != "")
                        options.databases.push_back(names);

                    names = "";
                } else
                    names += argv[i + 1][j];
            }

            if (names.size() > 0)
                options.databases.push_back(names);

            i++;

            if (options.databases.size() > 0) gotDB = true;
        } else if(arg == "-mat") {
            const std::string val = argv[++i];

            options.sequenceType  = cudasw4::SequenceType::Protein; // by default

            #ifdef CAN_USE_FULL_BLOSUM
            if(val == "pam30") options.subMatrixType = cudasw4::SubMatrixType::PAM30;
            if(val == "pam70") options.subMatrixType = cudasw4::SubMatrixType::PAM70;
            if(val == "blosum45") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM45;
            if(val == "blosum50") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM50;
            if(val == "blosum62") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM62;
            if(val == "blosum80") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM80;
            #else
            if(val == "pam30") options.subMatrixType = cudasw4::SubMatrixType::PAM30_20;
            if(val == "pam70") options.subMatrixType = cudasw4::SubMatrixType::PAM70_20;
            if(val == "blosum45") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM45_20;
            if(val == "blosum50") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM50_20;
            if(val == "blosum62") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM62_20;
            if(val == "blosum80") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM80_20;
            #endif

            if(val == "dna") {
                options.subMatrixType = cudasw4::SubMatrixType::DNA;
                options.sequenceType  = cudasw4::SequenceType::Nucleotide;
            }

            if(val == "nuc44") {
                options.subMatrixType = cudasw4::SubMatrixType::NUC44;
                options.sequenceType  = cudasw4::SequenceType::Nucleotide;
            }

            if(val == "pam30_20") options.subMatrixType = cudasw4::SubMatrixType::PAM30_20;
            if(val == "pam70_20") options.subMatrixType = cudasw4::SubMatrixType::PAM70_20;
            if(val == "blosum45_20") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM45_20;
            if(val == "blosum50_20") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM50_20;
            if(val == "blosum62_20") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM62_20;
            if(val == "blosum80_20") options.subMatrixType = cudasw4::SubMatrixType::BLOSUM80_20;
        } else if(arg == "-use_single") {
            // use single GPU
        } else if(arg == "-reverse") {
            // calculate scores and alignments also for reverse complement
            options.reverseComplement = true;
        } else if(arg == "-dpx") {
            // use DPX instructions (Hopper GPU cards: H100)
            gotDPX = true;
        } else if(arg == "-plain") {
            options.outputMode = ProgramOptions::OutputMode::Plain;
        } else if(arg == "-out") {
            options.outputfile = argv[++i];
        } else if(arg == "-outfmt") {
            // qacc, qlen, sacc, slen, score, length, nident, gaps, qstart, qend, sstart, send, positive, btop, topline, middleline, bottomline, reversed
            // options.outputFormat = argv[++i];
        	std::stringstream ss(argv[++i]);
        	std::string columnName;
        	const int availableColumnNamesCount = 19;
        	std::string availableColumnNames[availableColumnNamesCount] = {"qacc", "qlen", "sacc", "slen", "score", "length", "nident",
																	  "gaps", "qstart", "qend", "sstart", "send", "positive", "btop",
																	  "topline", "middleline", "bottomline", "reversed", "qcovs"};
        	while (ss >> columnName) {
        		bool isSupported = false;
        		for (int i = 0; i < availableColumnNamesCount; ++i) if (columnName == availableColumnNames[i]) { isSupported = true; break; }
        		if (!isSupported) {
        			std::cout << "Unknown outfmt option " << columnName << std::endl;
        			return false;
        		}
        		options.csvColumns.push_back(columnName);
        	}        
        } else if(arg == "-progress_key") {
            // search identifier used in search progress notifications
            options.progressKey = std::string(argv[++i]);
        } else if(arg == "-progress_pipe") {
            // progress pipe file path
            options.progressPipePath = std::string(argv[++i]);
        } else {
            std::cout << "Unexpected arg " << arg << "\n";
        }
    }

    //set specific gop gex for substitution matrix if no gop gex was set
    if (options.subMatrixType == cudasw4::SubMatrixType::BLOSUM45 || options.subMatrixType == cudasw4::SubMatrixType::BLOSUM45_20) {
        if(!gotGop) options.gop = -13;
        if(!gotGex) options.gex = -2;
    }
    if (options.subMatrixType == cudasw4::SubMatrixType::BLOSUM50 || options.subMatrixType == cudasw4::SubMatrixType::BLOSUM50_20) {
        if(!gotGop) options.gop = -13;
        if(!gotGex) options.gex = -2;
    }
    if (options.subMatrixType == cudasw4::SubMatrixType::BLOSUM62 || options.subMatrixType == cudasw4::SubMatrixType::BLOSUM62_20) {
        if(!gotGop) options.gop = -11;
        if(!gotGex) options.gex = -1;
    }
    if (options.subMatrixType == cudasw4::SubMatrixType::BLOSUM80 || options.subMatrixType == cudasw4::SubMatrixType::BLOSUM80_20) {
        if(!gotGop) options.gop = -10;
        if(!gotGex) options.gex = -1;
    }

    if (gotDPX) {
        options.singlePassType = cudasw4::KernelType::DPXs16;
        options.manyPassType_small = cudasw4::KernelType::DPXs16;
        options.manyPassType_large = cudasw4::KernelType::DPXs32;
        options.overflowType = cudasw4::KernelType::DPXs32;
    }

    if (!gotQuery) {
        // std::cout << "Query is missing\n";
        return false;
    }

    if (!gotDB) {
        // std::cout << "DB prefix is missing\n";
        return false;
    }

    return true;
}

void printVersion() {
    std::cout << "Version: " << cudasw4::PROGRAM_VERSION << std::endl;
}

void printHelp(int /*argc*/, char** argv) {
    ProgramOptions defaultoptions;

    std::cout << "Usage: " << argv[0] << " [options]\n";
    std::cout << "Options: \n";

    std::cout << "   Mandatory\n";
    std::cout << "      -query queryfile : Mandatory. Fasta or Fastq. Can be gzip'ed. Repeat this option for multiple query files\n";
    std::cout << "      -db \"<string> <string> ... <string>\" : Mandatory. The DB to query against. The same <string> as used for makedb\n";
    std::cout << "\n";

    std::cout << "   Scoring\n";
    std::cout << "      -min_score val : Specify the minimum score reported. Default val = 1\n";

    std::cout << "      -top or -topscore_num val : Output the val best scores. Default val = " << defaultoptions.numTopOutputs << "\n";
    std::cout << "      -gop or -gapo val : Gap open score. Overwrites our matrix-dependent default score.\n";
    std::cout << "      -gex or -gape val : Gap extend score. Overwrites our matrix-dependent default score.\n";
    #ifdef CAN_USE_FULL_BLOSUM
    std::cout << "      -mat val: Set substitution matrix. Supported values: dna, nuc44, pam30, pam70, pam30_20, pam70_20,\n";
    std::cout << "                 blosum45, blosum50, blosum62, blosum80, blosum45_20, blosum50_20, blosum62_20, blosum80_20.\n";
    std::cout << "                 Default: " << "blosum62_20" << "\n";
    #else 
    std::cout << "      -mat val: Set substitution matrix. Supported values: dna, nuc44, pam30, pam70, blosum45, blosum50, blosum62, blosum80. "
                        "Default: " << "blosum62" << "\n";
    #endif
    std::cout << "      -reverse	: to calculate scores and alignments also for reverse complement\n";

    std::cout << "\n";

    std::cout << "   Output options\n";
    std::cout << "      -out : Result output file. Parent directory must exist. Default: console output (/dev/stdout)\n";
    std::cout << "      -outfmt \"<string> <string> ... <string>\" 	: Specify the output file columns\n";
    std::cout << "       supported column names: qacc, qlen, sacc, slen, score, length, nident, gaps, qstart, qend, sstart,\n"; 
    std::cout << "                               send, positive, btop, topline, middleline, bottomline, reversed\n";
    std::cout << "      -progress_key	: Search identifier used in search progress notifications\n";
    std::cout << "      -progress_pipe	: Progress pipe file path\n";
    std::cout << "\n";

    std::cout << "   Misc\n";
    std::cout << "      -subject_length_min val : Specify the minimum subject sequence length. Default val = 1\n";
    std::cout << "      -subject_length_max val : Specify the maximum subject sequence length. Default val = 1,000,000\n";

    std::cout << "      -verbose : More console output. Shows timings. \n";
    std::cout << "      -version : Print program version\n";
    std::cout << "      -help : Print this message\n";
    std::cout << "\n";

    std::cout << "   Performance and benchmarking\n";
    std::cout << "      -use_single <integer>	: force to use the single GPU with ID #integer\n";
    std::cout << "\n";

}

    // }
    // std::cout << "memory limit per gpu: " << (options.maxGpuMem == std::numeric_limits<size_t>::max() ? 
    //     "unlimited" : std::to_string(options.maxGpuMem)) << "\n"; 

    // std::cout << "singlePassType: " << to_string(options.singlePassType) << "\n";
    // std::cout << "manyPassType_small: " << to_string(options.manyPassType_small) << "\n";
    // std::cout << "manyPassType_large: " << to_string(options.manyPassType_large) << "\n";
    // std::cout << "overflowType: " << to_string(options.overflowType) << "\n";
    // if(options.usePseudoDB){
    //     std::cout << "Using built-in pseudo db with " << options.pseudoDBSize << " sequences of length " << options.pseudoDBLength << "\n";
    // }else{

    // std::cout << "maxBatchBytes: " << options.maxBatchBytes << "\n";
    // std::cout << "maxBatchSequences: " << options.maxBatchSequences << "\n";
    // std::cout << "maxTempBytes: " << options.maxTempBytes << "\n";

    // std::cout << "loadFullDBToGpu: " << options.loadFullDBToGpu << "\n";
    // std::cout << "prefetchDBFile: " << options.prefetchDBFile << "\n";

    // std::cout << "   Memory\n";
    // std::cout << "      --maxGpuMem val : Try not to use more than val bytes of gpu memory per gpu. Uses all available gpu memory by default\n";
    // std::cout << "      --maxTempBytes val : Size of temp storage in GPU memory. Can use suffix K,M,G. Default val = " << defaultoptions.maxTempBytes << "\n";
    // std::cout << "      --maxBatchBytes val : Process DB in batches of at most val bytes. Can use suffix K,M,G. Default val = " << defaultoptions.maxBatchBytes << "\n";
    // std::cout << "      --maxBatchSequences val : Process DB in batches of at most val sequences. Default val = " << defaultoptions.maxBatchSequences << "\n";
    // std::cout << "\n";

//    std::cout << "      -dpx : Use DPX instructions. Hardware support requires Hopper (sm_90) or newer. Older GPUs fall back to software emulation.\n";
//    std::cout << "      -tsv : Print results as tab-separated values instead of plain text. \n";
    // std::cout << "      --printLengthPartitions : Print number of sequences per length partition in db.\n";
    // std::cout << "      --prefetchDBFile : Load DB into RAM immediately at program start instead of waiting for the first access.\n";
    // std::cout << "      --uploadFull : If enough GPU memory is available to store full db, copy full DB to GPU before processing queries.\n";
    // std::cout << "      --pseudodb num length : Use a generated DB which contains `num` equal sequences of length `length`.\n";
    // std::cout << "      --singlePassType val, --manyPassType_small val, --manyPassType_large val, --overflowType val :\n";
    // std::cout << "           Select kernel types for different length partitions. Valid values: Half2, DPXs16, DPXs32, Float.\n";
    // std::cout << "           Misc option --dpx is equivalent to --singlePassType DPXs16 --manyPassType_small DPXs16 --manyPassType_large DPXs32 --overflowType DPXs32.\n";
    // std::cout << "           Default is --singlePassType Half2 --manyPassType_small Half2 --manyPassType_large Float --overflowType Float.\n";

    // fprintf(stderr,"Usage:\n");
    // fprintf(stderr,"./cudasw [options]\n");
    // fprintf(stderr,"Standard options:\n");
    // fprintf(stderr,"\t-mat <string>\t: specify the substitution matrix name (default blosum62)\n");
    // fprintf(stderr,"\t\tsupported matrix names: dna, pam30, pam70, blosum45, blosum50, blosum62 and blosum80\n");
    // fprintf(stderr,"\t-query <string>\t: specify the query sequence file (fasta format)\n");
    // fprintf(stderr,"\t-db <string>\t: specify the database sequence file (blastdb format)\n");
    // fprintf(stderr,"\t-gapo <integer>\t: specify the gap open panelty (0 ~ 255), (default %d)\n", DEFAULT_GAPO);
    // fprintf(stderr,"\t-gape <integer>\t: specify the gap extension panelty (0 ~ 255), (default %d)\n", DEFAULT_GAPE);
    // fprintf(stderr,"\t-min_score <integer>\t: specify the minimum score reported(default %d)\n", DEFAULT_MIN_SCORE);
    // fprintf(stderr,"\t-topscore_num <integer>\t: specify the number of top scores reported(default %d)\n", DEFAULT_TOPSCORE_NUM);
    // fprintf(stderr,"\t-subject_length_min <integer>\t: specify the minimum subject sequence length(default %d)\n", 0);
    // fprintf(stderr,"\t-subject_length_max <integer>\t: specify the maximum subject sequence length(default %d)\n", MAX_SEQUENCE_LENGTH);
    // fprintf(stderr,"\t-use_single <integer>\t: force to use the single GPU with ID #integer\n");
    // fprintf(stderr,"\t-out <string>\t: specify the output file, (default out.csv)\n");
    // fprintf(stderr,"\t-reverse\t: to calculate scores and alignments also for reverse complement\n");
    // fprintf(stderr,"\t-outfmt \"<string> <string> ... <string> \t: specify the output file columns\n");
    // fprintf(stderr,"\t\t supported column names: qacc, qlen, sacc, slen, score, length, nident, gaps, qstart, qend, sstart, send, positive, btop, topline, middleline, bottomline, reversed\n");
    // fprintf(stderr,"\t-progress_key\t: search identifier used in search progress notifications\n");
    // fprintf(stderr,"\t-progress_pipe\t: progress pipe file path\n");
    // fprintf(stderr,"\t-version\t: print out the version\n");
	
