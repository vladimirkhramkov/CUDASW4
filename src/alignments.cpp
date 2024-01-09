#include <algorithm>
#include <iostream>
#include <thread>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "options.hpp"
#include "two_seq_aligner_optimized.h"

namespace cudasw2
{
    using namespace std;

    #define MAX_AMINO_ACIDS     	23

    char amino_acids[MAX_AMINO_ACIDS] = { 'A', 'B', 'C', 'D', 'E',
		'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
		'W', 'X', 'Y', 'Z' };

    int amino_acids_trans[256] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 10, 11, 12, 13, 0,
		14, 15, 16, 17, 18, 0, 19, 20, 21, 22, 23, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4,
		5, 6, 7, 8, 9, 0, 10, 11, 12, 13, 0, 14, 15, 16, 17, 18, 0, 19, 20, 21,
		22, 23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    char amino_acids_dna[MAX_AMINO_ACIDS] = { 'A',  'C',  'G',  'T',  'U',  'R',  'Y',  'M',  'W',  'S',  'K',  'D',  'H',  'V',  'B',  'N',  'X' };

    int amino_acids_trans_dna[256] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 15, 2, 12, 0, 0, 3, 13, 0, 0, 11, 0, 8, 16, 0, 0, 0, 6, 10, 4, 5, 14, 9, 17, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, };

    struct ResSeq {
        int name;
        int score;
        bool operator<(const ResSeq& o) const {
            if (score == o.score) return name > o.name;
            return score > o.score;
        }
    };

    struct AnswerEntry {
        int score;
        int queryIndex;
        int queryCoverageLength;
        bool reversed;
        string sequence;
        string sequenceName;
        vector<pair<int, int> > alignment;
        string getSacc() {
            return sequenceName;
        }
        int getSlen() {
            return sequence.length();
        }
        int getScore() {
            return score;
        }
        int getQstart() {
            return alignment[0].first;
        }
        int getQend() {
            return alignment.back().first;
        }
        int getSstart() {
            return alignment[0].second;
        }
        int getSend() {
            return alignment.back().second;
        }
    };

    struct QueryAnswer {
        string queryName;
        vector<AnswerEntry> entries;
    };

    struct AlignmentJob {
        string query;
        string sequence;
        int queryIndex;
        int entryIndex;
        AlignmentJob(const string& query, const string& sequence, int queryIndex,
                int entryIndex) :
                query(query), sequence(sequence), queryIndex(queryIndex), entryIndex(
                        entryIndex) {
        }
    };

    struct AlignmentThreadParameters {
        vector<AlignmentJob>* alignmentJobs;
        vector<cudasw4::HitResult>* totalResults;
        int numThreads;
        int remainder;
        int gape;
        int gapo;
        int scoringMatrix[32][32];
        int acid_trans[256];

        AlignmentThreadParameters() {
        }

        AlignmentThreadParameters(
            vector<AlignmentJob>* alignmentJobs,
            vector<cudasw4::HitResult>* totalResults,
            int numThreads, 
            int remainder,
            int gapo, 
            int gape, 
            int scoringMatrix[32][32], 
            int acid_trans[256]
        ) :
            alignmentJobs(alignmentJobs), 
            totalResults(totalResults), 
            numThreads(numThreads), 
            remainder(remainder), 
            gapo(gapo), 
            gape(gape) 
        {
            for (int i = 0; i < 32; ++i)
                for (int j = 0; j < 32; ++j)
                    this->scoringMatrix[i][j] = scoringMatrix[i][j];

            for (int i = 0; i < 256; ++i)
                this->acid_trans[i] = acid_trans[i];
        }
    };

    void* alignmentThreadFunction(void* ptr) {
        AlignmentThreadParameters parameters = *((AlignmentThreadParameters*) ptr);
        TwoSeqAlignerOptimized twoSeqAligner;

        for (int jobIndex = parameters.remainder; jobIndex < parameters.alignmentJobs->size(); jobIndex += parameters.numThreads) {
            const AlignmentJob& alignmentJob = (*parameters.alignmentJobs)[jobIndex];
            
            int alignmentScore;
            
            vector<pair<int, int> > alignment = twoSeqAligner.calcAlignment(
                alignmentJob.query, 
                alignmentJob.sequence, 
                parameters.gapo,
                parameters.gape, 
                parameters.acid_trans,
                parameters.scoringMatrix,
                alignmentScore
            );

            (*parameters.totalResults)[alignmentJob.entryIndex].alignment = alignment;

            // int previousScore = (*parameters.totalResults)[alignmentJob.entryIndex].score;

            // 32767 - overflowed score, calculated on cpu, 65535 - overflowed score calculated on gpu
            // if (previousScore != 32767 && previousScore != 65535 && alignmentScore != previousScore) {
            //     fprintf(stderr, "\n\nERROR: scores differ (job index: %d, %d != %d)\n\n", jobIndex, previousScore, alignmentScore);
            // }

            (*parameters.totalResults)[alignmentJob.entryIndex].score = alignmentScore;
        }

        return 0;
    }	

    void calculateAlignments(
        ostream& os, 
        vector<cudasw4::HitResult>& totalResults, 
        const ProgramOptions& options,
        const char* hostSubMatrix, 
        int hostSubMatrixDim
    ) {
        // options.numTopOutputs,
        // options.sequenceType,
        // options.subMatrixType, 
        // options.verbose,
        // options.progressKey,

        bool nucleotideSearch = (options.sequenceType == cudasw4::SequenceType::Nucleotide);

        const int n = min((int)totalResults.size(), options.numTopOutputs);

        int scoringMatrix[32][32] = {0};

        for (int i = 0; i < hostSubMatrixDim; ++i) {
            for (int j = 0; j < hostSubMatrixDim; ++j) {
                scoringMatrix[i][j] = static_cast<int>(hostSubMatrix[i*hostSubMatrixDim + j]);
            }
        }

       	fprintf(stderr, "prepairing answers\n");

        vector<AlignmentJob> alignmentJobs;
        // TwoSeqAlignerOptimized twoSeqAligner;

        for(int i = 0; i < n; i++){
            //     os << totalResults[i].queryHeader || totalResults[i].queryId;
            // } else if (parameter == "qlen") {
            //     os << totalResults[i].queryLength;
            // } else if (parameter == "sacc") {
            //     os << totalResults[i].subjectHeader;
            // } else if (parameter == "slen") {
            //     os << totalResults[i].subjectLength;
            // } else if (parameter == "score") {
            //     os << totalResults[i].score;


            // struct AlignmentJob {
            //     string query;
            //     string sequence;
            //     int queryIndex;
            //     int entryIndex;
            // };

            alignmentJobs.push_back (
                AlignmentJob (
                    totalResults[i].querySequence, 
                    totalResults[i].subjectSequence, 
                    totalResults[i].queryId, 
                    i
                )
            );
        }

        fprintf(stderr, "Calculation of alignments\n");

        // Get the number of supported concurrent threads
        unsigned int num_threads = std::thread::hardware_concurrency();

        vector<pthread_t> alignmentThreads(num_threads);

        vector<AlignmentThreadParameters> threadParameters(alignmentThreads.size());

        for (int threadIndex = 0; threadIndex < alignmentThreads.size(); ++threadIndex) {

            threadParameters[threadIndex] = AlignmentThreadParameters(
                &alignmentJobs, 
                &totalResults, 
                alignmentThreads.size(),
                threadIndex, 
                abs(options.gop), 
                abs(options.gex),
                scoringMatrix, 
                amino_acids_trans
            );
            
            pthread_create(&alignmentThreads[threadIndex], 0, alignmentThreadFunction, &threadParameters[threadIndex]);
        }

        for (int threadIndex = 0; threadIndex < alignmentThreads.size(); ++threadIndex) {
            pthread_join(alignmentThreads[threadIndex], NULL);
        }

        for(int i = 0; i < n; i++) {
            int coverageLength = totalResults[i].alignment.size();

            vector<string> stringAlignment = TwoSeqAlignerOptimized::getStringResult(
                totalResults[i].querySequence,
                totalResults[i].subjectSequence,
                totalResults[i].alignment,
                options.isDnaDatabase(),
                amino_acids_trans, 
                scoringMatrix
            );

            // adding pipe for T<->U matches
            if (nucleotideSearch)
                for (int i = 0; i < stringAlignment[0].length(); ++i)
                    if (stringAlignment[0][i] == 'T' && stringAlignment[2][i] == 'U' || stringAlignment[0][i] == 'U' && stringAlignment[2][i] == 'T')
                        stringAlignment[1][i] = '|';

            string BTOP = TwoSeqAlignerOptimized::getBTOPResult(
                totalResults[i].querySequence,
                totalResults[i].subjectSequence, 
                totalResults[i].alignment,
                options.isDnaDatabase(), 
                amino_acids_trans,
                scoringMatrix
            );

            for (int paramIndex = 0; paramIndex < options.csvColumns.size(); ++paramIndex) {
                string parameter = options.csvColumns[paramIndex];
                if (paramIndex)
                    os << ",";
                if (parameter == "qacc")
                    os << totalResults[i].queryHeader || totalResults[i].queryId;
                else if (parameter == "qlen")
                    os << totalResults[i].queryLength;
                else if (parameter == "sacc")
                    os << totalResults[i].subjectHeader;
                else if (parameter == "slen")
                    os << totalResults[i].subjectLength;
                else if (parameter == "score")
                    os << totalResults[i].score;
                else if (parameter == "length")
                    os << stringAlignment[0].length();
                else if (parameter == "nident") {
                    int matches = 0;
                    for (int i = 0; i < stringAlignment[0].length(); ++i)
                        if (stringAlignment[0][i] != '-' && (stringAlignment[0][i] == stringAlignment[2][i] || stringAlignment[1][i] == '|' && nucleotideSearch))
                            ++matches;
                    os << matches;
                } else if (parameter == "gaps") {
                    int gaps = 0;
                    for (int i = 0; i < stringAlignment[0].length(); ++i)
                        gaps += (stringAlignment[0][i] == '-')
                                + (stringAlignment[2][i] == '-');
                    os << gaps;
                } else if (parameter == "qstart")
                    os << totalResults[i].getQstart() + 1;
                else if (parameter == "qend")
                    os << totalResults[i].getQend() + 1;
                else if (parameter == "sstart")
                    os << totalResults[i].getSstart() + 1;
                else if (parameter == "send")
                    os << totalResults[i].getSend() + 1;
                else if (parameter == "positive") {
                    int positiveMatches = 0;
                    for (int i = 0; i < stringAlignment[0].length(); ++i) {
                        if (stringAlignment[0][i] != '-' && stringAlignment[2][i] != '-') {
                            positiveMatches += scoringMatrix[amino_acids_trans[stringAlignment[0][i]]][amino_acids_trans[stringAlignment[2][i]]] > 0;
                        }
                    }
                    os << positiveMatches;
                } else if (parameter == "btop") {
                    os << BTOP;
                } else if (parameter == "topline")
                    os << stringAlignment[0];
                else if (parameter == "middleline")
                    os << stringAlignment[1];
                else if (parameter == "bottomline")
                    os << stringAlignment[2];
                else if (parameter == "reversed")
                    os << totalResults[i].reversed;
                else if (parameter == "qcovs") {
                    os << fixed << setprecision(2) << 100.0 * coverageLength / totalResults[i].queryLength;
                }
            }
            os << "\n";
        }
    }

    void printScanResultCSV(
        std::ostream& os, 
        const std::vector<cudasw4::HitResult>& totalResults, 
        const ProgramOptions& options
    ) {
        constexpr char sep = ',';

        const int n = std::min((int)totalResults.size(), options.numTopOutputs);

        for(int i = 0; i < n; i++){
            for (unsigned int paramIndex = 0; paramIndex < options.csvColumns.size(); ++paramIndex) {
                std::string parameter = options.csvColumns[paramIndex];
                if (paramIndex)
                    os << sep;

                if (parameter == "qacc") {
                    os << totalResults[i].queryHeader || totalResults[i].queryId;
                } else if (parameter == "qlen") {
                    os << totalResults[i].queryLength;
                } else if (parameter == "sacc") {
                    os << totalResults[i].subjectHeader;
                } else if (parameter == "slen") {
                    os << totalResults[i].subjectLength;
                } else if (parameter == "score") {
                    os << totalResults[i].score;
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

} // namespace cudasw2
