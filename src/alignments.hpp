#ifndef ALIGNMENTS_HPP
#define ALIGNMENTS_HPP

#include <vector>
#include "types.hpp"
#include "options.hpp"

namespace cudasw2 {

    using namespace std;

    void calculateAlignments(
        ostream& os, 
        vector<cudasw4::HitResult>& totalResults, 
        const ProgramOptions& options,
        const char* hostSubMatrix, 
        int hostSubMatrixDim
    );

    void printScanResultCSV(
        std::ostream& os, 
        const std::vector<cudasw4::HitResult>& totalResults, 
        const ProgramOptions& options
    );
}

#endif
