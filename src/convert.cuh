#ifndef CONVERT_CUH
#define CONVERT_CUH

namespace cudasw4{

struct ConvertNA{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& NA) {
        // ORDER of Nucleotides (NCBI): A  C  G  T  U  R  Y  M  W  S  K  D  H  V  B  N  X
        if (NA == 'A') return 0;
        if (NA == 'C') return 1;
        if (NA == 'G') return 2;
        if (NA == 'T') return 3;
        if (NA == 'U') return 4;
        if (NA == 'R') return 5;
        if (NA == 'Y') return 6;
        if (NA == 'M') return 7;
        if (NA == 'W') return 8;
        if (NA == 'S') return 9;
        if (NA == 'K') return 10;
        if (NA == 'D') return 11;
        if (NA == 'H') return 12;
        if (NA == 'V') return 13;
        if (NA == 'B') return 14;
        if (NA == 'N') return 15;
        if (NA == 'X') return 16;
        return 17; //  else
    }
};

struct ConvertAA{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& AA) {
        // ORDER of AminoAcids (NCBI): A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X
        if (AA == 'A') return 0;
        if (AA == 'R') return 1;
        if (AA == 'N') return 2;
        if (AA == 'D') return 3;
        if (AA == 'C') return 4;
        if (AA == 'Q') return 5;
        if (AA == 'E') return 6;
        if (AA == 'G') return 7;
        if (AA == 'H') return 8;
        if (AA == 'I') return 9;
        if (AA == 'L') return 10;
        if (AA == 'K') return 11;
        if (AA == 'M') return 12;
        if (AA == 'F') return 13;
        if (AA == 'P') return 14;
        if (AA == 'S') return 15;
        if (AA == 'T') return 16;
        if (AA == 'W') return 17;
        if (AA == 'Y') return 18;
        if (AA == 'V') return 19;
        if (AA == 'B') return 20;
        if (AA == 'J') return 21;
        if (AA == 'Z') return 22;
        if (AA == 'X') return 23;
        return 24; //  else
    }
};

struct ConvertAA_20{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& AA) {
        // ORDER of AminoAcids (NCBI): A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
        if (AA == 'A') return 0;
        if (AA == 'R') return 1;
        if (AA == 'N') return 2;
        if (AA == 'D') return 3;
        if (AA == 'C') return 4;
        if (AA == 'Q') return 5;
        if (AA == 'E') return 6;
        if (AA == 'G') return 7;
        if (AA == 'H') return 8;
        if (AA == 'I') return 9;
        if (AA == 'L') return 10;
        if (AA == 'K') return 11;
        if (AA == 'M') return 12;
        if (AA == 'F') return 13;
        if (AA == 'P') return 14;
        if (AA == 'S') return 15;
        if (AA == 'T') return 16;
        if (AA == 'W') return 17;
        if (AA == 'Y') return 18;
        if (AA == 'V') return 19;
        return 20; //  else
    }
};

struct InverseConvertNA{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& NA) {
        // ORDER of Nucleotides (NCBI): A  C  G  T  U  R  Y  M  W  S  K  D  H  V  B  N  X
        if (NA == 0) return 'A';
        if (NA == 1) return 'C';
        if (NA == 2) return 'G';
        if (NA == 3) return 'T';
        if (NA == 4) return 'U';
        if (NA == 5) return 'R';
        if (NA == 6) return 'Y';
        if (NA == 7) return 'M';
        if (NA == 8) return 'W';
        if (NA == 9) return 'S';
        if (NA == 10) return 'K';
        if (NA == 11) return 'D';
        if (NA == 12) return 'H';
        if (NA == 13) return 'V';
        if (NA == 14) return 'B';
        if (NA == 15) return 'N';
        if (NA == 16) return 'X';
        return '-'; //  else
    }
};

struct InverseConvertAA_20{
    #ifdef __CUDACC__
    __host__ __device__
    #endif
    char operator()(const char& AA) {
        // ORDER of AminoAcids (NCBI): A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
        if (AA == 0) return 'A';
        if (AA == 1) return 'R';
        if (AA == 2) return 'N';
        if (AA == 3) return 'D';
        if (AA == 4) return 'C';
        if (AA == 5) return 'Q';
        if (AA == 6) return 'E';
        if (AA == 7) return 'G';
        if (AA == 8) return 'H';
        if (AA == 9) return 'I';
        if (AA == 10) return 'L';
        if (AA == 11) return 'K';
        if (AA == 12) return 'M';
        if (AA == 13) return 'F';
        if (AA == 14) return 'P';
        if (AA == 15) return 'S';
        if (AA == 16) return 'T';
        if (AA == 17) return 'W';
        if (AA == 18) return 'Y';
        if (AA == 19) return 'V';
        return '-'; //  else
    }
};

} //namespace cudasw4

#endif