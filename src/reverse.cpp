#include "reverse.hpp"
#include <iostream>
#include <algorithm>
#include <string>

// Function to get the complement of a nucleotide
char getComplement(char nucleotide) {
    switch (nucleotide) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'U': return 'A'; // In RNA, U pairs with A
        case 'R': return 'Y'; // Purine to Pyrimidine
        case 'Y': return 'R'; // Pyrimidine to Purine
        case 'M': return 'K'; // Amino to Keto
        case 'W': return 'W'; // Weak to Weak
        case 'S': return 'S'; // Strong to Strong
        case 'K': return 'M'; // Keto to Amino
        case 'D': return 'H'; // Not C to Not G
        case 'H': return 'D'; // Not G to Not C
        case 'V': return 'B'; // Not T (U) to Not A
        case 'B': return 'V'; // Not A to Not T (U)
        case 'N': // Any Nucleotide
        case 'X': // Any Nucleotide
            return 'N';
        default:  return nucleotide; // For non-standard characters
    }
}

// Function to get the reverse complement of a DNA sequence
std::string getReverseComplement(const std::string& sequence) {
    std::string reverseComplement = sequence;

    // Reverse the string
    std::reverse(reverseComplement.begin(), reverseComplement.end());

    // Replace each nucleotide with its complement
    std::transform(reverseComplement.begin(), reverseComplement.end(), reverseComplement.begin(), getComplement);

    return reverseComplement;
}
