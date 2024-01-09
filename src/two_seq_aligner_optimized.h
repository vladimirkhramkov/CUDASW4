#ifndef TWO_SEQ_ALIGNER_CPP
#define TWO_SEQ_ALIGNER_CPP

#include <iostream>
#include <vector>
#include <string>

namespace cudasw2
{
	using namespace std;

	const int TWO_SEQ_ALIGNER_THRESHOLD = 2000000000;

	struct SWCell {
		int E, F, H;
		SWCell(): E(0), F(0), H(0) {}
		SWCell(int E, int F, int H): E(E), F(F), H(H) {}
	};

	struct PrevCell {
		int i, j;
		char valueIndex;
		PrevCell():i(0), j(0), valueIndex('-') {}
		PrevCell(int i, int j, char valueIndex):i(i), j(j), valueIndex(valueIndex) {}
	};


	class TwoSeqAlignerOptimized {
	public:
		vector<pair<int, int> > calcAlignment(const string& s1, const string& s2, int gapo, int gape, const int amino_acids_trans[256], const int matr[32][32], int& score);
		static string getBTOPResult(const string& s1, const string& s2, const vector<pair<int, int> >& result, bool dna, int acid_trans[256], int scoringMatrix[32][32]);
		static vector<string> getStringResult(const string & s1, const string & s2, const vector<pair<int,int> >& result, bool dna, int acid_trans[256], int scoringMatrix[32][32]);
	private:
		vector<pair<int, int> > align(int fromI, int toI, int fromJ, int toJ, const SWCell& initialValue, char valueIndex, SWCell* resultCell);
		vector<pair<int, int> > matches;
		int gapo;
		int gape;
		int trans[256];
		int matr[32][32];
		string s1, s2;
	};

} // namespace cudasw2

#endif //TWO_SEQ_ALIGNER_CPP
