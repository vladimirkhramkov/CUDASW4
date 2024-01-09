#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "two_seq_aligner_optimized.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>

namespace cudasw2
{
	using namespace std;

	vector<pair<int,int> > TwoSeqAlignerOptimized::calcAlignment(const string & s1, const string & s2, int gapo, int gape, const int amino_acids_trans[256], const int matr[32][32], int& score) {
		for (int i = 0; i < 256; ++i) trans[i] = amino_acids_trans[i];
		for (int i = 0; i < 32; ++i) for (int j = 0; j < 32; ++j) this->matr[i][j] = matr[i][j];
		this->s1 = s1, this->s2 = s2;
		this->gapo = gapo;
		this->gape = gape;

		int best = 0 ;
		int bestI = 0;
		int bestJ = 0;

		vector<int> tr1(s1.length()), tr2(s2.length());
		for (int i = 0; i < s1.length(); ++i) tr1[i] = trans[s1[i]];
		for (int i = 0; i < s2.length(); ++i) tr2[i] = trans[s2[i]];
		vector<int> vH(s2.length() + 1);
		vector<int> vF(s2.length() + 1);

		for (int i = 1; i <= s1.length(); ++i) {
			int prevH = 0;
			int prevE = 0;
			const int* mm = matr[tr1[i-1]];
			for (int j = 1; j <= s2.length(); ++j) {
				if (prevE < vH[j-1] - gapo) prevE = vH[j-1] - gapo;
				prevE -= gape;

				if (vF[j] < vH[j] - gapo) vF[j] = vH[j] - gapo;
				vF[j] -= gape;

				prevH += mm[tr2[j-1]];
				if (prevH < prevE) prevH = prevE;
				if (prevH < vF[j]) prevH = vF[j];
				if (prevH < 0) prevH = 0;

				if (prevH > best) {
					best = prevH;
					bestI = i-1;
					bestJ = j-1;
				}
				swap(prevH, vH[j]);
			}
		}

		score = best;
		SWCell resultCell;
		vector<pair<int, int> > answer = align(0, bestI, 0, bestJ, SWCell(), 'H', &resultCell);

		sort(answer.begin(), answer.end());
		answer.erase(unique(answer.begin(), answer.end()), answer.end());

		return answer;
	}

	vector<pair<int,int> > TwoSeqAlignerOptimized::align(int fromI, int toI, int fromJ, int toJ, const SWCell& initialValue, char valueIndex, SWCell* resultCell) {
		if (fromI > toI || fromJ > toJ) return vector<pair<int, int> >();

		fromI = max(fromI, 0);
		fromJ = max(fromJ, 0);

		int height = toI - fromI + 1;
		int width = toJ - fromJ + 1;

		if (height < 3  || (height + 1) * (long long)(width + 1) < TWO_SEQ_ALIGNER_THRESHOLD) {
			vector<int> tr1(height + 1), tr2(width + 1);
			for (int i = 0; i <= height; ++i) tr1[i] = trans[s1[i + fromI]];
			for (int i = 0; i <= width; ++i) tr2[i] = trans[s2[i + fromJ]];

			vector<int> vH(width + 1);
			vector<int> vF(width + 1);

			vector<vector<unsigned int> > rE(height + 1, vector<unsigned int>(1+(width + 1) / 32));
			vector<vector<unsigned int> > rF(height + 1, vector<unsigned int>(1+(width + 1) / 32));
			vector<vector<unsigned int> > rH(height + 1, vector<unsigned int>(1+(width + 1) / 16));

			int prevE = 0;

			for (int i = 1; i <= height; ++i) {
				int prevH = 0;
				prevE = 0;
				int b12 = 0;
				int* mm = matr[tr1[i - 1]];
				unsigned int bE = 0;
				unsigned int bF = 0;
				unsigned int bH = 0;
				for (int j = 1; j <= width; ++j) {
					if ((fromI || fromJ) && i == 1 && j == 1) {
						prevE = initialValue.E;
						vF[j] = initialValue.F;
						vH[j] = initialValue.H;
						b12 = 0;
						if (initialValue.H <= 0) b12 = 3;
						else if (initialValue.H == initialValue.E) b12 = 1;
						else if (initialValue.H == initialValue.F) b12 = 2;
						bH |= b12 << ((j&15)<<1);
					} else {
						if (prevE < vH[j-1] - gapo) {
							prevE = vH[j-1] - gapo;
							bE |= 1u << (j&31);
						}
						prevE -= gape;

						if (vF[j] < vH[j] - gapo) {
							vF[j] = vH[j] - gapo;
							bF |= 1u << (j&31);
						}
						vF[j] -= gape;

						prevH += mm[tr2[j - 1]];
						b12 = 0;

						if (prevH <= vF[j]) {
							prevH = vF[j];
							b12 = 2;
						}
						if (prevH <= prevE) {
							prevH = prevE;
							b12 = 1;
						}

						if (prevH <= 0) {
							prevH = 0;
							b12 = 3;
						}
						bH |= b12 << ((j&15)<<1);
						swap(prevH, vH[j]);
					}
					if ((j&15) == 15 || j == width) {
						rH[i][j>>4] = bH;
						bH = 0;
						if ((j&31) == 31 || j == width) {
							rE[i][j>>5] = bE;
							rF[i][j>>5] = bF;
							bE = bF = 0;
						}
					}
				}
			}

			resultCell->E = prevE;
			resultCell->F = vF[width];
			resultCell->H = vH[width];
			int curI = height;
			int curJ = width;
			char curChar = valueIndex;
			vector<pair<int, int> > answer;
			
			while (curI && curJ) {
				if (curChar == 'E') {
					if (!((rE[curI][curJ>>5]>>(curJ&31))&1)) {
						curJ--;
					} else {
						curJ--;
						curChar = 'H';
					}
				} else if (curChar == 'F') {
					if (!((rF[curI][curJ>>5]>>(curJ&31))&1)) {
						curI--;
					} else {
						curI--;
						curChar = 'H';
					}
				} else if (curChar == 'H') {
					int val = (rH[curI][curJ>>4]>>((curJ&15)<<1))&3;
					if (val == 3) {
						break;
					} else if (val == 1) {
						curChar = 'E';
					} else if (val == 2) {
						curChar = 'F';
					} else {
						answer.push_back(make_pair(fromI + curI - 1, fromJ + curJ - 1));
						curChar = 'H';
						curI--;
						curJ--;
					}
				}
			}

			return answer;
		} else {
			int middle = height / 2 + 1;

			vector<SWCell> dp(width + 1);
			vector<PrevCell> HPrev(width + 1);
			vector<PrevCell> EPrev(width + 1);
			vector<PrevCell> FPrev(width + 1);

			for (int i = 0; i <= width; ++i) {
				HPrev[i] = PrevCell(middle, i, 'H');
				EPrev[i] = PrevCell(middle, i, 'E');
				FPrev[i] = PrevCell(middle, i, 'F');
			}

			for (int i = 1; i <= height; ++i) {
				char c1 = trans[s1[i-1+fromI]];
				int prevH = dp[0].H;
				PrevCell prevHCell = HPrev[0];
				for (int j = 1; j <= width; ++j) {
					if ((fromI || fromJ) && i == 1 && j == 1) {
						dp[j] = initialValue;
						continue;
					}
					char c2 = trans[s2[j-1+fromJ]];
					dp[j].E = dp[j-1].E - gape;
					if (i > middle) EPrev[j] = EPrev[j-1];

					if (dp[j].E < dp[j-1].H - gapo - gape) {
						dp[j].E = dp[j-1].H - gapo - gape;
						if (i > middle) EPrev[j] = HPrev[j-1];
					}

					dp[j].F = dp[j].F - gape;

					if (dp[j].F < dp[j].H - gapo - gape) {
						dp[j].F = dp[j].H - gapo - gape;
						if (i > middle) FPrev[j] = HPrev[j];
					}

					int tmp = prevH + matr[c1][c2];
					PrevCell tmpPrev = prevHCell;
					prevH = dp[j].E, prevHCell = EPrev[j];
					if (prevH < dp[j].F) prevH = dp[j].F, prevHCell = FPrev[j];
					if (prevH < tmp) prevH = tmp, prevHCell = tmpPrev;
					if (prevH <= 0) {
						prevH = 0;
						if (i > middle) prevHCell = PrevCell(middle, 0, 'H');
					}
					if (i <= middle) prevHCell = HPrev[j];
					swap(prevH, dp[j].H);
					swap(prevHCell, HPrev[j]);
				}
			}


			PrevCell middleCell = HPrev[width];
			
			if (valueIndex == 'F') middleCell = FPrev[width];
			else if (valueIndex == 'E') middleCell = EPrev[width];

			SWCell middleResult;
			vector<pair<int, int> > firstResultPart = align(fromI, fromI + middleCell.i-1, fromJ, fromJ + middleCell.j-1, initialValue, middleCell.valueIndex, &middleResult);
			vector<pair<int, int> > lastResultPart = align(fromI + middleCell.i-1, toI, fromJ + middleCell.j-1, toJ, middleResult, valueIndex, resultCell);
			vector<pair<int, int> > totalResult = firstResultPart;
			
			totalResult.insert(totalResult.end(), lastResultPart.begin(), lastResultPart.end());
			
			return totalResult;
		}
	}

	vector<string> TwoSeqAlignerOptimized::getStringResult(const string & s1, const string & s2, const vector<pair<int,int> > & result, bool dna, int acid_trans[256], int scoringMatrix[32][32]) {
		string r1 = "";
		string r2 = "";
		string r3 = "";
		int mptr = 0;
		if (result.size()) {
			int from1 = result[0].first;
			int from2 = result[0].second;

			int ptr1 = from1;
			int ptr2 = from2;
			while (mptr < result.size()) {
				if (ptr1 == result[mptr].first && ptr2 == result[mptr].second) {
					if (s1[ptr1] == s2[ptr2]) {
						if (dna)
							r3 += "|";
						else
							r3 += s1[ptr1];
					} else {
						if (dna)
							r3 += " ";
						else {
							if (scoringMatrix[acid_trans[s1[ptr1]]][acid_trans[s2[ptr2]]] > 0) {
								r3 += "+";
							} else {
								r3 += " ";
							}
						}
					}
					r1 += s1[ptr1++];
					r2 += s2[ptr2++];

					++mptr;
				} else if (ptr1 < result[mptr].first) {
					r1 += s1[ptr1++];
					r2 += "-";
					r3 += " ";
				} else if (ptr2 < result[mptr].second) {
					r1 += "-";
					r2 += s2[ptr2++];
					r3 += " ";
				}
			}
		}
		vector<string> res;
		res.push_back(r1);
		res.push_back(r3);
		res.push_back(r2);
		return res;
	}

	string intToString(int x) {
		string res = "";
		do {
			res += (char)('0' + x % 10);
			x /= 10;
		} while (x);
		reverse(res.begin(), res.end());
		return res;
	}

	string TwoSeqAlignerOptimized::getBTOPResult(const string& s1, const string& s2,
			const vector<pair<int, int> >& result, bool dna, int acid_trans[256], int scoringMatrix[32][32]) {
		vector<string> stringVector = getStringResult(s1, s2, result, dna, acid_trans, scoringMatrix);
		string res = "";
		int ptr = 0;
		while (ptr < stringVector[0].length()) {
			if (stringVector[0][ptr] != stringVector[2][ptr]) res += stringVector[0][ptr], res += stringVector[2][ptr], ++ptr;
			else {
				int matchCount = 0;
				while (ptr < stringVector[0].length() && stringVector[0][ptr] == stringVector[2][ptr]) ++matchCount, ++ptr;
				res += intToString(matchCount);
			}
		}
		return res;
	}

} // namespace cudasw2
