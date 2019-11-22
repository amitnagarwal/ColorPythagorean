#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <vector>
#include <assert.h>
#include <sstream>
#include <string>
#include "PythagoreanTriple.h"
#include "Color.h"
#include "FlipOption.h"
#include <set>

using namespace std;

long long maxGlobalElementIndex = -1;
long long maxUsedIndx = 0;

PytList allPTriples;
vector<pPytList> pTriplesOrig;
vector<pPytList> pTriplesMod;

void Color::equateNegatively(const Color& that) 
{
	_ndx = that._ndx; 
	_state = !that._state; 
	if (_id > maxUsedIndx)
	{
		maxUsedIndx = _id;
	}
}

int numPTriplesToBeFixed(FlipOption& fo, long long indx)
{
	vector<Color>& nColors = fo._nodeColors;

	long long numTriples = pTriplesOrig[indx].size();
	int numConflicts = 0;

	for (long long iTr = 0; iTr < numTriples; iTr++)
	{
		long long gTripleIndx = pTriplesOrig[indx][iTr];
		if (gTripleIndx >= maxGlobalElementIndex) break;

		if (gTripleIndx == 91)
		{
			int asa = 5;
		}

		const PythogoreanTriple& t = allPTriples[gTripleIndx];

		// cout << "      " << gTripleIndx << "] " << t.print(nColors) << endl;

		if (t.isColoredValidly(nColors)) continue;

		if (nColors[t._leg[0]]._ndx == -1 || nColors[t._leg[1]]._ndx == -1 || nColors[t._leg[2]]._ndx == -1)
		{
			t.handleSingleUnassigned(nColors);
			// cout << "      " << gTripleIndx << "} " << t.print(nColors) << endl;
		}
		else
		{
			int n1 = t.getLegIndex(indx);
			int n2 = (n1 + 1) % 3;
			int n3 = (n1 + 2) % 3;

			if (nColors[t._leg[n2]]._ndx != nColors[t._leg[n3]]._ndx)
			{
				mergeIndices(nColors, nColors[t._leg[n2]], nColors[t._leg[n3]]);
				// cout << "      " << gTripleIndx << "| " << t.print(nColors) << endl;
				continue;
			}

			fo._triplesToFix.push_back(gTripleIndx);
			numConflicts++;
			cout << " [" << t.print(nColors) << "]";
		}

	} // iTr

	return numConflicts;
}

int checkTriplesWithGivenNode(long long indx, vector<Color>& colors, bool recurse = false)
{
	cout << "Switchin " << indx << endl;
	if (recurse) cout << "              recursively." << endl;
	vector<Color> inColors = colors;

	long long numTriples = pTriplesOrig[indx].size();
	int numConflicts = 0;

	for (long long iTr = 0; iTr < numTriples; iTr++)
	{
		long long gTripleIndx = pTriplesOrig[indx][iTr];
		if (gTripleIndx >= maxGlobalElementIndex) break;

		const PythogoreanTriple& t = allPTriples[gTripleIndx];

		// cout << "      " << gTripleIndx << "] " << t.print(inColors) << endl;

		if (t.isColoredValidly(inColors)) continue;

		if (inColors[t._leg[0]]._ndx == -1 || inColors[t._leg[1]]._ndx == -1 || inColors[t._leg[2]]._ndx == -1)
		{
			t.handleSingleUnassigned(inColors);
			// cout << "      " << gTripleIndx << "} " << t.print(inColors) << endl;
		}
		else
		{
			int n1; 
			if (indx == t._leg[0])
				n1 = 0;
			else if (indx == t._leg[1])
				n1 = 1;
			else if (indx == t._leg[2])
				n1 = 2;
			int n2 = (n1 + 1) % 3;
			int n3 = (n1 + 2) % 3;

			if (inColors[t._leg[n2]]._ndx != inColors[t._leg[n3]]._ndx)
			{
				mergeIndices(inColors, inColors[t._leg[n2]], inColors[t._leg[n3]]);
				// cout << "      " << gTripleIndx << "| " << t.print(inColors) << endl;
				continue;
			}

			if (recurse)
				numConflicts += t.allColorsSameHandled(inColors);
			else
			{
				numConflicts++;
				cout << "TBF:     " << t.print(inColors) << endl;
			}

			if (numConflicts >= 1000) break;
		}

	} // iTr

	if (numConflicts == 0)
	{
		cout << "Switched " << indx << endl;
		colors = inColors;
	}
	else
	{
		cout << "     conflicts = " << numConflicts << endl;
	}

	return numConflicts;
}

void mergeIndices(vector<Color>& inpColors, Color& color1, Color& color2)
{
	long long discardIndx;
	long long keepIndx;
	if (color1._ndx < color2._ndx)
	{
		discardIndx = color2._ndx;	
		keepIndx = color1._ndx;
	}
	else
	{
		discardIndx = color1._ndx;
		keepIndx = color2._ndx;
	}
	bool equals = (color1._state != color2._state);
	inpColors[discardIndx]._ndx = keepIndx;
	inpColors[discardIndx]._state = !equals;

	for (long k = 0; k <= maxUsedIndx; k++)
	{
		if (inpColors[k]._ndx == discardIndx)
		{
			inpColors[k]._ndx = keepIndx;
			inpColors[k]._state = (inpColors[k]._state == equals);
		}
	} // k
}

//long long lim = 4230;
long long lim = 5200;
int lim10th = lim / 10;
void colorAPythagoreanTriple(long long i, vector<Color>& colors);

void findAllPythagoreanTriples()
{
	pTriplesOrig.resize(lim);
	long long ptCount = 0;
	//vector<Color> colors(lim);

	for (long long c = 5; c < lim; c++)
	{
		if (c % lim10th == 0) cout << c << endl;
		long long aStart = sqrt(2 * c - 1);
		int count = 0;
		for (long long a = aStart; a < c-1; a++)
		{
			long long bsq = c*c - a*a;
			double b = sqrt(bsq);
			long long lb = long(b);
			if (b == lb && lb > a)
			{
				//            cout << b << " , " << a << " , " << c << endl;
				PythogoreanTriple t;
				t._leg[0] = a;
				t._leg[1] = lb;
				t._leg[2] = c;
				//cout << t.print() << endl;

				//if (!t.isPrimitive()) continue;

				count++;
				allPTriples.push_back(t);
				pTriplesOrig[t._leg[0]].push_back(ptCount);
				pTriplesOrig[t._leg[1]].push_back(ptCount);
				pTriplesOrig[t._leg[2]].push_back(ptCount);
				//colorAPythagoreanTriple(ptCount, colors);
				ptCount++;
			}
		} // a
		if (count >= 1000)
		{
			int cc = pTriplesOrig[c].size();
			for (int j = 0; j < cc; j++)
			{
				int jj = pTriplesOrig[c][j];
				const PythogoreanTriple& t = allPTriples[jj];
				if (t._leg[2] == c)
				{
					cout << "     " << t._leg[0] << ", " << t._leg[1] << ", " << t._leg[2] << endl;
				}
			}
			cout << endl;
		}
	} // c

	pTriplesMod = pTriplesOrig;

//	reverse(allPTriples.begin(), allPTriples.end());
	/*
	vector<PythogoreanTriple>::iterator ptIter;
	for (ptIter = allPTriples.begin(); ptIter != allPTriples.end(); ptIter++)
	{
		pTriples[ptIter->_leg[0]].push_back(ptCount);
		pTriples[ptIter->_leg[1]].push_back(ptCount);
		pTriples[ptIter->_leg[2]].push_back(ptCount);
		ptCount++;
	} // ptIter
	*/
}

void printAllTriplesToFile()
{
	ofstream outFile("PythagoreanTriples.txt");
	for (long i = 0; i < allPTriples.size(); i++)
	{
		const PythogoreanTriple& t = allPTriples[i];
		int total = pTriplesOrig[t._leg[0]].size() + pTriplesOrig[t._leg[1]].size() + pTriplesOrig[t._leg[2]].size();
		outFile << "     " << t._leg[0] << ", " << t._leg[1] << ", " << t._leg[2] << " - " << total << "  ";

		if (pTriplesOrig[t._leg[0]].size() == 1 && pTriplesOrig[t._leg[1]].size() == 1 && pTriplesOrig[t._leg[2]].size() == 1)
		{
			outFile << " *** ";
		}
		outFile << endl;
	}

	for (long i = 0; i < lim; i++)
	{
		if (pTriplesOrig[i].size() > 0)
		{
			outFile << i << "   " << pTriplesOrig[i].size() << endl;
			for (unsigned long j = 0; j < pTriplesOrig[i].size(); j++)
			{
				const PythogoreanTriple& t = allPTriples[pTriplesOrig[i][j]];
				outFile << "     " << t._leg[0] << ", " << t._leg[1] << ", " << t._leg[2] << endl;
			} // j
		}
	}
}

void clearPreviousIndices();

FlipOption::FlipOption(int c, long long inNode1, long long inNode2)
{
	_cnf = c;
	_fN1 = inNode1;
	_fN2 = inNode2;
}

bool FlipOption::operator<(const FlipOption& that) const 
{
	if (_cnf == that._cnf)
	{
		return (_triplesToFix.size() < that._triplesToFix.size());
	}
	return (_cnf < that._cnf);
}

bool isFlippable(const set<long long>& inNodes, long long nodeIndx)
{
	set<long long>::iterator kIter = find(inNodes.begin(), inNodes.end(), nodeIndx);
	return (kIter != inNodes.end());
}

list<FlipOption> analyzePTriple(
	const PythogoreanTriple& pt,   // In
	const vector<Color>& inColors, // In
	set<long long>& fixedNodes   // Modified
)
{
	cout << "************************************************" << endl;
	cout << "TBF1:     " << pt.print(inColors) << endl;
	list<FlipOption> flipOptions;

	bool flipAllowed[3];
	
	flipAllowed[0] = !isFlippable(fixedNodes, pt._leg[0]);
	flipAllowed[1] = !isFlippable(fixedNodes, pt._leg[1]);
	flipAllowed[2] = !isFlippable(fixedNodes, pt._leg[2]);

	if (!flipAllowed[0] && !flipAllowed[1] && !flipAllowed[2])
	{
		return flipOptions;
	}

	fixedNodes.insert(pt._leg[0]);
	fixedNodes.insert(pt._leg[1]);
	fixedNodes.insert(pt._leg[2]);

	if (flipAllowed[0]) flipOptions.push_back(FlipOption(pTriplesMod[pt._leg[0]].size(), pt._leg[0]));
	if (flipAllowed[1]) flipOptions.push_back(FlipOption(pTriplesMod[pt._leg[1]].size(), pt._leg[1]));
	if (flipAllowed[2]) flipOptions.push_back(FlipOption(pTriplesMod[pt._leg[2]].size(), pt._leg[2]));

	if (flipAllowed[2] && flipAllowed[1])
	{
		//flipOptions.push_back(FlipOption(pTriples[pt._leg[2]].size()+ pTriples[pt._leg[1]].size(), pt._leg[2], pt._leg[1]));
	}

	if (flipAllowed[2] && flipAllowed[0])
	{
		//flipOptions.push_back(FlipOption(pTriples[pt._leg[2]].size() + pTriples[pt._leg[0]].size(), pt._leg[2], pt._leg[0]));
	}

	if (flipAllowed[1] && flipAllowed[0])
	{
		//flipOptions.push_back(FlipOption(pTriples[pt._leg[1]].size() + pTriples[pt._leg[0]].size(), pt._leg[1], pt._leg[0]));
	}

	flipOptions.sort();

	list<FlipOption>::iterator foIter;
	int opt = 1;

	for (foIter = flipOptions.begin(); foIter != flipOptions.end(); foIter++)
	{
		foIter->_fixedNodes = fixedNodes;
		foIter->_nodeColors = inColors;
		foIter->_nodeColors[foIter->_fN1]._state = !foIter->_nodeColors[foIter->_fN1]._state;
		cout << "     O" << opt++ << ": " << foIter->_fN1;
		if (foIter->_fN2 != -1)
		{
			foIter->_nodeColors[foIter->_fN2]._state = !foIter->_nodeColors[foIter->_fN2]._state;
			cout << " & " << foIter->_fN2;
		}
		else
		{
			cout << "      ";
		}

		int conf = numPTriplesToBeFixed(*foIter, foIter->_fN1);
		if (foIter->_fN2 != -1)
		{
			conf += numPTriplesToBeFixed(*foIter, foIter->_fN2);
		}
		foIter->_cnf = conf;
		cout << " => C = " << conf << endl;

	}

	flipOptions.sort();
	if (flipOptions.size() > 10)
	{
		flipOptions.resize(10);
	}

	return flipOptions;
}

bool allColorsSameHandledBFS(const PythogoreanTriple& pt, vector<Color>& inColors, set<long long>& fixedNodesHL, int maxPops=10000)
{
	if (pt.isColoredValidly(inColors)) return true;
	int numPops = 0;
	list<FlipOption> flipOptions = analyzePTriple(pt, inColors, fixedNodesHL);
	if (flipOptions.size() == 0)
	{
		int aabv = 4;
	}

	list<FlipOption>::iterator foIt, foIt2;

	set<long long> fixedNodes;

	while (flipOptions.size() > 0)
	{
		flipOptions.sort();
		if (flipOptions.size() > 100)
		{
			flipOptions.resize(100);
		}
		FlipOption fo = flipOptions.front();
		numPops++;
		if (numPops > maxPops) return false;
		flipOptions.pop_front();
		vector<long long>::iterator pPTIter = fo._triplesToFix.begin();
		if (pPTIter == fo._triplesToFix.end())
		{
			if (fo._fN2 == -1)
			{
				cout << "Flipping  " << fo._fN1 << endl;
			}
			else
			{
				cout << "Flipping  " << fo._fN1 << " & " << fo._fN2 << endl;
			}
			inColors = fo._nodeColors;
			return true;
		}
		PythogoreanTriple pt = allPTriples[*pPTIter];
		list<FlipOption> fosLower = analyzePTriple(pt, fo._nodeColors, fo._fixedNodes);

		for (pPTIter++; pPTIter != fo._triplesToFix.end(); pPTIter++)
		{
			PythogoreanTriple pt = allPTriples[*pPTIter];
			list<FlipOption> fOps;

			for (foIt = fosLower.begin(); foIt != fosLower.end(); foIt++)
			{
				list<FlipOption> fosCurr = analyzePTriple(pt, foIt->_nodeColors, foIt->_fixedNodes);
				for (foIt2 = fosCurr.begin(); foIt2 != fosCurr.end(); foIt2++)
				{
					// Merge the two triplesToFix
					vector<long long>::iterator qPTIt = foIt->_triplesToFix.begin();
					for (; qPTIt != foIt->_triplesToFix.end(); qPTIt++)
					{
						foIt2->_triplesToFix.push_back(*qPTIt);
					}
					foIt2->_cnf += foIt->_cnf;
					fOps.push_back(*foIt2);
				}
			}
			fosLower.clear();
			if (fOps.size() > 20)
			{
				fOps.sort();
				fOps.resize(20);
			}
			fosLower = fOps;
			fOps.clear();
		}

		fosLower.sort();

		for (foIt = fosLower.begin(); foIt != fosLower.end(); foIt++)
		{
			foIt->_cnf += fo._cnf;
			flipOptions.push_back(*foIt);
		}
		fosLower.clear();
	}

	return false;
}

void matchOppositely(vector<Color>& colors, long long a, long long b)
{
	if (colors[a]._ndx == -1 && colors[b]._ndx == -1)
	{
		if (a < b)
		{
			colors[a]._ndx = a;
			colors[b].equateNegatively(colors[a]);
		}
		else
		{
			colors[b]._ndx = b;
			colors[a].equateNegatively(colors[b]);
		}

		return;
	}

	if (colors[a]._ndx == -1)
	{
		colors[a].equateNegatively(colors[b]);
		return;
	}

	if (colors[b]._ndx == -1)
	{
		colors[b].equateNegatively(colors[a]);
		return;
	}

	if (colors[a]._ndx == colors[b]._ndx) return; 

	mergeIndices(colors, colors[a], colors[b]);
}

void updateTriples(const PythogoreanTriple& t, long long e)
{
	pPytList& ppl2 = pTriplesMod[t._leg[2]];
	pPytListIt pIt2 = find(ppl2.begin(), ppl2.end(), e);
	ppl2.erase(pIt2);
	pPytList& ppl1 = pTriplesMod[t._leg[1]];
	pPytListIt pIt1 = find(ppl1.begin(), ppl1.end(), e);
	ppl1.erase(pIt1);
	pPytList& ppl0 = pTriplesMod[t._leg[0]];
	pPytListIt pIt0 = find(ppl0.begin(), ppl0.end(), e);
	ppl0.erase(pIt0);
}

void colorAPythagoreanTriple(long long i, vector<Color>& colors)
{
	maxGlobalElementIndex = i;
	PythogoreanTriple& t = allPTriples[i];
	cout << i << ": " << t.print(pTriplesMod) << endl;
	if (t.isColoredValidly(colors)) goto printColors;

	maxUsedIndx = max(maxUsedIndx, t._leg[2]);

	long maxCount = pTriplesMod[t._leg[0]].size();
	t.maxCountIndx = 0;
	if (pTriplesMod[t._leg[1]].size() >= maxCount)
	{
		maxCount = pTriplesMod[t._leg[1]].size();
		t.maxCountIndx = 1;
	}
	if (pTriplesMod[t._leg[2]].size() >= maxCount)
	{
		maxCount = pTriplesMod[t._leg[2]].size();
		t.maxCountIndx = 2;
	}
	int n1 = t.maxCountIndx;
	int n2 = (n1 + 1) % 3;
	int n3 = (n1 + 2) % 3;
	matchOppositely(colors, t._leg[n2], t._leg[n3]);

	short numUnassigned = 0;
	if (colors[t._leg[0]].isUnassigned()) numUnassigned++;
	if (colors[t._leg[1]].isUnassigned()) numUnassigned++;
	if (colors[t._leg[2]].isUnassigned()) numUnassigned++;

	switch (numUnassigned)
	{
	case 3:
	{
		// If all 3 are unassigned, then pick the two with the lowest counts.
		if (t._leg[n2] < t._leg[n3])
		{
			colors[t._leg[n2]]._ndx = t._leg[n2];
			colors[t._leg[n3]].equateNegatively(colors[t._leg[n2]]);
		}
		else
		{
			colors[t._leg[n3]]._ndx = t._leg[n3];
			colors[t._leg[n2]].equateNegatively(colors[t._leg[n3]]);
		}

		break;
	}
	case 2:
		if (colors[t._leg[0]].isAssigned())
		{
			if (pTriplesMod[t._leg[1]].size() <= pTriplesMod[t._leg[2]].size())
			{
				colors[t._leg[1]].equateNegatively(colors[t._leg[0]]);
			}
			else
			{
				colors[t._leg[2]].equateNegatively(colors[t._leg[0]]);
			}
		}
		else if (colors[t._leg[1]].isAssigned())
		{
			if (pTriplesMod[t._leg[0]].size() <= pTriplesMod[t._leg[2]].size())
			{
				colors[t._leg[0]].equateNegatively(colors[t._leg[1]]);
			}
			else
			{
				colors[t._leg[2]].equateNegatively(colors[t._leg[1]]);
			}
		}
		else if (colors[t._leg[2]].isAssigned())
		{
			if (pTriplesMod[t._leg[0]].size() <= pTriplesMod[t._leg[1]].size())
			{
				colors[t._leg[0]].equateNegatively(colors[t._leg[2]]);
			}
			else
			{
				colors[t._leg[1]].equateNegatively(colors[t._leg[2]]);
			}
		}
		break;

	case 1:

		t.handleSingleUnassigned(colors);
		break;
	case 0: // All three indices are set
	{

		if (t.isColoredValidly(colors)) goto printColors;

		if (t.allColorsDifferentHandled(colors)) goto printColors;

		for (int n1 = 0; n1 <= 2; n1++)
		{
			int n2 = (n1 + 1) % 3;
			if (colors[t._leg[n1]]._ndx != colors[t._leg[n2]]._ndx)
			{
				mergeIndices(colors, colors[t._leg[n1]], colors[t._leg[n2]]);
				goto printColors;
			}

		} // n1

		set<long long> fixedNodes;
		bool passed = allColorsSameHandledBFS(t, colors, fixedNodes);

		if (!passed)
		{
			int nc = t.allColorsSameHandled(colors);
			int a = 3;
		}

		break;
	}
	default:
		assert(0);
	} // switch (numUnassigned)

printColors:

	clearPreviousIndices();

	cout << i << ": " << t.print(colors) << endl;
	updateTriples(t, i);

	for (long ik = 0; ik <= i; ik++)
	{
		const PythogoreanTriple& t = allPTriples[ik];
		Color c0 = colors[t._leg[0]];
		Color c1 = colors[t._leg[1]];
		Color c2 = colors[t._leg[2]];

		if (!t.isColoredValidly(colors))
		{
			if (colors[t._leg[0]]._ndx == -1)
			{
				matchOppositely(colors, t._leg[0], t._leg[2]);
			}
			else if (colors[t._leg[1]]._ndx == -1)
			{
				matchOppositely(colors, t._leg[1], t._leg[2]);
			}
			else
			{
				set<long long> fixNodes;
				bool passed = allColorsSameHandledBFS(t, colors, fixNodes);
			}
		}
	}

}

void readColorsFromFile(long long& startElem, long long& startNode, string fileName, vector<Color>& inColors);

void saveColorsToFile(long long ptEnd, long long nodeEnd, vector<Color>& inColors, string fileName);

int getConflict(vector<Color>& colors, long long a, long long b, vector<long long>& nodeToFlip)
{
	vector<Color> tColors = colors;
	tColors[a]._state = !tColors[a]._state;
	int cnf0 = checkTriplesWithGivenNode(a, tColors);
	tColors = colors;
	tColors[b]._state = !tColors[b]._state;
	int cnf1 = checkTriplesWithGivenNode(b, tColors);

	if (cnf0 <= cnf1)
		nodeToFlip.push_back(a);
	else
		nodeToFlip.push_back(b);

	return min(cnf0, cnf1);
	//return cnf0+cnf1;
}

bool tryMultipleTriples(vector<Color>& colors, long long& k)
{
	PythogoreanTriple& t = allPTriples[k];
	cout << "-------------------------------------------------------------" << endl;
	cout << k << ": All " << t.print(colors) << endl;

	long long tc = t._leg[2];
	long long i;

	for (i = 1; i < 1000; i++)
	{
		if (k + i >= allPTriples.size()) break;
		PythogoreanTriple& t2 = allPTriples[k+i];
		if (t2._leg[2] != tc) break;
		cout << k+i << ": All " << t2.print(colors) << endl;
	}

	if (i == 1) return false;

	cout << "==================================================" << endl;
	long long kEnd = k + i;
	int numConflicts = 0;
	long long conflictElem;
	int numValids = 0;
	vector<long long> trueVotes;
	vector<long long> falseVotes;

	for (i = k; i < kEnd; i++)
	{
		PythogoreanTriple& t = allPTriples[i];
		if (t.isColoredValidly(colors))
		{
			numValids++;
			continue;
		}
		Color c0 = colors[t._leg[0]];
		Color c1 = colors[t._leg[1]];
		Color c2 = colors[t._leg[2]];
		if (c0._ndx != -1 && c0._ndx == c1._ndx && c0._state == c1._state)
		{
			numConflicts++;
			conflictElem = i;
			cout << i << ": TBF " << t.print(colors) << endl;

			if (c2._ndx == -1)
			{
				if (c0._state == false)
				{
					trueVotes.push_back(i);
				}
				else
				{
					falseVotes.push_back(i);
				}
			}
			else
			{
				int asas = 4;
			}

		}
	}

	if (numValids == (kEnd - k))
	{
		k = kEnd - 1;
		return true;
	}

	if (numConflicts < 1) return false;

	if (trueVotes.size() == 0 && falseVotes.size() == 0)
	{
		int abc = 5;
	}

	if (trueVotes.size() > 0 && falseVotes.size() == 0)
	{
		PythogoreanTriple& t = allPTriples[conflictElem];
		colors[t._leg[2]]._ndx = colors[t._leg[0]]._ndx;
		colors[t._leg[2]]._state = true;
		cout << "Setting " << t._leg[2] << " to " << colors[t._leg[0]]._ndx << " T" << endl;

		return false;
	}

	if (falseVotes.size() > 0 && trueVotes.size() == 0)
	{
		PythogoreanTriple& t = allPTriples[conflictElem];
		colors[t._leg[2]]._ndx = colors[t._leg[0]]._ndx;
		colors[t._leg[2]]._state = false;
		cout << "Setting " << t._leg[2] << " to " << colors[t._leg[0]]._ndx << " F" << endl;

		return false;
	}

	vector<long long> tNodesToFlip;
	vector<long long> tMaxElems;
	int trueConflicts = 0;
	for (int j = 0; j < trueVotes.size(); j++)
	{
		long long k = trueVotes[j];
		PythogoreanTriple& t = allPTriples[k];
		maxGlobalElementIndex = k;
		tMaxElems.push_back(k);
		trueConflicts += getConflict(colors, t._leg[0], t._leg[1], tNodesToFlip);
	}

	vector<long long> fNodesToFlip;
	vector<long long> fMaxElems;
	int falseConflicts = 0;
	for (int j = 0; j < falseVotes.size(); j++)
	{
		long long k = falseVotes[j];
		PythogoreanTriple& t = allPTriples[k];
		maxGlobalElementIndex = k;
		fMaxElems.push_back(k);
		falseConflicts += getConflict(colors, t._leg[0], t._leg[1], fNodesToFlip);
	}

	PythogoreanTriple& pt = allPTriples[conflictElem];
	set<long long> fixedNodes;
	fixedNodes.insert(pt._leg[2]);
	bool stateOption;

	if (falseConflicts == trueConflicts)
	{
		stateOption = (trueVotes.size() >= falseVotes.size());
	}
	else
	{
		stateOption = trueConflicts > falseConflicts;
	}
	int numPopsMax = 10;
redo:

	if (stateOption)
	{
		vector<Color> tColors = colors;
		set<long long> tfn = fixedNodes;
		for (int j = 0; j < falseVotes.size(); j++)
		{
			long long k = falseVotes[j];
			PythogoreanTriple& t = allPTriples[k];
			maxGlobalElementIndex = k;
			bool passed = allColorsSameHandledBFS(t, tColors, tfn, numPopsMax);
			if (!passed)
			{
				stateOption = !stateOption;
				numPopsMax += 10;
				goto redo;
			}
			colors = tColors;
		}
	}
	else
	{
		vector<Color> tColors = colors;
		set<long long> tfn = fixedNodes;
		for (int j = 0; j < trueVotes.size(); j++)
		{
			long long k = trueVotes[j];
			PythogoreanTriple& t = allPTriples[k];
			maxGlobalElementIndex = k;
			bool passed = allColorsSameHandledBFS(t, tColors, tfn, numPopsMax);
			if (!passed)
			{
				stateOption = !stateOption;
				numPopsMax += 10;
				goto redo;
			}
			colors = tColors;
		}
	}

	colors[pt._leg[2]]._ndx = colors[pt._leg[0]]._ndx;
	colors[pt._leg[2]]._state = stateOption;
	cout << "Setting " << pt._leg[2] << " to ";
	if (colors[pt._leg[2]]._state) cout << "!";
	cout << colors[pt._leg[0]]._ndx << endl;

	// Check all others
	for (i = k; i < kEnd; i++)
	{
		PythogoreanTriple& t = allPTriples[i];
		if (!t.isColoredValidly(colors))
		{
			if (colors[t._leg[1]]._ndx == -1)
			{
				matchOppositely(colors, t._leg[1], t._leg[2]);
			}
			else if (colors[t._leg[0]]._ndx == -1)
			{
				matchOppositely(colors, t._leg[0], t._leg[2]);
			}
			else if (colors[t._leg[1]]._ndx != colors[t._leg[2]]._ndx)
			{
				mergeIndices(colors, colors[t._leg[1]], colors[t._leg[2]]);
			}
			else if (colors[t._leg[0]]._ndx != colors[t._leg[2]]._ndx)
			{
				mergeIndices(colors, colors[t._leg[0]], colors[t._leg[2]]);
			}
			else
			{
				set<long long> fixNodes;
				fixNodes.insert(t._leg[2]);
				maxGlobalElementIndex = kEnd;
				bool passed = allColorsSameHandledBFS(t, colors, fixNodes);
			}
		}
		updateTriples(t, i);
	}

	k = kEnd - 1;
	return true;
}

void colorPythagoreanNumbers()
{
	vector<Color> colors(lim);
	long long iPTStart = -1;
	readColorsFromFile(iPTStart, maxUsedIndx, "SavedColors.txt", colors);
	long long i;

	for (i = iPTStart+1; i < allPTriples.size(); i++)
	{
		if (!tryMultipleTriples(colors, i))
		{
			colorAPythagoreanTriple(i, colors);
		}
	}

	saveColorsToFile(i-1, maxUsedIndx, colors, "SavedColors.txt");
}

long main()
{
	findAllPythagoreanTriples();

	//printMaxFreqTriples();
	//printMinFreqTriples();

	//printAllTriplesToFile();

	colorPythagoreanNumbers(); // such that all three numbers of a triple are not the same color.

	return 0;
}

void printMaxFreqTriples()
{
	ofstream outFile("PythagoreanTriplesMax.txt");
	long maxF = 1;
	for (int k = 5; k < pTriplesOrig.size(); k++)
	{
		if (pTriplesOrig[k].size() > maxF)
		{
			maxF = pTriplesOrig[k].size();
			outFile << k << "    " << maxF << endl;
			for (int j = 0; j < maxF; j++)
			{
				int jj = pTriplesOrig[k][j];
				const PythogoreanTriple& t = allPTriples[jj];
				outFile << "     " << t._leg[0] << ", " << t._leg[1] << ", " << t._leg[2] << endl;
			}
			outFile << endl;
		}
	}
}

void printMinFreqTriples()
{
	ofstream outFile("PythagoreanTriplesMin.txt");
	long minF = 1;
	for (int k = 3; k < pTriplesOrig.size(); k++)
	{
		if (pTriplesOrig[k].size() == minF)
		{
			outFile << k << endl;
			for (int j = 0; j < minF; j++)
			{
				int jj = pTriplesOrig[k][j];
				const PythogoreanTriple& t = allPTriples[jj];
				outFile << "     " << t._leg[0] << ", " << t._leg[1] << ", " << t._leg[2] << endl;
			}
			outFile << endl;
		}
	}
}


void readColorsFromFile(long long& a, long long& b, string fileName, vector<Color>& inColors)
{
	ifstream inFile;
	inFile.open(fileName);

	if (!inFile.is_open()) return;

	inFile >> a;
	inFile >> b;
	long long j, n;
	bool s;
	for (long long k = 0; k <= b; k++)
	{
		inFile >> j >> n >> s;
		assert(inColors[j]._id == j);
		inColors[j]._ndx = n;
		inColors[j]._state = s;
	}
	inFile.close();
}

void saveColorsToFile(long long ptEnd, long long nodeEnd, vector<Color>& inColors, string fileName)
{
	ofstream myFile;
	myFile.open(fileName);
	myFile << ptEnd << endl;
	myFile << nodeEnd << endl;
	for (long long j = 0; j <= nodeEnd; j++)
	{
		assert(j == inColors[j]._id);
		myFile << j << "   " << inColors[j]._ndx << "   " << inColors[j]._state << endl;
	}

	myFile.close();
	return;
}
