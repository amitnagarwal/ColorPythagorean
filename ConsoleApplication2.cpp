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
vector<pPytList> pTriples;

void Color::equateNegatively(const Color& that) 
{
	_ndx = that._ndx; 
	_state = !that._state; 
	if (_id > maxUsedIndx) maxUsedIndx = _id;
}

int checkTriplesWithGivenNode(FlipOption& fo, long long indx)
{
	cout << "Switchin " << indx << endl;
	vector<Color>& tmpColors = fo._nodeColors;

	long long numTriples = pTriples[indx].size();
	int numConflicts = 0;

	for (long long iTr = 0; iTr < numTriples; iTr++)
	{
		long long gTripleIndx = pTriples[indx][iTr];
		if (gTripleIndx >= maxGlobalElementIndex) break;

		if (gTripleIndx == 91)
		{
			int asa = 5;
		}

		const PythogoreanTriple& t = allPTriples[gTripleIndx];

		// cout << "      " << gTripleIndx << "] " << t.print(tmpColors) << endl;

		if (t.isColoredValidly(tmpColors)) continue;

		if (tmpColors[t._leg[0]]._ndx == -1 || tmpColors[t._leg[1]]._ndx == -1 || tmpColors[t._leg[2]]._ndx == -1)
		{
			t.handleSingleUnassigned(tmpColors);
			// cout << "      " << gTripleIndx << "} " << t.print(tmpColors) << endl;
		}
		else
		{
			int n1 = t.getLegIndex(indx);
			int n2 = (n1 + 1) % 3;
			int n3 = (n1 + 2) % 3;

			if (tmpColors[t._leg[n2]]._ndx != tmpColors[t._leg[n3]]._ndx)
			{
				mergeIndices(tmpColors, tmpColors[t._leg[n2]], tmpColors[t._leg[n3]]);
				// cout << "      " << gTripleIndx << "| " << t.print(tmpColors) << endl;
				continue;
			}

			fo._triplesToFix.push_back(gTripleIndx);
			numConflicts++;
			cout << "TBF2:     " << t.print(tmpColors) << endl;
		}

	} // iTr

	if (numConflicts == 0)
	{
		cout << "Switched " << indx << endl;
	}
	else
	{
		cout << "     conflicts = " << numConflicts << endl;
	}

	return numConflicts;
}


int checkTriplesWithGivenNode(long long indx, vector<Color>& colors, bool recurse = false)
{
	cout << "Switchin " << indx << endl;
	if (recurse) cout << "              recursively." << endl;
	vector<Color> inColors = colors;

	long long numTriples = pTriples[indx].size();
	int numConflicts = 0;

	for (long long iTr = 0; iTr < numTriples; iTr++)
	{
		long long gTripleIndx = pTriples[indx][iTr];
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
			int n1 = indx;
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

long long lim = 3000;
int lim10th = lim / 10;
void colorAPythagoreanTriple(long i, vector<Color>& colors);

void findAllPythagoreanTriples()
{
	pTriples.resize(lim);
	long long ptCount = 0;
	vector<Color> colors(lim);

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
				pTriples[t._leg[0]].push_back(ptCount);
				pTriples[t._leg[1]].push_back(ptCount);
				pTriples[t._leg[2]].push_back(ptCount);
				colorAPythagoreanTriple(ptCount, colors);
				ptCount++;
			}
		} // a
		if (count >= 1000)
		{
			int cc = pTriples[c].size();
			for (int j = 0; j < cc; j++)
			{
				int jj = pTriples[c][j];
				const PythogoreanTriple& t = allPTriples[jj];
				if (t._leg[2] == c)
				{
					cout << "     " << t._leg[0] << ", " << t._leg[1] << ", " << t._leg[2] << endl;
				}
			}
			cout << endl;
		}
	} // c

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
		int total = pTriples[t._leg[0]].size() + pTriples[t._leg[1]].size() + pTriples[t._leg[2]].size();
		outFile << "     " << t._leg[0] << ", " << t._leg[1] << ", " << t._leg[2] << " - " << total << "  ";

		if (pTriples[t._leg[0]].size() == 1 && pTriples[t._leg[1]].size() == 1 && pTriples[t._leg[2]].size() == 1)
		{
			outFile << " *** ";
		}
		outFile << endl;
	}

	for (long i = 0; i < lim; i++)
	{
		if (pTriples[i].size() > 0)
		{
			outFile << i << "   " << pTriples[i].size() << endl;
			for (unsigned long j = 0; j < pTriples[i].size(); j++)
			{
				const PythogoreanTriple& t = allPTriples[pTriples[i][j]];
				outFile << "     " << t._leg[0] << ", " << t._leg[1] << ", " << t._leg[2] << endl;
			} // j
		}
	}
}

void clearPreviousIndices();

FlipOption::FlipOption(long long inNode, int c)
{
	_flipNodes.insert(inNode);
	_conflict = c;
}

bool isFlippable(const set<long long>& inNodes, long long nodeIndx)
{
	set<long long>::iterator kIter = find(inNodes.begin(), inNodes.end(), nodeIndx);
	return (kIter != inNodes.end());
}

int process(
	const PythogoreanTriple& pt, 
	vector<Color>& inColors, 
	set<long long>& fixedNodes, 
	list<FlipOption>& flipOptions
)
{
	cout << "************************************************" << endl;
	cout << "TBF1:     " << pt.print(inColors) << endl;

	bool flipAllowed[3];
	
	flipAllowed[0] = !isFlippable(fixedNodes, pt._leg[0]);
	flipAllowed[1] = !isFlippable(fixedNodes, pt._leg[1]);
	flipAllowed[2] = !isFlippable(fixedNodes, pt._leg[2]);

	if (!flipAllowed[0] && !flipAllowed[1] && !flipAllowed[2])
	{
		return 1000;
	}

	fixedNodes.insert(pt._leg[0]);
	fixedNodes.insert(pt._leg[1]);
	fixedNodes.insert(pt._leg[2]);

	if (flipAllowed[0]) flipOptions.push_back(FlipOption(pt._leg[0], pTriples[pt._leg[0]].size()));
	if (flipAllowed[1]) flipOptions.push_back(FlipOption(pt._leg[1], pTriples[pt._leg[1]].size()));
	if (flipAllowed[2]) flipOptions.push_back(FlipOption(pt._leg[2], pTriples[pt._leg[2]].size()));

	flipOptions.sort();

	list<FlipOption>::iterator foIter;
	int minConflict = 10000;

	for (foIter = flipOptions.begin(); foIter != flipOptions.end(); foIter++)
	{
		long long node = *(foIter->_flipNodes.begin());
		foIter->_nodeColors = inColors;
		foIter->_nodeColors[node]._state = !foIter->_nodeColors[node]._state;
		int conf = checkTriplesWithGivenNode(*foIter, node);
		minConflict = min(minConflict, conf);
		foIter->_conflict = conf;
		if (conf == 0)
		{
			inColors = foIter->_nodeColors;
			return conf;
		}
	}

	return minConflict;
}

int allColorsSameHandled2(const PythogoreanTriple& pt, vector<Color>& inColors)
{
	list<FlipOption> flipOptions;
	set<long long> fixedNodes;
	int conf = process(pt, inColors, fixedNodes, flipOptions);
	if (conf == 0) return 0;

	flipOptions.sort();
	list<FlipOption> fosLower;
	list<FlipOption> fosLowerAll;
	list<FlipOption>::iterator foIt;

	if (flipOptions.size() == 0)
	{
		int aabv = 4;
	}

	while (1)
	{
		FlipOption fo = flipOptions.front();
		flipOptions.pop_front();
		vector<long long>::iterator pPTIter;
		int totalConflict = 0;
		if (fo._triplesToFix.size() >= 2)
		{
			int and = 5;
		}
		for (pPTIter = fo._triplesToFix.begin(); pPTIter != fo._triplesToFix.end(); pPTIter++)
		{
			PythogoreanTriple pt = allPTriples[*pPTIter];
			int conf2 = process(pt, fo._nodeColors, fixedNodes, fosLower);
			totalConflict += conf2;			
			for (foIt = fosLower.begin(); foIt != fosLower.end(); foIt++)
			{
				fosLowerAll.push_back(*foIt);
			}
			fosLower.clear();
		}

		if (totalConflict == 0)
		{
			inColors = fo._nodeColors;
			return 0;
		}

		if (flipOptions.size() == 0)
		{
			flipOptions = fosLowerAll;
			fosLowerAll.clear();
			flipOptions.sort();
		}
	}
	return flipOptions.begin()->_conflict;
}

void colorAPythagoreanTriple(long i, vector<Color>& colors)
{
	maxGlobalElementIndex = i;
	PythogoreanTriple& t = allPTriples[i];

	if (i == 412)
	{
		int a = 5;
	}

	long maxCount = pTriples[t._leg[0]].size();
	t.maxCountIndx = 0;
	if (pTriples[t._leg[1]].size() >= maxCount)
	{
		maxCount = pTriples[t._leg[1]].size();
		t.maxCountIndx = 1;
	}
	if (pTriples[t._leg[2]].size() >= maxCount)
	{
		maxCount = pTriples[t._leg[2]].size();
		t.maxCountIndx = 2;
	}

	short numUnassigned = 0;
	if (colors[t._leg[0]].isUnassigned()) numUnassigned++;
	if (colors[t._leg[1]].isUnassigned()) numUnassigned++;
	if (colors[t._leg[2]].isUnassigned()) numUnassigned++;

	switch (numUnassigned)
	{
	case 3:
	{
		// If all 3 are unassigned, then pick the two with the lowest counts.
		int n1 = t.maxCountIndx;
		int n2 = (n1 + 1) % 3;
		int n3 = (n1 + 2) % 3;

		colors[t._leg[n2]]._ndx = t._leg[n2];
		colors[t._leg[n3]].equateNegatively(colors[t._leg[n2]]);

		break;
	}
	case 2:
		if (colors[t._leg[0]].isAssigned())
		{
			if (pTriples[t._leg[1]].size() <= pTriples[t._leg[2]].size())
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
			if (pTriples[t._leg[0]].size() <= pTriples[t._leg[2]].size())
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
			if (pTriples[t._leg[0]].size() <= pTriples[t._leg[1]].size())
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

		int nc = allColorsSameHandled2(t, colors);
		//int nc = t.allColorsSameHandled(colors);
		if (nc == 0) goto printColors;
		else
		{
			int a = 4;
		}

		break;
	} // switch (numUnassigned)

	int abc = 2;

printColors:

	clearPreviousIndices();

	cout << i << ": " << t.print(colors) << endl;

	for (long ik = 0; ik <= i; ik++)
	{
		const PythogoreanTriple& pt = allPTriples[ik];
		Color c0 = colors[pt._leg[0]];
		Color c1 = colors[pt._leg[1]];
		Color c2 = colors[pt._leg[2]];

		if (ik == 245)
		{
			int a = 5;
		}

		if (pt.isColoredValidly(colors)) continue;

		int a = 3;
	}

}

void colorPythagoreanNumbers()
{
	vector<Color> colors(lim);

	for (long i = 0; i < allPTriples.size(); i++)
	{
		colorAPythagoreanTriple(i, colors);
	}
}

long main() {
	findAllPythagoreanTriples();

	//printMaxFreqTriples();
	//printMinFreqTriples();

	printAllTriplesToFile();

	//colorPythagoreanNumbers(); // such that all three numbers of a triple are not the same color.

	return 0;
}

void printMaxFreqTriples()
{
	ofstream outFile("PythagoreanTriplesMax.txt");
	long maxF = 1;
	for (int k = 5; k < pTriples.size(); k++)
	{
		if (pTriples[k].size() > maxF)
		{
			maxF = pTriples[k].size();
			outFile << k << "    " << maxF << endl;
			for (int j = 0; j < maxF; j++)
			{
				int jj = pTriples[k][j];
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
	for (int k = 3; k < pTriples.size(); k++)
	{
		if (pTriples[k].size() == minF)
		{
			outFile << k << endl;
			for (int j = 0; j < minF; j++)
			{
				int jj = pTriples[k][j];
				const PythogoreanTriple& t = allPTriples[jj];
				outFile << "     " << t._leg[0] << ", " << t._leg[1] << ", " << t._leg[2] << endl;
			}
			outFile << endl;
		}
	}
}


