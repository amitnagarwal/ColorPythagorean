#include <sstream>
#include "PythagoreanTriple.h"
#include "Color.h"
#include <assert.h>
#include <iostream>
#include <list>
#include <set>
#include "FlipOption.h"

int PythogoreanTriple::getLegIndex(long long nodeIndex) const
{
	if (_leg[0] == nodeIndex) return 0;
	if (_leg[1] == nodeIndex) return 1;
	if (_leg[2] == nodeIndex) return 2;
	return -1;
}

bool isPrimeNumber(long long inp)
{
	long lim = sqrt(inp) + 1;
	for (long i = 2; i < lim; i++)
	{
		if (inp%i == 0) return false;
	}
	return true;
}

bool PythogoreanTriple::isPrimitive() const
{
	if (isPrimeNumber(_leg[0])) return true;
	if (isPrimeNumber(_leg[2])) return true;

	// Check if there is a common factor
	long lim = sqrt(_leg[0]) + 1;
	for (long i = 2; i < lim; i++)
	{
		if (_leg[0] % i == 0)
		{
			if (_leg[1] % i == 0 && _leg[2] % i == 0) return false;
			long f = _leg[0] / i;
			if (_leg[1] % f == 0 && _leg[2] % f == 0) return false;
		}
	}

	return true;
}

string PythogoreanTriple::print() const
{
	stringstream buffer;
	buffer << "(" << _leg[0] << ", " << _leg[1] << ", " << _leg[2] << ")";
	return buffer.str();
}

string PythogoreanTriple::print(const vector<Color>& inColors) const
{
	stringstream buffer;

	buffer << "(" << _leg[0] << ", " << _leg[1] << ", " << _leg[2] << ") = (";
	if (inColors[_leg[0]]._state) buffer << "!";
	buffer << inColors[_leg[0]]._ndx << ", ";
	if (inColors[_leg[1]]._state) buffer << "!";
	buffer << inColors[_leg[1]]._ndx << ", ";
	if (inColors[_leg[2]]._state) buffer << "!";
	buffer << inColors[_leg[2]]._ndx << ")";

	return buffer.str();
}

string PythogoreanTriple::print(const vector<pPytList>& pTriples) const
{
	stringstream buffer;

	buffer << "(" << _leg[0] << ", " << _leg[1] << ", " << _leg[2] << ") = (";
	buffer << pTriples[_leg[0]].size() << ", ";
	buffer << pTriples[_leg[1]].size() << ", ";
	buffer << pTriples[_leg[2]].size() << ")";

	return buffer.str();
}

bool PythogoreanTriple::operator==(const PythogoreanTriple& pt) const
{
	return ((_leg[0] == pt._leg[0]) && (_leg[1] == pt._leg[1]) && (_leg[2] == pt._leg[2]));
}

bool PythogoreanTriple::isColoredValidly(const vector<Color>& colors) const
{
	if (colors[_leg[1]]._ndx == colors[_leg[0]]._ndx &&
		colors[_leg[1]]._state != colors[_leg[0]]._state)
	{
		return true;
	}
	if (colors[_leg[2]]._ndx == colors[_leg[0]]._ndx &&
		colors[_leg[2]]._state != colors[_leg[0]]._state)
	{
		return true;
	}
	if (colors[_leg[2]]._ndx == colors[_leg[1]]._ndx &&
		colors[_leg[2]]._state != colors[_leg[1]]._state)
	{
		return true;
	}

	return false;
}

bool PythogoreanTriple::isColoredValidlyOnFirstTwo(const vector<Color>& colors) const
{
	if (colors[_leg[1]]._ndx == colors[_leg[0]]._ndx &&
		colors[_leg[1]]._state != colors[_leg[0]]._state)
	{
		return true;
	}

	return false;
}

void PythogoreanTriple::handleSingleUnassigned(vector<Color>& inColors) const
{
	long long a = _leg[0];
	long long b = _leg[1];
	long long c = _leg[2];

	if (inColors[c].isUnassigned())
	{
		if (inColors[b]._ndx == inColors[a]._ndx)
		{
			if (inColors[b]._state != inColors[a]._state)
			{
				return;
			}
			else
			{
				inColors[c].equateNegatively(inColors[a]);
			}
		}
		else
		{
			mergeIndices(inColors, inColors[b], inColors[a]);
		}
	}
	else if (inColors[b].isUnassigned())
	{
		if (inColors[a]._ndx == inColors[c]._ndx)
		{
			if (inColors[a]._state != inColors[c]._state)
			{
				return;
			}
			else
			{
				inColors[b].equateNegatively(inColors[a]);
			}
		}
		else
		{
			mergeIndices(inColors, inColors[c], inColors[a]);
		}
	}
	else if (inColors[a].isUnassigned())
	{
		if (inColors[b]._ndx == inColors[c]._ndx)
		{
			if (inColors[b]._state != inColors[c]._state)
			{
				return;
			}
			else
			{
				inColors[a].equateNegatively(inColors[b]);
			}
		}
		else
		{
			mergeIndices(inColors, inColors[c], inColors[b]);
		}
	}
}

bool PythogoreanTriple::allColorsDifferentHandled(vector<Color>& inColors) const
{
	long long a = _leg[0];
	long long b = _leg[1];
	long long c = _leg[2];

	if ((inColors[b]._ndx != inColors[a]._ndx) &&
		(inColors[b]._ndx != inColors[c]._ndx) &&
		(inColors[a]._ndx != inColors[c]._ndx))
	{
		switch (maxCountIndx)
		{
		case 0:
			mergeIndices(inColors, inColors[_leg[1]], inColors[_leg[2]]);
			return true;
		case 1:
			mergeIndices(inColors, inColors[_leg[0]], inColors[_leg[2]]);
			return true;
		case 2:
			mergeIndices(inColors, inColors[_leg[0]], inColors[_leg[1]]);
			return true;
		default:
			assert(0);
		}
	}
	return false;
}

set<long long> previousIndices;

void clearPreviousIndices()
{
	previousIndices.clear();
}

bool previouslySwitched(long long nodeIndx)
{
	set<long long>::iterator kIter = find(previousIndices.begin(), previousIndices.end(), nodeIndx);
	return (kIter != previousIndices.end());
}

class KVPair
{
public:
	KVPair(short k, int v) :_key(k), _val(v) {}

	bool operator<(const KVPair& that) const { return (_val < that._val); }

	short _key;
	int   _val;
};

int checkTriplesWithGivenNode(long long indx, vector<Color>& colors, bool recurse = false);

int PythogoreanTriple::allColorsSameHandled(vector<Color>& inColors) const
{
	cout << "TBF:     " << this->print(inColors) << endl;
	bool switchAllowed[3];

	switchAllowed[0] = !previouslySwitched(_leg[0]);
	switchAllowed[1] = !previouslySwitched(_leg[1]);
	switchAllowed[2] = !previouslySwitched(_leg[2]);

	if (!switchAllowed[0] && !switchAllowed[1] && !switchAllowed[2])
	{
		return 1000;
	}
	long long a = _leg[0];
	long long b = _leg[1];
	long long c = _leg[2];

	// None of the three nodes are allowed to be switched in the downstream recursive calls
	previousIndices.insert(_leg[2]);
	previousIndices.insert(_leg[1]);
	previousIndices.insert(_leg[0]);

	set<long long> previousIndicesSaved = previousIndices;
	vector<Color>  tmpColors;
	int numConflicts[3];

	list<KVPair> conflictMap;
	list<KVPair>::iterator mapIt;

	if (switchAllowed[2])
	{
		// Flip c only and not a & b
		tmpColors = inColors;
		tmpColors[c]._state = !tmpColors[c]._state;
		numConflicts[2] = checkTriplesWithGivenNode(c, tmpColors);
		if (numConflicts[2] == 0) goto completed;
		conflictMap.push_back(KVPair(2, numConflicts[2]));
	}

	if (switchAllowed[1])
	{
		// Flip b only and not a & c
		tmpColors = inColors;
		tmpColors[b]._state = !tmpColors[b]._state;
		numConflicts[1] = checkTriplesWithGivenNode(b, tmpColors);
		if (numConflicts[1] == 0) goto completed;
		conflictMap.push_back(KVPair(1, numConflicts[1]));

		if (switchAllowed[2])
		{
			conflictMap.push_back(KVPair(5, numConflicts[1] + numConflicts[2]));
		}
	}

	if (switchAllowed[0])
	{
		// Flip a only and not b & c
		tmpColors = inColors;
		tmpColors[a]._state = !tmpColors[a]._state;
		numConflicts[0] = checkTriplesWithGivenNode(a, tmpColors);
		if (numConflicts[0] == 0) goto completed;
		conflictMap.push_back(KVPair(0, numConflicts[0]));

		if (switchAllowed[2])
		{
			conflictMap.push_back(KVPair(4, numConflicts[0] + numConflicts[2]));
		}
		if (switchAllowed[1])
		{
			conflictMap.push_back(KVPair(3, numConflicts[0] + numConflicts[1]));
		}
	}
	conflictMap.sort();

	int nConflicts;
	for (mapIt = conflictMap.begin(); mapIt != conflictMap.end(); mapIt++)
	{
		int mk = mapIt->_key;
		previousIndices = previousIndicesSaved;
		tmpColors = inColors;

		if (mk < 3)
		{
			tmpColors[_leg[mk]]._state = !tmpColors[_leg[mk]]._state;
			nConflicts = checkTriplesWithGivenNode(_leg[mk], tmpColors, true);
		}
		else
		{
			switch (mk)
			{
			case 3:
				tmpColors[_leg[0]]._state = !tmpColors[_leg[0]]._state;
				tmpColors[_leg[1]]._state = !tmpColors[_leg[1]]._state;
				nConflicts = checkTriplesWithGivenNode(_leg[0], tmpColors, true);
				nConflicts += checkTriplesWithGivenNode(_leg[1], tmpColors, true);
				break;
			case 4:
				tmpColors[_leg[0]]._state = !tmpColors[_leg[0]]._state;
				tmpColors[_leg[2]]._state = !tmpColors[_leg[2]]._state;
				nConflicts = checkTriplesWithGivenNode(_leg[0], tmpColors, true);
				nConflicts += checkTriplesWithGivenNode(_leg[2], tmpColors, true);
				break;
			case 5:
				tmpColors[_leg[1]]._state = !tmpColors[_leg[1]]._state;
				tmpColors[_leg[2]]._state = !tmpColors[_leg[2]]._state;
				nConflicts = checkTriplesWithGivenNode(_leg[1], tmpColors, true);
				nConflicts += checkTriplesWithGivenNode(_leg[2], tmpColors, true);
				break;
			default:
				assert(0);
			} // switch

		} // if

		if (nConflicts == 0)
			goto completed;
		else
		{
			int my = 5;
		}
	}

	return nConflicts;

completed:

	inColors = tmpColors;
//	previousIndices.clear();
	return 0;
}


