#pragma once
#include <string>
#include <vector>

using namespace std; 

class Color;
typedef vector<long long> pPytList;

class PythogoreanTriple
{
public:
	bool isPrimitive() const;
	string print() const;
	string print(const vector<Color>& inColors) const;
	string print(const vector<pPytList>& pTriples) const;
	bool operator==(const PythogoreanTriple& pt) const;

	int getLegIndex(long long nodeIndex) const;
	bool isSingleColored(const vector<Color>& colors) const;
	bool isColoredValidly(const vector<Color>& colors) const;
	bool isColoredValidlyOnFirstTwo(const vector<Color>& colors) const;
	short getNumUnassigned(const vector<Color>& colors) const;
	void handleSingleUnassigned(vector<Color>& inColors) const;
	void handleDoubleUnassigned(vector<Color>& inColors) const;
	bool allColorsDifferentHandled(vector<Color>& inColors) const;
	int allColorsSameHandled(vector<Color>& inColors) const;

	long long _leg[3];
	short maxCountIndx;
};

typedef vector<PythogoreanTriple> PytList;
typedef vector<PythogoreanTriple>::iterator PytListIt;
typedef vector<long long>::iterator pPytListIt;
