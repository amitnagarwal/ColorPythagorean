#pragma once
#include <string>
#include <vector>

using namespace std; 

class Color;

class PythogoreanTriple
{
public:
	bool isPrimitive() const;
	string print(const vector<Color>& inColors) const;
	bool isColoredValidly(const vector<Color>& colors) const;
	void handleSingleUnassigned(vector<Color>& inColors) const;
	bool allColorsDifferentHandled(vector<Color>& inColors) const;
	int allColorsSameHandled(vector<Color>& inColors) const;

	long long _leg[3];
	short maxCountIndx;
};
