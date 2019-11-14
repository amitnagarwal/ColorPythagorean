#pragma once
#include <set>
#include <vector>
#include "Color.h"
using namespace std;

class FlipOption
{
public:
	FlipOption() : _cnf(-1) {}
	FlipOption(int c, long long inNode1, long long inNode2=-1);

	bool operator<(const FlipOption& that) const;

	long long       _fN1;
	long long       _fN2;
	int             _cnf;
	vector<Color>   _nodeColors;
	pPytList        _triplesToFix;
	set<long long>  _fixedNodes;
};