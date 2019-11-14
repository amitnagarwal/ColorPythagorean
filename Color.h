#pragma once
#include <vector>

static long startID = 0;


class Color
{
public:
	Color() : _ndx(-1), _state(false), _id(startID++) {}

	bool isAssigned() const { return (_ndx != -1); }
	bool isUnassigned() const { return (_ndx == -1); }
	bool operator==(const Color& that) const { return ((_ndx == that._ndx) && (_state == that._state)); }
	void equate(const Color& that) { _ndx = that._ndx; _state = that._state; }
	void equateNegatively(const Color& that);

	long _id;
	long long _ndx;
	bool _state;
};

void mergeIndices(vector<Color>& inpColors, Color& color1, Color& color2);
