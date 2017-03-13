#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "ShapeObject.h"

enum ShadeType {
	flat,
	smooth
};

class Triangle : public ShapeObject {

public:

	std::vector<double*> v;		// Vector of length 3, contents array of length 3 (x, y, z)
	std::vector<double*> vt;	// Vector of length 3, contents array of length 2 (u, v)
	std::vector<double*> vn;	// Vector of length 3, contents array of length 3 (x, y, z)
	ShadeType shading;

	Triangle();

};

#endif
