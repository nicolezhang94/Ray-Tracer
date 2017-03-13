#ifndef SPHERE_H
#define SPHERE_H

#include "ShapeObject.h"

class Sphere : public ShapeObject {

public:

	double location[3];		// [x, y, z]
	double r;				// radius

	Sphere();

};

#endif
