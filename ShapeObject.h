#ifndef SHAPEOBJECT_H
#define SHAPEOBJECT_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

enum ShapeType {
	sphere,
	triangle
};

class ShapeObject {

public:

	ShapeObject();

	ShapeType type;			// Object shape type
	double mtld[3];			// [odr, odg, odb]
	double mtls[3];			// [osr, osg, osb]
	double mtlk[3];			// [ka, kd, ks]
	double n;				// Specular falloff exponent
	double opacity;			// alpha
	double refraction;		// eta
	int ID;					// Object ID
	int texture;			// Texture file index

};

#endif
