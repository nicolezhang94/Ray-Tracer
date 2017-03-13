#include "Triangle.h"

using namespace std;

// Default constructor
Triangle::Triangle() {
	type = triangle;
	mtld = {-1.0, -1.0, -1.0};
	mtls = {-1.0, -1.0, -1.0};
	mtlk = {-1.0, -1.0, -1.0};
	n = -1.0;
	opacity = -1.0;
	refraction = -1.0;
	ID = -1;
	texture = -1;
	
	shading = flat;
}
