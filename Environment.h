#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include "Light.h"
#include "Sphere.h"
#include "Triangle.h"
#include <math.h>
#include <algorithm>

struct Texture {
	int width;
	int height;
	std::vector<double*> colors;
};

class Environment {

public:

	double eye[3];		// [x, y, z]
	double vdir[3];		// [x, y, z]
	double up[3];		// [x, y, z]
	double fovh;		// in degrees
	int width;			// in pixel units
	int height;			// in pixel units
	double bkg[3];		// [r, g, b]
	double viewDist;	// view distance
	int parallel;		// parallel projection flag
	std::vector<Light*> lightVector;				// Vector containing lights
	std::vector<Texture> textureVector;				// Vector containing textures
	std::vector<ShapeObject*> objectVector;			// Vector containing shape objects
	std::vector<double*> vertexCoordVector;			// Vector containing vertex coordinates
	std::vector<double*> textureCoordVector;		// Vector containing texture coordinates
	std::vector<double*> surfaceNormCoordVector;	// Vector containing surface normal coordinates
	std::vector<double*> outputPixels;				// Vector containing the output image data
	std::vector<double*> viewWindow;				// Viewing window coordinates, contains the scene "rays"

	Environment();

	void addLight(Light* light);
	void addTexture(Texture texture);
	void addObject(ShapeObject* object);
	void addVertexCoord(double* vertexCoord);
	void addTextureCoord(double* textureCoord);
	void addSurfaceNormCoord(double* surfaceNormCoord);
	void generateScene();

};

#endif
