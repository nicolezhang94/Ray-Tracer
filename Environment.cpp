#include "Environment.h"

#define PI 3.14159265358979323846264338327950288

using namespace std;

int limit = 0;

// Default constructor
Environment::Environment() {
	eye = {0.0, 0.0, 0.0};
	vdir = {0.0, 0.0, 0.0};
	up = {0.0, 0.0, 0.0};
	fovh = 200.0;
	width = 0;
	height = 0;
	bkg = {-1.0, -1.0, -1.0};
	viewDist = 5;
	parallel = 0;
}

// Cross product of two vectors of length 3
double* cross(double *u, double *v) {
	double* solution = new double[3];

	solution[0] = (u[1] * v[2]) - (u[2] * v[1]);
	solution[1] = (u[2] * v[0]) - (u[0] * v[2]);
	solution[2] = (u[0] * v[1]) - (u[1] * v[0]);

	return solution;
}

// Dot product of two vectors of length 3
double dot(double *u, double *v) {
	return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]);
}

// Magnitude of vector of length 3
double vectorLength(double *v) {
	return sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
}

// Multiply a vector of length 3 by a scalar
double* multScalar(double *v, double s) {
	double* solution = new double[3];

	solution[0] = v[0] * s;
	solution[1] = v[1] * s;
	solution[2] = v[2] * s;

	return solution;
}

// Multiply two vectors of length 3 together
double* mult(double *u, double* v) {
	double* solution = new double[3];

	solution[0] = u[0] * v[0];
	solution[1] = u[1] * v[1];
	solution[2] = u[2] * v[2];

	return solution;
}

// Divide a vector of length 3 by a scalar
double* divideScalar(double *v, double s) {
	double* solution = new double[3];

	solution[0] = v[0] / s;
	solution[1] = v[1] / s;
	solution[2] = v[2] / s;

	return solution;
}

// Add two vectors of length 3
double* add(double *u, double *v) {
	double* solution = new double[3];

	solution[0] = u[0] + v[0];
	solution[1] = u[1] + v[1];
	solution[2] = u[2] + v[2];

	return solution;
}

// Subtract two vectors of length 3
double* subtract(double *u, double *v) {
	double* solution = new double[3];

	solution[0] = u[0] - v[0];
	solution[1] = u[1] - v[1];
	solution[2] = u[2] - v[2];

	return solution;
}

// Makes a vector unit length
double* normalize(double* u) {
	return divideScalar(u, vectorLength(u));
}

// Determines ray/sphere intersection; returns array of the two solutions for t
double* raySphereIntersect(double* rayDir, double* origin, Sphere* currentSphere) {
	double A = 1;
	double B = 2 * ((rayDir[0] * (origin[0] - currentSphere->location[0])) +
					(rayDir[1] * (origin[1] - currentSphere->location[1])) +
					(rayDir[2] * (origin[2] - currentSphere->location[2])));
	double C =  pow((origin[0] - currentSphere->location[0]), 2) +
				pow((origin[1] - currentSphere->location[1]), 2) +
				pow((origin[2] - currentSphere->location[2]), 2) -
				pow(currentSphere->r, 2);
	double t1, t2;
	t1 = ((-1 * B) + sqrt(pow(B, 2) - (4 * A * C))) / (2 * A);
	t2 = ((-1 * B) - sqrt(pow(B, 2) - (4 * A * C))) / (2 * A);
	
	double* ts = new double[2]{t1, t2};
	return ts;
}

// Determines ray/triangle intersection; returns t, alpha, beta, and gamma
double* rayTriangleIntersect(double* rayDir, double* origin, Triangle* currentTriangle) {
	double* p0 = currentTriangle->v.at(0);
	double* p1 = currentTriangle->v.at(1);
	double* p2 = currentTriangle->v.at(2);
	
	double a = p0[0] - p1[0];
	double b = p0[1] - p1[1];
	double c = p0[2] - p1[2];
	double d = p0[0] - p2[0];
	double e = p0[1] - p2[1];
	double f = p0[2] - p2[2];
	double g = rayDir[0];
	double h = rayDir[1];
	double i = rayDir[2];
	double j = p0[0] - origin[0];
	double k = p0[1] - origin[1];
	double l = p0[2] - origin[2];
	
	double m = e*i - h*f;
	double n = g*f - d*i;
	double o = d*h - e*g;
	double p = a*k - j*b;
	double q = j*c - a*l;
	double r = b*l - k*c;
	
	double s = a*m + b*n + c*o;
	
	double beta = (j*m + k*n + l*o) / s;
	double gamma = (i*p + h*q + g*r) / s;
	double t = -1.0 * (f*p + e*q + d*r) / s;
	double alpha = 1.0 - (beta + gamma);
	
	double* solution = new double[4]{t, alpha, beta, gamma};
	
	if (t < 0) {
		solution[0] = -1;
		return solution;
	}
	if (gamma < 0 || 1 < gamma) {
		solution[0] = -1;
		return solution;
	}
	if (beta < 0 || 1-gamma < beta) {
		solution[0] = -1;
		return solution;
	}
	return solution;
}

// Determines a textured object's diffuse color at the current pixel
double* calculateTexture(ShapeObject* currentObject, double* rayIntersection, Texture currentTexture, double* barycentric) {
	double* texturedColor = new double[3];
	int textureI, textureJ;
	double u, v;
	if (currentObject->type == sphere) {
		Sphere* currentSphere = (Sphere*)currentObject;
		double Nx = (rayIntersection[0] - currentSphere->location[0]) / currentSphere->r;
		double Ny = (rayIntersection[1] - currentSphere->location[1]) / currentSphere->r;
		double Nz = (rayIntersection[2] - currentSphere->location[2]) / currentSphere->r;
		double phi = acos(Nz);
		double theta = atan2(Ny, Nx);
		
		u = (theta + PI) / (2 * PI);
		v = phi / PI;
		
		textureI = (int)(0.5 + u * ((double)currentTexture.width - 1.0));
		textureJ = (int)(0.5 + v * ((double)currentTexture.height - 1.0));
	} else if (currentObject->type = triangle) {
		Triangle* currentTriangle = (Triangle*)currentObject;
		double alpha = barycentric[1];
		double beta = barycentric[2];
		double gamma = barycentric[3];
		double u0 = currentTriangle->vt.at(0)[0];
		double u1 = currentTriangle->vt.at(1)[0];
		double u2 = currentTriangle->vt.at(2)[0];
		double v0 = currentTriangle->vt.at(0)[1];
		double v1 = currentTriangle->vt.at(1)[1];
		double v2 = currentTriangle->vt.at(2)[1];
		
		u = alpha*u0 + beta*u1 + gamma*u2;
		v = alpha*v0 + beta*v1 + gamma*v2;
		
		textureI = (int)(u * ((double)currentTexture.width - 1.0));
		textureJ = (int)(v * ((double)currentTexture.height - 1.0));
	}
	texturedColor = currentTexture.colors.at(textureI + textureJ * currentTexture.height);
	
	return texturedColor;
}

// Calculates pixel output color using the Phong illumination model
double* calculateIllumination(ShapeObject* currentObject, double* origin, double* rayIntersection, double* rayDir, vector<Light*> lightVector, vector<ShapeObject*> objectVector, vector<Texture> textureVector, double* barycentric, double* bkg, ShapeObject* oldObject) {
	double* illumination;
	double* ambientTerm = multScalar(currentObject->mtld, currentObject->mtlk[0]);
	double* diffuseTerm = multScalar(currentObject->mtld, currentObject->mtlk[1]);
	double* specularTerm = multScalar(currentObject->mtls, currentObject->mtlk[2]);
	double* sumOfLights = new double[3]{0, 0, 0};
	double* N;
	double* L;
	double* H;
	
	// Determines whether to use smooth or flat shading
	if (currentObject->type == sphere) {
		Sphere* currentSphere = (Sphere*)currentObject;
		N = divideScalar(subtract(rayIntersection, currentSphere->location), currentSphere->r);
	} else if (currentObject->type == triangle) {
		Triangle* currentTriangle = (Triangle*)currentObject;
		if (currentTriangle->shading == flat) {
			double* p0 = currentTriangle->v.at(0);
			double* p1 = currentTriangle->v.at(1);
			double* p2 = currentTriangle->v.at(2);
			double* e1 = subtract(p1, p0);
			double* e2 = subtract(p2, p0);
			N = normalize(cross(e1, e2));
		} else if (currentTriangle->shading == smooth) {
			double alpha = barycentric[1];
			double beta = barycentric[2];
			double gamma = barycentric[3];
			double* n0 = currentTriangle->vn.at(0);
			double* n1 = currentTriangle->vn.at(1);
			double* n2 = currentTriangle->vn.at(2);
			double* sum = add(add(multScalar(n0, alpha), multScalar(n1, beta)), multScalar(n2, gamma));
			N = normalize(sum);
		}
	}
	
	double totalLights = lightVector.size();
	double totalObjects = objectVector.size();
	int currentObjectID = currentObject->ID;
	
	for (int i = 0; i < totalLights; i++) {
		Light* currentLight = lightVector.at(i);
		
		double S = 1;
		
		double* shadowRayDir;
		double t;
		double* shadowRay;
		
		if (currentLight->type == 0) {
			// Directional lights
			shadowRayDir = normalize(multScalar(currentLight->location, -1));
			
			// Check to see if object is in a shadow
			for (int j = 0; j < totalObjects; j++) {
				if (currentObjectID != objectVector.at(j)->ID) {
					if (objectVector.at(j)->type == sphere) {
						Sphere* shadowSphere = (Sphere*)objectVector.at(j);
						double* ts = raySphereIntersect(shadowRayDir, rayIntersection, shadowSphere);
						if (ts[0] > 0 || ts[1] > 0) {
							S *= (1-shadowSphere->opacity);
						}
					} else if (objectVector.at(j)->type == triangle) {
						Triangle* shadowTriangle = (Triangle*)objectVector.at(j);
						double* ts = rayTriangleIntersect(shadowRayDir, rayIntersection, shadowTriangle);
						if (ts[0] > 0) {
							S *= (1-shadowTriangle->opacity);
						}
					}
				}
			}
			L = normalize(multScalar(currentLight->location, -1));
		} else if (currentLight->type == 1) {
			// Positional lights
			shadowRayDir = normalize(subtract(currentLight->location, rayIntersection));
			t = vectorLength(subtract(currentLight->location, rayIntersection));
			
			// Check to see if object is in a shadow
			for (int j = 0; j < totalObjects; j++) {
				if (currentObjectID != objectVector.at(j)->ID) {
					if (objectVector.at(j)->type == sphere) {
						Sphere* shadowSphere = (Sphere*)objectVector.at(j);
						double* ts = raySphereIntersect(shadowRayDir, rayIntersection, shadowSphere);
						if (ts[0] < t && ts[0] > 0) {
							S *= (1-shadowSphere->opacity);
						}
						if (ts[1] < t && ts[1] > 0) {
							S *= (1-shadowSphere->opacity);
						}
					} else if (objectVector.at(j)->type == triangle) {
						Triangle* shadowTriangle = (Triangle*)objectVector.at(j);
						double* ts = rayTriangleIntersect(shadowRayDir, rayIntersection, shadowTriangle);
						if (ts[0] < t && ts[0] > 0) {
							S *= (1-shadowTriangle->opacity);
						}
					}
				}
			}
			L = normalize(subtract(currentLight->location, rayIntersection));
		}
		
		// Changes a textured object's diffuse color based on the specified texture color at that pixel
		if (currentObject->texture > -1) {
			Texture currentTexture = textureVector.at(currentObject->texture);
			double* texturedDiffuse = calculateTexture(currentObject, rayIntersection, currentTexture, barycentric);
			ambientTerm = multScalar(texturedDiffuse, currentObject->mtlk[0]);
			diffuseTerm = multScalar(texturedDiffuse, currentObject->mtlk[1]);
		}
		
		double* diffuseTotal;
		double* specularTotal;
		
		// Final Phong illumination calculation
		diffuseTotal = multScalar(diffuseTerm, max(0.0, dot(N, L)));
		double* rayDir = normalize(subtract(origin, rayIntersection));
		H = normalize(add(L, rayDir));
		specularTotal = multScalar(specularTerm, pow(max(0.0, dot(N, H)), currentObject->n)); 
		sumOfLights = add(sumOfLights, multScalar(mult(currentLight->color, add(diffuseTotal, specularTotal)), S));
	}
	
	double* reflection = bkg;
	
	// Recursive call to calculate reflections, transparency, and refractions
	if (limit < 5) {
		limit++;
		
		double ni = 1.0;
		double nt = currentObject->refraction;
		
		if (limit == 1) {
			ni = 1.0;
			nt = currentObject->refraction;
		} else {
			ni = oldObject->refraction;
			nt = currentObject->refraction;
		}
		double opacity = currentObject->opacity;
		double* I = multScalar(normalize(rayDir), -1);
		if (dot(I, N) < 0) {
			N = multScalar(N, -1);
		}
		double angleIncident = acos(dot(I, N));
		//double angleTransmitted;
		double a = cos(angleIncident);
		double F0 = pow((nt-ni)/(nt+ni), 2);
		double Fr = F0 + (1-F0) * pow((1-cos(angleIncident)), 5);
		double* R = normalize(subtract(multScalar(N, 2*a), I)); 
		double* Rlambda = bkg;
		double* T = R;
		double* Tlambda = Rlambda;
		bool TIR = false;
		
		if (sin(angleIncident) > (nt/ni)) {
			TIR = true;
		} else {
			//angleTransmitted = asin((ni/nt) * sin(angleIncident));
			T = add(multScalar(multScalar(N, -1), sqrt(1-(pow(ni/nt, 2) * (1-pow(a, 2))))), multScalar(subtract(multScalar(N, a), I), ni/nt));
			Tlambda = bkg;
		}
		
		double rt = 10000;
		double tt = 10000;
		for (int j = 0; j < totalObjects; j++) {
			ShapeObject* newObject = objectVector.at(j);
			if (currentObjectID != newObject->ID) {
				if (newObject->type == sphere) {
					Sphere* newSphere = (Sphere*)newObject;
					double* rts = raySphereIntersect(R, rayIntersection, newSphere);
					
					if (rts[0] < rts[1]) {
						if (rts[0] > 0 && rts[0] < rt) {
							rt = rts[0];
							double* newRayIntersection = add(rayIntersection, multScalar(R, rt));
							Rlambda = calculateIllumination(newObject, rayIntersection, newRayIntersection, R, lightVector, objectVector, textureVector, newRayIntersection, bkg, currentObject);
						}
					} else {
						if (rts[1] > 0 && rts[1] < rt) {
							rt = rts[1];
							double* newRayIntersection = add(rayIntersection, multScalar(R, rt));
							Rlambda = calculateIllumination(newObject, rayIntersection, newRayIntersection, R, lightVector, objectVector, textureVector, newRayIntersection, bkg, currentObject);
						}
					}
					if (TIR) {
						Tlambda = Rlambda;
					}
					if (!TIR) {
						double* tts = raySphereIntersect(T, rayIntersection, newSphere);
						
						if (tts[0] < tts[1]) {
							if (tts[0] > 0 && tts[0] < tt) {
								tt = tts[0];
								double* newRayIntersection = add(rayIntersection, multScalar(T, tt));
								Tlambda = calculateIllumination(newObject, rayIntersection, newRayIntersection, T, lightVector, objectVector, textureVector, newRayIntersection, bkg, currentObject);
							}
						} else {
							if (tts[1] > 0 && tts[1] < tt) {
								tt = tts[1];
								double* newRayIntersection = add(rayIntersection, multScalar(T, tt));
								Tlambda = calculateIllumination(newObject, rayIntersection, newRayIntersection, T, lightVector, objectVector, textureVector, newRayIntersection, bkg, currentObject);
							}
						}
					}
				} else if (newObject->type == triangle) {
					Triangle* newTriangle = (Triangle*)newObject;
					double* newBarycentricR = rayTriangleIntersect(R, rayIntersection, newTriangle);
					if (newBarycentricR[0] > 0 && newBarycentricR[0] < rt) {
						rt = newBarycentricR[0];
						double* newRayIntersection = add(rayIntersection, multScalar(R, rt));
						Rlambda = calculateIllumination(newObject, rayIntersection, newRayIntersection, R, lightVector, objectVector, textureVector, newBarycentricR, bkg, currentObject);
					}
					if (TIR) {
						Tlambda = Rlambda;
					}
					if (!TIR) {
						double* newBarycentricT = rayTriangleIntersect(T, rayIntersection, newTriangle);
						if (newBarycentricT[0] > 0 && newBarycentricT[0] < tt) {
							tt = newBarycentricT[0];
							double* newRayIntersection = add(rayIntersection, multScalar(T, tt));
							Tlambda = calculateIllumination(newObject, rayIntersection, newRayIntersection, T, lightVector, objectVector, textureVector, newBarycentricT, bkg, currentObject);
						}
					}
				}
			}
		}
		
		reflection = add(multScalar(Rlambda, Fr), multScalar(Tlambda, (1-Fr)*(1-opacity)));
		limit = 0;
	}
	
	illumination = add(add(ambientTerm, sumOfLights), reflection);
	
	// Final check to avoid overflow
	for (int i = 0; i < 3; i++) {
		if (illumination[i] < 0) {
			illumination[i] = 0;
		}
		if (illumination[i] > 1) {
			illumination[i] = 1;
		}
	}
	
	return illumination;
}

// Add light from input file to lightVector
void Environment::addLight(Light* light) {
	lightVector.push_back(light);
}

// Add texture from input file to textureVector
void Environment::addTexture(Texture texture) {
	textureVector.push_back(texture);
}

// Add shape from input file to objectVector
void Environment::addObject(ShapeObject* object) {
	objectVector.push_back(object);
}

// Add vertex coordinates from input file to vertexCoordVector
void Environment::addVertexCoord(double* vertexCoord) {
	vertexCoordVector.push_back(vertexCoord);
}

// Add texture coordinates from input file to textureCoordVector
void Environment::addTextureCoord(double* textureCoord) {
	textureCoordVector.push_back(textureCoord);
}

// Add surface normal coordinates from input file to surfaceNormCoordVector
void Environment::addSurfaceNormCoord(double* surfaceNormCoord) {
	surfaceNormCoordVector.push_back(surfaceNormCoord);
}

// Computes all the calculations required to generate the correct scene
void Environment::generateScene() {
	double* u1 = cross(this->vdir, this->up);
	double u1Length = vectorLength(u1);
	double* u = new double[3];
	for (int i = 0; i < 3; i++) {
		u[i] = u1[i] / u1Length;
	}

	double* v1 = cross(u, this->vdir);
	double v1Length = vectorLength(v1);
	double* v = new double[3];
	for (int i = 0; i < 3; i++) {
		v[i] = v1[i] / v1Length;
	}

	double vdirLength = vectorLength(this->vdir);

	double* n = new double[3];
	for (int i = 0; i < 3; i++) {
		n[i] = this->vdir[i] / vdirLength;
	}

	// Define viewing window
	int d = this->viewDist;
	double aspect = (double)this->width / (double)this->height;
	double viewW = 2 * d * tan((this->fovh / 2) * PI / 180);
	double viewH = viewW / aspect;

	double* ul = new double[3];
	double* ur = new double[3];
	double* ll = new double[3];
	double* lr = new double[3];

	// Find viewing window corners
	if (this->parallel == 1) {
		// Parallel projection
		for (int i = 0; i < 3; i++) {
			ul[i] = this->eye[i] + (viewH / 2 * v[i]) - (viewW / 2 * u[i]);
			ur[i] = this->eye[i] + (viewH / 2 * v[i]) + (viewW / 2 * u[i]);
			ll[i] = this->eye[i] - (viewH / 2 * v[i]) - (viewW / 2 * u[i]);
			lr[i] = this->eye[i] - (viewH / 2 * v[i]) + (viewW / 2 * u[i]);
		}
	} else {
		// Perspective projection
		for (int i = 0; i < 3; i++) {
			ul[i] = this->eye[i] + (d * n[i]) + (viewH / 2 * v[i]) - (viewW / 2 * u[i]);
			ur[i] = this->eye[i] + (d * n[i]) + (viewH / 2 * v[i]) + (viewW / 2 * u[i]);
			ll[i] = this->eye[i] + (d * n[i]) - (viewH / 2 * v[i]) - (viewW / 2 * u[i]);
			lr[i] = this->eye[i] + (d * n[i]) - (viewH / 2 * v[i]) + (viewW / 2 * u[i]);
		}
	}

	// Calculate pixel offset in viewing window
	double* offsetH = new double[3];
	double* offsetV = new double[3];
	
	for (int i = 0; i < 3; i++) {
		offsetH[i] = (ur[i] - ul[i]) / (this->width - 1);
		offsetV[i] = (ll[i] - ul[i]) / (this->height - 1);
	}
	
	// Add viewing window coordinates to vector
	for (int i = 0; i < this->height; i++) {
		for (int j = 0; j < this->width; j++) {
			double* viewWindowCoordinate = add(multScalar(offsetV, i), multScalar(offsetH, j));
			viewWindow.push_back(add(ul, viewWindowCoordinate));
		}
	}
	
	double totalPixels = viewWindow.size();
	double totalObjects = objectVector.size();
	
	// Ray tracing and final pixel color output calculations
	for (int i = 0; i < totalPixels; i++) {
		double* pixelColor = this->bkg;
		double* rayDir;
		double t = 10000;
		if (this->parallel == 1) {
			// Parallel projection
			rayDir = this->vdir;
			this->eye[0] = viewWindow.at(i)[0];
			this->eye[1] = viewWindow.at(i)[1];
			this->eye[2] = viewWindow.at(i)[2];
		} else {
			// Perspective projection
			rayDir = normalize(subtract(viewWindow.at(i), this->eye));
		}
		
		for (int j = 0; j < totalObjects; j++) {
			ShapeObject* currentObject = objectVector.at(j);
			if (currentObject->type == sphere) {
				Sphere* currentSphere = (Sphere*)currentObject;
				double* ts = raySphereIntersect(rayDir, this->eye, currentSphere);
				
				if (ts[0] < ts[1]) {
					if (ts[0] > 0 && ts[0] < t) {
						t = ts[0];
						double* rayIntersection = add(this->eye, multScalar(rayDir, t));
						pixelColor = calculateIllumination(currentObject, this->eye, rayIntersection, rayDir, this->lightVector, this->objectVector, this->textureVector, rayIntersection, this->bkg, currentObject);
					}
				} else {
					if (ts[1] > 0 && ts[1] < t) {
						t = ts[1];
						double* rayIntersection = add(this->eye, multScalar(rayDir, t));
						pixelColor = calculateIllumination(currentObject, this->eye, rayIntersection, rayDir, this->lightVector, this->objectVector, this->textureVector, rayIntersection, this->bkg, currentObject);
					}
				}
			} else if (currentObject->type == triangle) {
				Triangle* currentTriangle = (Triangle*)currentObject;
				double* barycentric = rayTriangleIntersect(rayDir, this->eye, currentTriangle);
				if (barycentric[0] > 0 && barycentric[0] < t) {
					t = barycentric[0];
					double* rayIntersection = add(this->eye, multScalar(rayDir, t));
					pixelColor = calculateIllumination(currentObject, this->eye, rayIntersection, rayDir, this->lightVector, this->objectVector, this->textureVector, barycentric, this->bkg, currentObject);
				}
			}
		}
		outputPixels.push_back(pixelColor);
	}
}


