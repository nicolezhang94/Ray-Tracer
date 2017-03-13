/* Nicole Zhang 
 * CSCI 5607 Assignment 1d
 * 11/05/2015
 */

#include "Environment.h"
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char** argv) {
	if (argc <= 1) {
		printf("Please specify input text file.\n");
		return 1;
	}

	// Opens input file for reading
	ifstream inFile;
	inFile.open(argv[1]);

	string line;

	Environment* env = new Environment();

	int objectID = 0;
	int lightID = 0;
	int currentTriangle = 0;
	int textureIndex = -1;
	int materialIndex = -1;

	double* currentColor = new double[12]{-1.0};

	// Reads the input file, line by line, and makes sure the input is valid
	while (getline(inFile, line)) {
		istringstream lineStream(line);
		string word;
		lineStream >> word;
		if (word.compare("eye") == 0) {
			int index = 0;
			while (lineStream >> word) {
				if (index > 2) {
					printf("Invalid view origin: too many arguments\n");
					return 1;
				}
				env->eye[index] = stod(word);
				index++;
			}
			if (index < 2) {
				printf("Invalid view origin: too few arguments\n");
				return 1;
			}
		} else if (word.compare("viewdir") == 0) {
			int index = 0;
			while (lineStream >> word) {
				if (index > 2) {
					printf("Invalid viewing direction: too many arguments\n");
					return 1;
				}
				env->vdir[index] = stod(word);
				index++;
			}
			if (index < 2) {
				printf("Invalid viewing direction: too few arguments\n");
				return 1;
			}
		} else if (word.compare("updir") == 0) {
			int index = 0;
			while (lineStream >> word) {
				if (index > 2) {
					printf("Invalid up direction: too many arguments\n");
					return 1;
				}
				env->up[index] = stod(word);
				index++;
			}
			if (index < 2) {
				printf("Invalid up direction: too few arguments\n");
				return 1;
			}
		} else if (word.compare("fovh") == 0) {
			int index = 0;
			while (lineStream >> word) {
				if (index > 1) {
					printf("Invalid field of view: too many arguments\n");
					return 1;
				}
				env->fovh = stod(word);
				index++;
			}
			if (env->fovh >= 180 || env->fovh <= 0) {
				printf("Invalid field of view: value out of range\n");
				return 1;
			}
			if (index < 1) {
				printf("Invalid field of view: too few arguments\n");
				return 1;
			}
		} else if (word.compare("imsize") == 0) {
			int index = 0;
			while (lineStream >> word) {
				if (index > 1) {
					printf("Invalid image size: too many arguments\n");
					return 1;
				}
				if (index == 0) {
					env->width = stoi(word);
				} else if (index == 1) {
					env->height = stoi(word);
				}
				index++;
			}
			if (env->width <= 0 || env->height <= 0) {
				printf("Invalid image size: values must be positive\n");
			}
			if (index < 1) {
				printf("Invalid image size: too few arguments\n");
				return 1;
			}
		} else if (word.compare("bkgcolor") == 0) {
			int index = 0;
			while (lineStream >> word) {
				if (index > 2) {
					printf("Invalid background color: too many arguments\n");
					return 1;
				}
				env->bkg[index] = stod(word);
				if (env->bkg[index] < 0 || env->bkg[index] > 1) {
					printf("Invalid background color: value out of range\n");
					return 1;
				}
				index++;
			}
			if (index < 2) {
				printf("Invalid background color: too few arguments\n");
				return 1;
			}
		} else if (word.compare("light") == 0) {
			Light* newLight = new Light();
			int index = 0;
			while (lineStream >> word) {
				if (index > 6) {
					printf("Invalid light: too many arguments\n");
					return 1;
				}
				if (index < 3) {
					newLight->location[index] = stod(word);
				} else if (index == 3) {
					if (stoi(word) != 1 && stoi(word) != 0) {
						printf("Invalid light type: w value must be 1 or 0\n");
						return 1;
					}
					newLight->type = stoi(word);
				} else {
					newLight->color[index-4] = stod(word);
					if (newLight->color[index-4] < 0 || newLight->color[index-4] > 1) {
						printf("Invalid light color: value out of range\n");
						return 1;
					}
				}
				index++;
			}
			if (index < 6) {
				printf("Invalid light: too few arguments\n");
				return 1;
			}
			newLight->ID = lightID;
			env->addLight(newLight);
			lightID++;
		} else if (word.compare("mtlcolor") == 0) {
			materialIndex++;
			int index = 0;
			while (lineStream >> word) {
				if (index > 11) {
					printf("Invalid material color: too many arguments\n");
					return 1;
				}
				currentColor[index] = stod(word);
				if (index < 9 || index == 10) {
					if (currentColor[index] < 0 || currentColor[index] > 1) {
						printf("Invalid material color: value out of range\n");
						return 1;
					}
				} else if (index == 9 || index == 11) {
					if (currentColor[index] < 0) {
						printf("Invalid material color: value out of range\n");
						return 1;
					}
				}
				index++;
			}
			if (index < 11) {
				printf("Invalid material color: too few arguments\n");
				return 1;
			}
		} else if (word.compare("texture") == 0) {
			textureIndex++;
			
			ifstream inTexture;
			lineStream >> word;
			inTexture.open(word);
			string textureLine;
			
			getline(inTexture, textureLine);
			istringstream textureLineStream(textureLine);
			string header;
			textureLineStream >> header;
			textureLineStream >> header;
			int textureWidth = stoi(header);
			textureLineStream >> header;
			int textureHeight = stoi(header);
			
			Texture newTexture;
			newTexture.width = textureWidth;
			newTexture.height = textureHeight;
			
			int texIndex = 0;
			double* currentTextureColor = new double[3];
			
			while (getline(inTexture, textureLine)) {
				currentTextureColor[texIndex] = stod(textureLine) / 255.0;
				texIndex++;
				if (texIndex == 3) {
					newTexture.colors.push_back(currentTextureColor);
					texIndex = 0;
					currentTextureColor = new double[3];
				}
			}
			inTexture.close();
			
			env->addTexture(newTexture);
		} else if (word.compare("v") == 0) {
			double* newVertexCoord = new double[3];
			int index = 0;
			while (lineStream >> word) {
				if (index > 3) {
					printf("Invalid triangle vertex: too many arguments\n");
					return 1;
				}
				newVertexCoord[index] = stod(word);
				index++;
			}
			if (index < 3) {
				printf("Invalid triangle vertex: too few arguments\n");
				return 1;
			}
			env->addVertexCoord(newVertexCoord);
		} else if (word.compare("vt") == 0) {
			double* newTextureCoord = new double[2];
			int index = 0;
			while (lineStream >> word) {
				if (index > 2) {
					printf("Invalid triangle texture coordinates: too many arguments\n");
					return 1;
				}
				newTextureCoord[index] = stod(word);
				index++;
			}
			if (index < 2) {
				printf("Invalid triangle texture coordinates: too few arguments\n");
				return 1;
			}
			env->addTextureCoord(newTextureCoord);
		} else if (word.compare("vn") == 0) {
			double* newSurfaceNormCoord = new double[3];
			int index = 0;
			while (lineStream >> word) {
				if (index > 3) {
					printf("Invalid triangle surface normal directions: too many arguments\n");
					return 1;
				}
				newSurfaceNormCoord[index] = stod(word);
				index++;
			}
			if (index < 3) {
				printf("Invalid triangle surface normal directions: too few arguments\n");
				return 1;
			}
			env->addSurfaceNormCoord(newSurfaceNormCoord);
		} else if (word.compare("f") == 0) {
			Triangle* newTriangle = new Triangle();
			string v1, v2, v3;
			lineStream >> word;
			v1 = word;
			lineStream >> word;
			v2 = word;
			lineStream >> word;
			v3 = word;
			
			string delimiter1 = "/";
			string delimiter2 = "//";
			int faceIndicator = v1.find(delimiter1);
			int faceIndicator2 = v1.find(delimiter2);
			if (faceIndicator == string::npos) {
				// For face input like v1
				newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v1)-1));
				newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v2)-1));
				newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v3)-1));
				newTriangle->shading = flat;
				newTriangle->texture = -1;
			} else if (faceIndicator2 == string::npos) {
				int inputIndex1 = v1.find(delimiter1);
				int inputIndex2 = v2.find(delimiter1);
				int inputIndex3 = v3.find(delimiter1);
				if (inputIndex1 == v1.find_last_of(delimiter1)) {
					// For face input like v1/vt1
					newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v1.substr(0, inputIndex1))-1));
					newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v2.substr(0, inputIndex2))-1));
					newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v3.substr(0, inputIndex3))-1));
					newTriangle->vt.push_back(env->textureCoordVector.at(stoi(v1.substr(inputIndex1+1, v1.length()))-1));
					newTriangle->vt.push_back(env->textureCoordVector.at(stoi(v2.substr(inputIndex2+1, v2.length()))-1));
					newTriangle->vt.push_back(env->textureCoordVector.at(stoi(v3.substr(inputIndex3+1, v3.length()))-1));
					newTriangle->shading = flat;
					newTriangle->texture = textureIndex;
				} else {
					// For face input like v1/vt1/vn1
					newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v1.substr(0, inputIndex1))-1));
					newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v2.substr(0, inputIndex2))-1));
					newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v3.substr(0, inputIndex3))-1));
					newTriangle->vt.push_back(env->textureCoordVector.at(stoi(v1.substr(inputIndex1+1, v1.find(delimiter1, inputIndex1+1)))-1));
					newTriangle->vt.push_back(env->textureCoordVector.at(stoi(v2.substr(inputIndex2+1, v2.find(delimiter1, inputIndex2+1)))-1));
					newTriangle->vt.push_back(env->textureCoordVector.at(stoi(v3.substr(inputIndex3+1, v3.find(delimiter1, inputIndex3+1)))-1));
					newTriangle->vn.push_back(env->surfaceNormCoordVector.at(stoi(v1.substr(v1.find(delimiter1, inputIndex1+1)+1, v1.length()))-1));
					newTriangle->vn.push_back(env->surfaceNormCoordVector.at(stoi(v2.substr(v2.find(delimiter1, inputIndex2+1)+1, v2.length()))-1));
					newTriangle->vn.push_back(env->surfaceNormCoordVector.at(stoi(v3.substr(v3.find(delimiter1, inputIndex3+1)+1, v3.length()))-1));
					newTriangle->shading = smooth;
					newTriangle->texture = textureIndex;
				}
			} else {
				// For face input like v1//vn1
				newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v1.substr(0, v1.find(delimiter2)))-1));
				newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v2.substr(0, v2.find(delimiter2)))-1));
				newTriangle->v.push_back(env->vertexCoordVector.at(stoi(v3.substr(0, v3.find(delimiter2)))-1));
				newTriangle->vn.push_back(env->surfaceNormCoordVector.at(stoi(v1.substr(v1.find(delimiter2)+2, v1.length()-1))-1));
				newTriangle->vn.push_back(env->surfaceNormCoordVector.at(stoi(v2.substr(v2.find(delimiter2)+2, v2.length()-1))-1));
				newTriangle->vn.push_back(env->surfaceNormCoordVector.at(stoi(v3.substr(v3.find(delimiter2)+2, v3.length()-1))-1));
				newTriangle->shading = smooth;
				newTriangle->texture = -1;
			}
			if (materialIndex > -1) {
				newTriangle->mtld[0] = currentColor[0];
				newTriangle->mtld[1] = currentColor[1];
				newTriangle->mtld[2] = currentColor[2];
				newTriangle->mtls[0] = currentColor[3];
				newTriangle->mtls[1] = currentColor[4];
				newTriangle->mtls[2] = currentColor[5];
				newTriangle->mtlk[0] = currentColor[6];
				newTriangle->mtlk[1] = currentColor[7];
				newTriangle->mtlk[2] = currentColor[8];
				newTriangle->n = currentColor[9];
				newTriangle->opacity = currentColor[10];
				newTriangle->refraction = currentColor[11];
			} else {
				printf("A material color must be specified before creating objects.\n");
				return 1;
			}
			newTriangle->ID = objectID;
			newTriangle->type = triangle;
			env->addObject(newTriangle);
			objectID++;
		} else if (word.compare("sphere") == 0) {
			Sphere* newSphere = new Sphere();
			int index = 0;
			while (lineStream >> word) {
				if (index > 3) {
					printf("Invalid sphere: too many arguments\n");
					return 1;
				}
				if (index < 3) {
					newSphere->location[index] = stod(word);
				} else if (index == 3) {
					newSphere->r = stod(word);
				}
				index++;
			}
			if (newSphere->location[3] < 0) {
				printf("Invalid sphere: radius must be non-negative\n");
			}
			if (index < 3) {
				printf("Invalid sphere: too few arguments\n");
				return 1;
			}
			if (materialIndex > -1) {
				newSphere->mtld[0] = currentColor[0];
				newSphere->mtld[1] = currentColor[1];
				newSphere->mtld[2] = currentColor[2];
				newSphere->mtls[0] = currentColor[3];
				newSphere->mtls[1] = currentColor[4];
				newSphere->mtls[2] = currentColor[5];
				newSphere->mtlk[0] = currentColor[6];
				newSphere->mtlk[1] = currentColor[7];
				newSphere->mtlk[2] = currentColor[8];
				newSphere->n = currentColor[9];
				newSphere->opacity = currentColor[10];
				newSphere->refraction = currentColor[11];
			} else {
				printf("A material color must be specified before creating objects.\n");
				return 1;
			}
			newSphere->ID = objectID;
			newSphere->texture = textureIndex;
			newSphere->type = sphere;
			env->addObject(newSphere);
			objectID++;
		} else if (word.compare("parallel") == 0) {
			env->parallel = 1;
		} else if (word.compare("viewdist") == 0) {
			//lineStream >> word;
			//env->viewDist = stod(word);
		} else if (line.compare("") == 0 || line.compare(" ") == 0 || word.compare("#") == 0) {
			// Do nothing on empty lines or lines with comments
		} else {
			printf("Invalid keyword in input file.\n");
			return 1;
		}
	}
	
	// Makes sure the up direction and viewing direction are not equal
	int parallelCheck = 0;
	for (int i = 0; i < 3; i++) {
		if (env->up[i] == env->vdir[i]) {
			parallelCheck++;
		}
	}
	if (parallelCheck > 1) {
		printf("Invalid input: up direction and viewing direction cannot be parallel (or corrupt file).\n");
		return 1;
	}

	// Makes call to create scene and make calculations
	env->generateScene();

	vector<double*> pixels = env->outputPixels;
	int totalPixels = pixels.size();
 
	// Gets the raw input file name without the extension
	string inFileString = argv[1];
	int inFileNameEnd = inFileString.find_first_of(".");
	string inFileName = inFileString.substr(0, inFileNameEnd);

	inFile.close();

	// Creates output ppm file with the same name of the input file
	ofstream outFile;
	string outFileName = inFileName.append(".ppm");
	outFile.open(outFileName);

	// Writes the ppm image data into the output file
	outFile << "P3 " << env->width << " " << env->height << " 255" << endl;
	for (int i = 0; i < totalPixels; i++) {
		for (int j = 0; j < 3; j++) {
			double pixel = pixels.at(i)[j] * 255;
			int pixelScaled = (int)pixel;
			outFile << pixelScaled << endl;
		}
	}
	outFile.close();

	delete[] currentColor;
	delete env;

	return 0;
}
