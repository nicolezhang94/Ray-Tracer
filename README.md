# Ray Tracer

### I. PURPOSE ###

The purpose of this program is to create a scene including textured and non-textured
spheres and triangle-based polygons with flat or smooth shading based on a given input
text file using the Phong illumination model. Transparent and reflective objects can
also be created using recursive ray tracing. The output image will be in ASCII PPM format.


### II. COMPILATION ###

To compile this program, simply type "make" into a Linux machine terminal. Note: This
program uses C++11, so a compatible compiler is required. However, if an older compiler
is used, the program will still be able to run, as the makefile prepares for such a case
by using the -std=c++0x flag.


### III. RUNNING THE PROGRAM ###

To run the program, type "./raytracer" with a text file of correct syntax into the terminal
in the directory of the executable. For example, to run the provided sample text file, type
"./raytracer image1.txt" into the terminal.


### IV. WHAT THE PROGRAM DOES ###

This program takes in an text file and creates a scene based on the specified parameters.
The scene contains spheres, triangles, and lights that are specified by the input file. The
program uses ray casting to determine whether each ray originating from the viewing origin
that corresponds to a pixel in the output image is intersected by an object or not. If a ray
is not intersected by anything, the corresponding pixel will be the background color. If a
ray is intersected, then the corresponding pixel's color will be calculated using the Phong
illumination model based on the specified scene lights and other global objects. The program
uses the Phong illumination model and recursive ray tracing to calculate the pixel's color
based on the intersected object's material color, the global reflections and transparency
from other objects in the environment, and the properties of lights in the scene. The program
implements reflections and transparency by calculating the direction of reflection and
direction of transmission and recursively tracing a ray in those directions. The recursively
traced rays return a correct-computed intensity contribution using the Phong illumination
model. Schlick's approximation of the Fresnel reflectance is used to determine the relative
contributions of the colors returned by the reflected and transmitted rays to the computed
illumination at the ray/intersection point on a transparent surface. The program also checks
if the object is in the shadow of other objects in the scene, in which case the amount of
light blocked is calculated based on the opacity of the blocking objects. Triangle-based
polygons may also be added to the scene and is defined by its vertices on a plane. Objects
in the scene can be textured or non-textured with either flat or smooth shading. The program
takes in an ASCII PPM image in correct format and applies that to the specified textured
objects. Parallel projection may also be used instead of perspective projection.


### V. SPECIFICATIONS ###

The input file must include the minimum information to create a scene: the view origin,
viewing direction, horizontal field of view, the output image size, and a background color,
indicated by their respective keywords as follows: "eye", "viewdir", "updir", "fovh", 
"imsize", and "bkgcolor". Spheres may be added to the scene using the keyword "sphere".
Polygon vertices, texture coordinates, and surface normal coordinates may be added to the
scene using the keywords "v", "vt", and "vn", respectively. Triangles may be added to scene
after at least three vertices have been defined. Objects can be textured if a texture is
defined before it. The material color must specify the object's diffuse color, reflectivity,
specular highlight color, opacity, and index of refraction. Lights may be added to the scene
using the keyword "light", followed by specified light position/direction, light type
(positional or directional), and light color. All coordinates and vectors are in 3D space and
require an x, y, and z value. The input field of view is measured in degrees. Width and height
are measured in pixel units. Colors use the RGB color scale ranging from 0-1. Spheres require
x, y, and z coordinates along with a radius size. Non-unit length vectors will be normalized.
The "parallel" keyword turns on parallel projection. The following details the correct format
for each keyword in the input file:

```
eye        eyex eyey eyez
viewdir    vdirx vdiry vdirz
updir      upx upy upz
fovh       fovh
imsize     width height
bkgcolor   r g b
light      x y z w r g b
mtlcolor   odr odg odb osr osg osb ka kd ks n alpha eta
texture    texturefilename.ppm
sphere     cx cy cz r
v          xi yi zi
vt         ui vi
vn         nxi nyi nzi
f          v1 v2 v3
f          v1/vt1 v2/vt2 v3/vt3
f          v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3
f          v1//vn1 v2//vn2 v3//vn3
parallel
```

### VI. ASSUMPTIONS ###

1. The input file must be a text file.
2. The input file is in the same directory as the executable file, or a path must be provided.
3. Texture files must have one-line headers with no comments using the following ASCII PPM format:
   * P3 width height 255
   * r
   * g
   * b
   * r
   * g
   * b
   * ...
4. Triangle face input must not point to invalid/out-of-bounds vector indices.
5. Triangles may not be defined before sufficient appropriate vertices have been defined.
6. A material color must be specified before any objects are defined.
7. Texture files must exist and be in the same directory as the executable, or a path must be
   provided.


### VII. ERROR HANDLING ###

Any of the following circumstances will end the program:
1. Zero-length input vectors
2. Parallel up and viewing directions
3. Fields of view greater than or equal to 180 degrees
4. Non-positive image width or height
5. Incomplete or invalid scene or object definitions
6. Creating objects without specifying a material color first
7. RGB values not within 0-1 range


### VIII. NOTES ###

1. Occasionally, files with certain invalid artifacts will cause the program to end if there
   are spaces or empty lines in the file.
2. Total internal reflection does not always work.
3. Bug fixed: Changing image size dimensions does not shrink image.
4. Bug fixed: Surface normals are correct now.
5. Bug fixed: Adding lights will no longer make the scene darker.
