CC = g++
FLAGS = -std=c++0x
SOURCES = main.cpp Environment.cpp Light.cpp ShapeObject.cpp Sphere.cpp Triangle.cpp
HEADERS = Environment.h Light.h ShapeObject.h Sphere.h Triangle.h

all:
	$(CC) $(FLAGS) $(SOURCES) -o raytracer

clean: 
	rm raytracer
