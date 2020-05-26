#ifndef EXTRACTISOSURFACE
#define EXTRACTISOSURFACE

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
//#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <GL/glut.h>
#include <cstdlib>
#include "Isosurface.h"
#include "Array3D.h"
#include <Eigen/Dense>
#include <Eigen/LU>
using namespace std;
//
// decimate
//
// Draws a bunch of GL_TRIANGLES (with correct normals and vertex order)
// in order to render the given surface at the given isolevel.
//

void solid_construction_Bspline(const Isosurface& surface, float xMin, float xMax, float yMin, float yMax, float zMin, float zMax, float isolevel, size_t resolution);
void solid_construction_Bspline_save(const Isosurface& surface, float xMin, float xMax, float yMin, float yMax, float zMin, float zMax, float isolevel, size_t resolution,string s);
void solid_construction_Tspline(const Isosurface& surface, float isolevel, size_t resolution);
#endif
