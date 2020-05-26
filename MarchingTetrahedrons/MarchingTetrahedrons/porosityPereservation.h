#ifndef POROSITYPERESERVATION
#define POROSITYPERESERVATION

#define  MINN_MINN 1e-9

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
#include <time.h>
using namespace std;

void read_bspline_solid(string & filename);
void setPorosity(vector<Point3D>& porosity, vector<Point3D>& Cvalue);
void constructKnotVector(const int x_point_num, const int y_point_num, const int z_point_num, vector<vector<double >>&kont_vector);
void tdf_construction();
template<typename T>
void allocMem(const int x_points, const int y_points, const int z_points, T ***&para_cube, T values);
double getNi4(const std::vector<double>&knot_vector, double t, int i, int seg, int zone, std::string  &deb);
void calculateBaseFunctionOfSolidBSpline(Vector3D para, std::vector<std::vector<double> >&kont_vector, double(*result)[4][4], int &Bi_start_index, int &Bj_start_index, int &Bk_start_index);
void LSPIAPorosity(vector<Point3D> required_porosity, int iterations);
void calculate_cvalue_TDF_quick(Vector3D para, double &Cvalue);
void calculate_posistion_TBSS_quick(Vector3D para, Vector3D &phy_pos);

void cal_tetVolume(Point3D v[4], double& volume);
double cal_volumeGrid(Array3D<Point3D>& phyGrid, Array3D<double>& volumeGrid, size_t dim);
void save_tdf_file(string filename, int type_TPMS, int type_structure);
void cal_all_base_fun(vector<Point3D>& required_porosity, int resolution, double step[3]);
double cal_porosity(int type_TPMS, int type_structure, const Isosurface& surface, float isolevel, vector<Point3D>& required_porosity, vector<Array3D<double>>& phy_volume, vector<double>& volume, int id, int resolution, double step[3]);
void cal_partial_energy(int type_TPMS, int type_structure, const Isosurface& surface, float isolevel, vector<Point3D>& required_porosity, vector<Array3D<double>>& phy_volume, vector<double>& volume, int resolution, double step[3]);
void porosity_pereservation(const Isosurface& surface, float isolevel, int type_TPMS, int type_structure);
#endif