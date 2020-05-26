#ifndef CURVATURE2THRESHOLD
#define CURVATURE2THRESHOLD

#define  MINN_MINN 1e-9
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
//#include <sstream>
#include <vector>
#include <map>
#include <iostream> 
#include "math3D.h"
#include "porosityPereservation.h"
using namespace std;
// -------------------- OpenMesh  
#include <OpenMesh/Core/IO/MeshIO.hh>  
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>  
// ----------------------------------------------------------------------------  
typedef OpenMesh::PolyMesh_ArrayKernelT<>  MyMesh;

struct para_and_value
{
	OpenMesh::Vec3d para;//参数坐标
	double value;//曲率值
	para_and_value(){
		para = OpenMesh::Vec3d(0, 0, 0);
		value = 0;
	}
};

void curvature2threshold(int typeTPMS, int typeStructure);
void read_bsolid(const std::string & filename);
void create_triangular_mesh(const int num_x, const int num_y, const int num_z);
void cal_mean_curvature(MyMesh*& mesh);
void discrete_curvature();
void zoom_TDF(int type_TPMS, int type_structure);
void scale_stretch(double begin_minC, double begin_maxC, double end_minC, double end_maxC);
void laplacianSmoothInner(para_and_value ***&cube, const int *seg, const int iter_num);
void laplacianSmoothOuter(para_and_value ***&cube, const int *seg, const int iter_num);
void LSPIACurvature(para_and_value ***&cube, int iterations);
void construct_TDF();
void save_TDF_file(int type_TPMS, int type_structure, string filename);

void calculate_physcial_point(Vector3D para, Vector3D &phy_pos);

#endif