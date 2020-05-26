#ifndef SCAFFOLDMESH_H
#define SCAFFOLDMESH_H
#include <cassert>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <string>
#include "Isosurface.h"
#include "Array3D.h"
using namespace std;

#define  MINN_MINN 1e-9

class ScaffoldMesh
{
public:
	ScaffoldMesh(string Filename);
	~ScaffoldMesh();
	int type_TPMS;//TPMS类型（P,G,D,IWP）
	int type_structure;//结构类型（Rod,Pore,Sheet）
	vector<Vector3D> Vertex;
	vector<Vector3D> VertexNormal;
	vector<Triangle> Face;
	vector<Vector3D> FaceNormal;
	Vector3D center;
	double diagonal_length;
	size_t resolution;

	void scaffold_generation();

private:
	int omega_x, omega_y, omega_z;//三个方向上的周期
	int u_points, v_points, w_points;//TDF控制网格数
	vector<vector<vector<double>>> TDF_control_grid;
	vector<vector<double>> TDF_knot_vector;
	int x_points, y_points, z_points;//TBSS控制网格数
	vector<vector<vector<Vector3D>>> TBSS_control_grid;
	vector<vector<double>> TBSS_knot_vector;

	map<int, vector<int>> vec_face_map;


	void get_model_size_center(double& l,Vector3D& v);
	void choose_scaffold_type();
	void structure_generation(const Isosurface& surface, float xMin, float xMax, float yMin, float yMax, float zMin, float zMax, float isolevel);
	void structure_mapping_TBSS();
	void cal_mesh_normal();
	void clear();

	void zoom_TDF();
	void scale_stretch(double begin_minC, double begin_maxC, double end_minC, double end_maxC);
	void discreate_TDF(double*** &C);
	void calculate_cvalue_TDF(Vector3D para, double &Cvalue);
	void calculate_posistion_TBSS(Vector3D para, Vector3D &phy_pos);
	void calculate_cvalue_TDF_quick(Vector3D para, double &Cvalue);
	void calculate_posistion_TBSS_quick(Vector3D para, Vector3D &phy_pos);
	int findSpan(int n, double u, std::vector<double> U);
	double* basisFuns(int i, double u, std::vector<double> U);
	void calculateBaseFunctionOfSolidBSpline(Vector3D para, std::vector<std::vector<double> >&kont_vector, double(*result)[4][4], int &Bi_start_index, int &Bj_start_index, int &Bk_start_index);
	double getNi4(const std::vector<double>&knot_vector, double t, int i, int seg, int zone, std::string  &deb);
	void constructEdge(Array3D<Point3D>& PointGrid, vector<edge>& Edge);

	void cutVertex(Vector3D& _point);
	void calInterPoint(const Isosurface& surface, const Point3D& p1, const Point3D& p2, float isolevel, Point3D& interP);
	void calInterVert(const Isosurface& surface, const Point3D& p1, const Point3D& p2, float isolevel, Vector3D& interP);
	void calInterVert2(const Isosurface& surface, const int edge_id, float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, int& index, int sur_id);
	void extract_isosurface(const Isosurface& surface, float isolevel, int type, vector<edge>& Edge, vector<Point3D>& pointset, Array3D<Point3D>& PointGrid, int sur_id);
	void extractTriangleFromTetrahedron(const Isosurface& surface, const Point3D p[4], const int edges[6], float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, int sur_id);
	void close_isosurface_network(const Isosurface& surface, float isolevel, int type, vector<edge>& Edge, vector<Point3D>& pointset, Array3D<Point3D>& PointGrid, int sur_id);
	void extractTriangleOuter2(const Isosurface& surface, const Point3D p[3], const int edges[3], const int Pindex[3], float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, int sur_id);
	void extractTriangleInner2(const Isosurface& surface, const Point3D p[3], const int edges[3], const int Pindex[3], float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, int sur_id);
	void getPointinsurface(const int id, vector<Point3D>& pointset, int& index);
	void close_isosurface_sheet(const Isosurface& surface, float maxisolevel, float minisolevel, int type, vector<edge>& Edge, vector<Point3D>& pointset, Array3D<Point3D>& PointGrid);
	void extractTriangleFromSurface(const Isosurface& surface, const Point3D p[3], const int edges[3], const int Pindex[3], float maxisolevel, float minisolevel, vector<Point3D>& pointset, vector<edge>& Edge);
	void addFace(int f1, int f2, int f3);

	void alter_face_direction();

	void predictPorosity();
	void outPutPorosity(const Isosurface& surface, Array3D<Point3D>& phyPorGrid, Array3D<Point3D>& errorPorGrid, int dim);
	void getPorosity(const Isosurface& surface, Vector3D& pos, Point3D& point_por);
	void calValidVolume(Point3D v[4], double& volume, double& validvolume);

	double** dersBasisFuns(int i, double u, std::vector<double> nodeU);
	Vector3D* solidDerivisOneOrder(Vector3D para);
	void detJacobian(Array3D<double>& JacobianGrid);
	void phyPosVolume(Array3D<double>& VolumeGrid);
	void calVolume(Point3D v[4], double& volume);
	void setPhyPorosity(Array3D<Point3D>& phyPorGrid, Array3D<Point3D>& errorPorGrid, int dim);
	void LSPIA(Array3D<Point3D>& DisTDFGrid, int dim);
	void phyPorosity2TDF(Array3D<double>& VolumeGrid, Array3D<double>& JacobianGrid, Array3D<Point3D>& phyPorGrid, Array3D<Point3D>& errorPorGrid, Array3D<Point3D>& DisTDFGrid, int dim);


	template<typename T>
	static void allocMem(const int x_points, const int y_points, const int z_points, T ***&para_cube, T values);


	double randomUniform(double dMinValue, double dMaxValue);
	void laplacianSmoothInner(double ***&cube, const int iter_num);
	void laplacianSmoothOuter(double ***&cube, const int iter_num);
	void save_porosity(double*** &C, string filename);
};

#endif // OPENGLWIDGET_H