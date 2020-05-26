#include "ScaffoldMesh.h" 
#include "Psurface.h"
#include "Gsurface.h"
#include "Dsurface.h"
#include "I_WPsurface.h"
#include <QDebug>

ScaffoldMesh::ScaffoldMesh(string Fname)
{
	type_TPMS = 0;
	type_structure = 0;
	//resolution = 50;
	resolution = 80;
	//resolution = 200;

	clear();

	ifstream meshfile(Fname);
	string sline;
	getline(meshfile, sline);
	getline(meshfile, sline);
	istringstream sin(sline);
	sin >> type_TPMS;

	getline(meshfile, sline);
	getline(meshfile, sline);
	sin.clear();
	sin.str(sline);
	sin >> type_structure;
	
	getline(meshfile, sline);
	getline(meshfile, sline);
	sin.clear();
	sin.str(sline);
	sin >> omega_x >> omega_y >> omega_z;
	//qDebug("debug=%d", omega_x);

	getline(meshfile, sline);
	getline(meshfile, sline);
	sin.clear();
	sin.str(sline);
	sin >> u_points >> v_points >> w_points;

	TDF_control_grid.resize(u_points);
	for (int i = 0; i < u_points; i++)
	{
		TDF_control_grid[i].resize(v_points);
		for (int j = 0; j < v_points; j++)
		{
			TDF_control_grid[i][j].resize(w_points);
		}
	}
	TDF_knot_vector.resize(3);
	TDF_knot_vector[0].resize(u_points + 4);
	TDF_knot_vector[1].resize(v_points + 4);
	TDF_knot_vector[2].resize(w_points + 4);

	getline(meshfile, sline);
	for (int i = 0; i < u_points; i++)
	{
		for (int j = 0; j < v_points; j++)
		{
			for (int k = 0; k < w_points; k++)
			{
				getline(meshfile, sline);
				sin.clear();
				sin.str(sline);
				sin >> TDF_control_grid[i][j][k];
			}
		}
	}

	getline(meshfile, sline);
	getline(meshfile, sline);
	sin.clear();
	sin.str(sline);
	for (int i = 0; i < u_points + 4; i++)
	{
		sin >> TDF_knot_vector[0][i];
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	sin.clear();
	sin.str(sline);
	for (int i = 0; i < v_points + 4; i++)
	{
		sin >> TDF_knot_vector[1][i];
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	sin.clear();
	sin.str(sline);
	for (int i = 0; i < w_points + 4; i++)
	{
		sin >> TDF_knot_vector[2][i];
	}

	getline(meshfile, sline);
	getline(meshfile, sline);
	sin.clear();
	sin.str(sline);
	sin >> x_points >> y_points >> z_points;
	TBSS_control_grid.resize(x_points);
	for (int i = 0; i < x_points; i++)
	{
		TBSS_control_grid[i].resize(y_points);
		for (int j = 0; j < y_points; j++)
		{
			TBSS_control_grid[i][j].resize(z_points);
		}
	}
	TBSS_knot_vector.resize(3);
	TBSS_knot_vector[0].resize(x_points + 4);
	TBSS_knot_vector[1].resize(y_points + 4);
	TBSS_knot_vector[2].resize(z_points + 4);

	getline(meshfile, sline);
	for (int i = 0; i < x_points; i++)
	{
		for (int j = 0; j < y_points; j++)
		{
			for (int k = 0; k < z_points; k++)
			{
				getline(meshfile, sline);
				sin.clear();
				sin.str(sline);
				sin >> TBSS_control_grid[i][j][k].x >> TBSS_control_grid[i][j][k].y >> TBSS_control_grid[i][j][k].z;
			}
		}
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	sin.clear();
	sin.str(sline);
	for (int i = 0; i < x_points + 4; i++)
	{
		sin >> TBSS_knot_vector[0][i];
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	sin.clear();
	sin.str(sline);
	for (int i = 0; i < y_points + 4; i++)
	{
		sin >> TBSS_knot_vector[1][i];
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	sin.clear();
	sin.str(sline);
	for (int i = 0; i < z_points + 4; i++)
	{
		sin >> TBSS_knot_vector[2][i];
	}

}

ScaffoldMesh::~ScaffoldMesh()
{
	clear();
}

void ScaffoldMesh::scaffold_generation()
{
	get_model_size_center(diagonal_length, center);
	//zoom_TDF();//放缩TDF至合理区间

	//predictPorosity();

	//clock_t time_start, time_end,time_end1;//计算时间
	//time_start = clock();

	choose_scaffold_type();

	//time_end = clock();
	//alter_face_direction();
	structure_mapping_TBSS();

	//time_end1 = clock();

	//double duration = (static_cast<double>(time_end - time_start));//返回毫秒
	//double duration1 = (static_cast<double>(time_end1 - time_end));//返回毫秒
	//std::ofstream ss_time("time.txt", std::ios::app);
	//ss_time << "structure: " << duration << " ms" << "\n";
	//ss_time << "scaffold: " << duration1 << " ms" << "\n";
	//ss_time.close();

	cal_mesh_normal();
}

void ScaffoldMesh::get_model_size_center(double& l,Vector3D& v)
{
	v = Vector3D(0.0, 0.0, 0.0);
	double min_x, min_y, min_z, max_x, max_y, max_z;
	min_x = max_x = TBSS_control_grid[0][0][0].x;
	min_y = max_y = TBSS_control_grid[0][0][0].y;
	min_z = max_z = TBSS_control_grid[0][0][0].z;
	for (int i = 0; i < x_points; i++)
	{
		for (int j = 0; j < y_points; j++)
		{
			for (int k = 0; k < z_points; k++)
			{
				v += TBSS_control_grid[i][j][k];
				min_x = min_x < TBSS_control_grid[i][j][k].x ? min_x : TBSS_control_grid[i][j][k].x;
				max_x = max_x > TBSS_control_grid[i][j][k].x ? max_x : TBSS_control_grid[i][j][k].x;
				min_y = min_y < TBSS_control_grid[i][j][k].y ? min_y : TBSS_control_grid[i][j][k].y;
				max_y = max_y > TBSS_control_grid[i][j][k].y ? max_y : TBSS_control_grid[i][j][k].y;
				min_z = min_z < TBSS_control_grid[i][j][k].z ? min_z : TBSS_control_grid[i][j][k].z;
				max_z = max_z > TBSS_control_grid[i][j][k].z ? max_z : TBSS_control_grid[i][j][k].z;
			}
		}
	}
	v /= (x_points*y_points*z_points);
	l = sqrt(pow(max_x - min_x, 2) + pow(max_y - min_y, 2) + pow(max_z - min_z, 2));
}

void ScaffoldMesh::choose_scaffold_type()
{
	if (type_TPMS == 0)
	{
		Psurface surface;
		structure_generation(surface, 0, 1, 0, 1, 0, 1, 0);
	}
	else if (type_TPMS == 1)
	{
		Gsurface surface;
		structure_generation(surface, 0, 1, 0, 1, 0, 1, 0);
	}
	else if (type_TPMS == 2)
	{
		Dsurface surface;
		structure_generation(surface, 0, 1, 0, 1, 0, 1, 0);
	}
	else if (type_TPMS == 3)
	{
		I_WPsurface surface;
		structure_generation(surface, 0, 1, 0, 1, 0, 1, 0);
	}
}

void ScaffoldMesh::structure_generation(const Isosurface& surface, float xMin, float xMax, float yMin, float yMax, float zMin, float zMax, float isolevel)
{
	vector<edge> Edge;
	vector<Point3D> pointset;
	size_t pointRes = resolution + 1; // indicates the # of points per side
	Array3D<Point3D> PointGrid(pointRes, pointRes, pointRes);
	float xrange = xMax - xMin;
	float yrange = yMax - yMin;
	float zrange = zMax - zMin;
	double*** C;
	discreate_TDF(C);
	for (size_t i = 0; i <= resolution; ++i) {
		float x = (float)i / resolution * xrange + xMin;
		for (size_t j = 0; j <= resolution; ++j) {
			float y = (float)j / resolution * yrange + yMin;
			for (size_t k = 0; k <= resolution; ++k) {
				float z = (float)k / resolution * zrange + zMin;
				float value = surface.valueAt(omega_x * x, omega_y * y, omega_z * z) - C[i][j][k];
				PointGrid.set(i, j, k, Point3D(x, y, z, value));
				pointset.push_back(Point3D(x, y, z, value));
			}
		}
	}
	constructEdge(PointGrid, Edge);
	if (type_structure != 2 && type_TPMS != 3)
	{
		int type = type_structure;
		extract_isosurface(surface, isolevel, type, Edge, pointset, PointGrid, 0);
		close_isosurface_network(surface, isolevel,type, Edge, pointset, PointGrid, 0);
	}
	else if (type_structure != 2 && type_TPMS == 3)
	{
		int type = 1 - type_structure;
		extract_isosurface(surface, isolevel, type, Edge, pointset, PointGrid, 0);
		close_isosurface_network(surface, isolevel,type, Edge, pointset, PointGrid, 0);
	}
	else if (type_structure == 2)
	{
		int type = 0;
		if (type_TPMS == 3)
		{
			type = 1;
		}

		extract_isosurface(surface, isolevel, type, Edge, pointset, PointGrid, 0);
		extract_isosurface(surface, isolevel - 0.2, type, Edge, pointset, PointGrid, 1);
		close_isosurface_sheet(surface, isolevel, isolevel - 0.2, type,Edge, pointset, PointGrid);
	}
}

void ScaffoldMesh::structure_mapping_TBSS()
{
	Vector3D phy_pos;
	for (int i = 0; i < Vertex.size(); i++)
	{
		//calculate_posistion_TBSS(Vertex[i], phy_pos);
		calculate_posistion_TBSS_quick(Vertex[i], phy_pos);
		Vertex[i] = phy_pos;
	}
}

void ScaffoldMesh::cal_mesh_normal()
{
	Vector3D result;
	double x, y, z, sum;
	double v1v2x, v1v2y, v1v2z, v1v3x, v1v3y, v1v3z;
	if ((type_TPMS != 3 && type_structure == 1) || (type_TPMS == 3 && type_structure != 1))
	{
		for (int i = 0; i < Face.size(); i++)
		{
			vec_face_map[Face[i].f1].push_back(i);
			vec_face_map[Face[i].f2].push_back(i);
			vec_face_map[Face[i].f3].push_back(i);
		}

		for (int i = 0; i < Face.size(); i++)
		{
			v1v2x = Vertex[Face[i].f2].x - Vertex[Face[i].f1].x;
			v1v2y = Vertex[Face[i].f2].y - Vertex[Face[i].f1].y;
			v1v2z = Vertex[Face[i].f2].z - Vertex[Face[i].f1].z;
			v1v3x = Vertex[Face[i].f3].x - Vertex[Face[i].f1].x;
			v1v3y = Vertex[Face[i].f3].y - Vertex[Face[i].f1].y;
			v1v3z = Vertex[Face[i].f3].z - Vertex[Face[i].f1].z;
			x = v1v2y * v1v3z - v1v2z * v1v3y;
			y = v1v2z * v1v3x - v1v2x * v1v3z;
			z = v1v2x * v1v3y - v1v2y * v1v3x;
			//归一化
			sum = sqrt(x*x + y*y + z*z);
			result.x = x / sum;
			result.y = y / sum;
			result.z = z / sum;
			FaceNormal.push_back(result);
		}
	}
	else
	{
		for (int i = 0; i < Face.size(); i++)
		{
			vec_face_map[Face[i].f1].push_back(i);
			vec_face_map[Face[i].f2].push_back(i);
			vec_face_map[Face[i].f3].push_back(i);
		}

		for (int i = 0; i < Face.size(); i++)
		{
			v1v2x = Vertex[Face[i].f3].x - Vertex[Face[i].f1].x;
			v1v2y = Vertex[Face[i].f3].y - Vertex[Face[i].f1].y;
			v1v2z = Vertex[Face[i].f3].z - Vertex[Face[i].f1].z;
			v1v3x = Vertex[Face[i].f2].x - Vertex[Face[i].f1].x;
			v1v3y = Vertex[Face[i].f2].y - Vertex[Face[i].f1].y;
			v1v3z = Vertex[Face[i].f2].z - Vertex[Face[i].f1].z;
			x = v1v2y * v1v3z - v1v2z * v1v3y;
			y = v1v2z * v1v3x - v1v2x * v1v3z;
			z = v1v2x * v1v3y - v1v2y * v1v3x;
			//归一化
			sum = sqrt(x*x + y*y + z*z);
			result.x = x / sum;
			result.y = y / sum;
			result.z = z / sum;
			FaceNormal.push_back(result);
		}
	}
	for (int i = 0; i < vec_face_map.size(); i++)
	{
		x = 0;
		y = 0;
		z = 0;
		for (int j = 0; j < vec_face_map[i].size(); j++)
		{
			x += FaceNormal[vec_face_map[i][j]].x;
			y += FaceNormal[vec_face_map[i][j]].y;
			z += FaceNormal[vec_face_map[i][j]].z;
		}
		//归一化
		sum = sqrt(x*x + y*y + z*z);
		result.x = x / sum;
		result.y = y / sum;
		result.z = z / sum;
		VertexNormal.push_back(result);
	}
}

void ScaffoldMesh::clear()
{
	TDF_control_grid.clear();
	TDF_knot_vector.clear();
	TBSS_control_grid.clear();
	TBSS_knot_vector.clear();

	vec_face_map.clear();
	Vertex.clear();
	VertexNormal.clear();
	Face.clear();
	FaceNormal.clear();
}

void ScaffoldMesh::zoom_TDF()
{
	double maxC, minC;
	maxC = minC = TDF_control_grid[0][0][0];
	for (int i = 0; i < u_points; i++)
	{
		for (int j = 0; j < v_points; j++)
		{
			for (int k = 0; k < w_points; k++)
			{
				minC = minC<TDF_control_grid[i][j][k] ? minC : TDF_control_grid[i][j][k];
				maxC = maxC>TDF_control_grid[i][j][k] ? maxC : TDF_control_grid[i][j][k];
			}
		}
	}
	if (type_TPMS == 0 || type_TPMS == 1)
	{
		if (type_structure == 0 || type_structure == 1)
		{
			scale_stretch(minC, maxC, -0.8, 0.8);
		}
		else
		{
			scale_stretch(minC, maxC, -0.6, 0.8);
		}
	}
	else if (type_TPMS == 2)
	{
		if (type_structure == 0 || type_structure == 1)
		{
			scale_stretch(minC, maxC, -0.6, 0.6);
		}
		else
		{
			scale_stretch(minC, maxC, -0.4, 0.6);
		}
	}
	else
	{
		if (type_structure == 0 || type_structure == 1)
		{
			scale_stretch(minC, maxC, -2.0, 2.0);
		}
		else
		{
			scale_stretch(minC, maxC, -1.8, 2.0);
		}
	}
}

void ScaffoldMesh::scale_stretch(double begin_minC, double begin_maxC, double end_minC, double end_maxC)
{
	double begin_range = begin_maxC - begin_minC;
	double end_range = end_maxC - end_minC;
	for (int i = 0; i < u_points; i++)
	{
		for (int j = 0; j < v_points; j++)
		{
			for (int k = 0; k < w_points; k++)
			{
				TDF_control_grid[i][j][k] = end_range*(TDF_control_grid[i][j][k] - begin_minC) / begin_range + end_minC;
			}
		}
	}
}


void ScaffoldMesh::discreate_TDF(double*** &C)
{
	//ofstream ss("dis_tdf.txt");
	double min_z, max_z;
	min_z = max_z = TBSS_control_grid[0][0][0].z;
	for (int i = 0; i < x_points; i++)
	{
		for (int j = 0; j < y_points; j++)
		{
			for (int k = 0; k < z_points; k++)
			{
				min_z = min_z < TBSS_control_grid[i][j][k].z ? min_z : TBSS_control_grid[i][j][k].z;
				max_z = max_z > TBSS_control_grid[i][j][k].z ? max_z : TBSS_control_grid[i][j][k].z;
			}
		}
	}
	
	size_t pointRes = resolution + 1;
	C = new double **[pointRes];
	for (int i = 0; i < pointRes; i++)
	{
		C[i] = new double *[pointRes];
		for (int j = 0; j < pointRes; j++)
		{
			C[i][j] = new double[pointRes];
		}
	}
	double step = 1.0 / resolution;
	for (int i = 0; i < pointRes; i++)
	{
		double u = i*step;
		for (int j = 0; j < pointRes; j++)
		{
			double v = j*step;
			for (int k = 0; k < pointRes; k++)
			{
				double w = k*step;
				//calculate_cvalue_TDF(Vector3D(u, v, w), C[i][j][k]);
				//C[i][j][k] = randomUniform(-0.8,0.8);
				calculate_cvalue_TDF_quick(Vector3D(u, v, w), C[i][j][k]);
				//C[i][j][k] = 0.3*(1 - (w - min_z) / (max_z - min_z)) + 0.7*(w - min_z) / (max_z - min_z);
				//ss << C[i][j][k] << endl;
			}
		}
	}
	//laplacianSmoothOuter(C, 2);
	//laplacianSmoothInner(C, 2);
	//save_porosity(C, "random");
	//ss.close();
}

//利用控制网格求离散TDF
void ScaffoldMesh::calculate_cvalue_TDF(Vector3D para, double &Cvalue)
{
	Cvalue = 0;
	double Mat[4][4];
	double vec[4];
	int uspan = findSpan(u_points, para.x, TDF_knot_vector[0]);
	double* Nu = basisFuns(uspan, para.x, TDF_knot_vector[0]);
	int vspan = findSpan(v_points, para.y, TDF_knot_vector[1]);
	double* Nv = basisFuns(vspan, para.y, TDF_knot_vector[1]);
	int wspan = findSpan(w_points, para.z, TDF_knot_vector[2]);
	double* Nw = basisFuns(wspan, para.z, TDF_knot_vector[2]);

	int uind = uspan - 3;
	for (int k = 0; k <= 3; k++)
	{
		int wind = wspan - 3 + k;
		for (int j = 0; j <= 3; j++)
		{
			int vind = vspan - 3 + j;
			Mat[k][j] = 0;
			for (int i = 0; i <= 3; i++)
			{
				Mat[k][j] += TDF_control_grid[uind + i][vind][wind] * Nu[i];
			}
		}
	}
	for (int k = 0; k <= 3; k++)
	{
		vec[k] = 0;
		for (int j = 0; j <= 3; j++)
		{
			vec[k] += Mat[k][j] * Nv[j];
		}
	}
	for (int k = 0; k <= 3; k++)
	{
		Cvalue += vec[k] * Nw[k];
	}
	delete[] Nu;
	delete[] Nv;
	delete[] Nw;

}
//利用控制网格拟合样条体
void ScaffoldMesh::calculate_posistion_TBSS(Vector3D para, Vector3D &phy_pos)
{
	phy_pos = Vector3D(0, 0, 0);
	Vector3D Mat[4][4];
	Vector3D vec[4];
	int uspan = findSpan(x_points, para.x, TBSS_knot_vector[0]);
	double* Nu = basisFuns(uspan, para.x, TBSS_knot_vector[0]);
	int vspan = findSpan(y_points, para.y, TBSS_knot_vector[1]);
	double* Nv = basisFuns(vspan, para.y, TBSS_knot_vector[1]);
	int wspan = findSpan(z_points, para.z, TBSS_knot_vector[2]);
	double* Nw = basisFuns(wspan, para.z, TBSS_knot_vector[2]);

	int uind = uspan - 3;
	for (int k = 0; k <= 3; k++)
	{
		int wind = wspan - 3 + k;
		for (int j = 0; j <= 3; j++)
		{
			int vind = vspan - 3 + j;
			Mat[k][j] = Vector3D(0, 0, 0);
			for (int i = 0; i <= 3; i++)
			{
				Mat[k][j] += TBSS_control_grid[uind + i][vind][wind] * Nu[i];
			}
		}
	}
	for (int k = 0; k <= 3; k++)
	{
		vec[k] = Vector3D(0, 0, 0);
		for (int j = 0; j <= 3; j++)
		{
			vec[k] += Mat[k][j] * Nv[j];
		}
	}
	for (int k = 0; k <= 3; k++)
	{
		phy_pos += vec[k] * Nw[k];
	}
	delete[] Nu;
	delete[] Nv;
	delete[] Nw;
}

//二分法确定u所属区间下标(非周期节点矢量）
int ScaffoldMesh::findSpan(int n, double u, std::vector<double> U)
{
	//n为控制点数，p为基函数次数，U为节点划分
	if (u == U[n])
		return(n - 1);
	int low = 3;
	int high = n;
	int mid = (low + high) / 2;
	while (u < U[mid] || u >= U[mid + 1])
	{
		if (u < U[mid])
			high = mid;
		else
			low = mid;
		mid = (low + high) / 2;
	}
	return(mid);
}

//计算所有非零B样条基函数的值
double* ScaffoldMesh::basisFuns(int i, double u, std::vector<double> U)
{
	//i为对应u所属的区间下标，p为基函数的次数，U为对应的节点矢量
	double* N = new double[4];
	double* left = new double[4];
	double* right = new double[4];
	double saved, temp;
	N[0] = 1.0;
	for (int j = 1; j <= 3; j++)
	{
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		saved = 0.0;
		for (int r = 0; r < j; r++)
		{
			temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;

	}
	delete[] left;
	delete[] right;
	return N;
}

//利用控制网格求离散TDF
void ScaffoldMesh::calculate_cvalue_TDF_quick(Vector3D para, double &Cvalue)
{
	double matri[4][4][4];//这个表示每一个基函数Bi*Bj*Bk
	int Bi_start_index, Bj_start_index, Bk_start_index;
	calculateBaseFunctionOfSolidBSpline(para, TDF_knot_vector, matri, Bi_start_index, Bj_start_index, Bk_start_index);
	Cvalue = 0;
	for (int ii = Bi_start_index; ii<Bi_start_index + 4; ii++)
	for (int jj = Bj_start_index; jj<Bj_start_index + 4; jj++)
	for (int kk = Bk_start_index; kk<Bk_start_index + 4; kk++)
		Cvalue = Cvalue + TDF_control_grid[ii][jj][kk] * matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];
}

//利用控制网格拟合样条体
void ScaffoldMesh::calculate_posistion_TBSS_quick(Vector3D para, Vector3D &phy_pos)
{
	double matri[4][4][4];//这个表示每一个基函数Bi*Bj*Bk
	int Bi_start_index, Bj_start_index, Bk_start_index;
	calculateBaseFunctionOfSolidBSpline(para, TBSS_knot_vector, matri, Bi_start_index, Bj_start_index, Bk_start_index);
	phy_pos = Vector3D(0, 0, 0);
	for (int ii = Bi_start_index; ii<Bi_start_index + 4; ii++)
	for (int jj = Bj_start_index; jj<Bj_start_index + 4; jj++)
	for (int kk = Bk_start_index; kk<Bk_start_index + 4; kk++)
		phy_pos = phy_pos + TBSS_control_grid[ii][jj][kk] * matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];
}

/****************************************
计算一个参数点返回Bi(u)Bj(v)Bk(w)
para为计算点的参数
kont_vector节点向量，大小为3，表示xyz三个轴
result表示4*4*4的结果矩阵，每一维表示Bi*Bj*Bk
Bi_start_index表示返回值Bi的向量的基函数的起始节点向量的下标
****************************************/
void ScaffoldMesh::calculateBaseFunctionOfSolidBSpline(Vector3D para, std::vector<std::vector<double> >&kont_vector, double(*result)[4][4], int &Bi_start_index, int &Bj_start_index, int &Bk_start_index)
{
	std::vector<double>Bi, Bj, Bk;
	Bi.resize(4);
	Bj.resize(4);
	Bk.resize(4);
	int seg_x = kont_vector[0].size() - 7;//表示x方向的节点向量被分成了几段
	int seg_y = kont_vector[1].size() - 7;
	int seg_z = kont_vector[2].size() - 7;
	Bi_start_index = static_cast<int>(para.x * seg_x);//每一个方向基函数在节点向量的起始下标
	int symbol = 0;
	if (Bi_start_index == seg_x)//这个表示节点向量的最大值不能取到seg，因为取到后计算Bi_start_index+3这段基函数会节点向量越界
	{
		Bi_start_index -= 1;
		symbol++;
	}
	Bj_start_index = static_cast<int>(para.y * seg_y);
	if (Bj_start_index == seg_y)
	{
		Bj_start_index -= 1;
		symbol++;
	}
	Bk_start_index = static_cast<int>(para.z * seg_z);
	if (Bk_start_index == seg_z)
	{
		Bk_start_index -= 1;
		symbol++;
	}
	if (symbol == 1)
	{//取到1的时候，那么只有最后一个控制点的B为1，其余B都为0
		for (int i = 0; i<4; i++)
		for (int j = 0; j<4; j++)
		for (int k = 0; k<4; k++)
			result[i][j][k] = 0;
		//return ;
	}
	std::string s = "-1";
	std::stringstream ss;
	//计算Bi方向的四个基函数Bi_start_index,Bi_start_index+1,Bi_start_index+2,Bi_start_index+3
	for (int i = 0; i<4; i++)
	{
		Bi[i] = getNi4(kont_vector[0], para.x, Bi_start_index + i, seg_x, Bi_start_index, s);
		if (Bi[i]<0)
			ss << "x  " << Bi[i] << " " << para.x << " " << Bi_start_index << " " << i << " " << s << "\n";
		if (Bi[i]<0)
			Bi[i] = 0;
	}
	for (int i = 0; i<4; i++)
	{
		Bj[i] = getNi4(kont_vector[1], para.y, Bj_start_index + i, seg_y, Bj_start_index, s);
		if (Bj[i]<0)
			ss << "y  " << Bj[i] << " " << para.y << " " << Bj_start_index << " " << i << " " << s << "\n";
		if (Bj[i]<0)
			Bj[i] = 0;
	}
	for (int i = 0; i<4; i++)
	{
		Bk[i] = getNi4(kont_vector[2], para.z, Bk_start_index + i, seg_z, Bk_start_index, s);
		if (Bk[i]<0)
			ss << "z  " << Bk[i] << " " << para.z << " " << Bk_start_index << " " << i << " " << s << "\n";
		if (Bk[i]<0)
			Bk[i] = 0;
	}
	//把Bi*Bj*Bk放入结果
	for (int i = 0; i<4; i++)
	for (int j = 0; j<4; j++)
	for (int k = 0; k<4; k++)
		result[i][j][k] = Bi[i] * Bj[j] * Bk[k];
}

/****************************************
计算基函数Ni,4
t为参数值
seg 表示被0-1被划分了 seg段
zone 代表t所在节点向量的起始下标减去3
deb这个参数只是我调试用的
****************************************/
double ScaffoldMesh::getNi4(const std::vector<double>&knot_vector, double t, int i, int seg, int zone, std::string  &deb)
{
	double tt[6] = { 0.0 };
	for (int j = 1; j<6; j++)
		tt[j] = knot_vector[i + j - 1];//把节点向量先赋值到t1----t5
	int left_most = 0, right_most = 0;//这个标记表示t是否取到某个区间的边界
	//首先计算t应该在节点向量的下标起点
	int index_t = static_cast<int>(t*seg);

	if (index_t == seg)//表示这个t是1.00，那么它需要特殊处理，
	{
		//right_most=1;//取到了右边界
		index_t -= 1;
	}
	index_t = zone;
	//因为一个t对应4个N值，分别为Nindex_t,4  Nindex_t+1,4  Nindex_t+2,4  Nindex_t+3,4
	int index_t_N = index_t;//表示t的相对的上面4个基函数的起点
	index_t += 3;//节点向量前面有4个0的缘故,表示下标开始的区间Ni不为0
	if (fabs(knot_vector[index_t] - t)<MINN_MINN)//
		left_most = 1;//取到了左边界
	if (fabs(knot_vector[index_t + 1] - t)<MINN_MINN)
		right_most = 1;//取到右边界
	//那么N(3-(i-index_t_N),1)这个下标是1
	int domain = index_t - i;
	double ans = 0;
	if (domain == 0)
	{
		if (0 == left_most && 0 == right_most)
		{
			ans = (t - tt[1])*(t - tt[1])*(t - tt[1]) / ((tt[4] - tt[1])*(tt[3] - tt[1])*(tt[2] - tt[1]));
			deb = "0_0";
			return ans;
		}
		else if (left_most)
		{
			ans = 0;
			deb = "0_1";
			return ans;
		}
		else
		{
			ans = (t - tt[1])*(t - tt[1]) / ((tt[4] - tt[1])*(tt[3] - tt[1]));
			deb = "0_2";
			return ans;
		}
	}//domain=0;
	else if (domain == 1)
	{
		double first = 0, second = 0, third = 0;//分别代表第一项，第二项
		if (0 == left_most && 0 == right_most)
		{
			first = (t - tt[2])*(t - tt[2])*(tt[5] - t) / ((tt[5] - tt[2])*(tt[4] - tt[2])*(tt[3] - tt[2]));
			second = (t - tt[2])*(t - tt[1])*(tt[4] - t) / ((tt[4] - tt[1])*(tt[4] - tt[2])*(tt[3] - tt[2]));
			third = (tt[3] - t)*(t - tt[1])*(t - tt[1]) / ((tt[4] - tt[1])*(tt[3] - tt[1])*(tt[3] - tt[2]));
			ans = first + second + third;
			deb = "1_0";
			return ans;
		}
		else if (left_most)
		{
			third = (t - tt[1])*(t - tt[1]) / ((tt[4] - tt[1])*(tt[3] - tt[1]));
			ans = third;
			deb = "1_1";
			return ans;
		}
		else
		{
			first = (t - tt[2])*(tt[5] - t) / ((tt[5] - tt[2])*(tt[4] - tt[2]));
			second = (t - tt[1])*(tt[4] - t) / ((tt[4] - tt[1])*(tt[4] - tt[2]));
			ans = first + second;
			deb = "1_2";
			return ans;
		}
	}//domain 1
	else if (domain == 2)
	{
		double first = 0, second = 0, third = 0;//分别代表第一项，第二项
		if (0 == left_most && 0 == right_most)
		{
			first = (tt[4] - t)*(tt[4] - t)*(t - tt[1]) / ((tt[4] - tt[1])*(tt[4] - tt[2])*(tt[4] - tt[3]));
			second = (tt[4] - t)*(tt[5] - t)*(t - tt[2]) / ((tt[5] - tt[2])*(tt[4] - tt[2])*(tt[4] - tt[3]));
			third = (t - tt[3])*(tt[5] - t)*(tt[5] - t) / ((tt[5] - tt[2])*(tt[5] - tt[3])*(tt[4] - tt[3]));
			ans = first + second + third;
			deb = "2_0";
			return ans;
		}
		else if (left_most)
		{
			first = (tt[4] - t)*(t - tt[1]) / ((tt[4] - tt[1])*(tt[4] - tt[2]));
			second = (tt[5] - t)*(t - tt[2]) / ((tt[5] - tt[2])*(tt[4] - tt[2]));
			ans = first + second;
			deb = "2_1";
			return ans;
		}
		else
		{
			third = (tt[5] - t)*(tt[5] - t) / ((tt[5] - tt[2])*(tt[5] - tt[3]));
			ans = third;
			deb = "2_2";
			return ans;
		}
	}//doamin2
	else
	{
		if (0 == left_most && 0 == right_most)
		{
			ans = (tt[5] - t)*(tt[5] - t)*(tt[5] - t) / ((tt[5] - tt[2])*(tt[5] - tt[3])*(tt[5] - tt[4]));
			deb = "3_0";
			return ans;
		}
		else if (left_most)
		{
			ans = (tt[5] - t)*(tt[5] - t) / ((tt[5] - tt[2])*(tt[5] - tt[3]));
			deb = "3_1";
			return ans;
		}
		else
		{
			ans = 0;
			deb = "3_2";
			return ans;
		}
	}//domain3
}


void ScaffoldMesh::constructEdge(Array3D<Point3D>& PointGrid, vector<edge>& Edge)
{
	edge e;
	e.index[0] = -1; e.index[1] = -1;
	for (size_t k = 0; k <= resolution; ++k)
	{
		for (size_t j = 0; j <= resolution; ++j)
		{
			for (size_t i = 0; i < resolution; ++i)
			{
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				Edge.push_back(e);
			}
		}
	}//index: 0 - res*(res+1)*(res+1)-1
	for (size_t k = 0; k <= resolution; ++k)
	{
		for (size_t i = 0; i <= resolution; ++i)
		{
			for (size_t j = 0; j < resolution; ++j)
			{
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k;
				Edge.push_back(e);
			}
		}
	}//index: res*(res+1)*(res+1) - 2*res*(res+1)*(res+1)-1
	for (size_t i = 0; i <= resolution; ++i)
	{
		for (size_t j = 0; j <= resolution; ++j)
		{
			for (size_t k = 0; k < resolution; ++k)
			{
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k + 1;
				Edge.push_back(e);
			}
		}
	}//index: 2*res*(res+1)*(res+1) - 3*res*(res+1)*(res+1)-1
	for (size_t k = 0; k <= resolution; ++k)//ij
	{
		for (size_t j = 0; j < resolution; ++j)
		{
			for (size_t i = 0; i < resolution; ++i)
			{
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k;
				Edge.push_back(e);
			}
		}
	}//index: 3*res*(res+1)*(res+1) - 3*res*(res+1)*(res+1)+res*res*(res+1)-1
	for (size_t j = 0; j <= resolution; ++j)//ik
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			for (size_t i = 0; i < resolution; ++i)
			{
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k + 1;
				e.id2 = (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				Edge.push_back(e);
			}
		}
	}//index: 3*res*(res+1)*(res+1)+res*res*(res+1) - 3*res*(res+1)*(res+1)+2*res*res*(res+1)-1
	for (size_t i = 0; i <= resolution; ++i)//jk
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			for (size_t j = 0; j < resolution; ++j)
			{
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1;
				Edge.push_back(e);
			}
		}
	}//index: 3*res*(res+1)*(res+1)+2*res*res*(res+1) - 3*res*(res+1)*(res+1)+3*res*res*(res+1)-1
	for (size_t k = 0; k < resolution; ++k)
	{
		for (size_t j = 0; j < resolution; ++j)
		{
			for (size_t i = 0; i < resolution; ++i)
			{
				e.id1 = i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k;
				e.id2 = (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + (k + 1);
				Edge.push_back(e);
			}
		}
	}//index: 3*res*(res+1)*(res+1)+3*res*res*(res+1) - 3*res*(res+1)*(res+1)+3*res*res*(res+1)+res*res*res-1
}


void ScaffoldMesh::cutVertex(Vector3D& _point)
{//避免计算误差
	_point.x = 1.0*(int(_point.x * 1000000 + 0.5)) / 1000000;
	_point.y = 1.0*(int(_point.y * 1000000 + 0.5)) / 1000000;
	_point.z = 1.0*(int(_point.z * 1000000 + 0.5)) / 1000000;
}

void ScaffoldMesh::calInterPoint(const Isosurface& surface, const Point3D& p1, const Point3D& p2, float isolevel, Point3D& interP)
{

	float v1 = p1.value;
	float v2 = p2.value;

	float x, y, z;

	if (v2 == v1)
	{
		x = (p1.x + p2.x) / 2.0f;
		y = (p1.y + p2.y) / 2.0f;
		z = (p1.z + p2.z) / 2.0f;
	}
	else
	{

		/*

		<----+-----+---+----->
		v1    |   v2
		isolevel


		<----+-----+---+----->
		0     |   1
		interp

		*/


		// interp == 0: vert should be at p1
		// interp == 1: vert should be at p2
		float interp = (isolevel - v1) / (v2 - v1);
		float oneMinusInterp = 1 - interp;

		x = p1.x * oneMinusInterp + p2.x * interp;
		y = p1.y * oneMinusInterp + p2.y * interp;
		z = p1.z * oneMinusInterp + p2.z * interp;
	}
	interP.x = x;
	interP.y = y;
	interP.z = z;
	interP.value = isolevel;
}

void ScaffoldMesh::calInterVert(const Isosurface& surface, const Point3D& p1, const Point3D& p2, float isolevel, Vector3D& interP)
{

	float v1 = p1.value;
	float v2 = p2.value;

	float x, y, z;

	if (v2 == v1)
	{
		x = (p1.x + p2.x) / 2.0f;
		y = (p1.y + p2.y) / 2.0f;
		z = (p1.z + p2.z) / 2.0f;
	}
	else
	{

		/*

		<----+-----+---+----->
		v1    |   v2
		isolevel


		<----+-----+---+----->
		0     |   1
		interp

		*/


		// interp == 0: vert should be at p1
		// interp == 1: vert should be at p2
		float interp = (isolevel - v1) / (v2 - v1);
		float oneMinusInterp = 1 - interp;

		x = p1.x * oneMinusInterp + p2.x * interp;
		y = p1.y * oneMinusInterp + p2.y * interp;
		z = p1.z * oneMinusInterp + p2.z * interp;
	}
	interP.x = x;
	interP.y = y;
	interP.z = z;
	cutVertex(interP);
}

void ScaffoldMesh::calInterVert2(const Isosurface& surface, const int edge_id, float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, int& index, int sur_id)
{
	if (Edge[edge_id].index[sur_id] == -1)
	{
		float value1 = pointset[Edge[edge_id].id1].value;
		float value2 = pointset[Edge[edge_id].id2].value;

		float x, y, z;

		if (value2 == value1)
		{
			x = (pointset[Edge[edge_id].id1].x + pointset[Edge[edge_id].id2].x) / 2.0f;
			y = (pointset[Edge[edge_id].id1].y + pointset[Edge[edge_id].id2].y) / 2.0f;
			z = (pointset[Edge[edge_id].id1].z + pointset[Edge[edge_id].id2].z) / 2.0f;
		}
		else
		{
			float interp = (isolevel - value1) / (value2 - value1);
			interp = 1.0*(int(interp * 100000 + 0.5)) / 100000;
			float oneMinusInterp = 1 - interp;
			if (interp == 1)
			{
				if (pointset[Edge[edge_id].id2].index == -1)
				{
					x = pointset[Edge[edge_id].id2].x;
					y = pointset[Edge[edge_id].id2].y;
					z = pointset[Edge[edge_id].id2].z;
					Vertex.push_back(Vector3D(x, y, z));
					pointset[Edge[edge_id].id2].index = Vertex.size() - 1;
				}
				Edge[edge_id].index[sur_id] = pointset[Edge[edge_id].id2].index;
			}
			else if (interp == 0)
			{
				if (pointset[Edge[edge_id].id1].index == -1)
				{
					x = pointset[Edge[edge_id].id1].x;
					y = pointset[Edge[edge_id].id1].y;
					z = pointset[Edge[edge_id].id1].z;
					Vertex.push_back(Vector3D(x, y, z));
					pointset[Edge[edge_id].id1].index = Vertex.size() - 1;
				}
				Edge[edge_id].index[sur_id] = pointset[Edge[edge_id].id1].index;
			}
			else
			{
				x = pointset[Edge[edge_id].id1].x* oneMinusInterp + pointset[Edge[edge_id].id2].x* interp;
				y = pointset[Edge[edge_id].id1].y* oneMinusInterp + pointset[Edge[edge_id].id2].y* interp;
				z = pointset[Edge[edge_id].id1].z* oneMinusInterp + pointset[Edge[edge_id].id2].z* interp;
				Vertex.push_back(Vector3D(x, y, z));
				Edge[edge_id].index[sur_id] = Vertex.size() - 1;
			}
		}
		//cutVertex(p);
	}
	index = Edge[edge_id].index[sur_id];
}

void ScaffoldMesh::extract_isosurface(const Isosurface& surface, float isolevel, int type, vector<edge>& Edge, vector<Point3D>& pointset, Array3D<Point3D>& PointGrid, int sur_id) // resolution indicates # of cubes
{
	for (size_t i = 0; i < resolution; ++i)
	{
		for (size_t j = 0; j < resolution; ++j)
		{
			for (size_t k = 0; k < resolution; ++k)
			{
				/*
				Coordinates:
				z
				|
				|___ y
				/
				/
				x

				Cube layout:
				4-------7
				/|      /|
				/ |     / |
				5-------6  |
				|  0----|--3
				| /     | /
				|/      |/
				1-------2

				Tetrahedrons are:
				0, 7, 3, 2
				0, 7, 2, 6
				0, 4, 6, 7
				0, 6, 1, 2
				0, 6, 1, 4
				5, 6, 1, 4
				*/

				const Point3D v[8] = {
					{ PointGrid.get(i, j, k) },
					{ PointGrid.get(i + 1, j, k) },
					{ PointGrid.get(i + 1, j + 1, k) },
					{ PointGrid.get(i, j + 1, k) },
					{ PointGrid.get(i, j, k + 1) },
					{ PointGrid.get(i + 1, j, k + 1) },
					{ PointGrid.get(i + 1, j + 1, k + 1) },
					{ PointGrid.get(i, j + 1, k + 1) }
				};
				const Point3D tetrahedra[6][4] = {
					{ v[0], v[7], v[3], v[2] },
					{ v[0], v[7], v[2], v[6] },
					{ v[0], v[4], v[7], v[6] },
					{ v[0], v[1], v[6], v[2] },
					{ v[0], v[4], v[6], v[1] },
					{ v[5], v[1], v[6], v[4] }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tet_edge[6][6] = {
					{ 3 * len1 + 2 * len2 + i*resolution*resolution + k*resolution + j, //v[0]v[7]
					len1 + k*(resolution + 1)*resolution + i*resolution + j, //v[0]v[3]
					3 * len1 + k*resolution*resolution + j*resolution + i, //v[0]v[2]
					2 * len1 + i*(resolution + 1)*resolution + (j + 1)*resolution + k, //v[7]v[3]
					3 * len1 + len2 + (j + 1)*resolution*resolution + k*resolution + i,//v[7]v[2]
					k*(resolution + 1)*resolution + (j + 1)*resolution + i//v[3]v[2]
					},
					{ 3 * len1 + 2 * len2 + i*resolution*resolution + k*resolution + j,//v[0]v[7]
					3 * len1 + k*resolution*resolution + j*resolution + i,//v[0]v[2]
					3 * len1 + 3 * len2 + k*resolution*resolution + j*resolution + i,//v[0]v[6]
					3 * len1 + len2 + (j + 1)*resolution*resolution + k*resolution + i,//v[7]v[2]
					(k + 1)*(resolution + 1)*resolution + (j + 1)*resolution + i,//v[7]v[6]
					2 * len1 + (i + 1)*(resolution + 1)*resolution + (j + 1)*resolution + k //v[2]v[6]
					},
					{ 2 * len1 + i*(resolution + 1)*resolution + j*resolution + k, //v[0]v[4]
					3 * len1 + 2 * len2 + i*resolution*resolution + k*resolution + j,//v[0]v[7]
					3 * len1 + 3 * len2 + k*resolution*resolution + j*resolution + i,//v[0]v[6]
					len1 + (k + 1)*(resolution + 1)*resolution + i*resolution + j, //v[4]v[7]
					3 * len1 + (k + 1)*resolution*resolution + j*resolution + i, //v[4]v[6]
					(k + 1)*(resolution + 1)*resolution + (j + 1)*resolution + i//v[7]v[6]
					},
					{ k*(resolution + 1)*resolution + j*resolution + i,//v[0]v[1]
					3 * len1 + 3 * len2 + k*resolution*resolution + j*resolution + i,//v[0]v[6]
					3 * len1 + k*resolution*resolution + j*resolution + i,//v[0]v[2]
					3 * len1 + 2 * len2 + (i + 1)*resolution*resolution + k*resolution + j, //v[1]v[6]
					len1 + k*(resolution + 1)*resolution + (i + 1)*resolution + j, //v[1]v[2]
					2 * len1 + (i + 1)*(resolution + 1)*resolution + (j + 1)*resolution + k //v[6]v[2]
					},
					{ 2 * len1 + i*(resolution + 1)*resolution + j*resolution + k, //v[0]v[4]
					3 * len1 + 3 * len2 + k*resolution*resolution + j*resolution + i,//v[0]v[6]
					k*(resolution + 1)*resolution + j*resolution + i,//v[0]v[1]
					3 * len1 + (k + 1)*resolution*resolution + j*resolution + i, //v[4]v[6]
					3 * len1 + len2 + j*resolution*resolution + k*resolution + i,//v[4]v[1]
					3 * len1 + 2 * len2 + (i + 1)*resolution*resolution + k*resolution + j //v[1]v[6]
					},
					{ 2 * len1 + (i + 1)*(resolution + 1)*resolution + j*resolution + k, //v[5]v[1]
					len1 + (k + 1)*(resolution + 1)*resolution + (i + 1)*resolution + j, //v[5]v[6]
					(k + 1)*(resolution + 1)*resolution + j*resolution + i,//v[5]v[4]
					3 * len1 + 2 * len2 + (i + 1)*resolution*resolution + k*resolution + j, //v[1]v[6]
					3 * len1 + len2 + j*resolution*resolution + k*resolution + i,//v[1]v[4]
					3 * len1 + (k + 1)*resolution*resolution + j*resolution + i, //v[6]v[4]
					}
				};
				for (int t = 0; t < 6; ++t)
				{
					extractTriangleFromTetrahedron(surface, tetrahedra[t], tet_edge[t], isolevel, pointset, Edge, sur_id);
				}
			}
		}
	}
}

void ScaffoldMesh::extractTriangleFromTetrahedron(const Isosurface& surface, const Point3D p[4], const int edges[6], float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, int sur_id)
{
	Triangle face;
	int index1 = -1, index2 = -1, index3 = -1, index4 = -1;
	unsigned char index = 0;
	for (int i = 0; i < 4; ++i)
	if (p[i].value < isolevel)
		index |= (1 << i);

	switch (index) {

		// we don't do anything if everyone is inside or outside
	case 0x00:
	case 0x0F:
		break;

		// only vert 0 is inside
	case 0x01:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index2, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index3, sur_id);//p[0]p[2]
		addFace(index1, index2, index3);
		break;

		// only vert 1 is inside
	case 0x02:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, index2, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, index3, sur_id);//p[1]p[3]
		addFace(index1, index2, index3);
		break;

		// only vert 2 is inside
	case 0x04:
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index1, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, index2, sur_id);//p[2]p[3]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, index3, sur_id);//p[1]p[2]
		addFace(index1, index2, index3);
		break;

		// only vert 3 is inside
	case 0x08:
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, index1, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, index2, sur_id);//p[2]p[3]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index3, sur_id);//p[0]p[3]
		addFace(index1, index2, index3);
		break;

		// verts 0, 1 are inside
	case 0x03:
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index1, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, index3, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, index4, sur_id);//p[1]p[2]
		addFace(index1, index2, index3);
		addFace(index3, index2, index4);
		break;

		// verts 0, 2 are inside
	case 0x05:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index2, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, index3, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, index4, sur_id);//p[2]p[3]
		addFace(index1, index2, index3);
		addFace(index3, index2, index4);
		break;

		// verts 0, 3 are inside
	case 0x09:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, index2, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index3, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, index4, sur_id);//p[2]p[3]
		addFace(index1, index2, index3); 
		addFace(index3, index2, index4);
		break;

		// verts 1, 2 are inside
	case 0x06:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, index3, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, index4, sur_id);//p[2]p[3]
		addFace(index1, index2, index3);
		addFace(index3, index2, index4);
		break;

		// verts 2, 3 are inside
	case 0x0C:
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index1, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, index2, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index3, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, index4, sur_id);//p[1]p[2]
		addFace(index1, index2, index3);
		addFace(index3, index2, index4);
		break;

		// verts 1, 3 are inside
	case 0x0A:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, index2, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index3, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, index4, sur_id);//p[2]p[3]
		addFace(index1, index2, index3);
		addFace(index3, index2, index4);
		break;

		// verts 0, 1, 2 are inside
	case 0x07:
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index1, sur_id);//p[0]p[3]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, index2, sur_id);//p[2]p[3]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, index3, sur_id);//p[1]p[3]
		addFace(index1, index2, index3);
		break;

		// verts 0, 1, 3 are inside
	case 0x0B:
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, index1, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[5], isolevel, pointset, Edge, index2, sur_id);//p[2]p[3]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index3, sur_id);//p[0]p[2]
		addFace(index1, index2, index3);
		break;

		// verts 0, 2, 3 are inside
	case 0x0D:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[4], isolevel, pointset, Edge, index2, sur_id);//p[1]p[3]
		calInterVert2(surface, edges[3], isolevel, pointset, Edge, index3, sur_id);//p[1]p[2]
		addFace(index1, index2, index3);
		break;

		// verts 1, 2, 3 are inside
	case 0x0E:
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index3, sur_id);//p[0]p[3]
		addFace(index1, index2, index3);
		break;

		// what is this I don't even
	default:
		assert(false);
	}
}

void ScaffoldMesh::close_isosurface_network(const Isosurface& surface, float isolevel, int type, vector<edge>& Edge, vector<Point3D>& pointset, Array3D<Point3D>& PointGrid, int sur_id)
{
	//提取六张表面生成封闭模型
	float x1, y1, z1, x2, y2, z2;
	for (size_t i = 0; i < resolution; ++i)
	{
		for (size_t j = 0; j < resolution; ++j)
		{
			const Point3D v[8] = {
				{ PointGrid.get(i, j, 0) },
				{ PointGrid.get(i + 1, j, 0) },
				{ PointGrid.get(i + 1, j + 1, 0) },
				{ PointGrid.get(i, j + 1, 0) },
				{ PointGrid.get(i, j, resolution) },
				{ PointGrid.get(i + 1, j, resolution) },
				{ PointGrid.get(i + 1, j + 1, resolution) },
				{ PointGrid.get(i, j + 1, resolution) },
			};
			if (type == 0)
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[2], v[1] },
					{ v[0], v[3], v[2] },
					{ v[4], v[5], v[6] },
					{ v[4], v[6], v[7] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution, i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ 3 * len1 + j*resolution + i,//v[0]v[2]
					j*resolution + i,//v[0]v[1]
					len1 + (i + 1)*resolution + j//v[2]v[1]
					},
					{ len1 + i*resolution + j, //v[0]v[3]
					3 * len1 + j*resolution + i,//v[0]v[2]
					(j + 1)*resolution + i//v[3]v[2]
					},
					{ resolution*(resolution + 1)*resolution + j*resolution + i,//v[4]v[5]
					3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					len1 + resolution*(resolution + 1)*resolution + (i + 1)*resolution + j//v[5]v[6]
					},
					{ 3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					len1 + resolution*(resolution + 1)*resolution + i*resolution + j, //v[4]v[7]
					resolution*(resolution + 1)*resolution + (j + 1)*resolution + i//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleInner2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, sur_id);
				}
			}
			else
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[1], v[2] },
					{ v[0], v[2], v[3] },
					{ v[4], v[6], v[5] },
					{ v[4], v[7], v[6] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ j*resolution + i,//v[0]v[1]
					3 * len1 + j*resolution + i,//v[0]v[2]
					len1 + (i + 1)*resolution + j//v[2]v[1]
					},
					{ 3 * len1 + j*resolution + i,//v[0]v[2]
					len1 + i*resolution + j, //v[0]v[3]
					(j + 1)*resolution + i//v[3]v[2]
					},
					{ 3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					resolution*(resolution + 1)*resolution + j*resolution + i,//v[4]v[5]
					len1 + resolution*(resolution + 1)*resolution + (i + 1)*resolution + j//v[5]v[6]
					},
					{ len1 + resolution*(resolution + 1)*resolution + i*resolution + j, //v[4]v[7]
					3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					resolution*(resolution + 1)*resolution + (j + 1)*resolution + i//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleOuter2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, sur_id);
				}
			}
		}
	}
	for (size_t i = 0; i < resolution; ++i)
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			const Point3D v[8] = {
				{ PointGrid.get(i, 0, k) },
				{ PointGrid.get(i + 1, 0, k) },
				{ PointGrid.get(i + 1, 0, k + 1) },
				{ PointGrid.get(i, 0, k + 1) },
				{ PointGrid.get(i, resolution, k) },
				{ PointGrid.get(i + 1, resolution, k) },
				{ PointGrid.get(i + 1, resolution, k + 1) },
				{ PointGrid.get(i, resolution, k + 1) },
			};
			if (type == 0)
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[1], v[3] },
					{ v[1], v[2], v[3] },
					{ v[4], v[7], v[5] },
					{ v[5], v[7], v[6] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + k + 1 },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + k + 1, i*(resolution + 1)*(resolution + 1) + k + 1 },
					{ i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ k*(resolution + 1)*resolution + i,//v[0]v[1]
					2 * len1 + i*(resolution + 1)*resolution + k,//v[0]v[3]
					3 * len1 + len2 + k*resolution + i//v[1]v[3]
					},
					{ 2 * len1 + (i + 1)*(resolution + 1)*resolution + k, //v[1]v[2]
					3 * len1 + len2 + k*resolution + i,//v[1]v[3]
					(k + 1)*(resolution + 1)*resolution + i//v[2]v[3]
					},
					{ 2 * len1 + i*(resolution + 1)*resolution + resolution*resolution + k,//v[4]v[7]
					k*(resolution + 1)*resolution + resolution*resolution + i,//v[4]v[5]
					3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i//v[7]v[5]
					},
					{ 3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i,//v[5]v[7]
					2 * len1 + (i + 1)*(resolution + 1)*resolution + resolution*resolution + k, //v[5]v[6]
					(k + 1)*(resolution + 1)*resolution + resolution*resolution + i//v[7]v[6]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleInner2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, sur_id);
				}
			}
			else
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[3], v[1] },
					{ v[1], v[3], v[2] },
					{ v[4], v[5], v[7] },
					{ v[5], v[6], v[7] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + k },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + k + 1 },
					{ i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1 },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ 2 * len1 + i*(resolution + 1)*resolution + k,//v[0]v[3]
					k*(resolution + 1)*resolution + i,//v[0]v[1]
					3 * len1 + len2 + k*resolution + i//v[1]v[3]
					},
					{ 3 * len1 + len2 + k*resolution + i,//v[1]v[3]
					2 * len1 + (i + 1)*(resolution + 1)*resolution + k, //v[1]v[2]
					(k + 1)*(resolution + 1)*resolution + i//v[2]v[3]
					},
					{ k*(resolution + 1)*resolution + resolution*resolution + i,//v[4]v[5]
					2 * len1 + i*(resolution + 1)*resolution + resolution*resolution + k,//v[4]v[7]
					3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i//v[7]v[5]
					},
					{ 2 * len1 + (i + 1)*(resolution + 1)*resolution + resolution*resolution + k, //v[5]v[6]
					3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i,//v[5]v[7]
					(k + 1)*(resolution + 1)*resolution + resolution*resolution + i//v[7]v[6]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleOuter2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, sur_id);
				}
			}
		}
	}
	for (size_t j = 0; j < resolution; ++j)
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			const Point3D v[8] = {
				{ PointGrid.get(0, j, k) },
				{ PointGrid.get(0, j + 1, k) },
				{ PointGrid.get(0, j + 1, k + 1) },
				{ PointGrid.get(0, j, k + 1) },
				{ PointGrid.get(resolution, j, k) },
				{ PointGrid.get(resolution, j + 1, k) },
				{ PointGrid.get(resolution, j + 1, k + 1) },
				{ PointGrid.get(resolution, j, k + 1) },
			};
			if (type == 0)
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[2], v[1] },
					{ v[0], v[3], v[2] },
					{ v[4], v[5], v[6] },
					{ v[4], v[6], v[7] },
				};
				const int Pindex[4][3] = {
					{ j*(resolution + 1) + k, (j + 1)*(resolution + 1) + k + 1, (j + 1)*(resolution + 1) + k },
					{ j*(resolution + 1) + k, j*(resolution + 1) + k + 1, (j + 1)*(resolution + 1) + k + 1 },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1 },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1, resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ 3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					len1 + k*(resolution + 1)*resolution + j,//v[0]v[1]
					2 * len1 + (j + 1)*resolution + k//v[2]v[1]
					},
					{ 2 * len1 + j*resolution + k, //v[0]v[3]
					3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					len1 + (k + 1)*(resolution + 1)*resolution + j//v[3]v[2]
					},
					{ len1 + k*(resolution + 1)*resolution + resolution*resolution + j,//v[4]v[5]
					3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					2 * len1 + resolution*(resolution + 1)*resolution + (j + 1)*resolution + k//v[5]v[6]
					},
					{ 3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					2 * len1 + resolution*(resolution + 1)*resolution + j*resolution + k, //v[4]v[7]
					len1 + (k + 1)*(resolution + 1)*resolution + resolution*resolution + j//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleInner2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, sur_id);
				}
			}
			else
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[1], v[2] },
					{ v[0], v[2], v[3] },
					{ v[4], v[6], v[5] },
					{ v[4], v[7], v[6] },
				};
				const int Pindex[4][3] = {
					{ j*(resolution + 1) + k, (j + 1)*(resolution + 1) + k, (j + 1)*(resolution + 1) + k + 1 },
					{ j*(resolution + 1) + k, (j + 1)*(resolution + 1) + k + 1, j*(resolution + 1) + k + 1 },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k + 1, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ len1 + k*(resolution + 1)*resolution + j,//v[0]v[1]
					3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					2 * len1 + (j + 1)*resolution + k//v[2]v[1]
					},
					{ 3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					2 * len1 + j*resolution + k, //v[0]v[3]
					len1 + (k + 1)*(resolution + 1)*resolution + j//v[3]v[2]
					},
					{ 3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					len1 + k*(resolution + 1)*resolution + resolution*resolution + j,//v[4]v[5]
					2 * len1 + resolution*(resolution + 1)*resolution + (j + 1)*resolution + k//v[5]v[6]
					},
					{ 2 * len1 + resolution*(resolution + 1)*resolution + j*resolution + k, //v[4]v[7]
					3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					len1 + (k + 1)*(resolution + 1)*resolution + resolution*resolution + j//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleOuter2(surface, triangle[t], tri_edge[t], Pindex[t], isolevel, pointset, Edge, sur_id);
				}
			}
		}
	}
}

void ScaffoldMesh::extractTriangleOuter2(const Isosurface& surface, const Point3D p[3], const int edges[3], const int Pindex[3], float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, int sur_id)
{
	int index1 = -1, index2 = -1, index3 = -1, index4 = -1;
	Triangle face;
	if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//0,1,2 outer
		getPointinsurface(Pindex[0], pointset, index1);
		getPointinsurface(Pindex[1], pointset, index2);
		getPointinsurface(Pindex[2], pointset, index3);
		addFace(index1, index2, index3);
	}
	else if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//0,1 outer 2 inner
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index1, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[0], pointset, index3);
		getPointinsurface(Pindex[1], pointset, index4);
		addFace(index3, index4, index1);
		addFace(index3, index1, index2);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//0,2 outer 1 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index2, sur_id);//p[1]p[2]
		getPointinsurface(Pindex[0], pointset, index3);
		getPointinsurface(Pindex[2], pointset, index4);
		addFace(index3, index1, index2);
		addFace(index3, index2, index4);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0 outer 1,2 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[0], pointset, index3);
		addFace(index3, index1, index2);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//1,2 outer 0 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[1], pointset, index3);
		getPointinsurface(Pindex[2], pointset, index4);
		addFace(index1, index3, index4);
		addFace(index1, index4, index2);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//1 outer 0,2 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index2, sur_id);//p[1]p[2]
		getPointinsurface(Pindex[1], pointset, index3);
		addFace(index1, index3, index2);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//2 outer 0,1 inner
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index1, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[2], pointset, index3);
		addFace(index1, index3, index2);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0,1,2 inner
	}
}

void ScaffoldMesh::extractTriangleInner2(const Isosurface& surface, const Point3D p[3], const int edges[3], const int Pindex[3], float isolevel, vector<Point3D>& pointset, vector<edge>& Edge, int sur_id)
{
	int index1 = -1, index2 = -1, index3 = -1, index4 = -1;
	Triangle face;
	if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//0,1,2 outer
	}
	else if (p[0].value >= isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//0,1 outer 2 inner
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index1, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[2], pointset, index3);
		addFace(index1, index3, index2);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//0,2 outer 1 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index2, sur_id);//p[1]p[2]
		getPointinsurface(Pindex[1], pointset, index3);
		addFace(index1, index3, index2);
	}
	else if (p[0].value >= isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0 outer 1,2 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[1], pointset, index3);
		getPointinsurface(Pindex[2], pointset, index4);
		addFace(index1, index3, index4);
		addFace(index1, index4, index2);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value >= isolevel)
	{//1,2 outer 0 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[0], pointset, index3);
		addFace(index3, index1, index2);
	}
	else if (p[0].value < isolevel && p[1].value >= isolevel && p[2].value < isolevel)
	{//1 outer 0,2 inner
		calInterVert2(surface, edges[0], isolevel, pointset, Edge, index1, sur_id);//p[0]p[1]
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index2, sur_id);//p[1]p[2]
		getPointinsurface(Pindex[0], pointset, index3);
		getPointinsurface(Pindex[2], pointset, index4);
		addFace(index3, index1, index2);
		addFace(index3, index2, index4);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value >= isolevel)
	{//2 outer 0,1 inner
		calInterVert2(surface, edges[2], isolevel, pointset, Edge, index1, sur_id);//p[1]p[2]
		calInterVert2(surface, edges[1], isolevel, pointset, Edge, index2, sur_id);//p[0]p[2]
		getPointinsurface(Pindex[0], pointset, index3);
		getPointinsurface(Pindex[1], pointset, index4);
		addFace(index3, index4, index1);
		addFace(index3, index1, index2);
	}
	else if (p[0].value < isolevel && p[1].value < isolevel && p[2].value < isolevel)
	{//0,1,2 inner
		getPointinsurface(Pindex[0], pointset, index1);
		getPointinsurface(Pindex[1], pointset, index2);
		getPointinsurface(Pindex[2], pointset, index3);
		addFace(index1, index2, index3);
	}
}

void ScaffoldMesh::getPointinsurface(const int id, vector<Point3D>& pointset, int& index)
{
	if (pointset[id].index == -1)
	{
		Vertex.push_back(Vector3D(pointset[id].x, pointset[id].y, pointset[id].z));
		pointset[id].index = Vertex.size() - 1;
	}
	index = pointset[id].index;
}

void ScaffoldMesh::close_isosurface_sheet(const Isosurface& surface, float maxisolevel, float minisolevel, int type, vector<edge>& Edge, vector<Point3D>& pointset, Array3D<Point3D>& PointGrid)
{
	//提取六张表面生成封闭模型
	float x1, y1, z1, x2, y2, z2;
	for (size_t i = 0; i < resolution; ++i)
	{
		for (size_t j = 0; j < resolution; ++j)
		{
			const Point3D v[8] = {
				{ PointGrid.get(i, j, 0) },
				{ PointGrid.get(i + 1, j, 0) },
				{ PointGrid.get(i + 1, j + 1, 0) },
				{ PointGrid.get(i, j + 1, 0) },
				{ PointGrid.get(i, j, resolution) },
				{ PointGrid.get(i + 1, j, resolution) },
				{ PointGrid.get(i + 1, j + 1, resolution) },
				{ PointGrid.get(i, j + 1, resolution) },
			};
			if (type == 0)
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[2], v[1] },
					{ v[0], v[3], v[2] },
					{ v[4], v[5], v[6] },
					{ v[4], v[6], v[7] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution, i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ 3 * len1 + j*resolution + i,//v[0]v[2]
					j*resolution + i,//v[0]v[1]
					len1 + (i + 1)*resolution + j//v[2]v[1]
					},
					{ len1 + i*resolution + j, //v[0]v[3]
					3 * len1 + j*resolution + i,//v[0]v[2]
					(j + 1)*resolution + i//v[3]v[2]
					},
					{ resolution*(resolution + 1)*resolution + j*resolution + i,//v[4]v[5]
					3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					len1 + resolution*(resolution + 1)*resolution + (i + 1)*resolution + j//v[5]v[6]
					},
					{ 3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					len1 + resolution*(resolution + 1)*resolution + i*resolution + j, //v[4]v[7]
					resolution*(resolution + 1)*resolution + (j + 1)*resolution + i//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleFromSurface(surface, triangle[t], tri_edge[t], Pindex[t], maxisolevel, minisolevel, pointset, Edge);
				}
			}
			else{
				const Point3D triangle[4][3] = {
					{ v[0], v[1], v[2] },
					{ v[0], v[2], v[3] },
					{ v[4], v[6], v[5] },
					{ v[4], v[7], v[6] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1), (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1), i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution },
					{ i*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + resolution, i*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution, (i + 1)*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + resolution }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ j*resolution + i,//v[0]v[1]
					3 * len1 + j*resolution + i,//v[0]v[2]
					len1 + (i + 1)*resolution + j//v[2]v[1]
					},
					{ 3 * len1 + j*resolution + i,//v[0]v[2]
					len1 + i*resolution + j, //v[0]v[3]
					(j + 1)*resolution + i//v[3]v[2]
					},
					{ 3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					resolution*(resolution + 1)*resolution + j*resolution + i,//v[4]v[5]
					len1 + resolution*(resolution + 1)*resolution + (i + 1)*resolution + j//v[5]v[6]
					},
					{ len1 + resolution*(resolution + 1)*resolution + i*resolution + j, //v[4]v[7]
					3 * len1 + resolution*resolution*resolution + j*resolution + i,//v[4]v[6]
					resolution*(resolution + 1)*resolution + (j + 1)*resolution + i//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleFromSurface(surface, triangle[t], tri_edge[t], Pindex[t], maxisolevel, minisolevel, pointset, Edge);
				}
			}
		}
	}
	for (size_t i = 0; i < resolution; ++i)
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			const Point3D v[8] = {
				{ PointGrid.get(i, 0, k) },
				{ PointGrid.get(i + 1, 0, k) },
				{ PointGrid.get(i + 1, 0, k + 1) },
				{ PointGrid.get(i, 0, k + 1) },
				{ PointGrid.get(i, resolution, k) },
				{ PointGrid.get(i + 1, resolution, k) },
				{ PointGrid.get(i + 1, resolution, k + 1) },
				{ PointGrid.get(i, resolution, k + 1) },
			};
			if (type == 0){
				const Point3D triangle[4][3] = {
					{ v[0], v[1], v[3] },
					{ v[1], v[2], v[3] },
					{ v[4], v[7], v[5] },
					{ v[5], v[7], v[6] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + k + 1 },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + k + 1, i*(resolution + 1)*(resolution + 1) + k + 1 },
					{ i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ k*(resolution + 1)*resolution + i,//v[0]v[1]
					2 * len1 + i*(resolution + 1)*resolution + k,//v[0]v[3]
					3 * len1 + len2 + k*resolution + i//v[1]v[3]
					},
					{ 2 * len1 + (i + 1)*(resolution + 1)*resolution + k, //v[1]v[2]
					3 * len1 + len2 + k*resolution + i,//v[1]v[3]
					(k + 1)*(resolution + 1)*resolution + i//v[2]v[3]
					},
					{ 2 * len1 + i*(resolution + 1)*resolution + resolution*resolution + k,//v[4]v[7]
					k*(resolution + 1)*resolution + resolution*resolution + i,//v[4]v[5]
					3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i//v[7]v[5]
					},
					{ 3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i,//v[5]v[7]
					2 * len1 + (i + 1)*(resolution + 1)*resolution + resolution*resolution + k, //v[5]v[6]
					(k + 1)*(resolution + 1)*resolution + resolution*resolution + i//v[7]v[6]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleFromSurface(surface, triangle[t], tri_edge[t], Pindex[t], maxisolevel, minisolevel, pointset, Edge);
				}
			}
			else
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[3], v[1] },
					{ v[1], v[3], v[2] },
					{ v[4], v[5], v[7] },
					{ v[5], v[6], v[7] },
				};
				const int Pindex[4][3] = {
					{ i*(resolution + 1)*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + k },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + k + 1, (i + 1)*(resolution + 1)*(resolution + 1) + k + 1 },
					{ i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1 },
					{ (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k, (i + 1)*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1, i*(resolution + 1)*(resolution + 1) + resolution*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ 2 * len1 + i*(resolution + 1)*resolution + k,//v[0]v[3]
					k*(resolution + 1)*resolution + i,//v[0]v[1]
					3 * len1 + len2 + k*resolution + i//v[1]v[3]
					},
					{ 3 * len1 + len2 + k*resolution + i,//v[1]v[3]
					2 * len1 + (i + 1)*(resolution + 1)*resolution + k, //v[1]v[2]
					(k + 1)*(resolution + 1)*resolution + i//v[2]v[3]
					},
					{ k*(resolution + 1)*resolution + resolution*resolution + i,//v[4]v[5]
					2 * len1 + i*(resolution + 1)*resolution + resolution*resolution + k,//v[4]v[7]
					3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i//v[7]v[5]
					},
					{ 2 * len1 + (i + 1)*(resolution + 1)*resolution + resolution*resolution + k, //v[5]v[6]
					3 * len1 + len2 + resolution*resolution*resolution + k*resolution + i,//v[5]v[7]
					(k + 1)*(resolution + 1)*resolution + resolution*resolution + i//v[7]v[6]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleFromSurface(surface, triangle[t], tri_edge[t], Pindex[t], maxisolevel, minisolevel, pointset, Edge);
				}
			}
		}
	}
	for (size_t j = 0; j < resolution; ++j)
	{
		for (size_t k = 0; k < resolution; ++k)
		{
			const Point3D v[8] = {
				{ PointGrid.get(0, j, k) },
				{ PointGrid.get(0, j + 1, k) },
				{ PointGrid.get(0, j + 1, k + 1) },
				{ PointGrid.get(0, j, k + 1) },
				{ PointGrid.get(resolution, j, k) },
				{ PointGrid.get(resolution, j + 1, k) },
				{ PointGrid.get(resolution, j + 1, k + 1) },
				{ PointGrid.get(resolution, j, k + 1) },
			};
			if (type == 0){
				const Point3D triangle[4][3] = {
					{ v[0], v[2], v[1] },
					{ v[0], v[3], v[2] },
					{ v[4], v[5], v[6] },
					{ v[4], v[6], v[7] },
				};
				const int Pindex[4][3] = {
					{ j*(resolution + 1) + k, (j + 1)*(resolution + 1) + k + 1, (j + 1)*(resolution + 1) + k },
					{ j*(resolution + 1) + k, j*(resolution + 1) + k + 1, (j + 1)*(resolution + 1) + k + 1 },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1 },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1, resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ 3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					len1 + k*(resolution + 1)*resolution + j,//v[0]v[1]
					2 * len1 + (j + 1)*resolution + k//v[2]v[1]
					},
					{ 2 * len1 + j*resolution + k, //v[0]v[3]
					3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					len1 + (k + 1)*(resolution + 1)*resolution + j//v[3]v[2]
					},
					{ len1 + k*(resolution + 1)*resolution + resolution*resolution + j,//v[4]v[5]
					3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					2 * len1 + resolution*(resolution + 1)*resolution + (j + 1)*resolution + k//v[5]v[6]
					},
					{ 3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					2 * len1 + resolution*(resolution + 1)*resolution + j*resolution + k, //v[4]v[7]
					len1 + (k + 1)*(resolution + 1)*resolution + resolution*resolution + j//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleFromSurface(surface, triangle[t], tri_edge[t], Pindex[t], maxisolevel, minisolevel, pointset, Edge);
				}
			}
			else
			{
				const Point3D triangle[4][3] = {
					{ v[0], v[1], v[2] },
					{ v[0], v[2], v[3] },
					{ v[4], v[6], v[5] },
					{ v[4], v[7], v[6] },
				};
				const int Pindex[4][3] = {
					{ j*(resolution + 1) + k, (j + 1)*(resolution + 1) + k, (j + 1)*(resolution + 1) + k + 1 },
					{ j*(resolution + 1) + k, (j + 1)*(resolution + 1) + k + 1, j*(resolution + 1) + k + 1 },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k },
					{ resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k, resolution*(resolution + 1)*(resolution + 1) + j*(resolution + 1) + k + 1, resolution*(resolution + 1)*(resolution + 1) + (j + 1)*(resolution + 1) + k + 1 }
				};
				int len1 = resolution*(resolution + 1)*(resolution + 1);
				int len2 = resolution*resolution*(resolution + 1);
				const int tri_edge[4][3] = {
					{ len1 + k*(resolution + 1)*resolution + j,//v[0]v[1]
					3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					2 * len1 + (j + 1)*resolution + k//v[2]v[1]
					},
					{ 3 * len1 + 2 * len2 + k*resolution + j,//v[0]v[2]
					2 * len1 + j*resolution + k, //v[0]v[3]
					len1 + (k + 1)*(resolution + 1)*resolution + j//v[3]v[2]
					},
					{ 3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					len1 + k*(resolution + 1)*resolution + resolution*resolution + j,//v[4]v[5]
					2 * len1 + resolution*(resolution + 1)*resolution + (j + 1)*resolution + k//v[5]v[6]
					},
					{ 2 * len1 + resolution*(resolution + 1)*resolution + j*resolution + k, //v[4]v[7]
					3 * len1 + 2 * len2 + resolution*resolution*resolution + k*resolution + j,//v[4]v[6]
					len1 + (k + 1)*(resolution + 1)*resolution + resolution*resolution + j//v[6]v[7]
					}
				};
				for (int t = 0; t < 4; ++t)
				{
					extractTriangleFromSurface(surface, triangle[t], tri_edge[t], Pindex[t], maxisolevel, minisolevel, pointset, Edge);
				}
			}
		}
	}
}

void ScaffoldMesh::extractTriangleFromSurface(const Isosurface& surface, const Point3D p[3], const int edges[3], const int Pindex[3], float maxisolevel, float minisolevel, vector<Point3D>& pointset, vector<edge>& Edge)
{
	int index1 = -1, index2 = -1, index3 = -1, index4 = -1, index5 = -1;
	if (p[0].value >= maxisolevel && p[1].value >= maxisolevel && p[2].value >= maxisolevel)
	{//0,1,2 outer
	}
	else if (p[0].value >= maxisolevel && p[1].value >= maxisolevel && p[2].value < maxisolevel)
	{
		if (p[2].value >= minisolevel)
		{
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index1, 0);
			getPointinsurface(Pindex[2], pointset, index2);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index3, 0);

			addFace(index1, index2, index3);
		}
		else if (p[2].value < minisolevel)
		{
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index1, 0);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index2, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index3, 1);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index4, 0);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
	}
	else if (p[0].value >= maxisolevel && p[1].value < maxisolevel && p[2].value >= maxisolevel)
	{
		if (p[1].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index1, 0);
			getPointinsurface(Pindex[1], pointset, index2);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index3, 0);

			addFace(index1, index2, index3);
		}
		else if (p[1].value < minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index1, 0);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index2, 1);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index3, 1);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index4, 0);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
	}
	else if (p[0].value < maxisolevel && p[1].value >= maxisolevel && p[2].value >= maxisolevel)
	{
		if (p[0].value >= minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, index1);
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index2, 0);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index3, 0);

			addFace(index1, index2, index3);
		}
		else if (p[0].value < minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index1, 1);
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index2, 0);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index3, 0);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index4, 1);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
	}
	else if (p[0].value >= maxisolevel && p[1].value < maxisolevel && p[2].value < maxisolevel)
	{
		if (p[1].value >= minisolevel && p[2].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index1, 0);
			getPointinsurface(Pindex[1], pointset, index2);
			getPointinsurface(Pindex[2], pointset, index3);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index4, 0);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
		else if (p[1].value >= minisolevel && p[2].value < minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index1, 0);
			getPointinsurface(Pindex[1], pointset, index2);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index3, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index4, 1);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index5, 0);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
			addFace(index1, index4, index5);
		}
		else if (p[1].value < minisolevel && p[2].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index1, 0);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index2, 1);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index3, 1);
			getPointinsurface(Pindex[2], pointset, index4);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index5, 0);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
			addFace(index1, index4, index5);
		}
		else if (p[1].value < minisolevel && p[2].value < minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index1, 0);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index2, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index3, 1);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index4, 0);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
	}
	else if (p[0].value < maxisolevel && p[1].value >= maxisolevel && p[2].value < maxisolevel)
	{
		if (p[0].value >= minisolevel && p[2].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index1, 0);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index2, 0);
			getPointinsurface(Pindex[2], pointset, index3);
			getPointinsurface(Pindex[0], pointset, index4);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
		else if (p[0].value >= minisolevel && p[2].value < minisolevel)
		{
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index1, 0);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index2, 0);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index3, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index4, 1);
			getPointinsurface(Pindex[0], pointset, index5);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
			addFace(index1, index4, index5);
		}
		else if (p[0].value < minisolevel && p[2].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index1, 1);
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index2, 0);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index3, 0);
			getPointinsurface(Pindex[2], pointset, index4);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index5, 1);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
			addFace(index1, index4, index5);
		}
		else if (p[0].value < minisolevel && p[2].value < minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index1, 1);
			calInterVert2(surface, edges[0], maxisolevel, pointset, Edge, index2, 0);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index3, 0);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index4, 1);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
	}
	else if (p[0].value < maxisolevel && p[1].value < maxisolevel && p[2].value >= maxisolevel)
	{
		if (p[0].value >= minisolevel && p[1].value >= minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, index1);
			getPointinsurface(Pindex[1], pointset, index2);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index3, 0);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index4, 0);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
		else if (p[0].value >= minisolevel && p[1].value < minisolevel)
		{
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index1, 0);
			getPointinsurface(Pindex[0], pointset, index2);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index3, 1);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index4, 1);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index5, 0);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
			addFace(index1, index4, index5);
		}
		else if (p[0].value < minisolevel && p[1].value >= minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index1, 1);
			getPointinsurface(Pindex[1], pointset, index2);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index3, 0);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index4, 0);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index5, 1);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
			addFace(index1, index4, index5);
		}
		else if (p[0].value < minisolevel && p[1].value < minisolevel)
		{
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index1, 1);
			calInterVert2(surface, edges[2], maxisolevel, pointset, Edge, index2, 0);
			calInterVert2(surface, edges[1], maxisolevel, pointset, Edge, index3, 0);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index4, 1);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
	}
	else if (p[0].value < maxisolevel && p[1].value < maxisolevel && p[2].value < maxisolevel)
	{
		if (p[0].value > minisolevel && p[1].value > minisolevel && p[2].value > minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, index1);
			getPointinsurface(Pindex[1], pointset, index2);
			getPointinsurface(Pindex[2], pointset, index3);

			addFace(index1, index2, index3);
		}
		else if (p[0].value > minisolevel && p[1].value > minisolevel && p[2].value <= minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, index1);
			getPointinsurface(Pindex[1], pointset, index2);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index3, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index4, 1);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
		else if (p[0].value > minisolevel && p[1].value <= minisolevel && p[2].value > minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, index1);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index2, 1);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index3, 1);
			getPointinsurface(Pindex[2], pointset, index4);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
		else if (p[0].value <= minisolevel && p[1].value > minisolevel && p[2].value > minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index1, 1);
			getPointinsurface(Pindex[1], pointset, index2);
			getPointinsurface(Pindex[2], pointset, index3);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index4, 1);

			addFace(index1, index2, index3);
			addFace(index1, index3, index4);
		}
		else if (p[0].value > minisolevel && p[1].value <= minisolevel && p[2].value <= minisolevel)
		{
			getPointinsurface(Pindex[0], pointset, index1);
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index2, 1);
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index3, 1);

			addFace(index1, index2, index3);
		}
		else if (p[0].value <= minisolevel && p[1].value > minisolevel && p[2].value <= minisolevel)
		{
			calInterVert2(surface, edges[0], minisolevel, pointset, Edge, index1, 1);
			getPointinsurface(Pindex[1], pointset, index2);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index3, 1);

			addFace(index1, index2, index3);
		}
		else if (p[0].value <= minisolevel && p[1].value <= minisolevel && p[2].value > minisolevel)
		{
			calInterVert2(surface, edges[1], minisolevel, pointset, Edge, index1, 1);
			calInterVert2(surface, edges[2], minisolevel, pointset, Edge, index2, 1);
			getPointinsurface(Pindex[2], pointset, index3);

			addFace(index1, index2, index3);
		}
		else if (p[0].value <= minisolevel && p[1].value <= minisolevel && p[2].value <= minisolevel)
		{
		}
	}
}

void ScaffoldMesh::addFace(int f1, int f2, int f3)
{
	Triangle face;
	if (f1 != f2&&f1 != f3&&f2 != f3){
		face.f1 = f1;
		face.f2 = f2;
		face.f3 = f3;
		Face.push_back(face);
	}
}

void ScaffoldMesh::alter_face_direction()
{
	if (type_structure != 1)
	{
		for (int i = 0; i < Face.size(); i++)
		{
			int index_temp;
			for (int i = 0; i < Face.size(); i++)
			{
				index_temp = Face[i].f2;
				Face[i].f2 = Face[i].f3;
				Face[i].f3 = index_temp;
			}
		}
	}
}

void ScaffoldMesh::predictPorosity()
{
	size_t pointRes = resolution + 1; // indicates the # of points per side
	Array3D<double> VolumeGrid(resolution, resolution, resolution);
	Array3D<double> JacobianGrid(resolution, resolution, resolution);

	detJacobian(JacobianGrid);
	phyPosVolume(VolumeGrid);

	int dimDisTDF = 24;
	Array3D<Point3D> phyPorGrid(dimDisTDF, dimDisTDF, dimDisTDF);
	Array3D<Point3D> errorPorGrid(dimDisTDF, dimDisTDF, dimDisTDF);
	Array3D<Point3D> DisTDFGrid(dimDisTDF, dimDisTDF, dimDisTDF);

	setPhyPorosity(phyPorGrid, errorPorGrid, dimDisTDF);//设置模型孔隙率
	if (type_TPMS == 0)
	{
		Psurface surface;
		//phyPorosity2TDF(VolumeGrid, JacobianGrid, phyPorGrid, errorPorGrid, DisTDFGrid, dimDisTDF);
		for (int iter = 0; iter < 5; iter++)
		{
			phyPorosity2TDF(VolumeGrid, JacobianGrid, phyPorGrid, errorPorGrid, DisTDFGrid, dimDisTDF);
			outPutPorosity(surface, phyPorGrid, errorPorGrid, dimDisTDF);
		}

		ofstream ss("new_por.txt");
		ss << u_points << " " << v_points << " " << w_points << "\n";
		for (int i = 0; i < u_points; i++)
		{
			for (int j = 0; j < v_points; j++)
			{
				for (int k = 0; k < w_points; k++)
				{
					ss << TDF_control_grid[i][j][k] << "\n";
				}
			}
		}
		ss.close();
	}
	else if (type_TPMS == 1)
	{
		Gsurface surface;
	}
	else if (type_TPMS == 2)
	{
		Dsurface surface;
	}
	else if (type_TPMS == 3)
	{
		I_WPsurface surface;
		//phyPorosity2TDF(VolumeGrid, JacobianGrid, phyPorGrid, errorPorGrid, DisTDFGrid, dimDisTDF);
		for (int iter = 0; iter < 5; iter++)
		{
			phyPorosity2TDF(VolumeGrid, JacobianGrid, phyPorGrid, errorPorGrid, DisTDFGrid, dimDisTDF);
			outPutPorosity(surface, phyPorGrid, errorPorGrid, dimDisTDF);
		}

		ofstream ss("new_por.txt");
		ss << u_points << " " << v_points << " " << w_points << "\n";
		for (int i = 0; i < u_points; i++)
		{
			for (int j = 0; j < v_points; j++)
			{
				for (int k = 0; k < w_points; k++)
				{
					ss << TDF_control_grid[i][j][k] << "\n";
				}
			}
		}
		ss.close();
	}
}

void ScaffoldMesh::outPutPorosity(const Isosurface& surface, Array3D<Point3D>& phyPorGrid, Array3D<Point3D>& errorPorGrid, int dim)
{
	ofstream ss("debug.txt");
	ofstream ss1("debug1.txt");

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
			{
				double u = phyPorGrid.get(i, j, k).x, v = phyPorGrid.get(i, j, k).y, w = phyPorGrid.get(i, j, k).z;
				double Por = phyPorGrid.get(i, j, k).value;
				double errorPor = errorPorGrid.get(i, j, k).value;

				Point3D point_por;
				getPorosity(surface, Vector3D(u, v, w), point_por);
				errorPorGrid.set(i, j, k, Point3D(u, v, w, errorPor + Por - point_por.value));
				//ss << point_por.x << " " << point_por.y << " " << point_por.z << " " << Por << " " << point_por.value << " " << fabs(Por - point_por.value) / Por << endl;
				//ss << point_por.x << " " << point_por.y << " " << point_por.z << " " << point_por.value << endl;
				ss << point_por.value << endl;
				ss1 << errorPor + Por << endl;
			}
		}
	}
	ss.close();
	ss1.close();
}

void ScaffoldMesh::getPorosity(const Isosurface& surface, Vector3D& pos, Point3D& point_por)
{
	Vector3D phy_pos;
	calculate_posistion_TBSS_quick(pos, phy_pos);
	point_por = Point3D(phy_pos.x, phy_pos.y, phy_pos.z);

	Vector3D origin = pos - Vector3D(0.05, 0.05, 0.05);
	origin.x = origin.x > 0 ? origin.x : 0;
	origin.x = origin.x < 0.9 ? origin.x : 0.9;
	origin.y = origin.y > 0 ? origin.y : 0;
	origin.y = origin.y < 0.9 ? origin.y : 0.9;
	origin.z = origin.z > 0 ? origin.z : 0;
	origin.z = origin.z < 0.9 ? origin.z : 0.9;


	//Vector3D origin = pos - Vector3D(0.04, 0.04, 0.04);
	int seg = 11;
	Array3D<Point3D> phy_PointGrid(seg, seg, seg);
	for (int i = 0; i < seg; i++)
	{
		for (int j = 0; j < seg; j++)
		{
			for (int k = 0; k < seg; k++)
			{
				double cvalue;
				Vector3D phy;
				Vector3D para = origin + Vector3D(0.01*i, 0.01*j, 0.01*k);
				calculate_cvalue_TDF_quick(para, cvalue);
				calculate_posistion_TBSS_quick(para, phy);
				float value = surface.valueAt(omega_x * para.x, omega_y * para.y, omega_z * para.z) - cvalue;

				phy_PointGrid.set(i, j, k, Point3D(phy.x, phy.y, phy.z, value));
			}
		}
	}
	double sumVolume = 0;
	double sumValidVolume = 0;
	for (int i = 0; i < seg-1; i++)
	{
		for (int j = 0; j < seg-1; j++)
		{
			for (int k = 0; k < seg-1; k++)
			{
				Point3D v[8] = {
					{ phy_PointGrid.get(i, j, k) },
					{ phy_PointGrid.get(i + 1, j, k) },
					{ phy_PointGrid.get(i + 1, j + 1, k) },
					{ phy_PointGrid.get(i, j + 1, k) },
					{ phy_PointGrid.get(i, j, k + 1) },
					{ phy_PointGrid.get(i + 1, j, k + 1) },
					{ phy_PointGrid.get(i + 1, j + 1, k + 1) },
					{ phy_PointGrid.get(i, j + 1, k + 1) }
				};
				Point3D tetrahedra[6][4] = {
					{ v[0], v[7], v[3], v[2] },
					{ v[0], v[7], v[2], v[6] },
					{ v[0], v[4], v[7], v[6] },
					{ v[0], v[1], v[6], v[2] },
					{ v[0], v[4], v[6], v[1] },
					{ v[5], v[1], v[6], v[4] }
				};

				for (int t = 0; t < 6; ++t)
				{
					double volume, validvolume;
					calValidVolume(tetrahedra[t], volume, validvolume);
					sumVolume += volume;
					sumValidVolume += validvolume;
				}
			}
		}
	}
	if (sumVolume == 0 || sumValidVolume == 0)
	{
		point_por.value = 1.0;
	}
	else
	{
		point_por.value = 1.0 - sumValidVolume / sumVolume;
	}
}

void ScaffoldMesh::calValidVolume(Point3D v[4], double& volume, double& validvolume)
{
	Vector3D v0v1 = Vector3D(v[1].x - v[0].x, v[1].y - v[0].y, v[1].z - v[0].z);
	Vector3D v0v2 = Vector3D(v[2].x - v[0].x, v[2].y - v[0].y, v[2].z - v[0].z);
	Vector3D v0v3 = Vector3D(v[3].x - v[0].x, v[3].y - v[0].y, v[3].z - v[0].z);
	volume = fabs((v0v1.cross(v0v2)).dot(v0v3) / 6);
	//volume = fabs(v0v3.x*(v0v1.y*v0v2.z - v0v2.y*v0v1.z) + v0v3.y*(v0v1.x*v0v2.z - v0v2.x*v0v1.z) + v0v3.z*(v0v1.x*v0v2.y - v0v2.x*v0v1.y)) / 6;
	if (v[0].value <= 0 && v[1].value <= 0 && v[2].value <= 0 && v[3].value <= 0)
		validvolume = 0;
	else if (v[0].value > 0 && v[1].value > 0 && v[2].value > 0 && v[3].value > 0)
		validvolume = volume;
	else
		validvolume = volume / 2;
}

//计算非零B样条基函数及一阶二阶B样条基函数导矢
double** ScaffoldMesh::dersBasisFuns(int i, double u, std::vector<double> nodeU)
{
	// Input: i, u, p, nodeU
	// Output: ders
	double saved, temp;
	double* left = new double[4];
	double* right = new double[4];

	double ** ndu;
	double ** a;
	double ** ders;
	ndu = new double*[4];
	a = new double*[2];
	ders = new double*[3];
	for (int j = 0; j <= 3; j++)
	{
		ndu[j] = new double[4];
	}
	for (int j = 0; j <= 1; j++)
	{
		a[j] = new double[4];
	}
	for (int j = 0; j <= 2; j++)
	{
		ders[j] = new double[4];
	}

	ndu[0][0] = 1.0;
	for (int j = 1; j <= 3; j++)
	{
		left[j] = u - nodeU[i + 1 - j];
		right[j] = nodeU[i + j] - u;
		saved = 0.0;
		for (int r = 0; r < j; r++)
		{
			// Lower triangle
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];
			// Upper triangle
			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}
	for (int j = 0; j <= 3; j++)
	{
		// Load the basis functions
		ders[0][j] = ndu[j][3];
	}
	int s1, s2, rk, pk, j1, j2;
	double d;
	// This section computes the derivatives (Eq.[2.9])
	for (int r = 0; r <= 3; r++)
	{
		s1 = 0; s2 = 1; // Alternate rows in array a
		a[0][0] = 1.0;
		// Loop to compute k'th derivative
		for (int k = 1; k <= 2; k++)
		{
			d = 0.0;
			rk = r - k;
			pk = 3 - k;
			if (r >= k)
			{
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}
			if (rk >= -1)
				j1 = 1;
			else
				j1 = -rk;
			if (r - 1 <= pk)
				j2 = k - 1;
			else
				j2 = 3 - r;
			for (int j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
				d += a[s2][j] * ndu[rk + j][pk];
			}
			if (r <= pk)
			{
				a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
				d += a[s2][k] * ndu[r][pk];
			}
			ders[k][r] = d;
			// Switch rows
			int j = s1;
			s1 = s2;
			s2 = j;
		}
	}
	// Multiply through by the correct factors
	// (Eq.2.9)
	int r = 3;
	for (int k = 1; k <= 2; k++)
	{
		for (int j = 0; j <= 3; j++)
		{
			ders[k][j] *= r;
		}
		r *= (3 - k);
	}

	//free the memory
	for (int j = 0; j <= 3; j++)
	{
		delete[] ndu[j];
	}
	delete[] ndu;
	for (int j = 0; j <= 1; j++)
	{
		delete[] a[j];
	}
	delete[] a;
	delete[] left;
	delete[] right;

	return ders;
}

//计算离散网格上的Jacobian value
void ScaffoldMesh::detJacobian(Array3D<double>& JacobianGrid)
{
	//size_t pointRes = resolution + 1; // indicates the # of points per side
	//Array3D<double> JacobianGrid(pointRes, pointRes, pointRes);
	for (size_t i = 0; i < resolution; ++i) {
		double u = (double)(i+0.5) / resolution;
		for (size_t j = 0; j < resolution; ++j) {
			double v = (double)(j+0.5) / resolution;
			for (size_t k = 0; k < resolution; ++k) {
				double w = (double)(k+0.5) / resolution;
				Vector3D* Ders = solidDerivisOneOrder(Vector3D(u, v, w));
				double Jac = Ders[0].x * Ders[1].y * Ders[2].z + Ders[0].y * Ders[1].z * Ders[2].x + Ders[0].z * Ders[1].x * Ders[2].y
					- Ders[0].z * Ders[1].y * Ders[2].x - Ders[0].y * Ders[1].x * Ders[2].z - Ders[0].x * Ders[1].z * Ders[2].y;
				JacobianGrid.set(i, j, k, fabs(Jac));
			}
		}
	}
}

//计算体样条在三个方向上的一阶偏导矢
Vector3D* ScaffoldMesh::solidDerivisOneOrder(Vector3D para)
{
	Vector3D* der = new Vector3D[3];
	for (int i = 0; i < 3; i++)
	{
		der[i] = Vector3D(0, 0, 0);
	}
	Vector3D Mat_u[4][4];
	Vector3D Mat_v[4][4];
	Vector3D Mat_w[4][4];
	Vector3D vec_u[4];
	Vector3D vec_v[4];
	Vector3D vec_w[4];

	int uspan = findSpan(x_points, para.x, TBSS_knot_vector[0]);
	double** Nu = dersBasisFuns(uspan, para.x, TBSS_knot_vector[0]);
	int vspan = findSpan(y_points, para.y, TBSS_knot_vector[1]);
	double** Nv = dersBasisFuns(vspan, para.y, TBSS_knot_vector[1]);
	int wspan = findSpan(z_points, para.z, TBSS_knot_vector[2]);
	double** Nw = dersBasisFuns(wspan, para.z, TBSS_knot_vector[2]);

	int uind = uspan - 3;
	for (int k = 0; k <= 3; k++)
	{
		int wind = wspan - 3 + k;
		for (int j = 0; j <= 3; j++)
		{
			int vind = vspan - 3 + j;
			Mat_u[k][j] = Vector3D(0, 0, 0);
			Mat_v[k][j] = Vector3D(0, 0, 0);
			Mat_w[k][j] = Vector3D(0, 0, 0);
			for (int i = 0; i <= 3; i++)
			{
				Mat_u[k][j] += TBSS_control_grid[uind + i][vind][wind] * Nu[1][i];
				Mat_v[k][j] += TBSS_control_grid[uind + i][vind][wind] * Nu[0][i];
				Mat_w[k][j] += TBSS_control_grid[uind + i][vind][wind] * Nu[0][i];
			}
		}
	}
	for (int k = 0; k <= 3; k++)
	{
		vec_u[k] = Vector3D(0, 0, 0);
		vec_v[k] = Vector3D(0, 0, 0);
		vec_w[k] = Vector3D(0, 0, 0);
		for (int j = 0; j <= 3; j++)
		{
			vec_u[k] += Mat_u[k][j] * Nv[0][j];
			vec_v[k] += Mat_v[k][j] * Nv[1][j];
			vec_w[k] += Mat_w[k][j] * Nv[0][j];
		}
	}
	for (int k = 0; k <= 3; k++)
	{
		der[0] += vec_u[k] * Nw[0][k];
		der[1] += vec_v[k] * Nw[0][k];
		der[2] += vec_w[k] * Nw[1][k];
	}

	for (int i = 0; i < 3; i++)
	{
		delete[] Nu[i];
		delete[] Nv[i];
		delete[] Nw[i];
	}
	delete[] Nu;
	delete[] Nv;
	delete[] Nw;
	return der;
}

//计算每一个六面体单元的体积
void ScaffoldMesh::phyPosVolume(Array3D<double>& VolumeGrid)
{
	size_t pointRes = resolution + 1; // indicates the # of points per side
	Array3D<Point3D> PhyPosGrid(pointRes, pointRes, pointRes);
	for (size_t i = 0; i <= resolution; ++i) {
		double u = (double)i / resolution;
		for (size_t j = 0; j <= resolution; ++j) {
			double v = (double)j / resolution;
			for (size_t k = 0; k <= resolution; ++k) {
				double w = (double)k / resolution;
				Vector3D phy_pos;
				calculate_posistion_TBSS_quick(Vector3D(u, v, w), phy_pos);
				PhyPosGrid.set(i, j, k, Point3D(phy_pos.x, phy_pos.y, phy_pos.z));
			}
		}
	}
	double sumVolume = 0;
	for (int i = 0; i < resolution; i++)
	{
		for (int j = 0; j < resolution; j++)
		{
			for (int k = 0; k < resolution; k++)
			{
				Point3D v[8] = {
					{ PhyPosGrid.get(i, j, k) },
					{ PhyPosGrid.get(i + 1, j, k) },
					{ PhyPosGrid.get(i + 1, j + 1, k) },
					{ PhyPosGrid.get(i, j + 1, k) },
					{ PhyPosGrid.get(i, j, k + 1) },
					{ PhyPosGrid.get(i + 1, j, k + 1) },
					{ PhyPosGrid.get(i + 1, j + 1, k + 1) },
					{ PhyPosGrid.get(i, j + 1, k + 1) }
				};
				Point3D tetrahedra[6][4] = {
					{ v[0], v[7], v[3], v[2] },
					{ v[0], v[7], v[2], v[6] },
					{ v[0], v[4], v[7], v[6] },
					{ v[0], v[1], v[6], v[2] },
					{ v[0], v[4], v[6], v[1] },
					{ v[5], v[1], v[6], v[4] }
				};

				for (int t = 0; t < 6; ++t)
				{
					double volume;
					calVolume(tetrahedra[t], volume);
					sumVolume += volume;
				}
				VolumeGrid.set(i, j, k, sumVolume);
				sumVolume = 0;
			}
		}
	}
}

void ScaffoldMesh::calVolume(Point3D v[4], double& volume)
{
	Vector3D v0v1 = Vector3D(v[1].x - v[0].x, v[1].y - v[0].y, v[1].z - v[0].z);
	Vector3D v0v2 = Vector3D(v[2].x - v[0].x, v[2].y - v[0].y, v[2].z - v[0].z);
	Vector3D v0v3 = Vector3D(v[3].x - v[0].x, v[3].y - v[0].y, v[3].z - v[0].z);
	volume = fabs((v0v1.cross(v0v2)).dot(v0v3) / 6);
}

//预设初始的三维多孔模型上各点的孔隙率值
void ScaffoldMesh::setPhyPorosity(Array3D<Point3D>& phyPorGrid, Array3D<Point3D>& errorPorGrid, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
			{
				//int ii = 5 * i + 2, jj = 5 * j + 2, kk = 5 * k + 2;
				//int ii = 2 * i + 2, jj = 2 * j + 2, kk = 2 * k + 2;
				int ii = 4 * i + 4, jj = 4 * j + 4, kk = 4 * k + 4;
				//int ii = 8 * i + 8, jj = 8 * j + 8, kk = 8 * k + 8;
				double u = (double)ii / resolution, v = (double)jj / resolution, w = (double)kk / resolution;

				//double Por = w*0.3 + 0.3;
				double Por = 0.5;

				phyPorGrid.set(i, j, k, Point3D(u, v, w, Por));
				errorPorGrid.set(i, j, k, Point3D(u, v, w, 0));
			}
		}
	}
}

//基于预设孔隙率值和多次累计误差孔隙率映射到参数域构建离散TDF分布
void ScaffoldMesh::phyPorosity2TDF(Array3D<double>& VolumeGrid, Array3D<double>& JacobianGrid, Array3D<Point3D>& phyPorGrid, Array3D<Point3D>& errorPorGrid, Array3D<Point3D>& DisTDFGrid, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
			{
				//int ii = 5 * i + 2, jj = 5 * j + 2, kk = 5 * k + 2;
				//int ii = 2 * i + 2, jj = 2 * j + 2, kk = 2 * k + 2;
				int ii = 4 * i + 4, jj = 4 * j + 4, kk = 4 * k + 4;
				//int ii = 8 * i + 8, jj = 8 * j + 8, kk = 8 * k + 8;
				double Jac = JacobianGrid.get(ii, jj, kk);
				double Vol = VolumeGrid.get(ii, jj, kk);
				double phyPor = phyPorGrid.get(i, j, k).value;
				double errorPor = errorPorGrid.get(i, j, k).value;
				double newPor = phyPor + errorPor;
				//double newPor = (phyPor + errorPor)*Vol / (Jac*(1.0 / (resolution*resolution*resolution)));//0.000001是参数域六面体单元体积，Vol是参数域六面体单元映射到实体空间的单元体积；
				//double disTDF = -3.4644*newPor + 1.7318; //P_rod
				double disTDF = 3.4644*newPor - 1.7326;//P_pore

				//double disTDF = -7.5062*newPor + 3.9353;

				DisTDFGrid.set(i, j, k, Point3D(phyPorGrid.get(i, j, k).x, phyPorGrid.get(i, j, k).y, phyPorGrid.get(i, j, k).z, disTDF));
				//DisTDFGrid.set(i, j, k, Point3D(phyPorGrid.get(i, j, k).x, phyPorGrid.get(i, j, k).y, phyPorGrid.get(i, j, k).z, newPor));
			}
		}
	}
	LSPIA(DisTDFGrid, dim);
}

void ScaffoldMesh::LSPIA(Array3D<Point3D>& DisTDFGrid, int dim)
{
	double matri[4][4][4];//这个表示每一个基函数Bi*Bj*Bk
	int Bi_start_index, Bj_start_index, Bk_start_index;
	double cur_Por, cur_delta;

	double ***cube_delta;
	double ***cube_delta_deno;
	allocMem<double>(u_points, v_points, w_points, cube_delta, 0.0);
	allocMem<double>(u_points, v_points, w_points, cube_delta_deno, 0.0);

	for (int iter_num = 0; iter_num < 2; iter_num++)
	{
		for (int i = 0; i < u_points; i++)
		for (int j = 0; j < v_points; j++)
		for (int k = 0; k < w_points; k++)
		{
			cube_delta[i][j][k] = 0.0;
			cube_delta_deno[i][j][k] = 0.0;
		}

		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				for (int k = 0; k < dim; k++)
				{
					double u = DisTDFGrid.get(i, j, k).x, v = DisTDFGrid.get(i, j, k).y, w = DisTDFGrid.get(i, j, k).z;

					calculateBaseFunctionOfSolidBSpline(Vector3D(u, v, w), TDF_knot_vector, matri, Bi_start_index, Bj_start_index, Bk_start_index);
					cur_Por = 0;
					for (int ii = Bi_start_index; ii < Bi_start_index + 4; ii++)
					for (int jj = Bj_start_index; jj < Bj_start_index + 4; jj++)
					for (int kk = Bk_start_index; kk < Bk_start_index + 4; kk++)
					{
						cur_Por += TDF_control_grid[ii][jj][kk] * matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];
					}
					cur_delta = DisTDFGrid.get(i, j, k).value - cur_Por;

					for (int ii = Bi_start_index; ii < Bi_start_index + 4; ii++)
					for (int jj = Bj_start_index; jj < Bj_start_index + 4; jj++)
					for (int kk = Bk_start_index; kk < Bk_start_index + 4; kk++)
					{
						cube_delta_deno[ii][jj][kk] += matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];//分母和
						cube_delta[ii][jj][kk] += cur_delta*matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];//差向量
					}
				}
			}
		}
		for (int i = 0; i < u_points; i++)
		{
			for (int j = 0; j < v_points; j++)
			{
				for (int k = 0; k<w_points; k++)
				{
					if (fabs(cube_delta_deno[i][j][k])>1e-20)//小心除0
					{
						TDF_control_grid[i][j][k] += cube_delta[i][j][k] / cube_delta_deno[i][j][k];
					}
				}
			}
		}
	}

}

template<typename T>
static void ScaffoldMesh::allocMem(const int x_points, const int y_points, const int z_points, T ***&para_cube, T values)
{
	para_cube = new T**[x_points];
	for (int i = 0; i<x_points; i++)
	{
		para_cube[i] = new T*[y_points];
		for (int j = 0; j<y_points; j++)
		{
			para_cube[i][j] = new T[z_points];
			for (int k = 0; k<z_points; k++)
			{
				para_cube[i][j][k] = values;
			}
		}
	}
}

double ScaffoldMesh::randomUniform(double dMinValue, double dMaxValue)
{
	double pRandomValue = (double)(rand() / (double)RAND_MAX);
	pRandomValue = pRandomValue*(dMaxValue - dMinValue) + dMinValue;
	return pRandomValue;
}

void ScaffoldMesh::laplacianSmoothInner(double ***&cube, const int iter_num)
{
	for (int times = 0; times < iter_num; times++)
	{
		for (int i = 1; i < resolution; i++)
		{
			for (int j = 1; j < resolution; j++)
			{
				for (int k = 1; k < resolution; k++)
				{
					cube[i][j][k] = (cube[i - 1][j - 1][k - 1] + cube[i - 1][j - 1][k] + cube[i - 1][j - 1][k + 1]
						+ cube[i - 1][j][k - 1] + cube[i - 1][j][k] + cube[i - 1][j][k + 1]
						+ cube[i - 1][j + 1][k - 1] + cube[i - 1][j + 1][k] + cube[i - 1][j + 1][k + 1]
						+ cube[i][j - 1][k - 1] + cube[i][j - 1][k] + cube[i][j - 1][k + 1]
						+ cube[i][j][k - 1] + cube[i][j][k + 1]
						+ cube[i][j + 1][k - 1] + cube[i][j + 1][k] + cube[i][j + 1][k + 1]
						+ cube[i + 1][j - 1][k - 1] + cube[i + 1][j - 1][k] + cube[i + 1][j - 1][k + 1]
						+ cube[i + 1][j][k - 1] + cube[i + 1][j][k] + cube[i + 1][j][k + 1]
						+ cube[i + 1][j + 1][k - 1] + cube[i + 1][j + 1][k] + cube[i + 1][j + 1][k + 1]) / 26;
				}
			}
		}
	}
}

void ScaffoldMesh::laplacianSmoothOuter(double ***&cube, const int iter_num)
{
	for (int times = 0; times < iter_num; times++)
	{
		cube[0][0][0] = (cube[0][0][1] + cube[0][1][1] + cube[0][1][0] + cube[1][0][0] + cube[1][0][1] + cube[1][1][0] + cube[1][1][1]) / 7;
		cube[0][0][resolution] = (cube[0][0][resolution - 1] + cube[0][1][resolution - 1] + cube[0][1][resolution] + cube[1][0][resolution] + cube[1][0][resolution - 1] + cube[1][1][resolution] + cube[1][1][resolution - 1]) / 7;
		cube[0][resolution][0] = (cube[0][resolution][1] + cube[0][resolution - 1][1] + cube[0][resolution - 1][0] + cube[1][resolution][0] + cube[1][resolution][1] + cube[1][resolution - 1][0] + cube[1][resolution - 1][1]) / 7;
		cube[0][resolution][resolution] = (cube[0][resolution][resolution - 1] + cube[0][resolution - 1][resolution - 1] + cube[0][resolution - 1][resolution] + cube[1][resolution][resolution] + cube[1][resolution][resolution - 1] + cube[1][resolution - 1][resolution] + cube[1][resolution - 1][resolution - 1]) / 7;
		cube[resolution][0][0] = (cube[resolution][0][1] + cube[resolution][1][1] + cube[resolution][1][0] + cube[resolution - 1][0][0] + cube[resolution - 1][0][1] + cube[resolution - 1][1][0] + cube[resolution - 1][1][1]) / 7;
		cube[resolution][0][resolution] = (cube[resolution][0][resolution - 1] + cube[resolution][1][resolution - 1] + cube[resolution][1][resolution] + cube[resolution - 1][0][resolution] + cube[resolution - 1][0][resolution - 1] + cube[resolution - 1][1][resolution] + cube[resolution - 1][1][resolution - 1]) / 7;
		cube[resolution][resolution][0] = (cube[resolution][resolution][1] + cube[resolution][resolution - 1][1] + cube[resolution][resolution - 1][0] + cube[resolution - 1][resolution][0] + cube[resolution - 1][resolution][1] + cube[resolution - 1][resolution - 1][0] + cube[resolution - 1][resolution - 1][1]) / 7;
		cube[resolution][resolution][resolution] = (cube[resolution][resolution][resolution - 1] + cube[resolution][resolution - 1][resolution - 1] + cube[resolution][resolution - 1][resolution] + cube[resolution - 1][resolution][resolution] + cube[resolution - 1][resolution][resolution - 1] + cube[resolution - 1][resolution - 1][resolution] + cube[resolution - 1][resolution - 1][resolution - 1]) / 7;

		for (int i = 1; i < resolution; i++)
		{
			cube[i][0][0] = (cube[i - 1][0][0] + cube[i + 1][0][0] + cube[i][1][0] + cube[i][0][1]) / 4;
			cube[i][0][resolution] = (cube[i - 1][0][resolution] + cube[i + 1][0][resolution] + cube[i][1][resolution] + cube[i][0][resolution - 1]) / 4;
			cube[i][resolution][0] = (cube[i - 1][resolution][0] + cube[i + 1][resolution][0] + cube[i][resolution - 1][0] + cube[i][resolution][1]) / 4;
			cube[i][resolution][resolution] = (cube[i - 1][resolution][resolution] + cube[i + 1][resolution][resolution] + cube[i][resolution - 1][resolution] + cube[i][resolution][resolution - 1]) / 4;
		}
		for (int j = 1; j < resolution; j++)
		{
			cube[0][j][0] = (cube[0][j - 1][0] + cube[0][j + 1][0] + cube[1][j][0] + cube[0][j][1]) / 4;
			cube[resolution][j][0] = (cube[resolution][j - 1][0] + cube[resolution][j + 1][0] + cube[resolution - 1][j][0] + cube[resolution][j][1]) / 4;
			cube[0][j][resolution] = (cube[0][j - 1][resolution] + cube[0][j + 1][resolution] + cube[1][j][resolution] + cube[0][j][resolution - 1]) / 4;
			cube[resolution][j][resolution] = (cube[resolution][j - 1][resolution] + cube[resolution][j + 1][resolution] + cube[resolution - 1][j][resolution] + cube[resolution][j][resolution - 1]) / 4;
		}
		for (int k = 1; k < resolution; k++)
		{
			cube[0][0][k] = (cube[0][0][k - 1] + cube[0][0][k + 1] + cube[1][0][k] + cube[0][1][k]) / 4;
			cube[resolution][0][k] = (cube[resolution][0][k - 1] + cube[resolution][0][k + 1] + cube[resolution - 1][0][k] + cube[resolution][1][k]) / 4;
			cube[0][resolution][k] = (cube[0][resolution][k - 1] + cube[0][resolution][k + 1] + cube[1][resolution][k] + cube[0][resolution - 1][k]) / 4;
			cube[resolution][resolution][k] = (cube[resolution][resolution][k - 1] + cube[resolution][resolution][k + 1] + cube[resolution - 1][resolution][k] + cube[resolution][resolution - 1][k]) / 4;
		}
		for (int i = 1; i < resolution; i++)
		{
			for (int j = 1; j < resolution; j++)
			{
				cube[i][j][0] = (cube[i - 1][j - 1][0] + cube[i][j - 1][0] + cube[i + 1][j - 1][0]
					+ cube[i - 1][j][0] + cube[i + 1][j][0]
					+ cube[i - 1][j + 1][0] + cube[i][j + 1][0] + cube[i + 1][j + 1][0]) / 8;
				cube[i][j][resolution] = (cube[i - 1][j - 1][resolution] + cube[i][j - 1][resolution] + cube[i + 1][j - 1][resolution]
					+ cube[i - 1][j][resolution] + cube[i + 1][j][resolution]
					+ cube[i - 1][j + 1][resolution] + cube[i][j + 1][resolution] + cube[i + 1][j + 1][resolution]) / 8;
			}
		}
		for (int i = 1; i < resolution; i++)
		{
			for (int k = 1; k < resolution; k++)
			{
				cube[i][0][k] = (cube[i - 1][0][k - 1] + cube[i - 1][0][k] + cube[i - 1][0][k + 1]
					+ cube[i][0][k - 1] + cube[i][0][k + 1]
					+ cube[i + 1][0][k - 1] + cube[i + 1][0][k] + cube[i + 1][0][k + 1]) / 8;
				cube[i][resolution][k] = (cube[i - 1][resolution][k - 1] + cube[i - 1][resolution][k] + cube[i - 1][resolution][k + 1]
					+ cube[i][resolution][k - 1] + cube[i][resolution][k + 1]
					+ cube[i + 1][resolution][k - 1] + cube[i + 1][resolution][k] + cube[i + 1][resolution][k + 1]) / 8;
			}
		}
		for (int j = 1; j < resolution; j++)
		{
			for (int k = 1; k < resolution; k++)
			{
				cube[0][j][k] = (cube[0][j - 1][k - 1] + cube[0][j - 1][k] + cube[0][j - 1][k + 1]
					+ cube[0][j][k - 1] + cube[0][j][k + 1]
					+ cube[0][j + 1][k - 1] + cube[0][j + 1][k] + cube[0][j + 1][k + 1]) / 8;
				cube[resolution][j][k] = (cube[resolution][j - 1][k - 1] + cube[resolution][j - 1][k] + cube[resolution][j - 1][k + 1]
					+ cube[resolution][j][k - 1] + cube[resolution][j][k + 1]
					+ cube[resolution][j + 1][k - 1] + cube[resolution][j + 1][k] + cube[resolution][j + 1][k + 1]) / 8;
			}
		}
	}
}

void ScaffoldMesh::save_porosity(double*** &C, string filename)
{
	ofstream fout(filename + "_porosity.txt");
	for (int i = 0; i <= resolution; i++)
	{
		for (int j = 0; j <= resolution; j++)
		{
			for (int k = 0; k <= resolution; k++)
			{
				fout << C[i][j][k] << endl;
			}
		}
	}
	fout.close();
}