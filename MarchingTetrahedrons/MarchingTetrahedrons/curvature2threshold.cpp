#include"curvature2threshold.h"
int x_num, y_num, z_num;//是针对parameter每一维个数
vector<vector<double> > mesh_knots;
vector<vector<vector<Vector3D>>> bsolid_cont_mesh;

int u_num, v_num, w_num;
vector<vector<double> > tdf_knots;
vector<vector<vector<double>>> tdf_cont_mesh;

int cube_cur_num_x, cube_cur_num_y, cube_cur_num_z;//是针对curvature离散的每一维个数
para_and_value ***cube_curvature;
std::vector<double> mean_curvature;//表面网格点的曲率


int rho_u, rho_v, rho_w;


//读取样条体控制网格和节点向量
void read_bsolid(const std::string & filename)
{
	ifstream meshfile(filename);
	string sline;
	getline(meshfile, sline);
	getline(meshfile, sline);
	istringstream sin(sline);
	sin >> x_num >> y_num >> z_num;
	bsolid_cont_mesh.resize(x_num);
	for (int i = 0; i < x_num; i++)
	{
		bsolid_cont_mesh[i].resize(y_num);
		for (int j = 0; j < y_num; j++)
		{
			bsolid_cont_mesh[i][j].resize(z_num);
		}
	}
	mesh_knots.resize(3);
	mesh_knots[0].resize(x_num + 4);
	mesh_knots[1].resize(y_num + 4);
	mesh_knots[2].resize(z_num + 4);

	getline(meshfile, sline);
	for (int i = 0; i < x_num; i++)
	{
		for (int j = 0; j < y_num; j++)
		{
			for (int k = 0; k < z_num; k++)
			{
				getline(meshfile, sline);
				istringstream sin1(sline);
				sin1 >> bsolid_cont_mesh[i][j][k].x >> bsolid_cont_mesh[i][j][k].y >> bsolid_cont_mesh[i][j][k].z;
			}
		}
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	istringstream sin2(sline);
	for (int i = 0; i < x_num + 4; i++)
	{
		sin2 >> mesh_knots[0][i];
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	istringstream sin3(sline);
	for (int i = 0; i < y_num + 4; i++)
	{
		sin3 >> mesh_knots[1][i];
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	istringstream sin4(sline);
	for (int i = 0; i < z_num + 4; i++)
	{
		sin4 >> mesh_knots[2][i];
	}

}

//绘制一个三角网格曲面
void create_triangular_mesh(const int num_x, const int num_y, const int num_z)
{
	std::map<int, int> surface_index_in_cube;
	int surface_index = 0;
	for (int i = 0; i < num_x; i++)
	{
		for (int j = 0; j < num_y; j++)
		{
			for (int k = 0; k < num_z; k++)
			{
				if (i == 0 || j == 0 || k == 0 || i == num_x - 1 || j == num_y - 1 || k == num_z - 1)
				{
					surface_index_in_cube[i*num_y*num_z + j*num_z + k] = surface_index;
					surface_index++;
				}
			}
		}
	}
	std::ofstream fout("surface.obj");
	double step_x = 1.0 / (num_x - 1), step_y = 1.0 / (num_y - 1), step_z = 1.0 / (num_z - 1);//点与点之间的间隔
	for (int i = 0; i < num_x; i++)
	{
		for (int j = 0; j < num_y; j++)
		{
			for (int k = 0; k < num_z; k++)
			{
				if (i == 0 || j == 0 || k == 0 || i == num_x - 1 || j == num_y - 1 || k == num_z - 1)
				{
					Vector3D para(step_x*i, step_y*j, step_z*k);
					Vector3D phy_pos;
					calculate_physcial_point(para, phy_pos);
					fout << "v " << phy_pos.x << " " << phy_pos.y << " " << phy_pos.z << "\n";
				}
			}
		}
	}

	//x-y plane
	for (int i = 0; i < (num_x - 1); i++)
	{
		for (int j = 0; j < (num_y - 1); j++)
		{
			//k=0
			fout << "f " << surface_index_in_cube[i * num_y * num_z + j * num_z] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + (j + 1) * num_z] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + j * num_z] + 1 << "\n";
			fout << "f " << surface_index_in_cube[i * num_y * num_z + j * num_z] + 1 << " " << surface_index_in_cube[i *num_y * num_z + (j + 1) * num_z] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + (j + 1) * num_z] + 1 << "\n";
			//k=num_z-1
			fout << "f " << surface_index_in_cube[i * num_y * num_z + j * num_z + (num_z - 1)] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + j * num_z + (num_z - 1)] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + (j + 1) * num_z + (num_z - 1)] + 1 << "\n";
			fout << "f " << surface_index_in_cube[i * num_y * num_z + j * num_z + (num_z - 1)] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + (j + 1) * num_z + (num_z - 1)] + 1 << " " << surface_index_in_cube[i * num_y * num_z + (j + 1) * num_z + (num_z - 1)] + 1 << "\n";
		}
	}
	//x-z plane
	for (int i = 0; i < (num_x - 1); i++)
	{
		for (int k = 0; k < (num_z - 1); k++)
		{
			//j=0
			fout << "f " << surface_index_in_cube[i * num_y * num_z + k] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + k] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + (k + 1)] + 1 << "\n";
			fout << "f " << surface_index_in_cube[i * num_y * num_z + k] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + (k + 1)] + 1 << " " << surface_index_in_cube[i *num_y * num_z + (k + 1)] + 1 << "\n";
			//j=num_y-1
			fout << "f " << surface_index_in_cube[i * num_y * num_z + k + (num_y - 1) * num_z] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + (k + 1) + (num_y - 1) * num_z] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + k + (num_y - 1) * num_z] + 1 << "\n";
			fout << "f " << surface_index_in_cube[i * num_y * num_z + k + (num_y - 1) * num_z] + 1 << " " << surface_index_in_cube[i *num_y * num_z + (k + 1) + (num_y - 1) * num_z] + 1 << " " << surface_index_in_cube[(i + 1) * num_y * num_z + (k + 1) + (num_y - 1) * num_z] + 1 << "\n";
		}
	}
	//y-z plane
	for (int j = 0; j < (num_y - 1); j++)
	{
		for (int k = 0; k < (num_z - 1); k++)
		{
			//i=0
			fout << "f " << surface_index_in_cube[j * num_z + k] + 1 << " " << surface_index_in_cube[(j + 1) * num_z + (k + 1)] + 1 << " " << surface_index_in_cube[(j + 1) * num_z + k] + 1 << "\n";
			fout << "f " << surface_index_in_cube[j * num_z + k] + 1 << " " << surface_index_in_cube[j * num_z + (k + 1)] + 1 << " " << surface_index_in_cube[(j + 1) * num_z + (k + 1)] + 1 << "\n";
			//i=num_x-1
			fout << "f " << surface_index_in_cube[j * num_z + k + (num_x - 1) * num_y * num_z] + 1 << " " << surface_index_in_cube[(j + 1) * num_z + k + (num_x - 1) * num_y * num_z] + 1 << " " << surface_index_in_cube[(j + 1) * num_z + (k + 1) + (num_x - 1) * num_y * num_z] + 1 << "\n";
			fout << "f " << surface_index_in_cube[j * num_z + k + (num_x - 1) * num_y * num_z] + 1 << " " << surface_index_in_cube[(j + 1) * num_z + (k + 1) + (num_x - 1) * num_y * num_z] + 1 << " " << surface_index_in_cube[j * num_z + (k + 1) + (num_x - 1) * num_y * num_z] + 1 << "\n";
		}
	}
	fout.close();

}

//计算曲率
void cal_mean_curvature(MyMesh*& mesh)
{
	double k_max, k_min;
	mean_curvature.clear();
	mesh->request_vertex_normals();
	mesh->update_normals();
	std::ofstream ss("curvature.txt");
	for (auto vit = mesh->vertices_begin(); vit != mesh->vertices_end(); ++vit)
	{
		k_max = -1000000;
		k_min = 1000000;
		for (auto vvit = mesh->vv_begin(*vit); vvit != mesh->vv_end(*vit); ++vvit)
		{
			MyMesh::Point v_vv = mesh->point(*vvit) - mesh->point(*vit);
			MyMesh::Point normal = mesh->normal(*vit);
			double k = -normal | v_vv *2.0 / pow(v_vv.norm(), 2);
			if (k > k_max) k_max = k;
			if (k < k_min) k_min = k;
		}
		double H = (k_max + k_min) / 2;
		ss << H << "\n";
		mean_curvature.push_back(H);
	}
	ss.close();
}

//样条曲面离散曲率
void discrete_curvature()
{
	cube_curvature = new para_and_value **[cube_cur_num_x];
	for (int i = 0; i < cube_cur_num_x; i++)
	{
		cube_curvature[i] = new para_and_value *[cube_cur_num_y];
		for (int j = 0; j < cube_cur_num_y; j++)
		{
			cube_curvature[i][j] = new para_and_value[cube_cur_num_z];
		}
	}
	int index = 0;
	for (int i = 0; i < cube_cur_num_x; i++)
	{
		double x_value = (i*1.0) / (cube_cur_num_x - 1);
		for (int j = 0; j < cube_cur_num_y; j++)
		{
			double y_value = (j*1.0) / (cube_cur_num_y - 1);
			for (int k = 0; k < cube_cur_num_z; k++)
			{
				double z_value = (k*1.0) / (cube_cur_num_z - 1);
				cube_curvature[i][j][k].para = OpenMesh::Vec3d(x_value, y_value, z_value);
				if (i == 0 || j == 0 || k == 0 || i == cube_cur_num_x - 1 || j == cube_cur_num_y - 1 || k == cube_cur_num_z - 1)
				{
					cube_curvature[i][j][k].value = mean_curvature[index];
					index++;
				}
				else
				{
					cube_curvature[i][j][k].value = 0;
				}
			}
		}
	}
	int seg[3] = { cube_cur_num_x - 1, cube_cur_num_y - 1, cube_cur_num_z - 1 };//这个表示cube被分成了几段
	laplacianSmoothOuter(cube_curvature, seg, 5);
	laplacianSmoothInner(cube_curvature, seg, 5);


	//std::ofstream ss("cube_curvature.txt");
	//ss << cube_cur_num_x << " " << cube_cur_num_y << " " << cube_cur_num_z << "\n";
	//for (int i = 0; i < cube_cur_num_x; i++)
	//{
	//	for (int j = 0; j < cube_cur_num_y; j++)
	//	{
	//		for (int k = 0; k < cube_cur_num_z; k++)
	//		{
	//			ss << cube_curvature[i][j][k] << "\n";
	//		}
	//	}
	//}
	//ss.close();
}

//曲率映射
void zoom_TDF(int type_TPMS, int type_structure)
{
	double maxC, minC;
	maxC = minC = fabs(cube_curvature[0][0][0].value);
	for (int i = 0; i < cube_cur_num_x; i++)
	{
		for (int j = 0; j < cube_cur_num_y; j++)
		{
			for (int k = 0; k < cube_cur_num_z; k++)
			{
				minC = minC<fabs(cube_curvature[i][j][k].value) ? minC : fabs(cube_curvature[i][j][k].value);
				maxC = maxC>fabs(cube_curvature[i][j][k].value) ? maxC : fabs(cube_curvature[i][j][k].value);
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

void scale_stretch(double begin_minC, double begin_maxC, double end_minC, double end_maxC)
{
	double begin_range = begin_maxC - begin_minC;
	double end_range = end_maxC - end_minC;
	for (int i = 0; i < cube_cur_num_x; i++)
	{
		for (int j = 0; j < cube_cur_num_y; j++)
		{
			for (int k = 0; k < cube_cur_num_z; k++)
			{
				cube_curvature[i][j][k].value = end_range*(fabs(cube_curvature[i][j][k].value) - begin_minC) / begin_range + end_minC;
			}
		}
	}
}


void laplacianSmoothInner(para_and_value ***&cube, const int *seg, const int iter_num)
{
	for (int times = 0; times < iter_num; times++)
	{
		for (int i = 1; i < seg[0]; i++)
		{
			for (int j = 1; j < seg[1]; j++)
			{
				for (int k = 1; k < seg[2]; k++)
				{
					cube[i][j][k].value = (cube[i - 1][j - 1][k - 1].value + cube[i - 1][j - 1][k].value + cube[i - 1][j - 1][k + 1].value
						+ cube[i - 1][j][k - 1].value + cube[i - 1][j][k].value + cube[i - 1][j][k + 1].value
						+ cube[i - 1][j + 1][k - 1].value + cube[i - 1][j + 1][k].value + cube[i - 1][j + 1][k + 1].value
						+ cube[i][j - 1][k - 1].value + cube[i][j - 1][k].value + cube[i][j - 1][k + 1].value
						+ cube[i][j][k - 1].value + cube[i][j][k + 1].value
						+ cube[i][j + 1][k - 1].value + cube[i][j + 1][k].value + cube[i][j + 1][k + 1].value
						+ cube[i + 1][j - 1][k - 1].value + cube[i + 1][j - 1][k].value + cube[i + 1][j - 1][k + 1].value
						+ cube[i + 1][j][k - 1].value + cube[i + 1][j][k].value + cube[i + 1][j][k + 1].value
						+ cube[i + 1][j + 1][k - 1].value + cube[i + 1][j + 1][k].value + cube[i + 1][j + 1][k + 1].value) / 26;
				}
			}
		}
	}
}
void laplacianSmoothOuter(para_and_value ***&cube, const int *seg, const int iter_num)
{
	for (int times = 0; times < iter_num; times++)
	{
		cube[0][0][0].value = (cube[0][0][1].value + cube[0][1][1].value + cube[0][1][0].value + cube[1][0][0].value + cube[1][0][1].value + cube[1][1][0].value + cube[1][1][1].value) / 7;
		cube[0][0][seg[2]].value = (cube[0][0][seg[2] - 1].value + cube[0][1][seg[2] - 1].value + cube[0][1][seg[2]].value + cube[1][0][seg[2]].value + cube[1][0][seg[2] - 1].value + cube[1][1][seg[2]].value + cube[1][1][seg[2] - 1].value) / 7;
		cube[0][seg[1]][0].value = (cube[0][seg[1]][1].value + cube[0][seg[1] - 1][1].value + cube[0][seg[1] - 1][0].value + cube[1][seg[1]][0].value + cube[1][seg[1]][1].value + cube[1][seg[1] - 1][0].value + cube[1][seg[1] - 1][1].value) / 7;
		cube[0][seg[1]][seg[2]].value = (cube[0][seg[1]][seg[2] - 1].value + cube[0][seg[1] - 1][seg[2] - 1].value + cube[0][seg[1] - 1][seg[2]].value + cube[1][seg[1]][seg[2]].value + cube[1][seg[1]][seg[2] - 1].value + cube[1][seg[1] - 1][seg[2]].value + cube[1][seg[1] - 1][seg[2] - 1].value) / 7;
		cube[seg[0]][0][0].value = (cube[seg[0]][0][1].value + cube[seg[0]][1][1].value + cube[seg[0]][1][0].value + cube[seg[0] - 1][0][0].value + cube[seg[0] - 1][0][1].value + cube[seg[0] - 1][1][0].value + cube[seg[0] - 1][1][1].value) / 7;
		cube[seg[0]][0][seg[2]].value = (cube[seg[0]][0][seg[2] - 1].value + cube[seg[0]][1][seg[2] - 1].value + cube[seg[0]][1][seg[2]].value + cube[seg[0] - 1][0][seg[2]].value + cube[seg[0] - 1][0][seg[2] - 1].value + cube[seg[0] - 1][1][seg[2]].value + cube[seg[0] - 1][1][seg[2] - 1].value) / 7;
		cube[seg[0]][seg[1]][0].value = (cube[seg[0]][seg[1]][1].value + cube[seg[0]][seg[1] - 1][1].value + cube[seg[0]][seg[1] - 1][0].value + cube[seg[0] - 1][seg[1]][0].value + cube[seg[0] - 1][seg[1]][1].value + cube[seg[0] - 1][seg[1] - 1][0].value + cube[seg[0] - 1][seg[1] - 1][1].value) / 7;
		cube[seg[0]][seg[1]][seg[2]].value = (cube[seg[0]][seg[1]][seg[2] - 1].value + cube[seg[0]][seg[1] - 1][seg[2] - 1].value + cube[seg[0]][seg[1] - 1][seg[2]].value + cube[seg[0] - 1][seg[1]][seg[2]].value + cube[seg[0] - 1][seg[1]][seg[2] - 1].value + cube[seg[0] - 1][seg[1] - 1][seg[2]].value + cube[seg[0] - 1][seg[1] - 1][seg[2] - 1].value) / 7;

		for (int i = 1; i < seg[0]; i++)
		{
			cube[i][0][0].value = (cube[i - 1][0][0].value + cube[i + 1][0][0].value + cube[i][1][0].value + cube[i][0][1].value) / 4;
			cube[i][0][seg[2]].value = (cube[i - 1][0][seg[2]].value + cube[i + 1][0][seg[2]].value + cube[i][1][seg[2]].value + cube[i][0][seg[2] - 1].value) / 4;
			cube[i][seg[1]][0].value = (cube[i - 1][seg[1]][0].value + cube[i + 1][seg[1]][0].value + cube[i][seg[1] - 1][0].value + cube[i][seg[1]][1].value) / 4;
			cube[i][seg[1]][seg[2]].value = (cube[i - 1][seg[1]][seg[2]].value + cube[i + 1][seg[1]][seg[2]].value + cube[i][seg[1] - 1][seg[2]].value + cube[i][seg[1]][seg[2] - 1].value) / 4;
		}
		for (int j = 1; j < seg[0]; j++)
		{
			cube[0][j][0].value = (cube[0][j - 1][0].value + cube[0][j + 1][0].value + cube[1][j][0].value + cube[0][j][1].value) / 4;
			cube[seg[0]][j][0].value = (cube[seg[0]][j - 1][0].value + cube[seg[0]][j + 1][0].value + cube[seg[0] - 1][j][0].value + cube[seg[0]][j][1].value) / 4;
			cube[0][j][seg[2]].value = (cube[0][j - 1][seg[2]].value + cube[0][j + 1][seg[2]].value + cube[1][j][seg[2]].value + cube[0][j][seg[2] - 1].value) / 4;
			cube[seg[0]][j][seg[2]].value = (cube[seg[0]][j - 1][seg[2]].value + cube[seg[0]][j + 1][seg[2]].value + cube[seg[0] - 1][j][seg[2]].value + cube[seg[0]][j][seg[2] - 1].value) / 4;
		}
		for (int k = 1; k < seg[0]; k++)
		{
			cube[0][0][k].value = (cube[0][0][k - 1].value + cube[0][0][k + 1].value + cube[1][0][k].value + cube[0][1][k].value) / 4;
			cube[seg[0]][0][k].value = (cube[seg[0]][0][k - 1].value + cube[seg[0]][0][k + 1].value + cube[seg[0] - 1][0][k].value + cube[seg[0]][1][k].value) / 4;
			cube[0][seg[2]][k].value = (cube[0][seg[2]][k - 1].value + cube[0][seg[2]][k + 1].value + cube[1][seg[2]][k].value + cube[0][seg[2] - 1][k].value) / 4;
			cube[seg[0]][seg[2]][k].value = (cube[seg[0]][seg[2]][k - 1].value + cube[seg[0]][seg[2]][k + 1].value + cube[seg[0] - 1][seg[2]][k].value + cube[seg[0]][seg[2] - 1][k].value) / 4;
		}
		for (int i = 1; i < seg[0]; i++)
		{
			for (int j = 1; j < seg[1]; j++)
			{
				cube[i][j][0].value = (cube[i - 1][j - 1][0].value + cube[i][j - 1][0].value + cube[i + 1][j - 1][0].value
					+ cube[i - 1][j][0].value + cube[i + 1][j][0].value
					+ cube[i - 1][j + 1][0].value + cube[i][j + 1][0].value + cube[i + 1][j + 1][0].value) / 8;
				cube[i][j][seg[2]].value = (cube[i - 1][j - 1][seg[2]].value + cube[i][j - 1][seg[2]].value + cube[i + 1][j - 1][seg[2]].value
					+ cube[i - 1][j][seg[2]].value + cube[i + 1][j][seg[2]].value
					+ cube[i - 1][j + 1][seg[2]].value + cube[i][j + 1][seg[2]].value + cube[i + 1][j + 1][seg[2]].value) / 8;
			}
		}
		for (int i = 1; i < seg[0]; i++)
		{
			for (int k = 1; k < seg[2]; k++)
			{
				cube[i][0][k].value = (cube[i - 1][0][k - 1].value + cube[i - 1][0][k].value + cube[i - 1][0][k + 1].value
					+ cube[i][0][k - 1].value + cube[i][0][k + 1].value
					+ cube[i + 1][0][k - 1].value + cube[i + 1][0][k].value + cube[i + 1][0][k + 1].value) / 8;
				cube[i][seg[1]][k].value = (cube[i - 1][seg[1]][k - 1].value + cube[i - 1][seg[1]][k].value + cube[i - 1][seg[1]][k + 1].value
					+ cube[i][seg[1]][k - 1].value + cube[i][seg[1]][k + 1].value
					+ cube[i + 1][seg[1]][k - 1].value + cube[i + 1][seg[1]][k].value + cube[i + 1][seg[1]][k + 1].value) / 8;
			}
		}
		for (int j = 1; j < seg[1]; j++)
		{
			for (int k = 1; k < seg[2]; k++)
			{
				cube[0][j][k].value = (cube[0][j - 1][k - 1].value + cube[0][j - 1][k].value + cube[0][j - 1][k + 1].value
					+ cube[0][j][k - 1].value + cube[0][j][k + 1].value
					+ cube[0][j + 1][k - 1].value + cube[0][j + 1][k].value + cube[0][j + 1][k + 1].value) / 8;
				cube[seg[0]][j][k].value = (cube[seg[0]][j - 1][k - 1].value + cube[seg[0]][j - 1][k].value + cube[seg[0]][j - 1][k + 1].value
					+ cube[seg[0]][j][k - 1].value + cube[seg[0]][j][k + 1].value
					+ cube[seg[0]][j + 1][k - 1].value + cube[seg[0]][j + 1][k].value + cube[seg[0]][j + 1][k + 1].value) / 8;
			}
		}
	}
}


//利用控制网格拟合样条体
void calculate_physcial_point(Vector3D para, Vector3D &phy_pos)
{
	double matri[4][4][4];//这个表示每一个基函数Bi*Bj*Bk
	int Bi_start_index, Bj_start_index, Bk_start_index;
	calculateBaseFunctionOfSolidBSpline(para, mesh_knots, matri, Bi_start_index, Bj_start_index, Bk_start_index);
	phy_pos = Vector3D(0, 0, 0);
	for (int ii = Bi_start_index; ii<Bi_start_index + 4; ii++)
	for (int jj = Bj_start_index; jj<Bj_start_index + 4; jj++)
	for (int kk = Bk_start_index; kk<Bk_start_index + 4; kk++)
		phy_pos = phy_pos + bsolid_cont_mesh[ii][jj][kk] * matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];
}



void construct_TDF()
{
	tdf_cont_mesh.resize(u_num);
	for (int i = 0; i < u_num; i++)
	{
		tdf_cont_mesh[i].resize(v_num);
		for (int j = 0; j < v_num; j++)
		{
			tdf_cont_mesh[i][j].resize(w_num);
			for (int k = 0; k < w_num; k++)
			{
				tdf_cont_mesh[i][j][k] = 0;
			}
		}
	}
	constructKnotVector(u_num, v_num, w_num, tdf_knots);
}


void LSPIACurvature(para_and_value ***&cube, int iterations)
{
	double matri[4][4][4];//这个表示每一个基函数Bi*Bj*Bk
	int Bi_start_index, Bj_start_index, Bk_start_index;
	double cur_value, cur_delta;

	double ***cube_delta;
	double ***cube_delta_deno;
	allocMem<double>(u_num, v_num, w_num, cube_delta, 0.0);
	allocMem<double>(u_num, v_num, w_num, cube_delta_deno, 0.0);

	for (int iter_num = 0; iter_num < iterations; iter_num++)
	{
		for (int i = 0; i < u_num; i++)
		for (int j = 0; j < v_num; j++)
		for (int k = 0; k < w_num; k++)
		{
			cube_delta[i][j][k] = 0.0;
			cube_delta_deno[i][j][k] = 0.0;
		}

		for (int i = 0; i < cube_cur_num_x; i++)
		for (int j = 0; j < cube_cur_num_y; j++)
		for (int k = 0; k < cube_cur_num_z; k++)
		{
			double u = cube_curvature[i][j][k].para[0];
			double v = cube_curvature[i][j][k].para[1];
			double w = cube_curvature[i][j][k].para[2];
			double value = cube_curvature[i][j][k].value;
			calculateBaseFunctionOfSolidBSpline(Vector3D(u, v, w), tdf_knots, matri, Bi_start_index, Bj_start_index, Bk_start_index);
			cur_value = 0;
			for (int ii = Bi_start_index; ii < Bi_start_index + 4; ii++)
			for (int jj = Bj_start_index; jj < Bj_start_index + 4; jj++)
			for (int kk = Bk_start_index; kk < Bk_start_index + 4; kk++)
			{
				cur_value += tdf_cont_mesh[ii][jj][kk] * matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];
			}
			cur_delta = value - cur_value;

			for (int ii = Bi_start_index; ii < Bi_start_index + 4; ii++)
			for (int jj = Bj_start_index; jj < Bj_start_index + 4; jj++)
			for (int kk = Bk_start_index; kk < Bk_start_index + 4; kk++)
			{
				cube_delta_deno[ii][jj][kk] += matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];//分母和
				cube_delta[ii][jj][kk] += cur_delta*matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];//差向量
			}
		}
		for (int i = 0; i < u_num; i++)
		{
			for (int j = 0; j < v_num; j++)
			{
				for (int k = 0; k< w_num; k++)
				{
					if (fabs(cube_delta_deno[i][j][k])>MINN_MINN)//小心除0
					{
						tdf_cont_mesh[i][j][k] += cube_delta[i][j][k] / cube_delta_deno[i][j][k];
					}
				}
			}
		}
	}

}


void save_TDF_file(int type_TPMS, int type_structure,string filename)
{
	//ofstream fout(filename + "_tdf_cont_mesh.txt");
	//for (int i = 0; i < u_points; i++)
	//{
	//	for (int j = 0; j < v_points; j++)
	//	{
	//		for (int k = 0; k < w_points; k++)
	//		{
	//			fout << tdf_control_mesh[i][j][k] << endl;
	//		}
	//	}
	//}
	//fout.close();
	ofstream fout(filename + "_porous_scaffold.tdf");
	fout << "#TPMS type, m=0: P-type; m=1: G-type; m=2: D-type; m=3: I-WP-type" << endl;
	fout << type_TPMS << endl;
	fout << "#structure type, n=0: rod structure; n=1: pore structure; n=2: sheet structure" << endl;
	fout << type_structure << endl;
	fout << "#period coefficients (rho_u,rho_v,rho_w)" << endl;
	fout << rho_u << " " << rho_v << " " << rho_w << endl;
	fout << "#resolution of the control grid of TDF" << endl;
	fout << u_num << " " << v_num << " " << w_num;
	fout << "#control points of TDF" << endl;
	for (int i = 0; i < u_num; i++)
	{
		for (int j = 0; j < v_num; j++)
		{
			for (int k = 0; k < w_num; k++)
			{

				fout << tdf_cont_mesh[i][j][k] << endl;
			}
		}
	}
	fout << "#knot vector in u-direction of TDF" << endl;
	for (int i = 0; i < tdf_knots[0].size(); i++){
		fout << tdf_knots[0][i] << " ";
	}
	fout << "#knot vector in v-direction of TDF" << endl;
	for (int i = 0; i < tdf_knots[1].size(); i++){
		fout << tdf_knots[1][i] << " ";
	}
	fout << "#knot vector in w-direction of TDF" << endl;
	for (int i = 0; i < tdf_knots[2].size(); i++){
		fout << tdf_knots[2][i] << " ";
	}
	fout << "#resolution of the control grid of TBSS" << endl;
	fout << x_num << " " << y_num << " " << z_num << endl;
	fout << "#control points of TBSS" << endl;
	for (int i = 0; i < x_num; i++)
	{
		for (int j = 0; j < y_num; j++)
		{
			for (int k = 0; k < z_num; k++)
			{

				fout << bsolid_cont_mesh[i][j][k].x << " " << bsolid_cont_mesh[i][j][k].y << " " << bsolid_cont_mesh[i][j][k].z << endl;
			}
		}
	}
	fout << "#knot vector in u-direction of TBSS" << endl;
	for (int i = 0; i < mesh_knots[0].size(); i++){
		fout << mesh_knots[0][i] << " ";
	}
	fout << "#knot vector in v-direction of TBSS" << endl;
	for (int i = 0; i < mesh_knots[1].size(); i++){
		fout << mesh_knots[1][i] << " ";
	}
	fout << "#knot vector in w-direction of TBSS" << endl;
	for (int i = 0; i < mesh_knots[2].size(); i++){
		fout << mesh_knots[2][i] << " ";
	}
	fout.close();
}

void curvature2threshold(int typeTPMS,int typeStructure){
	read_bsolid("balljoint_control_mesh.txt");

	cube_cur_num_x = 101;
	cube_cur_num_y = 101;
	cube_cur_num_z = 101;
	create_triangular_mesh(cube_cur_num_x, cube_cur_num_y, cube_cur_num_z);

	MyMesh* tri_mesh;
	tri_mesh = NULL;
	tri_mesh = new MyMesh;
	OpenMesh::IO::Options opt;
	OpenMesh::IO::read_mesh(*tri_mesh, "surface.obj",opt);
	cal_mean_curvature(tri_mesh);

	cube_curvature = NULL;
	discrete_curvature();
	zoom_TDF(typeTPMS, typeStructure);

	//新增TDF控制网格
	u_num = 20;
	v_num = 20;
	w_num = 20;
	construct_TDF();

	//利用孔隙率要求拟合
	LSPIACurvature(cube_curvature, 10);


	rho_u = 10;
	rho_v = 10;
	rho_w = 10;
	save_TDF_file(typeTPMS, typeStructure,"balljoint");
}