#include "porosityPereservation.h"

int omega_u, omega_v, omega_w;
int x_points, y_points, z_points;//是针对parameter每一维个数
vector<vector<double> > mesh_knot_vector;
vector<vector<vector<Vector3D>>> bspline_solid_control_mesh;

int u_points, v_points, w_points;
vector<vector<double> > tdf_knot_vector;
vector<vector<vector<double>>> tdf_control_mesh;
vector<vector<vector<double>>> partial_tdf;

//记录全局待求的TDF基函数值
vector<vector<vector<vector<vector<vector<vector<double>>>>>>> matrixBase;
vector<vector<vector<vector<vector<int>>>>> baseStartIndex;


//读取样条体控制网格和节点向量
void read_bspline_solid(const std::string & filename)
{
	ifstream meshfile(filename);
	string sline;
	getline(meshfile, sline);
	getline(meshfile, sline);
	istringstream sin(sline);
	sin >> x_points >> y_points >> z_points;
	bspline_solid_control_mesh.resize(x_points);
	for (int i = 0; i < x_points; i++)
	{
		bspline_solid_control_mesh[i].resize(y_points);
		for (int j = 0; j < y_points; j++)
		{
			bspline_solid_control_mesh[i][j].resize(z_points);
		}
	}
	mesh_knot_vector.resize(3);
	mesh_knot_vector[0].resize(x_points + 4);
	mesh_knot_vector[1].resize(y_points + 4);
	mesh_knot_vector[2].resize(z_points + 4);

	getline(meshfile, sline);
	for (int i = 0; i < x_points; i++)
	{
		for (int j = 0; j < y_points; j++)
		{
			for (int k = 0; k < z_points; k++)
			{
				getline(meshfile, sline);
				istringstream sin1(sline);
				sin1 >> bspline_solid_control_mesh[i][j][k].x >> bspline_solid_control_mesh[i][j][k].y >> bspline_solid_control_mesh[i][j][k].z;
			}
		}
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	istringstream sin2(sline);
	for (int i = 0; i < x_points + 4; i++)
	{
		sin2 >> mesh_knot_vector[0][i];
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	istringstream sin3(sline);
	for (int i = 0; i < y_points + 4; i++)
	{
		sin3 >> mesh_knot_vector[1][i];
	}
	getline(meshfile, sline);
	getline(meshfile, sline);
	istringstream sin4(sline);
	for (int i = 0; i < z_points + 4; i++)
	{
		sin4 >> mesh_knot_vector[2][i];
	}

}

//读取文件设置孔隙率
void setPorosity(vector<Point3D>& porosity, vector<Point3D>& Cvalue)
{
	//int num;
	//ifstream meshfile("required_porosity.txt");
	//string sline;
	//getline(meshfile, sline);
	//getline(meshfile, sline);
	//istringstream sin(sline);
	//sin >> num;
	//porosity.resize(num);
	//Cvalue.resize(num);
	//double u, v, w, por, c;
	//for (int i = 0; i < num; i++){
	//	getline(meshfile, sline);
	//	istringstream sin1(sline);
	//	sin1 >> u >> v >> w >> por;
	//	double c = 3.4644*por - 1.7322;//P_pore
	//	//double c = 2.996*por - 1.4974;//G_pore
	//	//double c = 1.6954*por - 0.8479;//D_pore
	//	//double c = 7.5062*por - 3.5709;//IWP_rod
	//	porosity.push_back(Point3D(u, v, w, por));
	//	Cvalue.push_back(Point3D(u, v, w, c));
	//}

	double x_step = 1.0 / omega_u;
	double y_step = 1.0 / omega_v;
	double z_step = 1.0 / omega_w;
	for (int i = 0; i < omega_u; i++)
	{
		double u = x_step*(i + 0.5);
		for (int j = 0; j < omega_v; j++)
		{
			double v = y_step*(j + 0.5);
			for (int k = 0; k < omega_w; k++)
			{
				double w = z_step*(k + 0.5);
				//double por = 0.4*w+0.3;
				double por = 0.5;
				double c = 3.4644*por - 1.7322;//P_pore
				//double c = 2.996*por - 1.4974;//G_pore
				//double c = 1.6954*por - 0.8479;//D_pore
				//double c = 7.5062*por - 3.5709;//IWP_rod

				if (v == 0.05){//只考虑一个面的孔隙率
					porosity.push_back(Point3D(u, v, w, por));
					Cvalue.push_back(Point3D(u, v, w, c));
				}
			}
		}
	}


}

void constructKnotVector(const int x_point_num, const int y_point_num, const int z_point_num, std::vector<std::vector<double > >&kont_vector)
{
	kont_vector.clear();
	kont_vector.resize(3);
	//首先计算节点向量
	for (int i = 0; i<3; i++)
	for (int j = 0; j<4; j++)
		kont_vector[i].push_back(0);//节点向量的前4个数都是0
	//x方向,均匀的节点向量
	for (int i = 1; i <= x_point_num - 4; i++)
		kont_vector[0].push_back((i*1.0) / (x_point_num - 3));
	//y方向,均匀的节点向量
	for (int i = 1; i <= y_point_num - 4; i++)
		kont_vector[1].push_back((i*1.0) / (y_point_num - 3));
	//x方向,均匀的节点向量
	for (int i = 1; i <= z_point_num - 4; i++)
		kont_vector[2].push_back((i*1.0) / (z_point_num - 3));
	for (int i = 0; i<3; i++)
	for (int j = 0; j<4; j++)
		kont_vector[i].push_back(1);//节点向量的后4个数都是1

}

void tdf_construction()
{
	tdf_control_mesh.resize(u_points);
	for (int i = 0; i < u_points; i++)
	{
		tdf_control_mesh[i].resize(v_points);
		for (int j = 0; j < v_points; j++)
		{
			tdf_control_mesh[i][j].resize(w_points);
			for (int k = 0; k < w_points; k++)
			{
				tdf_control_mesh[i][j][k] = 0;
			}
		}
	}
	constructKnotVector(u_points, v_points, w_points, tdf_knot_vector);

	partial_tdf.resize(u_points);
	for (int i = 0; i < u_points; i++)
	{
		partial_tdf[i].resize(v_points);
		for (int j = 0; j < v_points; j++)
		{
			partial_tdf[i][j].resize(w_points);
			for (int k = 0; k < w_points; k++)
			{
				partial_tdf[i][j][k] = 0;
			}
		}
	}
}

template<typename T>
void allocMem(const int points_x, const int points_y, const int points_z, T ***&para_cube, T values)
{
	para_cube = new T**[points_x];
	for (int i = 0; i<points_x; i++)
	{
		para_cube[i] = new T*[points_y];
		for (int j = 0; j<points_y; j++)
		{
			para_cube[i][j] = new T[points_z];
			for (int k = 0; k<points_z; k++)
			{
				para_cube[i][j][k] = values;
			}
		}
	}
}

double getNi4(const std::vector<double>&knot_vector, double t, int i, int seg, int zone, std::string  &deb)
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

void calculateBaseFunctionOfSolidBSpline(Vector3D para, std::vector<std::vector<double> >&kont_vector, double(*result)[4][4], int &Bi_start_index, int &Bj_start_index, int &Bk_start_index)
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

void LSPIAPorosity(vector<Point3D> required_porosity, int iterations)
{
	double matri[4][4][4];//这个表示每一个基函数Bi*Bj*Bk
	int Bi_start_index, Bj_start_index, Bk_start_index;
	double cur_Por, cur_delta;

	double ***cube_delta;
	double ***cube_delta_deno;
	allocMem<double>(u_points, v_points, w_points, cube_delta, 0.0);
	allocMem<double>(u_points, v_points, w_points, cube_delta_deno, 0.0);

	for (int iter_num = 0; iter_num < iterations; iter_num++)
	{
		for (int i = 0; i < u_points; i++)
		for (int j = 0; j < v_points; j++)
		for (int k = 0; k < w_points; k++)
		{
			cube_delta[i][j][k] = 0.0;
			cube_delta_deno[i][j][k] = 0.0;
		}

		for (int i = 0; i < required_porosity.size(); i++)
		{
			double u = required_porosity[i].x;
			double v = required_porosity[i].y;
			double w = required_porosity[i].z;
			double value = required_porosity[i].value;
			calculateBaseFunctionOfSolidBSpline(Vector3D(u, v, w), tdf_knot_vector, matri, Bi_start_index, Bj_start_index, Bk_start_index);
			cur_Por = 0;
			for (int ii = Bi_start_index; ii < Bi_start_index + 4; ii++)
			for (int jj = Bj_start_index; jj < Bj_start_index + 4; jj++)
			for (int kk = Bk_start_index; kk < Bk_start_index + 4; kk++)
			{
				cur_Por += tdf_control_mesh[ii][jj][kk] * matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];
			}
			cur_delta = value - cur_Por;

			for (int ii = Bi_start_index; ii < Bi_start_index + 4; ii++)
			for (int jj = Bj_start_index; jj < Bj_start_index + 4; jj++)
			for (int kk = Bk_start_index; kk < Bk_start_index + 4; kk++)
			{
				cube_delta_deno[ii][jj][kk] += matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];//分母和
				cube_delta[ii][jj][kk] += cur_delta*matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];//差向量
			}
		}
		for (int i = 0; i < u_points; i++)
		{
			for (int j = 0; j < v_points; j++)
			{
				for (int k = 0; k<w_points; k++)
				{
					if (fabs(cube_delta_deno[i][j][k])>MINN_MINN)//小心除0
					{
						tdf_control_mesh[i][j][k] += cube_delta[i][j][k] / cube_delta_deno[i][j][k];
					}
				}
			}
		}
	}

}

//利用控制网格求离散TDF
void calculate_cvalue_TDF_quick(Vector3D para, double &Cvalue)
{
	double matri[4][4][4];//这个表示每一个基函数Bi*Bj*Bk
	int Bi_start_index, Bj_start_index, Bk_start_index;
	calculateBaseFunctionOfSolidBSpline(para, tdf_knot_vector, matri, Bi_start_index, Bj_start_index, Bk_start_index);
	Cvalue = 0;
	for (int ii = Bi_start_index; ii<Bi_start_index + 4; ii++)
	for (int jj = Bj_start_index; jj<Bj_start_index + 4; jj++)
	for (int kk = Bk_start_index; kk<Bk_start_index + 4; kk++)
		Cvalue = Cvalue + tdf_control_mesh[ii][jj][kk] * matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];
}
//利用控制网格拟合样条体
void calculate_posistion_TBSS_quick(Vector3D para, Vector3D &phy_pos)
{
	double matri[4][4][4];//这个表示每一个基函数Bi*Bj*Bk
	int Bi_start_index, Bj_start_index, Bk_start_index;
	calculateBaseFunctionOfSolidBSpline(para, mesh_knot_vector, matri, Bi_start_index, Bj_start_index, Bk_start_index);
	phy_pos = Vector3D(0, 0, 0);
	for (int ii = Bi_start_index; ii<Bi_start_index + 4; ii++)
	for (int jj = Bj_start_index; jj<Bj_start_index + 4; jj++)
	for (int kk = Bk_start_index; kk<Bk_start_index + 4; kk++)
		phy_pos = phy_pos + bspline_solid_control_mesh[ii][jj][kk] * matri[ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];
}

void cal_tetVolume(Point3D v[4], double& volume)
{
	Vector3D v0v1 = Vector3D(v[1].x - v[0].x, v[1].y - v[0].y, v[1].z - v[0].z);
	Vector3D v0v2 = Vector3D(v[2].x - v[0].x, v[2].y - v[0].y, v[2].z - v[0].z);
	Vector3D v0v3 = Vector3D(v[3].x - v[0].x, v[3].y - v[0].y, v[3].z - v[0].z);
	volume = fabs((v0v1.cross(v0v2)).dot(v0v3) / 6);
}

double cal_volumeGrid(Array3D<Point3D>& phyGrid, Array3D<double>& volumeGrid, size_t dim)
{
	double cellVolume = 0;
	double voxelVolume = 0;
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			for (int k = 0; k < dim; k++)
			{
				Point3D v[8] = {
					{ phyGrid.get(i, j, k) },
					{ phyGrid.get(i + 1, j, k) },
					{ phyGrid.get(i + 1, j + 1, k) },
					{ phyGrid.get(i, j + 1, k) },
					{ phyGrid.get(i, j, k + 1) },
					{ phyGrid.get(i + 1, j, k + 1) },
					{ phyGrid.get(i + 1, j + 1, k + 1) },
					{ phyGrid.get(i, j + 1, k + 1) }
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
					cal_tetVolume(tetrahedra[t], volume);
					voxelVolume += volume;
				}
				volumeGrid.set(i, j, k, voxelVolume);
				cellVolume += voxelVolume;
				voxelVolume = 0;
			}
		}
	}
	return cellVolume;
}

void save_tdf_file(int type_TPMS, int type_structure, string filename)
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
	fout << omega_u << " " << omega_v << " " << omega_w << endl;
	fout << "#resolution of the control grid of TDF" << endl;
	fout << u_points << " " << v_points << " " << w_points;
	fout << "#control points of TDF" << endl;
	for (int i = 0; i < u_points; i++)
	{
		for (int j = 0; j < v_points; j++)
		{
			for (int k = 0; k < w_points; k++)
			{

				fout << tdf_control_mesh[i][j][k] << endl;
			}
		}
	}
	fout << "#knot vector in u-direction of TDF" << endl;
	for (int i = 0; i < tdf_knot_vector[0].size();i++){
		fout << tdf_knot_vector[0][i] << " ";
	}
	fout << "#knot vector in v-direction of TDF" << endl;
	for (int i = 0; i < tdf_knot_vector[1].size(); i++){
		fout << tdf_knot_vector[1][i] << " ";
	}
	fout << "#knot vector in w-direction of TDF" << endl;
	for (int i = 0; i < tdf_knot_vector[2].size(); i++){
		fout << tdf_knot_vector[2][i] << " ";
	}
	fout << "#resolution of the control grid of TBSS" << endl;
	fout << x_points << " " << y_points << " " << z_points << endl;
	fout << "#control points of TBSS" << endl;
	for (int i = 0; i < x_points; i++)
	{
		for (int j = 0; j < y_points; j++)
		{
			for (int k = 0; k < z_points; k++)
			{

				fout << bspline_solid_control_mesh[i][j][k].x << " " << bspline_solid_control_mesh[i][j][k].y << " " << bspline_solid_control_mesh[i][j][k].z << endl;
			}
		}
	}
	fout << "#knot vector in u-direction of TBSS" << endl;
	for (int i = 0; i < mesh_knot_vector[0].size(); i++){
		fout << mesh_knot_vector[0][i] << " ";
	}
	fout << "#knot vector in v-direction of TBSS" << endl;
	for (int i = 0; i < mesh_knot_vector[1].size(); i++){
		fout << mesh_knot_vector[1][i] << " ";
	}
	fout << "#knot vector in w-direction of TBSS" << endl;
	for (int i = 0; i < mesh_knot_vector[2].size(); i++){
		fout << mesh_knot_vector[2][i] << " ";
	}

	fout.close();
}

void cal_all_base_fun(vector<Point3D>& required_porosity, int resolution, double step[3]){
	matrixBase.resize(required_porosity.size());
	baseStartIndex.resize(required_porosity.size());
	for (int id = 0; id < required_porosity.size(); id++){
		matrixBase[id].resize(resolution);
		baseStartIndex[id].resize(resolution);
		for (int i = 0; i < resolution; i++){
			matrixBase[id][i].resize(resolution);
			baseStartIndex[id][i].resize(resolution);
			for (int j = 0; j < resolution; j++){
				matrixBase[id][i][j].resize(resolution);
				baseStartIndex[id][i][j].resize(resolution);
				for (int k = 0; k < resolution; k++){
					matrixBase[id][i][j][k].resize(4);
					baseStartIndex[id][i][j][k].resize(3);
					for (int ii = 0; ii < 4; ii++){
						matrixBase[id][i][j][k][ii].resize(4);
						for (int jj = 0; jj < 4; jj++){
							matrixBase[id][i][j][k][ii][jj].resize(4);
						}
					}
				}
			}
		}
	}

	double u_step = step[0];
	double v_step = step[1];
	double w_step = step[2];
	for (int id = 0; id < required_porosity.size(); id++)
	{
		Vector3D origin(required_porosity[id].x, required_porosity[id].y, required_porosity[id].z);
		for (int i = 0; i < resolution; i++)
		{
			for (int j = 0; j < resolution; j++)
			{
				for (int k = 0; k < resolution; k++)
				{
					Vector3D para = origin + Vector3D(u_step / resolution*(i + 0.5), v_step / resolution*(j + 0.5), w_step / resolution*(k + 0.5));	
					double matri[4][4][4];//这个表示每一个基函数Bi*Bj*Bk
					int Bi_start_index, Bj_start_index, Bk_start_index;
					calculateBaseFunctionOfSolidBSpline(para, tdf_knot_vector, matri, Bi_start_index, Bj_start_index, Bk_start_index);

					baseStartIndex[id][i][j][k][0] = Bi_start_index;
					baseStartIndex[id][i][j][k][1] = Bj_start_index;
					baseStartIndex[id][i][j][k][2] = Bk_start_index;
					for (int ii = 0; ii < 4; ii++){
						for (int jj = 0; jj < 4; jj++){
							for (int kk = 0; kk < 4; kk++){
								matrixBase[id][i][j][k][ii][jj][kk] = matri[ii][jj][kk];
							}
						}
					}
				}
			}
		}

	}

}


double cal_porosity(int type_TPMS, int type_structure, const Isosurface& surface, float isolevel, vector<Point3D>& required_porosity, vector<Array3D<double>>& phy_volume, vector<double>& volume, int id, int resolution, double step[3])
{
	double u_step = step[0];
	double v_step = step[1];
	double w_step = step[2];
	Vector3D origin(required_porosity[id].x, required_porosity[id].y, required_porosity[id].z);
	double poreVolume = 0;
	int count = 0;
	for (int i = 0; i < resolution; i++)
	{
		for (int j = 0; j < resolution; j++)
		{
			for (int k = 0; k < resolution; k++)
			{
				Vector3D para = origin + Vector3D(u_step / resolution*(i + 0.5), v_step / resolution*(j + 0.5), w_step / resolution*(k + 0.5));
				double Cvalue = 0;
				//calculate_cvalue_TDF_quick(para, cvalue);

				int Bi_start_index, Bj_start_index, Bk_start_index;
				Bi_start_index = baseStartIndex[id][i][j][k][0];
				Bj_start_index = baseStartIndex[id][i][j][k][1];
				Bk_start_index = baseStartIndex[id][i][j][k][2];
				for (int ii = Bi_start_index; ii < Bi_start_index + 4; ii++)
				for (int jj = Bj_start_index; jj < Bj_start_index + 4; jj++)
				for (int kk = Bk_start_index; kk < Bk_start_index + 4; kk++)
					Cvalue = Cvalue + tdf_control_mesh[ii][jj][kk] * matrixBase[id][i][j][k][ii - Bi_start_index][jj - Bj_start_index][kk - Bk_start_index];

				float value = surface.valueAt(omega_u * para.x, omega_v * para.y, omega_w * para.z) - Cvalue;
				if (value <= isolevel)
				{
					count++;
					poreVolume += phy_volume[id].get(i, j, k);
				}
			}
		}
	}
	if ((type_TPMS != 3 && type_structure == 1) || (type_TPMS == 3 && type_structure == 0)){
		//pore volume: value<=isolevel
		return poreVolume / volume[id];
	}
	else
	{
		//pore volume: value>=isolevel
		return 1 - poreVolume / volume[id];
	}
}

void cal_partial_energy(int type_TPMS, int type_structure, const Isosurface& surface, float isolevel, vector<Point3D>& required_porosity, vector<Array3D<double>>& phy_volume, vector<double>& volume, int resolution, double step[3])
{
	double energy = 0;
	for (int id = 0; id < required_porosity.size(); id++){
		double Por = cal_porosity(type_TPMS, type_structure, surface, isolevel, required_porosity, phy_volume, volume, id, resolution, step);
		energy += pow(Por - required_porosity[id].value, 2);
	}
	cout << "absolute error: " << energy / required_porosity.size() << endl;

	////all: calculate
	//double delta = 0.01;
	//for (int i = 0; i < u_points; i++){
	//	for (int j = 0; j < v_points; j++){
	//		for (int k = 0; k < w_points; k++){
	//			double E_delta = 0;
	//			tdf_control_mesh[i][j][k] += delta;
	//			for (int id = 0; id < required_porosity.size(); id++){
	//				double Por = cal_porosity(type_TPMS, type_structure, surface, isolevel, required_porosity, phy_volume, volume, id, resolution, step);
	//				E_delta += pow(Por - required_porosity[id].value, 2);
	//			}
	//			tdf_control_mesh[i][j][k] -= delta;
	//			partial_tdf[i][j][k] = (E_delta - energy) / delta;
	//		}
	//	}
	//}

	//part: calculate
	double delta = 0.02;
	for (int i = 0; i < u_points; i++){
		for (int j = 0; j < v_points; j++){
			for (int k = 0; k < w_points; k++){
				vector<int> index;
				double E=0, E_delta=0;
				for (int id = 0; id < required_porosity.size(); id++){
					Vector3D origin(required_porosity[id].x, required_porosity[id].y, required_porosity[id].z);
					if (origin.x*u_points>=(i - 3) && origin.x*u_points<=(i + 3) && origin.y*v_points>=(j - 3) && origin.y*v_points<=(j + 3) && origin.z*w_points >=(k - 3) && origin.z*w_points <=(k + 3)){
						double Por = cal_porosity(type_TPMS, type_structure, surface, isolevel, required_porosity, phy_volume, volume, id, resolution, step);
						E += pow(Por - required_porosity[id].value, 2);
						index.push_back(id);
					}
				}
				tdf_control_mesh[i][j][k] += delta;
				for (int in = 0; in < index.size(); in++){
					int id = index[in];
					double Por = cal_porosity(type_TPMS, type_structure, surface, isolevel, required_porosity, phy_volume, volume, id, resolution, step);
					E_delta += pow(Por - required_porosity[id].value, 2);
				}
				tdf_control_mesh[i][j][k] -= delta;

				partial_tdf[i][j][k] = (E_delta - E) / delta;
			}
		}
	}

}

void porosity_pereservation(const Isosurface& surface, float isolevel, int type_TPMS, int type_structure)
{
	//读取B样条控制网格
	read_bspline_solid("tooth_control_mesh.txt");
	
	//新增TDF控制网格
	u_points = 20;
	v_points = 20;
	w_points = 20;
	tdf_construction();

	omega_u = 10;
	omega_v = 10;
	omega_w = 10;
	vector<Point3D> required_porosity;
	vector<Point3D> modified_Cvalue;
	vector<Point3D> modified_porosity;
	setPorosity(required_porosity, modified_Cvalue);


	//利用孔隙率要求拟合
	LSPIAPorosity(modified_Cvalue, 10);


	vector<Array3D<Point3D>> phy_value;
	vector<Array3D<double>> phy_volume;
	vector<double> volume;

	//计算每一个孔隙率要求对应区域的所有体素体积（物理域上）
	size_t resolution = 15;
	size_t pointRes = resolution + 1;
	double u_step = 1.0 / omega_u, v_step = 1.0 / omega_v, w_step = 1.0 / omega_w;
	for (int id = 0; id < required_porosity.size(); id++){
		Vector3D pos(required_porosity[id].x, required_porosity[id].y, required_porosity[id].z);
		Vector3D origin = pos - Vector3D(0.5*u_step, 0.5*v_step, 0.5*w_step);
		origin.x = origin.x > 0 ? origin.x : 0;
		origin.x = origin.x < 1 - u_step ? origin.x : 1 - u_step;
		origin.y = origin.y > 0 ? origin.y : 0;
		origin.y = origin.y < 1 - v_step ? origin.y : 1 - v_step;
		origin.z = origin.z > 0 ? origin.z : 0;
		origin.z = origin.z < 1 - w_step ? origin.z : 1 - w_step;
		required_porosity[id].x = origin.x;
		required_porosity[id].y = origin.y;
		required_porosity[id].z = origin.z;
	}
	for (int id = 0; id < required_porosity.size(); id++)
	{
		Vector3D origin(required_porosity[id].x, required_porosity[id].y, required_porosity[id].z);
		Array3D<Point3D> phy_PointGrid(pointRes, pointRes, pointRes);
		for (int i = 0; i < pointRes; i++)
		{
			for (int j = 0; j < pointRes; j++)
			{
				for (int k = 0; k < pointRes; k++)
				{
					Vector3D phy;
					Vector3D para = origin + Vector3D(u_step / resolution*i, v_step / resolution*j, w_step / resolution*k);
					calculate_posistion_TBSS_quick(para, phy);
					phy_PointGrid.set(i, j, k, Point3D(phy.x, phy.y, phy.z, 0));
				}
			}
		}

		Array3D<double> volume_PointGrid(resolution, resolution, resolution);
		double cellVolume;
		cellVolume = cal_volumeGrid(phy_PointGrid, volume_PointGrid, resolution);
		phy_volume.push_back(volume_PointGrid);
		volume.push_back(cellVolume);
	}

	double step[3] = { u_step, v_step, w_step };
	cout << "iteration is starting" << endl;

	cal_all_base_fun(required_porosity, resolution, step);
	cout << "base function is calculated" << endl;

	double min_energy=0;
	for (int id = 0; id < required_porosity.size(); id++){
		double Por = cal_porosity(type_TPMS, type_structure, surface, isolevel, required_porosity, phy_volume, volume, id, resolution, step);
		min_energy += pow(Por - required_porosity[id].value, 2);
	}

	int iterations = 100;
	double max = 0.8, min = -0.8;//设置c值范围
	double ratio;
	int n = 20;
	clock_t time_start, time_end;
	time_start = clock();
	for (int iter = 0; iter < iterations; iter++){
		cal_partial_energy(type_TPMS, type_structure, surface, isolevel, required_porosity, phy_volume, volume, resolution, step);
		for (int t = 0; t < n; t++){
			for (int i = 0; i < u_points; i++){
				for (int j = 0; j < v_points; j++){
					for (int k = 0; k < w_points; k++){
						tdf_control_mesh[i][j][k] -= 1.0/n*partial_tdf[i][j][k];
						tdf_control_mesh[i][j][k] = tdf_control_mesh[i][j][k]>max ? max : tdf_control_mesh[i][j][k];
						tdf_control_mesh[i][j][k] = tdf_control_mesh[i][j][k] < min ? min : tdf_control_mesh[i][j][k];
					}
				}
			}
			double energy = 0;
			for (int id = 0; id < required_porosity.size(); id++){
				double Por = cal_porosity(type_TPMS, type_structure, surface, isolevel, required_porosity, phy_volume, volume, id, resolution, step);
				energy += pow(Por - required_porosity[id].value, 2);
			}
			if (energy<=min_energy)
			{
				ratio = 1.0 / n*(t + 1);
				min_energy = energy;
			}
		}

		for (int i = 0; i < u_points; i++){
			for (int j = 0; j < v_points; j++){
				for (int k = 0; k < w_points; k++){
					tdf_control_mesh[i][j][k] += (1-ratio)*partial_tdf[i][j][k];
					tdf_control_mesh[i][j][k] = tdf_control_mesh[i][j][k]>max ? max : tdf_control_mesh[i][j][k];
					tdf_control_mesh[i][j][k] = tdf_control_mesh[i][j][k] < min ? min : tdf_control_mesh[i][j][k];
				}
			}
		}
	}	
	time_end = clock();
	double duration = (static_cast<double>(time_end - time_start));//返回毫秒
	cout << duration << endl;
	save_tdf_file(type_TPMS, type_structure, "tooth_gradient");
}