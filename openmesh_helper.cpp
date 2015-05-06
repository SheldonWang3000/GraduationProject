#include "openmesh_helper.h"
#include "Eigen/Dense"
#include <iostream>
using namespace std;
//************************************
// Method:    OpenmeshHelper
// FullName:  OpenmeshHelper::OpenmeshHelper
// Access:    public 
// Returns:   
// Qualifier: 
// Parameter: string input_location input文件的位置
//************************************
OpenmeshHelper::OpenmeshHelper(string input_location)
{
	IplImage *image = cvLoadImage(input_location.c_str(), -1);
	this->image = image;
	cout << "Init complete" << endl;
}
OpenmeshHelper::~OpenmeshHelper()
{
	cvReleaseImage(&image);
	cout << "Destory complete" << endl;
}
//************************************
// Method:    InitPointData
// FullName:  OpenmeshHelper::InitPointData
// Access:    private 
// Returns:   void
// Qualifier:
// 初始化像素点数组，每个点代表一个原图像素点，包含颜色值，是否为撕裂点，以及openmesh的vertex对象
//************************************
void OpenmeshHelper::InitPointData()
{
	point_data = new MyMesh::VertexHandle *[image->height];
	for (int i = 0; i != image->height; ++i)
	{
		point_data[i] = new MyMesh::VertexHandle[image->width];
	}
	OpencvHelper helper;

	IplImage *edge = cvCreateImage(cvGetSize(image), image->depth, 1);
	
	helper.get_pic_edge(image, edge);
	uchar *data = (uchar*)edge->imageData;
	int step = edge->widthStep / sizeof(uchar);

	for (int i = 0; i < image->height; ++i)
	{
		for (int j = 0; j < image->width; ++j)
		{
			auto point = mesh.add_vertex(MyMesh::Point(i, j, 0));
			MyMesh::TexCoord3D t;
			t[0] = i;
			t[1] = j;
			t[2] = None;
			mesh.set_texcoord3D(point, t);
			auto color = cvGet2D(image, i, j);
			MyMesh::Color vertex_color;
			for (int k = 0; k != 3; ++k)
			{
				vertex_color[k] = (int)color.val[k];
			}
			mesh.set_color(point, vertex_color);
			if (data[i * step + j] == 255)
				mesh.status(point).set_feature(true);
			else
				mesh.status(point).set_feature(false);
			bool isEdge = mesh.status(point).feature();
			point_data[i][j] = point;
		}
	}
	ofstream out("D:/edge.txt");
	for (int i = 0; i < image->height; ++i)
	{
		for (int j = 0; j < image->width; ++j)
		{
			out << mesh.status(point_data[i][j]).feature() << ' ';
		}
		out << endl;
	}
	out.close();

	cvReleaseImage(&edge);
	cout << "Init point data complete" << endl;
}
//************************************
// Method:    ConnectMesh
// FullName:  OpenmeshHelper::ConnectMesh
// Access:    private 
// Returns:   void
// Qualifier:
// 连接网格，遵循之前的规则将点连接成网格
//************************************
void OpenmeshHelper::ConnectMesh(bool isContour)
{
	vector<MyMesh::VertexHandle> list;

	for (int i = 0; i != image->height - 1; ++ i)
	{
		for (int j = 0; j != image->width - 1; ++ j)
		{
			MyMesh::VertexHandle point[4];
			point[0] = point_data[i][j];
			point[1] = point_data[i + 1][j];
			point[2] = point_data[i + 1][j + 1];
			point[3] = point_data[i][j + 1];
			int num = 0;
			for (int k = 0; k < 4; ++k)
			{
				if (mesh.status(point[k]).feature())
					++num;
			}

			switch (num)
			{
			case 0:
				if (isContour)
				{
					list.clear();
					list.push_back(point[0]);
					list.push_back(point[1]);
					list.push_back(point[2]);
					list.push_back(point[3]);
					mesh.add_face(list);
				}
				break;
			case 1:
				list.clear();
				list.push_back(point[0]);
				list.push_back(point[1]);
				list.push_back(point[2]);
				list.push_back(point[3]);
				mesh.add_face(list);
				break;
			case 2:
				if (mesh.status(point[0]).feature() && mesh.status(point[2]).feature())
				{
					list.clear();
					list.push_back(point[0]);
					list.push_back(point[1]);
					list.push_back(point[2]);
					mesh.add_face(list);

					list.clear();
					list.push_back(point[2]);
					list.push_back(point[3]);
					list.push_back(point[0]);
					mesh.add_face(list);
				}
				else
				{
					if (mesh.status(point[1]).feature() && mesh.status(point[3]).feature())
					{
						list.clear();
						list.push_back(point[0]);
						list.push_back(point[1]);
						list.push_back(point[3]);
						mesh.add_face(list);

						list.clear();
						list.push_back(point[1]);
						list.push_back(point[2]);
						list.push_back(point[3]);
						mesh.add_face(list);
					}
					else
					{
						list.clear();
						list.push_back(point[0]);
						list.push_back(point[1]);
						list.push_back(point[2]);
						list.push_back(point[3]);
						mesh.add_face(list);
					}
				}
				break;
			case 3:
				list.clear();
				list.push_back(point[0]);
				list.push_back(point[1]);
				list.push_back(point[2]);
				list.push_back(point[3]);
				mesh.add_face(list);
				break;
			case 4:
				list.clear();
				list.push_back(point[0]);
				list.push_back(point[1]);
				list.push_back(point[2]);
				list.push_back(point[3]);
				mesh.add_face(list);
				break;
			default:
				break;
			}
		}
	}
	mesh.triangulate();
	cout << "Connect mesh complete" << endl;
}
//************************************
// Method:    CountVertices
// FullName:  OpenmeshHelper::CountVertices
// Access:    private 
// Returns:   void
// Qualifier:
// 计算openmesh中顶点数，用于停止简化
//************************************
void OpenmeshHelper::CountVertices()
{
	num_vertices = 0;
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++ i)
	{
		++ num_vertices;
	}
	num_all_vertices = num_vertices;
	cout << "Count vertices complete " << num_all_vertices << endl;
}
//************************************
// Method:    InitPairList
// FullName:  OpenmeshHelper::InitPairList
// Access:    private 
// Returns:   void
// Qualifier:
// 初始化顶点对集合，如果不是撕裂点与非撕裂点的连接就计算色差作为代价
// 此处会将所有的顶点对加入，为的是之后建立链表的时候可以按照编号找到对应边，之后再删除撕裂点与非撕裂点连接的边
//************************************
void OpenmeshHelper::InitPairList()
{
	for (auto j = mesh.halfedges_begin(); j != mesh.halfedges_end(); ++j, ++j)
	{
		auto from = mesh.from_vertex_handle(j);
		auto to = mesh.to_vertex_handle(j);

		if ((mesh.status(from).feature() && mesh.status(to).feature()) ||
			(!(mesh.status(from).feature() || mesh.status(to).feature())))
		{
			Pair temp;
			temp.edge = j;
			int difference = 0;
			for (int k = 0; k < 3; ++k)
			{
				int from_color = mesh.color(from).data()[k];
				int to_color = mesh.color(to).data()[k];
				difference += abs(from_color - to_color);
			}
			temp.value = difference;
			point_pair_list.push_back(temp);
		}
	}
	cout << "Init pair list complete" << endl;
}
void OpenmeshHelper::SortVertices()
{
	int id = 0;
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		//cout << id++ << endl;
		int crease = 0; 
		double x = mesh.point(i).data()[0];
		double y = mesh.point(i).data()[1];
		if ((x == 0 && y == 0) || (x == 0 & y == image->width - 1)
			|| (x == image->height - 1 && y == 0) || (x == image->height - 1 && y == image->width- 1))
		{
			auto t = mesh.texcoord3D(i);
			t[2] = Corner;
			mesh.set_texcoord3D(i, t);
			continue;
		}
		for (auto it = mesh.voh_iter(i); it; ++it)
		{
			auto position = find(crease_list.begin(), crease_list.end(), it);
			int id = it.handle().idx();
			if (position != crease_list.end())
			{
				++crease;
			}
		}
		
		if (crease == 0)
		{
			auto t = mesh.texcoord3D(i);
			t[2] = Smooth;
			mesh.set_texcoord3D(i, t);
			continue;
		}
		if (crease == 2)
		{
			auto t = mesh.texcoord3D(i);
			t[2] = Crease;
			mesh.set_texcoord3D(i, t);
			continue;
		}
		auto t = mesh.texcoord3D(i);
		t[2] = Corner;
		mesh.set_texcoord3D(i, t);
	}
	
	cout << "sort complete" << endl;
}
//************************************
// Method:    Output
// FullName:  OpenmeshHelper::Output
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: string output_location 输出的位置
//************************************
void OpenmeshHelper::Output(MyMesh output_mesh, string output_location)
{
	try
	{
		if (!OpenMesh::IO::write_mesh(output_mesh, output_location))
		{
			std::cerr << "Cannot write mesh to " + output_location << std::endl;
			return;
		}
	}
	catch (std::exception& x)
	{
		std::cerr << x.what() << std::endl;
		return;
	}
}
bool OpenmeshHelper::Compare(Pair a, Pair b)
{
	return a.value > b.value;
}

//************************************
// Method:    RegulatePosition
// FullName:  OpenmeshHelper::RegulatePosition
// Access:    private 
// Returns:   void
// Qualifier:
// Parameter: OpenMesh::Decimater::CollapseInfoT<MyMesh> info collapse之后的信息对象
// Parameter: int * x collapse之后的x位置
// Parameter: int * y collapse之后的y位置
// 计算点合并之后的位置，如果被消除的点是角落点，则角落点不变，如果2个点都是四边上的点，则如果在同一边则取中点
// 如果不是同一边，则位置不变，如果2点中有一个点是四边点，则点的位置计算到边上，如果都不是则取中点
//************************************
void OpenmeshHelper::RegulatePosition(OpenMesh::Decimater::CollapseInfoT<MyMesh> info, double *x, double *y)
{
	double remove_x = info.p0.data()[0];
	double remove_y = info.p0.data()[1];
	double remain_x = info.p1.data()[0];
	double remain_y = info.p1.data()[1];

	if ((remove_x == 0 && remove_y == 0) ||
		(remove_x == 0 && remove_y == image->width - 1) ||
		(remove_x == image->height - 1 && remove_y == 0) ||
		(remove_x == image->height - 1 && remove_y == image->width - 1))
	{
		*x = remove_x;
		*y = remove_y;
		return;
	}
	if ((remain_x == 0 && remain_y == 0) ||
		(remain_x == 0 && remain_y == image->width - 1) ||
		(remain_x == image->height - 1 && remain_y == 0) ||
		(remain_x == image->height - 1 && remain_y == image->width - 1))
	{
		*x = remain_x;
		*y = remain_y;
		return;
	}
	if (remove_x == 0 || remove_y == 0 || remove_x == image->height - 1 || remove_y == image->width - 1)
	{
		if (remain_x == 0 || remain_y == 0 || remain_x == image->height - 1 || remain_y == image->width - 1)
		{
			if ((remain_x == remove_x) || (remain_y == remove_y))
			{
				*x = ((info.p1.data()[0] + info.p0.data()[0]) / 2);
				*y = ((info.p1.data()[1] + info.p0.data()[1]) / 2);
				return;
			}
			else
			{
				*x = remain_x;
				*y = remain_y;
				return;
			}
		}
		else
		{
			*x = remove_x;
			*y = remove_y;
			return;
		}

	}
	else
	{
		if (!(remain_x == 0 || remain_y == 0 || remain_x == image->height - 1 || remain_y == image->width - 1))
		{
			*x = ((info.p1.data()[0] + info.p0.data()[0]) / 2);
			*y = ((info.p1.data()[1] + info.p0.data()[1]) / 2);
			return;
		}
		else
		{
			*x = remain_x;
			*y = remain_y;
			return;
		}
	}
}
//************************************
// Method:    CollapseEdge
// FullName:  OpenmeshHelper::CollapseEdge
// Access:    private 
// Returns:   void
// Qualifier:
// Parameter: MyMesh::HalfedgeHandle half
// 合并边，并在合并之后调整颜色值，调整位置，并删除合并之前的相关边，重新计算代价并加入相关边
//************************************
void OpenmeshHelper::CollapseEdge(MyMesh::HalfedgeHandle half)
{
	auto info = OpenMesh::Decimater::CollapseInfoT<MyMesh>::CollapseInfoT(mesh, half);
	auto remove = info.v0;
	auto remain = info.v1;
	for (auto i = mesh.voh_begin(remove); i != mesh.voh_end(remove); ++i)
	{
		Pair delete_pair;
		delete_pair.edge = i.handle();
		auto delete_it = find(point_pair_list.begin(), point_pair_list.end(), delete_pair);
		if (delete_it != point_pair_list.end())
			point_pair_list.erase(delete_it);
		delete_pair.edge = mesh.opposite_halfedge_handle(i);
		delete_it = find(point_pair_list.begin(), point_pair_list.end(), delete_pair);
		if (delete_it != point_pair_list.end())
			point_pair_list.erase(delete_it);
	}
	for (auto i = mesh.voh_begin(remain); i != mesh.voh_end(remain); ++i)
	{
		Pair delete_pair;
		delete_pair.edge = i.handle();
		auto delete_it = find(point_pair_list.begin(), point_pair_list.end(), delete_pair);
		if (delete_it != point_pair_list.end())
			point_pair_list.erase(delete_it);
		delete_pair.edge = mesh.opposite_halfedge_handle(i);
		delete_it = find(point_pair_list.begin(), point_pair_list.end(), delete_pair);
		if (delete_it != point_pair_list.end())
			point_pair_list.erase(delete_it);
	}
	mesh.collapse(half);
	//调整合并后点的位置
	double x, y;
	RegulatePosition(info, &x, &y);
	mesh.point(info.v1).data()[0] = x;
	mesh.point(info.v1).data()[1] = y;
	//调整点的颜色
	
	MyMesh::Color temp_color;
	for (int k = 0; k < 3; ++k)
	{
		int remove_color = mesh.color(remove)[k];
		int remain_color = mesh.color(remain)[k];
		temp_color[k] = (int)((remove_color + remain_color) / 2);
	}
	mesh.set_color(remain, temp_color);

	//重新计算相连的边的代价，并加入vector
	for (auto j = mesh.voh_begin(remain); j != mesh.voh_end(remain); ++j)
	{
		auto to = mesh.to_vertex_handle(j);

		if ((mesh.status(remain).feature() && mesh.status(to).feature()) ||
			(!(mesh.status(remain).feature() || mesh.status(to).feature())))
		{
			Pair temp;
			temp.edge = j;
			int difference = 0;
			for (int k = 0; k < 3; ++k)
			{
				int from_color = mesh.color(remain).data()[k];
				int to_color = mesh.color(to).data()[k];
				difference += abs(from_color - to_color);
			}
			temp.value = difference;
			point_pair_list.push_back(temp);
		}
	}
	make_heap(point_pair_list.begin(), point_pair_list.end(), Compare);
}

void OpenmeshHelper::InsertCreaseEdge()
{
	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
	{
		auto from = mesh.from_vertex_handle(i);
		auto to = mesh.to_vertex_handle(i);
		
		int from_x = mesh.point(from).data()[0];
		int from_y = mesh.point(from).data()[1];
		int to_x = mesh.point(to).data()[0];
		int to_y = mesh.point(to).data()[1];

		if (mesh.status(from).feature() && mesh.status(to).feature())
		{
			crease_list.push_back(i.handle());
		}
		if ((from_x == 0 && to_x == 0) || (from_x == image->height - 1 && to_x == image->height - 1) ||
			(from_y == 0 && to_y == 0) || (from_y == image->width - 1 && to_y == image->width - 1))
		{
			crease_list.push_back(i.handle());
		}
	}
}


//************************************
// Method:    ReduceVertices
// FullName:  OpenmeshHelper::ReduceVertices
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: double rate 简化率，表示简化后剩余边的比例
// Parameter: bool visual 标记，用于是否显示简化进度候
//************************************
void OpenmeshHelper::ReduceVertices(double rate, bool visual)
{
	mesh.request_vertex_colors();
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();
	mesh.request_halfedge_status();
	mesh.request_vertex_texcoords3D();
	InitPointData();
	ConnectMesh(true);
	InsertCreaseEdge();

	CountVertices();
	SortVertices();
	InitPairList();

	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		int x = mesh.point(i.handle()).data()[0];
		int y = mesh.point(i.handle()).data()[1];
		control_list.push_back(point_data[x][y]);
	}

	LoopReduce(rate, visual);
	vector<MyMesh::HalfedgeHandle> temp_crease_list;
	for (int i = 0; i != crease_list.size(); ++i)
	{
		auto temp = crease_list[i];
		if (mesh.from_vertex_handle(temp) != mesh.to_vertex_handle(temp))
		{
			temp_crease_list.push_back(crease_list[i]);
		}
	}
	crease_list = temp_crease_list;
	for (int i = 0; i != crease_list.size(); ++i)
	{
		Position from, to;
		from.x = mesh.point(mesh.from_vertex_handle(crease_list[i])).data()[0];
		from.y = mesh.point(mesh.from_vertex_handle(crease_list[i])).data()[1];
		to.x = mesh.point(mesh.to_vertex_handle(crease_list[i])).data()[0];
		to.y = mesh.point(mesh.to_vertex_handle(crease_list[i])).data()[1];
		crease_map.insert(make_pair(from, to));
	}

	ofstream crease_out("D:/crease.txt");
	for (auto i = crease_map.begin(); i != crease_map.end(); ++i)
	{
		crease_out << i->first.x << ' ' << i->first.y << ' ' << i->second.x << " " << i->second.y << endl;
	}
	crease_out.close();
	cout << "Crease Output Complete" << endl;
	mesh.garbage_collection();
	//OpenMesh::IO::read_mesh(mesh, "D:/before.off");
	//ifstream input("D:/crease.txt");
	//Position from, to;
	//while (input >> from.x >> from.y >> to.x >> to.y)
	//{
	//	crease_map.insert(make_pair(from, to));
	//}
	//input.close();
	//RebuildCreaseList(mesh);
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_vertex_normals();
	mesh.update_face_normals();
	mesh.release_face_normals();
	mesh.release_vertex_normals();

	Output(mesh, "D:/before.off");
	OptimizePosition();
	ReSampleColor();
	OptimizeColor();
	Output(mesh, "D:/origin.off");
	CountVertices();
	while (num_all_vertices < image->width * zoom * image->height * zoom)
	{
		SubDivision();
		SortVertices();
		CountVertices();
		cout << num_all_vertices << endl;
	}
	OutputImage();
	mesh.release_vertex_colors();
	mesh.release_vertex_texcoords3D();
	mesh.release_edge_status();
	mesh.release_vertex_status();
	mesh.release_halfedge_status();
	mesh.release_face_status();
	cout << "Reduce Vertrices complete" << endl;
}

//************************************
// Method:    IsCollapseOK
// FullName:  OpenmeshHelper::IsCollapseOK
// Access:    private 
// Returns:   bool // Qualifier:
// Parameter: MyMesh::HalfedgeHandle half
// 用于判断collapse之后的面是否会出现翻转的现象，在collapse之前执行
//************************************
bool OpenmeshHelper::IsCollapseOK(MyMesh::HalfedgeHandle half)
{
	if (!mesh.is_collapse_ok(half)) return false;
	auto info = OpenMesh::Decimater::CollapseInfoT<MyMesh>::CollapseInfoT(mesh, half); 
	double x = 0;
	double y = 0;
	RegulatePosition(info, &x, &y);
	double x0 = mesh.point(info.v0).data()[0];
	double y0 = mesh.point(info.v0).data()[1];
	double x1 = mesh.point(info.v1).data()[0];
	double y1 = mesh.point(info.v1).data()[1];

	
	set<MyMesh::VertexHandle> remain_set;
	for (auto i = mesh.vv_begin(info.v0); i != mesh.vv_end(info.v0); ++i)
	{
		remain_set.insert(i);
	}
	for (auto i = mesh.vv_begin(info.v1); i != mesh.vv_end(info.v1); ++i)
	{
		remain_set.insert(i);
	}
	remain_set.erase(info.v0);
	remain_set.erase(info.v1);
	bool remain_boundary = IsBoundary(info.v1);
	bool remove_boundary = IsBoundary(info.v0);
	if (!remove_boundary && !remain_boundary)
	{
		if (remain_set.size() < 4 || remain_set.size() > 8)
			return false;
	}
	else
	{
		if (remain_set.size() < 3 || remain_set.size() > 5)
			return false;
	}
	mesh.point(info.v0).data()[0] = x;
	mesh.point(info.v0).data()[1] = y;
	for (auto i = mesh.vf_begin(info.v0); i != mesh.vf_end(info.v0); ++i)
	{
		auto normals = mesh.calc_face_normal(i);
		if (normals.data()[2] < 0)
		{
			mesh.point(info.v0).data()[0] = x0;
			mesh.point(info.v0).data()[1] = y0;
			return false;
		}
		for (auto j = mesh.fh_begin(i); j != mesh.fh_end(i); ++j)
		{
			double angle = mesh.calc_sector_angle(j) * 180 / M_PI;
			if (abs(angle - 60) > degree)
			{
				mesh.point(info.v0).data()[0] = x0;
				mesh.point(info.v0).data()[1] = y0;
				return false;
			}
		}
	}
	mesh.point(info.v0).data()[0] = x0;
	mesh.point(info.v0).data()[1] = y0;

	mesh.point(info.v1).data()[0] = x;
	mesh.point(info.v1).data()[1] = y;
	
	for (auto i = mesh.vf_begin(info.v1); i != mesh.vf_end(info.v1); ++i)
	{
		auto normals = mesh.calc_face_normal(i.handle());
		if (normals.data()[2] < 0)
		{
			mesh.point(info.v1).data()[0] = x1;
			mesh.point(info.v1).data()[1] = y1;
			return false;
		}
		for (auto j = mesh.fh_begin(i); j != mesh.fh_end(i); ++j)
		{
			double angle = mesh.calc_sector_angle(j) * 180 / M_PI;
			if (abs(angle - 60) > degree)
			{
				mesh.point(info.v1).data()[0] = x1;
				mesh.point(info.v1).data()[1] = y1;
				return false;
			}
		}
	}

	mesh.point(info.v1).data()[0] = x1;
	mesh.point(info.v1).data()[1] = y1;



	return true;
}

//************************************
// Method:    LoopReduce
// FullName:  OpenmeshHelper::LoopReduce
// Access:    private 
// Returns:   void
// Qualifier:
// Parameter: double rate 简化率，表示简化后剩余边的比例
// Parameter: bool visual 标记，用于是否显示简化进度
// 简化边
//************************************
void OpenmeshHelper::LoopReduce(double rate, bool visual)
{
	make_heap(point_pair_list.begin(), point_pair_list.end(), Compare);
	int min = num_all_vertices * rate;
	int times = 0;
	int last_num = num_all_vertices;
	while (num_vertices > min)
	{
		while (point_pair_list.size() > 0)
		{
			pop_heap(point_pair_list.begin(), point_pair_list.end(), Compare);
			Pair pair = point_pair_list[point_pair_list.size() - 1];
			auto half = pair.edge;
			int id = half.idx();
			point_pair_list.pop_back();
			if (IsCollapseOK(half))
			{
				CollapseEdge(half);
				--num_vertices;
				if (num_vertices <= min)
					break;
			}
			else
			{
				half = mesh.opposite_halfedge_handle(half);
				if (IsCollapseOK(half))
				{
					cout << "opposite" << endl;
					--num_vertices;
					CollapseEdge(half);
					if (num_vertices <= min)
						break;
				}
				else
				{
					second_pair_list.push_back(pair);
				}
			}
			if (visual)
			{
			
				cout << num_vertices << endl;
			}
				//cout << num_vertices / num_all_vertices << endl;
		}
		point_pair_list = second_pair_list;
		second_pair_list.clear();
		make_heap(point_pair_list.begin(), point_pair_list.end(), Compare);
		cout << times++ << ":" << num_vertices / num_all_vertices << endl;
		if (num_vertices == last_num)
		{
			cout << num_vertices << "all" << num_all_vertices << endl;
			break;
		}
		last_num = num_vertices;
	}
	
	cout << "Loop_Reduce_Complete" << endl;
}

void OpenmeshHelper::OptimizePosition()
{
	SortVertices();
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		int id = i.handle().idx();
		double x = 0;
		double y = 0;
		int num = 0;
		auto t = mesh.texcoord3D(i);
		if (t[2] == Corner)
		{
			continue;
		}
		if (t[2] == Crease)
		{
			for (auto j = mesh.vv_begin(i); j != mesh.vv_end(i); ++j)
			{
				if (mesh.texcoord3D(j)[2] == Crease)
				{
					++num;
					x += abs(mesh.point(i).data()[0] - mesh.point(j).data()[0]);
					y += abs(mesh.point(i).data()[1] - mesh.point(j).data()[1]);
				}
			}
		}
		if (t[2] == Smooth)
		{
			for (auto j = mesh.vv_begin(i); j != mesh.vv_end(i); ++j)
			{
				++num;
				x += abs(mesh.point(i).data()[0] - mesh.point(j).data()[0]);
				y += abs(mesh.point(i).data()[1] - mesh.point(j).data()[1]);
			}
		}
		if (mesh.point(i).data()[0] == 0 || mesh.point(i).data()[0] == image->height - 1
			|| mesh.point(i).data()[1] == 0 || mesh.point(i).data()[1] == image->width - 1)
		//if (t[2] == Edge)
		{
			x = 0;
			y = 0;
			num = 0;
			for (auto j = mesh.vv_begin(i); j != mesh.vv_end(i); ++j)
			{
				double x1 = mesh.point(i).data()[0];
				double y1 = mesh.point(i).data()[1];
				double xx = mesh.point(j).data()[0];
				double yy = mesh.point(j).data()[1];
				if (mesh.point(i).data()[0] == 0)
				{
					if (mesh.point(j).data()[0] == 0)
					{
						x += abs(mesh.point(i).data()[0] - mesh.point(j).data()[0]);
						y += abs(mesh.point(i).data()[1] - mesh.point(j).data()[1]);
						continue;
					}
				}
				if (mesh.point(i).data()[1] == 0)
				{
					if (mesh.point(j).data()[1] == 0)
					{
						x += abs(mesh.point(i).data()[0] - mesh.point(j).data()[0]);
						y += abs(mesh.point(i).data()[1] - mesh.point(j).data()[1]);
						continue;
					}
				}
				if (mesh.point(i).data()[0] == image->height - 1)
				{
					if (mesh.point(j).data()[0] == image->height - 1)
					{
						x += abs(mesh.point(i).data()[0] - mesh.point(j).data()[0]);
						y += abs(mesh.point(i).data()[1] - mesh.point(j).data()[1]);
						continue;
					}
				}
				if (mesh.point(i).data()[1] == image->width - 1)
				{
					if (mesh.point(j).data()[1] == image->width - 1)
					{
						x += abs(mesh.point(i).data()[0] - mesh.point(j).data()[0]);
						y += abs(mesh.point(i).data()[1] - mesh.point(j).data()[1]);
						continue;
					}
				}

			}
		}
		
		double temp_x = mesh.point(i).data()[0];
		double temp_y = mesh.point(i).data()[1];
		if (x != 0 && y != 0)
		{
			mesh.point(i).data()[0] += lambda * x / num;
			mesh.point(i).data()[1] += lambda * y / num;
		}
		for (auto j = mesh.vf_begin(i); j != mesh.vf_end(i); ++j)
		{
			auto normals = mesh.calc_face_normal(j);
			if (normals.data()[2] < 0)
			{
				mesh.point(i).data()[0] = temp_x;
				mesh.point(i).data()[1] = temp_y;
				break;
			}
		}
	}
}
void OpenmeshHelper::OptimizeColor()
{
	control_list.clear();
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		control_list.push_back(i);
	}

	Eigen::MatrixXf X, Y;
	if (exact)
	{
		X = Eigen::MatrixXf::Zero(1, control_list.size() * 3);
	}
	else
	{
		X = Eigen::MatrixXf::Zero(3, control_list.size());
	}
	Y = Eigen::MatrixXf::Zero(control_list.size(), control_list.size());
	for (int i = 0; i != control_list.size(); ++ i)
	{
		//TODO 未考虑Tear点
		int x = mesh.point(control_list[i]).data()[0];
		int y = mesh.point(control_list[i]).data()[1];
		auto t = mesh.texcoord3D(control_list[i]);
		auto color = mesh.color(control_list[i]);
		X(0, i) = color[0];
		if (exact)
		{
			X(0, i + control_list.size()) = color[1];
			X(0, i + 2 * control_list.size()) = color[2];
		}
		else
		{
			X(1, i) = color[1];
			X(2, i) = color[2];
		}
		
		//处理Y
		if (t[2] == Corner)
		{
			Y(i, i) = 1;
			continue;
		}
		if (t[2] == Crease)
		{
			Y(i, i) = 4.0 / 6;
			for (auto j = mesh.voh_begin(control_list[i]); j != mesh.voh_end(control_list[i]); ++j)
			{
				auto position = find(crease_list.begin(), crease_list.end(), j.handle());
				if (position != crease_list.end())
				{
					Y(i, mesh.to_vertex_handle(j).idx()) = 1.0 / 6;
				}
			}
			continue;
		}
		auto point = control_list[i];
		int valence = 0;
		vector<MyMesh::VertexHandle> out_point_list;
		out_point_list.push_back(point);
		for (auto j = mesh.vv_begin(point); j != mesh.vv_end(point); ++j)
		{
			++valence;
			out_point_list.push_back(j);
		}
		bool is_boundary = false;
		int boundary[2];
		for (int j = 1; j != out_point_list.size(); ++j)
		{
			int times = 0;
			for (auto k = mesh.vv_begin(out_point_list[j]); k != mesh.vv_end(out_point_list[j]); ++k)
			{
				auto find_position = find(out_point_list.begin(), out_point_list.end(), k);
				if (find_position != out_point_list.end())
				{
					++times;
				}
			}
			if (times != 3)
			{
				if (!is_boundary)
				{
					is_boundary = true;
					boundary[0] = out_point_list[j].idx();
				}
				else
				{
					boundary[1] = out_point_list[j].idx();
					break;
				}
				
			}

		}
		if (is_boundary)
		{
			Y(i, boundary[0]) = 1.0 / 6;
			Y(i, boundary[1]) = 1.0 / 6;
			Y(i, i) = 4.0 / 6;
		}
		else
		{
			double alpha = 0.625 - pow(0.375 + 0.25 * cos(2 * M_PI / valence), 2);
			double omega = 0.375 * valence / alpha;
			double sum = omega + valence;
			Y(i, out_point_list[0].idx()) = omega / sum;
			for (int j = 1; j != valence + 1; ++j)
			{
				Y(i, out_point_list[j].idx()) = 1.0 / sum;
			}
		}
		
	}
	if (exact)
	{
		Eigen::MatrixXf YY;
		YY = Eigen::MatrixXf::Zero(3 * control_list.size(), 3 * control_list.size());
		for (int i = 0; i != control_list.size(); ++i)
		{
			for (int j = 0; j != control_list.size(); ++j)
			{
				if (Y(i, j) < 0)
				{
					cout << "--------" << endl;
					exit(0);
				}
				YY(i, j) = Y(i, j);
				YY(i + control_list.size(), j + control_list.size()) = Y(i, j);
				YY(i + 2 * control_list.size(), j + 2 * control_list.size()) = Y(i, j);
			}
		}
		Eigen::MatrixXf optimize = YY.colPivHouseholderQr().solve(X.transpose());
		for (int i = 0; i != control_list.size(); ++i)
		{
			auto color = mesh.color(control_list[i]);
			color[0] = optimize(i, 0);
			color[1] = optimize(i + control_list.size(), 0);
			color[2] = optimize(i + 2 * control_list.size(), 0);
			mesh.set_color(control_list[i], color);
		}
	}
	else
	{
		Eigen::MatrixXf optimize = Y.colPivHouseholderQr().solve(X.transpose());
		for (int i = 0; i != control_list.size(); ++i)
		{
			auto color = mesh.color(control_list[i]);
			color[0] = optimize(i, 0);
			color[1] = optimize(i, 1);
			color[2] = optimize(i, 2);
			mesh.set_color(control_list[i], color);
		}
	}
	

}

bool OpenmeshHelper::IsBoundary(MyMesh::VertexHandle point)
{
	vector<MyMesh::VertexHandle> out_point_list;
	for (auto i = mesh.vv_begin(point); i != mesh.vv_end(point); ++i)
	{
		out_point_list.push_back(i);
	}
	
	for (auto i = mesh.vv_begin(point); i != mesh.vv_end(point); ++i)
	{
		int times = 0;
		for (auto j = mesh.vv_begin(i); j != mesh.vv_end(i); ++j)
		{
			auto find_position = find(out_point_list.begin(), out_point_list.end(), j);
			if (find_position != out_point_list.end())
			{
				++times;
			}
	
		}
		if (times != 2)
		{
			return true;
		}
	}
	return false;
}

void OpenmeshHelper::ReSampleColor()
{
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		double real_x = mesh.point(i).data()[0];
		double real_y = mesh.point(i).data()[1];
		double x = real_x - (int)real_x;
		double y = real_y - (int)real_y;
		int xx = real_x;
		int yy = real_y;
		if (xx == image->height - 1)
		{
			if (yy == image->width - 1)
			{
				auto color_image = cvGet2D(image, xx, yy);
				MyMesh::Color color;
				for (int i = 0; i != 3; ++i)
				{
					color[i] = (int)color_image.val[i];
				}
				mesh.set_color(i, color);
				continue;
			}
			else
			{
				auto color1 = cvGet2D(image, xx, yy);
				auto color2 = cvGet2D(image, xx, yy + 1);
				MyMesh::Color color;
				for (int i = 0; i != 3; ++i)
				{
					color[i] = (int)((color1.val[i] + color2.val[i]) / 2);
				}
				mesh.set_color(i, color);
				continue;
			}
		}
		if (yy == image->width - 1)
		{
			auto color1 = cvGet2D(image, xx, yy);
			auto color2 = cvGet2D(image, xx + 1, yy);
			MyMesh::Color color;
			for (int i = 0; i != 3; ++i)
			{
				color[i] = (int)((color1.val[i] + color2.val[i]) / 2);
			}
			mesh.set_color(i, color);
			continue;
		}
		auto color00 = cvGet2D(image, xx, yy);
		auto color01 = cvGet2D(image, xx, yy + 1);
		auto color10 = cvGet2D(image, xx + 1, yy);
		auto color11 = cvGet2D(image, xx + 1, yy + 1);

		MyMesh::Color color;
		for (int j = 0; j != 3; ++j)
		{
			color[j] = (int)(color00.val[j] * (1 - x) * (1 - y) 
				+ color10.val[j] * x * (1 - y)
				+ color01.val[j] * y * (1 - x)
				+ color11.val[j] * x * y);
		}
		mesh.set_color(i, color);
	}
	ofstream out("D:/color.txt");
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		out << mesh.color(i).data()[0] << ' ' << mesh.color(i).data()[1] << ' ' << mesh.color(i).data()[2] << endl;
	}
	out.close();
	cout << "Color Complete" << endl;
}

void OpenmeshHelper::SubDivision()
{
	MyMesh subdivision_mesh;
	subdivision_mesh.request_vertex_colors();
	subdivision_mesh.request_vertex_status();
	subdivision_mesh.request_edge_status();
	subdivision_mesh.request_face_status();
	subdivision_mesh.request_halfedge_status();
	subdivision_mesh.request_vertex_texcoords3D();
	crease_map.clear();
	vertex_map.clear();
	int id = 0;
	for (auto i = mesh.faces_begin(); i != mesh.faces_end(); ++i)
	{
		//cout << id++ << endl;
		MyMesh::VertexHandle face_vertex[6];
		MyMesh::VertexHandle old_vertex[3];
		int j = 0;
		for (auto k = mesh.fv_begin(i); k != mesh.fv_end(i); ++k)
		{
			Position p;
			p.x = mesh.point(k).data()[0];
			p.y = mesh.point(k).data()[1];
			if (vertex_map[p].idx() == -1)
			{
				vertex_map[p] = subdivision_mesh.add_vertex(MyMesh::Point(p.x, p.y, 0));
				subdivision_mesh.set_color(vertex_map[p], mesh.color(k));
			}
			face_vertex[j++] = vertex_map[p];
			old_vertex[j - 1] = k.handle();
		}
		for (auto k = mesh.fh_begin(i); k != mesh.fh_end(i); ++k)
		{
			if ((mesh.from_vertex_handle(k) == old_vertex[0] &&
				mesh.to_vertex_handle(k) == old_vertex[1]) ||
				(mesh.to_vertex_handle(k) == old_vertex[0] &&
				mesh.from_vertex_handle(k) == old_vertex[1]))
			{
				face_vertex[3] = InsertVertex(subdivision_mesh, k);
			}
			if ((mesh.from_vertex_handle(k) == old_vertex[1] &&
				mesh.to_vertex_handle(k) == old_vertex[2]) ||
				(mesh.to_vertex_handle(k) == old_vertex[1] &&
				mesh.from_vertex_handle(k) == old_vertex[2]))
			{
				face_vertex[4] = InsertVertex(subdivision_mesh, k);
			}
			if ((mesh.from_vertex_handle(k) == old_vertex[2] &&
				mesh.to_vertex_handle(k) == old_vertex[0]) ||
				(mesh.to_vertex_handle(k) == old_vertex[2] &&
				mesh.from_vertex_handle(k) == old_vertex[0]))
			{
				face_vertex[5] = InsertVertex(subdivision_mesh, k);
			}
		}
		vector<MyMesh::VertexHandle> list;
		list.push_back(face_vertex[0]);
		list.push_back(face_vertex[3]);
		list.push_back(face_vertex[5]);
		subdivision_mesh.add_face(list);

		list.clear();
		list.push_back(face_vertex[3]);
		list.push_back(face_vertex[1]);
		list.push_back(face_vertex[4]);
		subdivision_mesh.add_face(list);

		list.clear();
		list.push_back(face_vertex[5]);
		list.push_back(face_vertex[4]);
		list.push_back(face_vertex[2]);
		subdivision_mesh.add_face(list);
		
		list.clear();
		list.push_back(face_vertex[3]);
		list.push_back(face_vertex[4]);
		list.push_back(face_vertex[5]);
		subdivision_mesh.add_face(list);
	}
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		int x = mesh.point(i).data()[0];
		int y = mesh.point(i).data()[1];
		auto t = mesh.texcoord3D(i);
		if (x == 0 || y == 0 || x == image->height - 1 || y == image->width - 1)
		{
			continue;
		}
		if (t[2] == Corner)
		{
			continue;
		}
		if (t[2] == Crease)
		{
			double temp_x = 0.75 * mesh.point(i).data()[0];
			double temp_y = 0.75 * mesh.point(i).data()[1];
			for (auto j = mesh.voh_begin(i); j != mesh.voh_end(i); ++j)
			{
				auto position = find(crease_list.begin(), crease_list.end(), j.handle());
				if (position != crease_list.end())
				{
					temp_x += 0.125 * mesh.point(mesh.to_vertex_handle(j)).data()[0];
					temp_y += 0.125 * mesh.point(mesh.to_vertex_handle(j)).data()[1];
				}
			}
			Position position;
			position.x = mesh.point(i).data()[0];
			position.y = mesh.point(i).data()[1];
			subdivision_mesh.point(vertex_map[position]).data()[0] = temp_x;
			subdivision_mesh.point(vertex_map[position]).data()[1] = temp_y;
			continue;
		}
		auto point = i;
		int valence = 0;
		vector<MyMesh::VertexHandle> out_point_list;
		out_point_list.push_back(point);
		for (auto j = mesh.vv_begin(point); j != mesh.vv_end(point); ++j)
		{
			++valence;
			out_point_list.push_back(j);
		}
		bool is_boundary = false;
		MyMesh::VertexHandle boundary[2];
		for (int j = 1; j != out_point_list.size(); ++j)
		{
			int times = 0;
			for (auto k = mesh.vv_begin(out_point_list[j]); k != mesh.vv_end(out_point_list[j]); ++k)
			{
				auto find_position = find(out_point_list.begin(), out_point_list.end(), k);
				if (find_position != out_point_list.end())
				{
					++times;
				}
			}
			if (times != 3)
			{
				if (!is_boundary)
				{
					is_boundary = true;
					boundary[0] = out_point_list[j];
				}
				else
				{
					boundary[1] = out_point_list[j];
					break;
				}

			}

		}
		if (is_boundary)
		{
			double temp_x = 0.75 * mesh.point(i).data()[0];
			double temp_y = 0.75 * mesh.point(i).data()[1];
			temp_x += 0.125 * mesh.point(boundary[0]).data()[0];
			temp_y += 0.125 * mesh.point(boundary[0]).data()[1];
			temp_x += 0.125 * mesh.point(boundary[1]).data()[0];
			temp_y += 0.125 * mesh.point(boundary[1]).data()[1];
			Position p;
			p.x = mesh.point(i).data()[0];
			p.y = mesh.point(i).data()[1];
			subdivision_mesh.point(vertex_map[p]).data()[0] = temp_x;
			subdivision_mesh.point(vertex_map[p]).data()[1] = temp_y;
		}
		else
		{
			double alpha = 0.625 - pow(0.375 + 0.25 * cos(2 * M_PI / valence), 2);
			double temp_x = (1 - alpha) * mesh.point(i).data()[0];
			double temp_y = (1 - alpha) * mesh.point(i).data()[1];
			
			for (int j = 1; j != valence + 1; ++j)
			{
				temp_x += alpha / valence * mesh.point(out_point_list[j]).data()[0];
				temp_y += alpha / valence * mesh.point(out_point_list[j]).data()[1];
			}
			Position p;
			p.x = mesh.point(i).data()[0];
			p.y = mesh.point(i).data()[1];
			subdivision_mesh.point(vertex_map[p]).data()[0] = temp_x;
			subdivision_mesh.point(vertex_map[p]).data()[1] = temp_y;
		}
	}
	RebuildCreaseList(subdivision_mesh);
	mesh = subdivision_mesh;
	Output(subdivision_mesh, "D:/subdivision.off");
	cout << "Subdivision Complete" << endl;
}

MyMesh::VertexHandle OpenmeshHelper::InsertVertex(MyMesh &subdivision_mesh, MyMesh::HalfedgeHandle half)
{
	auto from = mesh.from_vertex_handle(half);
	auto to = mesh.to_vertex_handle(half);
	auto info = OpenMesh::Decimater::CollapseInfoT<MyMesh>::CollapseInfoT(mesh, half);
	if (info.fl.idx() == -1 || info.fr.idx() == -1)
	{
		Position p;
		p.x = (mesh.point(from).data()[0] + mesh.point(to).data()[0]) / 2;
		p.y = (mesh.point(from).data()[1] + mesh.point(to).data()[1]) / 2;
		if (vertex_map[p].idx() == -1)
		{
			vertex_map[p] = subdivision_mesh.add_vertex(MyMesh::Point(p.x, p.y, 0)); 
			MyMesh::Color color;
			for (int k = 0; k != 3; ++k)
			{
				color[k] = (mesh.color(from).data()[k] + mesh.color(to).data()[k]) / 2;
			}
			subdivision_mesh.set_color(vertex_map[p], color);
		}
		return vertex_map[p];
	}
	auto position = find(crease_list.begin(), crease_list.end(), half);
	if (position != crease_list.end())
	{
		Position p;
		p.x = (mesh.point(from).data()[0] + mesh.point(to).data()[0]) / 2;
		p.y = (mesh.point(from).data()[1] + mesh.point(to).data()[1]) / 2;
		if (vertex_map[p].idx() == -1)
		{
			vertex_map[p] = subdivision_mesh.add_vertex(MyMesh::Point(p.x, p.y, 0)); 
			MyMesh::Color color;
			for (int k = 0; k != 3; ++k)
			{
				color[k] = (mesh.color(from).data()[k] + mesh.color(to).data()[k]) / 2;
			}
			subdivision_mesh.set_color(vertex_map[p], color);
		}
		Position from_position, to_position;
		from_position.x = mesh.point(from).data()[0];
		from_position.y = mesh.point(from).data()[1];
		to_position.x = mesh.point(to).data()[0];
		to_position.y = mesh.point(to).data()[1];
		crease_map.insert(make_pair(from_position, p));
		crease_map.insert(make_pair(p, to_position));
		return vertex_map[p];
	}

	MyMesh::VertexHandle left, right;
	for (auto i = mesh.fv_begin(info.fl); i != mesh.fv_end(info.fl); ++i)
	{
		if (i.handle() != from && i.handle() != to)
		{
			left = i.handle();
			break;
		}
	}
	for (auto i = mesh.fv_begin(info.fr); i != mesh.fv_end(info.fr); ++i)
	{
		if (i.handle() != from && i.handle() != to)
		{
			right = i.handle();
			break;
		}
	}
	Position p;
	p.x = mesh.point(from).data()[0] * 0.375 + mesh.point(to).data()[0] * 0.375 +
		mesh.point(right).data()[0] * 0.125 + mesh.point(left).data()[0] * 0.125;
	p.y = mesh.point(from).data()[1] * 0.375 + mesh.point(to).data()[1] * 0.375 +
		mesh.point(right).data()[1] * 0.125 + mesh.point(left).data()[1] * 0.125;
	if (vertex_map[p].idx() == -1)
	{
		vertex_map[p] = subdivision_mesh.add_vertex(MyMesh::Point(p.x, p.y, 0));
		MyMesh::Color color;
		for (int k = 0; k != 3; ++k)
		{
			color[k] = 0.375 * mesh.color(from).data()[k] + 0.375 * mesh.color(to).data()[k] +
				0.125 * mesh.color(left).data()[k] + 0.125 * mesh.color(right).data()[k];
		}
		subdivision_mesh.set_color(vertex_map[p], color);
	
	}
	return vertex_map[p];
}

void OpenmeshHelper::RebuildCreaseList(MyMesh mesh)
{
	crease_list.clear();
	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
	{
		Position from, to;
		from.x = mesh.point(mesh.from_vertex_handle(i.handle())).data()[0];
		from.y = mesh.point(mesh.from_vertex_handle(i.handle())).data()[1];
		to.x = mesh.point(mesh.to_vertex_handle(i.handle())).data()[0];
		to.y = mesh.point(mesh.to_vertex_handle(i.handle())).data()[1];

		for (auto k = crease_map.lower_bound(from); k != crease_map.upper_bound(from); ++k)
		{
			if (k->second == to)
			{
				auto position = find(crease_list.begin(), crease_list.end(), i.handle());
				if (position == crease_list.end())
					crease_list.push_back(i.handle());
				break;
			}
		}
	}
}

void OpenmeshHelper::OutputImage()
{
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		mesh.point(i).data()[0] = mesh.point(i).data()[0] * zoom;
		mesh.point(i).data()[1] = mesh.point(i).data()[1] * zoom;
	}
	IplImage *result = cvCreateImage(cvSize(image->height * zoom, image->width * zoom), image->depth, 3);
	uchar *data = (uchar *)result->imageData;
	int step = result->widthStep / sizeof(uchar);
	for (auto i = mesh.faces_begin(); i != mesh.faces_end(); ++i)
	{
		Eigen::Vector3f point[3];
		int k = 0;
		double x_min = image->height;
		double x_max = 0;
		double y_min = image->width;
		double y_max = 0;
		for (auto j = mesh.fv_begin(i); j != mesh.fv_end(i); ++j)
		{
			double x = mesh.point(j).data()[0];
			double y = mesh.point(j).data()[1];
			point[k++] << x, y, 0;
			if (x < x_min)
			{
				x_min = x;
			}
			if (x > x_max)
			{
				x_max = x;
			}
			if (y < y_min)
			{
				y_min = y;
			}
			if (y > y_max)
			{
				y_max = y;
			}
		}
		int x_end = (int)x_max;
		int y_end = (int)y_max;
		int x_begin, y_begin;
		double d = x_min - (int)x_min;
		if (d == 0)
			x_begin = (int)x_min;
		else
			x_begin = (int)x_min + 1;
		d = y_min - (int)y_min;
		if (d == 0)
			y_begin = (int)y_min;
		else
			y_begin = (int)y_min + 1;
		for (int x = x_begin; x <= x_end; ++x)
		{
			for (int y = y_begin; y <= y_end; ++y)
			{
				Eigen::Vector3f P(x, y, 0);
				auto AP = point[0] - P;
				auto BP = point[1] - P;
				auto CP = point[2] - P;
				double t1, t2, t3;
				t1 = AP.cross(BP)[2];
				t2 = BP.cross(CP)[2];
				t3 = CP.cross(AP)[2];
				if (t1 * t2 >= 0 && t1 * t3 >= 0)
				{
					int r = 0;
					int g = 0;
					int b = 0;
					for (auto tt = mesh.fv_begin(i); tt != mesh.fv_end(i); ++tt)
					{
						auto color = mesh.color(tt);
						b += (int)color.data()[0];
						g += (int)color.data()[1];
						r += (int)color.data()[2];
					}
					data[x * step + y * 3] = b / 3;
					data[x * step + y * 3 + 1] = g / 3;
					data[x * step + y * 3 + 2] = r / 3;
				}
			}
		}

	}
	
	cvSaveImage("D:/result.jpg", result);
}

