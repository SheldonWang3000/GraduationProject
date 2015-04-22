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
	point_data = new MyPoint*[image->height];
	for (int i = 0; i != image->height; ++i)
	{
		point_data[i] = new MyPoint[image->width];
	}

	OpencvHelper helper;

	//IplImage *contours = cvCreateImage(cvGetSize(image), image->depth, 1);
	//helper.get_pic_contours(image, contours);
	//uchar *data = (uchar *)contours->imageData;
	//int step = contours->widthStep / sizeof(uchar);
	//for (int i = 0; i < contours->height; ++i)
	//{
	//	for (int j = 0; j < contours->width; ++j)
	//	{
	//		point_data[i][j].point_type = None;
	//		//point_data[i][j].point = mesh.add_vertex(MyMesh::Point(j, image->height - i, 0));
	//		point_data[i][j].point = mesh.add_vertex(MyMesh::Point(i, j, 0));
	//		if (data[i * step + j] == 0)
	//		{
	//			point_data[i][j].isEdge = true;
	//		}
	//		else
	//		{
	//			point_data[i][j].isEdge = false;
	//		}
	//	}
	//}
	//ConnectMesh(false);

	//for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
	//{
	//	int from_x = mesh.point(mesh.from_vertex_handle(i.handle())).data()[0];
	//	int from_y = mesh.point(mesh.from_vertex_handle(i.handle())).data()[1];
	//	int to_x = mesh.point(mesh.to_vertex_handle(i.handle())).data()[0];
	//	int to_y = mesh.point(mesh.to_vertex_handle(i.handle())).data()[1];
	//	
	//	if (point_data[from_x][from_y].isEdge && point_data[to_x][to_y].isEdge)
	//	{
	//		//tear_list.push_back(i.handle());
	//		Position from, to;
	//		from.x = from_x;
	//		from.y = from_y;
	//		to.x = to_x;
	//		to.y = to_y;
	//		tear_map.insert(make_pair(from, to));
	//	}
	//}
	//mesh.clear();
	//mesh.garbage_collection();
	//for (int i = 0; i < image->height; ++i)
	//{
	//	for (int j = 0; j < image->width; ++j)
	//	{
	//		//point_data[i][j].point = mesh.add_vertex(MyMesh::Point(j, image->width - i - 1, 0));
	//		point_data[i][j].point = mesh.add_vertex(MyMesh::Point(i, j, 0));
	//	}
	//}

	//cvReleaseImage(&contours);
	//IplImage *area = cvCreateImage(cvGetSize(image), image->depth, 3);
	//helper.get_pic_area(image, area);
	//data = (uchar *)area->imageData;
	//step = area->widthStep / sizeof(uchar);
	//map<int, int> color_map;
	//int color_idx = 0;
	//for (int i = 0; i < area->height; ++ i)
	//{
	//	for (int j = 0; j < area->width; ++ j)
	//	{
	//		int color = 0;
	//		CvScalar temp_color = cvGet2D(area, i, j);
	//		
	//		for (int k = 0; k < 3; ++k)
	//		{
	//			color += temp_color.val[k];
	//			color *= 1000;
	//		}
	//		if (color_map.count(color) == 0)
	//		{
	//			color_map[color] = color_idx++;
	//		}
	//		point_data[i][j].area_idx = color_map[color];
	//		if (point_data[i][j].isEdge) point_data[i][j].area_idx = -1;
	//	}
	//}
	//cvReleaseImage(&area);
	IplImage *edge = cvCreateImage(cvGetSize(image), image->depth, 1);
	
	helper.get_pic_edge(image, edge);
	uchar *data = (uchar*)edge->imageData;
	int step = edge->widthStep / sizeof(uchar);
	for (int i = 0; i < image->height; ++i)
	{
		for (int j = 0; j < image->width; ++j)
		{
			point_data[i][j].x = i;
			point_data[i][j].y = j;
			point_data[i][j].color = cvGet2D(image, i, j);
			point_data[i][j].point_type = None;
			point_data[i][j].point = mesh.add_vertex(MyMesh::Point(i, j, 0));
			if (data[i * step + j] == 255)
				point_data[i][j].isEdge = true;
			else
				point_data[i][j].isEdge = false;
		}
	}
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
			MyPoint point[4];
			point[0] = point_data[i][j];
			point[1] = point_data[i + 1][j];
			point[2] = point_data[i + 1][j + 1];
			point[3] = point_data[i][j + 1];
			int num = 0;
			for (int i = 0; i < 4; ++i)
			{
				if (point[i].isEdge)
					++num;
			}

			switch (num)
			{
			case 0:
				if (isContour)
				{
					list.clear();
					list.push_back(point[0].point);
					list.push_back(point[1].point);
					list.push_back(point[2].point);
					list.push_back(point[3].point);
					mesh.add_face(list);
				}
				break;
			case 1:
				list.clear();
				list.push_back(point[0].point);
				list.push_back(point[1].point);
				list.push_back(point[2].point);
				list.push_back(point[3].point);
				mesh.add_face(list);
				break;
			case 2:
				if (point[0].isEdge && point[2].isEdge)
				{
					list.clear();
					list.push_back(point[0].point);
					list.push_back(point[1].point);
					list.push_back(point[2].point);
					mesh.add_face(list);

					list.clear();
					list.push_back(point[2].point);
					list.push_back(point[3].point);
					list.push_back(point[0].point);
					mesh.add_face(list);
				}
				else
				{
					if (point[1].isEdge && point[3].isEdge)
					{
						list.clear();
						list.push_back(point[0].point);
						list.push_back(point[1].point);
						list.push_back(point[3].point);
						mesh.add_face(list);

						list.clear();
						list.push_back(point[1].point);
						list.push_back(point[2].point);
						list.push_back(point[3].point);
						mesh.add_face(list);
					}
					else
					{
						list.clear();
						list.push_back(point[0].point);
						list.push_back(point[1].point);
						list.push_back(point[2].point);
						list.push_back(point[3].point);
						mesh.add_face(list);
					}
				}
				break;
			case 3:
				list.clear();
				list.push_back(point[0].point);
				list.push_back(point[1].point);
				list.push_back(point[2].point);
				list.push_back(point[3].point);
				mesh.add_face(list);
				break;
			case 4:
				list.clear();
				list.push_back(point[0].point);
				list.push_back(point[1].point);
				list.push_back(point[2].point);
				list.push_back(point[3].point);
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
	all_num = num_all_vertices;
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
	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++ i, ++ i)
	{
		int from_x = mesh.point(mesh.from_vertex_handle(i.handle())).data()[0];
		int from_y = mesh.point(mesh.from_vertex_handle(i.handle())).data()[1];
		int to_x = mesh.point(mesh.to_vertex_handle(i.handle())).data()[0];
		int to_y = mesh.point(mesh.to_vertex_handle(i.handle())).data()[1];
		int difference = 0;
		Pair newPair;
		newPair.edge = i.handle();
		if ((point_data[from_x][from_y].isEdge && point_data[to_x][to_y].isEdge) ||
			(!(point_data[from_x][from_y].isEdge || point_data[to_x][to_y].isEdge)))
		{
			CvScalar from_color = point_data[from_x][from_y].color;
			CvScalar to_color = point_data[to_x][to_y].color;
			for (int j = 0; j < 3; ++ j)
			{
				int temp = abs(from_color.val[j] - to_color.val[j]);
				difference += temp;
			}
		}
		newPair.value = difference;
		point_pair_list.push_back(newPair);
	}
	cout << "Init pair list complete" << endl;
}
//************************************
// Method:    LinkVertices
// FullName:  OpenmeshHelper::LinkVertices
// Access:    private
// Returns:   void
// Qualifier:
// 建立以点为核心的相关边链表，也就是找到与i-th点相联系的所有边，并连接成链表，相反的两天半边只计算一条
//************************************
void OpenmeshHelper::LinkVertices()
{
	ListLink temp_list_link;
	temp_list_link.head = nullptr;
	temp_list_link.last = nullptr;
	vector<ListLink> temp_list(num_all_vertices, temp_list_link);
	vertex_list = temp_list;

	EdgeLink *temp_link;
	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++ i, ++ i)
	{
		int from_x = mesh.point(mesh.from_vertex_handle(i.handle())).data()[0];
		int from_y = mesh.point(mesh.from_vertex_handle(i.handle())).data()[1];
		int to_x = mesh.point(mesh.to_vertex_handle(i.handle())).data()[0];
		int to_y = mesh.point(mesh.to_vertex_handle(i.handle())).data()[1];
		if ((point_data[from_x][from_y].isEdge && point_data[to_x][to_y].isEdge) ||
			(!(point_data[from_x][from_y].isEdge || point_data[to_x][to_y].isEdge)))	
		{
			int vertex_idx = mesh.to_vertex_handle(i.handle()).idx();
			if (vertex_list[vertex_idx].head == nullptr)
			{
				vertex_list[vertex_idx].head = new EdgeLink;
				vertex_list[vertex_idx].head->pair = point_pair_list[i.handle().idx() / 2];
				vertex_list[vertex_idx].head->next = nullptr;

				vertex_list[vertex_idx].last = vertex_list[vertex_idx].head;
			}
			else
			{
				temp_link = new EdgeLink;
				temp_link->pair = point_pair_list[i.handle().idx() / 2];
				temp_link->next = nullptr;
				vertex_list[vertex_idx].last->next = temp_link;
				vertex_list[vertex_idx].last = vertex_list[vertex_idx].last->next;
			}

			vertex_idx = mesh.from_vertex_handle(i.handle()).idx();
			if (vertex_list[vertex_idx].head == nullptr)
			{
				vertex_list[vertex_idx].head = new EdgeLink;
				vertex_list[vertex_idx].head->pair = point_pair_list[i.handle().idx() / 2];
				vertex_list[vertex_idx].head->next = nullptr;

				vertex_list[vertex_idx].last = vertex_list[vertex_idx].head;
			}
			else
			{
				temp_link = new EdgeLink;
				temp_link->pair = point_pair_list[i.handle().idx() / 2];
				temp_link->next = nullptr;
				vertex_list[vertex_idx].last->next = temp_link;
				vertex_list[vertex_idx].last = vertex_list[vertex_idx].last->next;
			}
		}
		
	}
	all_pair_list = point_pair_list;
	cout << "Link vertrices complete" << endl;
}
//************************************
// Method:    DeletePairList
// FullName:  OpenmeshHelper::DeletePairList
// Access:    private 
// Returns:   void
// Qualifier:
// 删除撕裂点与非撕裂点连接的边
//************************************
void OpenmeshHelper::DeletePairList(int flag)
{
	//flag == 1 代表是平滑点
	//flag == 2 代表是撕裂点
	vector<Pair> temp_list;
	for (auto i = all_pair_list.begin(); i != all_pair_list.end(); ++i)
	{
		/*int from_x = mesh.point(mesh.from_vertex_handle(i->edge)).data()[0];
		int from_y = mesh.point(mesh.from_vertex_handle(i->edge)).data()[1];
		int to_x = mesh.point(mesh.to_vertex_handle(i->edge)).data()[0];
		int to_y = mesh.point(mesh.to_vertex_handle(i->edge)).data()[1];*/
		int from_idx = mesh.from_vertex_handle(i->edge).idx();
		int to_idx = mesh.to_vertex_handle(i->edge).idx();
		if (flag == 1)
		{
			if (!(control_list[from_idx].isEdge || control_list[to_idx].isEdge))
			{
				temp_list.push_back(*i);
			}
			if (control_list[from_idx].isEdge && control_list[to_idx].isEdge)
			{
				temp_list.push_back(*i);
			}
		}
		else
		{
			if (control_list[from_idx].isEdge && control_list[to_idx].isEdge)
			{
				temp_list.push_back(*i);
			}
		}
		
	}
	point_pair_list = temp_list;
	cout << "Delete pair list complete" << endl;
}
void OpenmeshHelper::SortVertices()
{
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		int tear = 0;
		int crease = 0; 
		int x = mesh.point(i.handle()).data()[0];
		int y = mesh.point(i.handle()).data()[1];
		if ((x == 0 && y == 0) || (x == 0 & y == image->height - 1)
			|| (x == image->width - 1 && y == 0) || (x == image->width - 1 && y == image->height - 1))
		{
			point_data[x][y].point_type = Corner;
			continue;
		}
		for (auto it = mesh.voh_iter(i.handle()); it; ++it)
		{
			//auto position = find(tear_list.begin(), tear_list.end(), it.handle());
			//if (position != tear_list.end())
			//{
			//	++tear;
			//	continue;
			//}

			auto position = find(crease_list.begin(), crease_list.end(), it.handle());
			if (position != crease_list.end())
			{
				++crease;
			}
		}
		
		if (tear == 0 && crease == 0)
		{
			point_data[x][y].point_type = Smooth;
			continue;
		}
	/*	if (tear == 2)
		{
			point_data[x][y].point_type = Tear;
			continue;
		}*/
		if (crease == 2)
		{
			point_data[x][y].point_type = Crease;
			continue;
		}
		point_data[x][y].point_type = Corner;
	}
	
	ofstream outfile("D:/type.txt");
	for (int i = 0; i != image->height; ++i)
	{
		for (int j = 0; j != image->width; ++j)
		{
			outfile << point_data[i][j].point_type << ' ';
		}
		outfile << endl;
	}
	outfile.close();
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
void OpenmeshHelper::Output(string output_location)
{
	try
	{
		if (!OpenMesh::IO::write_mesh(mesh, output_location))
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
	EdgeLink *remove_head = vertex_list[info.v0.idx()].head;
	EdgeLink *remain_head = vertex_list[info.v1.idx()].head;

	mesh.collapse(half);
	//-- num_vertices;
	//调整合并后点的位置
	double x, y;
	RegulatePosition(info, &x, &y);
	mesh.point(info.v1).data()[0] = x;
	mesh.point(info.v1).data()[1] = y;
	//调整点的颜色
	int remove_idx = info.v0.idx();
	int remain_idx = info.v1.idx();

	CvScalar remove_color = control_list[remove_idx].color;
	CvScalar remain_color = control_list[remove_idx].color;
	for (int k = 0; k < 3; ++k)
	{
		control_list[remain_idx].color.val[k] = (remove_color.val[k] + remain_color.val[k]) / 2;
	}

	//删除相邻的边
	while (remove_head)
	{
		Pair temp_pair = remove_head->pair;
		int stop = temp_pair.edge.idx();
		remove_head = remove_head->next;
		auto delete_it = find(point_pair_list.begin(), point_pair_list.end(), temp_pair);
		if (delete_it != point_pair_list.end())
			point_pair_list.erase(delete_it);
		//remove(point_pair_list.begin(), point_pair_list.end(), temp_pair);
		//cout << "delete:" << temp_pair.edge.idx() << endl;
	}

	while (remain_head)
	{
		Pair temp_pair = remain_head->pair;
		remain_head = remain_head->next;
		auto delete_it = find(point_pair_list.begin(), point_pair_list.end(), temp_pair);
		if (delete_it != point_pair_list.end())
			point_pair_list.erase(delete_it);
		//remove(point_pair_list.begin(), point_pair_list.end(), temp_pair);
		//cout << "delete:" << temp_pair.edge.idx() << endl;
	}
	//重新计算相连的边的代价，并加入vector
	vertex_list[info.v1.idx()].head = nullptr;
	for (auto it = mesh.voh_iter(info.v1); it; ++it)
	{
		Pair temp_pair;
		if ((control_list[mesh.from_vertex_handle(it).idx()].isEdge && control_list[mesh.to_vertex_handle(it).idx()].isEdge) ||
			(!control_list[mesh.from_vertex_handle(it).idx()].isEdge && !control_list[mesh.to_vertex_handle(it).idx()].isEdge))
		{
			if (mesh.from_vertex_handle(it) != mesh.to_vertex_handle(it))
			{
				temp_pair.edge = it.handle();
				//cout << "input:" << it.handle().idx() << endl;
				int from_idx = mesh.from_vertex_handle(temp_pair.edge).idx();
				int to_idx = mesh.to_vertex_handle(temp_pair.edge).idx();
				int difference = 0;
				CvScalar from_color = control_list[from_idx].color;
				CvScalar to_color = control_list[to_idx].color;
				for (int j = 0; j < 3; ++j)
				{
					int temp = abs(from_color.val[j] - to_color.val[j]);
					difference += temp;
				}
				temp_pair.value = difference;
				point_pair_list.push_back(temp_pair);

				if (vertex_list[info.v1.idx()].head == nullptr)
				{
					vertex_list[info.v1.idx()].head = new EdgeLink;
					vertex_list[info.v1.idx()].head->pair = temp_pair;
					vertex_list[info.v1.idx()].head->next = nullptr;

					vertex_list[info.v1.idx()].last = vertex_list[info.v1.idx()].head;
				}
				else
				{
					EdgeLink *temp_link;
					temp_link = new EdgeLink;
					temp_link->pair = temp_pair;
					temp_link->next = nullptr;
					vertex_list[info.v1.idx()].last->next = temp_link;
					vertex_list[info.v1.idx()].last = vertex_list[info.v1.idx()].last->next;
				}
			}
		}
		
		
	}
	make_heap(point_pair_list.begin(), point_pair_list.end(), Compare);
}

void OpenmeshHelper::InsertCreaseEdge()
{
	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
	{
		int from_x = mesh.point(mesh.from_vertex_handle(i.handle())).data()[0];
		int from_y = mesh.point(mesh.from_vertex_handle(i.handle())).data()[1];
		int to_x = mesh.point(mesh.to_vertex_handle(i.handle())).data()[0];
		int to_y = mesh.point(mesh.to_vertex_handle(i.handle())).data()[1];
		if (point_data[from_x][from_y].isEdge && point_data[to_x][to_y].isEdge)
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
	InitPointData();
	ConnectMesh(true);
	InsertCreaseEdge();
		
	CountVertices();
	SortVertices();
	InitPairList();
	LinkVertices();
	
	mesh.request_edge_status();
	mesh.request_vertex_status();
	mesh.request_halfedge_status();
	mesh.request_face_status();
	
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		int x = mesh.point(i.handle()).data()[0];
		int y = mesh.point(i.handle()).data()[1];
		control_list.push_back(point_data[x][y]);
	}

	DeletePairList(1);
	LoopReduce(rate, visual);

	vector<MyPoint> temp_list;
	for (auto i = 0; i != control_list.size(); ++ i)
	{
		MyPoint temp_point = control_list[i];
		if (!mesh.status(temp_point.point).deleted())
		{
			temp_list.push_back(temp_point);
		}
	}
	control_list = temp_list;
	MyPoint testPoint;
	for (int i = 0; i != control_list.size(); ++i)
	{
		Position temp_position;
		temp_position.x = mesh.point(control_list[i].point).data()[0];
		temp_position.y = mesh.point(control_list[i].point).data()[1];
		control_map[temp_position] = control_list[i];
	}
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

	mesh.garbage_collection();
	crease_list.clear();
	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
	{
		Position from, to;
		from.x = mesh.point(mesh.from_vertex_handle(i)).data()[0];
		from.y = mesh.point(mesh.from_vertex_handle(i)).data()[1];
		to.x = mesh.point(mesh.to_vertex_handle(i)).data()[0];
		to.y = mesh.point(mesh.to_vertex_handle(i)).data()[1];

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
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_vertex_normals();
	mesh.update_face_normals();
	mesh.release_face_normals();
	mesh.release_vertex_normals();

	mesh.release_edge_status();
	mesh.release_vertex_status();
	mesh.release_halfedge_status();
	mesh.release_face_status();
	Output("D:/before.off");
	OptimizePosition();
	cout << "Reduce Vertrices complete" << endl;
}

//************************************
// Method:    IsCollapseOK
// FullName:  OpenmeshHelper::IsCollapseOK
// Access:    private 
// Returns:   bool
// Qualifier:
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

	mesh.point(info.v0).data()[0] = x;
	mesh.point(info.v0).data()[1] = y;
	
	for (auto i = mesh.vf_begin(info.v0); i != mesh.vf_end(info.v0); ++i)
	{
		auto normals = mesh.calc_face_normal(i.handle());
		if (normals.data()[2] < 0)
		{
			mesh.point(info.v0).data()[0] = x0;
			mesh.point(info.v0).data()[1] = y0;
			return false;
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
	//while (num_vertices / num_all_vertices > rate)
	int times = 0;
	int last_num = num_all_vertices;
	while (num_vertices / num_all_vertices > rate)
	{
		while (point_pair_list.size() > 0)
		{
			pop_heap(point_pair_list.begin(), point_pair_list.end(), Compare);
			Pair pair = point_pair_list[point_pair_list.size() - 1];
			auto half = pair.edge;
			point_pair_list.pop_back();
			if (IsCollapseOK(half))
			{
				CollapseEdge(half);
				--num_vertices;
			}
			else
			{
				half = mesh.opposite_halfedge_handle(half);
				if (IsCollapseOK(half))
				{
					cout << "opposite" << endl;
					--num_vertices;
					CollapseEdge(half);
				}
				else
				{
					second_pair_list.push_back(pair);
				}
			}
			if (visual)
				cout << point_pair_list.size() << endl;
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
	//这个部分是用于重组控制点的列表，因为在简化之后，删除点，进行garbage_collection，导致点类型的混乱，
	//因此在简化后通过记录点的位置，生成map，从而记录点的类型和颜色信息，在点garbage_collection之后
	//再通过点的位置重新找回这些点的颜色和类型信息
	//经过测试garbage_collection之后，并不只是改变点的序号，因此先前存储的点类型变量，需要进行更换
	control_list.clear();
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		Position temp_position;
		double x = mesh.point(i.handle()).data()[0];
		double y = mesh.point(i.handle()).data()[1];
		temp_position.x = x;
		temp_position.y = y;
		int tear = 0;
		int crease = 0;
		MyPoint temp_point;
		temp_point = control_map[temp_position];
		temp_point.point = i.handle();

		if ((x == 0 && y == 0) || (x == 0 & y == image->height - 1)
			|| (x == image->width - 1 && y == 0) || (x == image->width - 1 && y == image->height - 1))
		{
			temp_point.point_type = Corner;
			control_list.push_back(temp_point);
		}
		for (auto it = mesh.voh_iter(i.handle()); it; ++it)
		{
			/*	auto position = find(tear_list.begin(), tear_list.end(), it.handle());
			if (position != tear_list.end())
			{
			++tear;
			continue;
			}
			*/
			auto position = find(crease_list.begin(), crease_list.end(), it.handle());
			if (position != crease_list.end())
			{
				++crease;
			}
		}

		//TODO 点分类存在问题
		if (crease == 0)
		{
			temp_point.point_type = Smooth;
		}
		/*	if (tear == 2)
		{
		point_data[x][y].point_type = Tear;
		continue;
		}*/
		if (crease == 2)
		{
			temp_point.point_type = Crease;
		}
		if (crease != 2 && crease != 0)
		{
			temp_point.point_type = Corner;
		}
		

		control_list.push_back(temp_point);
	}

	Eigen::MatrixXf X, Y;
	if (exact)
	{
		X = Eigen::MatrixXf::Zero(1, control_list.size() * 2);
	}
	else
	{
		X = Eigen::MatrixXf::Zero(2, control_list.size());
	}
	Y = Eigen::MatrixXf::Zero(control_list.size(), control_list.size());
	for (int i = 0; i != control_list.size(); ++ i)
	{
		//TODO 未考虑Tear点
		int x = mesh.point(control_list[i].point).data()[0];
		int y = mesh.point(control_list[i].point).data()[1];
		X(0, i) = control_list[i].x;
		if (exact)
		{
			X(0, i + control_list.size()) = control_list[i].y;
		}
		else
		{
			X(1, i) = control_list[i].y;
		}
		
		//处理Y
		if (control_list[i].point_type == Corner)
		{
			Y(i, control_list[i].point.idx()) = 1;
			continue;
		}
		if (control_list[i].point_type == Crease)
		{
			Y(i, control_list[i].point.idx()) = (double)4 / 6;
			int times = 0;
			for (auto j = mesh.voh_begin(control_list[i].point); j != mesh.voh_end(control_list[i].point); ++j)
			{
				auto position = find(crease_list.begin(), crease_list.end(), j.handle());
				if (position != crease_list.end())
				{
					++times;
					Y(i, mesh.to_vertex_handle(j.handle()).idx()) = (double)1 / 6;
				}
			}
			if (times != 2)
			{
				cout << times << endl;
			}
			continue;
		}
		auto point = control_list[i].point;
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
			Y(i, out_point_list[0].idx()) = 4.0 / 6;
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
		YY = Eigen::MatrixXf::Zero(2 * control_list.size(), 2 * control_list.size());
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
			}
		}
		ofstream parameter("D:/parameter_matrix.txt");
		for (int i = 0; i != 2 * control_list.size(); ++i)
		{
			for (int j = 0; j != 2 * control_list.size(); ++j)
			{
				parameter << YY(i, j) << ' ';
			}
			parameter << endl;
		}
		cout << "parameter output" << endl;
		Eigen::MatrixXf optimize = YY.colPivHouseholderQr().solve(X.transpose());
		ofstream optimizeOut("D:/optimize.txt");
		for (int i = 0; i != control_list.size(); ++i)
		{
			optimizeOut << optimize(i, 0) << "," << optimize(i + control_list.size(), 0) << endl;
			mesh.point(control_list[i].point).data()[0] = optimize(i, 0);
			mesh.point(control_list[i].point).data()[1] = optimize(i + control_list.size(), 0);
		}
		optimizeOut.close();
	}
	else
	{
		ofstream parameter("D:/parameter_matrix.txt");
		for (int i = 0; i != control_list.size(); ++i)
		{
			for (int j = 0; j != control_list.size(); ++j)
			{
				parameter << Y(i, j) << ' ';
			}
			parameter << endl;
		}
		cout << "parameter output" << endl;
		Eigen::MatrixXf optimize = Y.colPivHouseholderQr().solve(X.transpose());
		ofstream optimizeOut("D:/optimize.txt");
		for (int i = 0; i != control_list.size(); ++i)
		{
			optimizeOut << optimize(i, 0) << "," << optimize(i, 1) << endl;
			mesh.point(control_list[i].point).data()[0] = optimize(i, 0);
			mesh.point(control_list[i].point).data()[1] = optimize(i, 1);
		}
		optimizeOut.close();
	}
	

}
