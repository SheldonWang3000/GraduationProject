#include "openmesh_helper.h"
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
	for (int i = 0; i != image->height; ++ i)
	{
		point_data[i] = new MyPoint[image->width];
	}

	IplImage *result = cvCreateImage(cvGetSize(image), image->depth, 1);
	OpencvHelper helper;
	helper.pic_edge(image, result);
	uchar *data = (uchar*)result->imageData;
	int step = result->widthStep / sizeof(uchar);

	for (int i = 0; i < image->height; ++ i)
	{
		for (int j = 0; j < image->width; ++ j)
		{
			point_data[i][j].color = cvGet2D(image, i, j);
			point_data[i][j].point = mesh.add_vertex(MyMesh::Point(i, j, 0)); 
			if (data[i * step + j] == 255)
				point_data[i][j].isEdge = true;
			else point_data[i][j].isEdge = false;
		}
	}
	cvReleaseImage(&result);
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
void OpenmeshHelper::ConnectMesh()
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
				list.clear();
				list.push_back(point[0].point);
				list.push_back(point[1].point);
				list.push_back(point[2].point);
				list.push_back(point[3].point);
				mesh.add_face(list);
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
void OpenmeshHelper::DeletePairList()
{
	vector<Pair> temp_list;
	for (auto i = point_pair_list.begin(); i != point_pair_list.end(); ++i)
	{
		int from_x = mesh.point(mesh.from_vertex_handle(i->edge)).data()[0];
		int from_y = mesh.point(mesh.from_vertex_handle(i->edge)).data()[1];
		int to_x = mesh.point(mesh.to_vertex_handle(i->edge)).data()[0];
		int to_y = mesh.point(mesh.to_vertex_handle(i->edge)).data()[1];
		if ((point_data[from_x][from_y].isEdge && point_data[to_x][to_y].isEdge) ||
			(!(point_data[from_x][from_y].isEdge || point_data[to_x][to_y].isEdge)))
		{
			temp_list.push_back(*i);
		}
	}
	point_pair_list = temp_list;
	cout << "Delete pair list complete" << endl;
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
// Parameter: OpenMesh::Decimater::CollapseInfoT<MyMesh> info
// 计算点合并之后的位置，如果被消除的点是角落点，则角落点不变，如果2个点都是四边上的点，则如果在同一边则取中点
// 如果不是同一边，则位置不变，如果2点中有一个点是四边点，则点的位置计算到边上，如果都不是则取中点
//************************************
void OpenmeshHelper::RegulatePosition(OpenMesh::Decimater::CollapseInfoT<MyMesh> info)
{
	int remove_x = info.p0.data()[0];
	int remove_y = info.p0.data()[1];
	int remain_x = info.p1.data()[0];
	int remain_y = info.p1.data()[1];

	if ((remove_x == 0 && remove_y == 0) ||
		(remove_x == 0 && remove_y == image->width - 1) ||
		(remove_x == image->height - 1 && remove_y == 0) ||
		(remove_x == image->height - 1 && remove_y == image->width - 1))
	{
		mesh.point(info.v1).data()[0] = remove_x;
		mesh.point(info.v1).data()[1] = remove_y;
		return;
	}
	if ((remain_x == 0 && remain_y == 0) ||
		(remain_x == 0 && remain_y == image->width - 1) ||
		(remain_x == image->height - 1 && remain_y == 0) ||
		(remain_x == image->height - 1 && remain_y == image->width - 1))
	{
		return;
	}
	if (remove_x == 0 || remove_y == 0 || remove_x == image->height - 1 || remove_y == image->width - 1)
	{
		if (remain_x == 0 || remain_y == 0 || remain_x == image->height - 1 || remain_y == image->width - 1)
		{
			if ((remain_x == remove_x) || (remain_y == remove_y))
			{
				mesh.point(info.v1).data()[0] = (int)((info.p1.data()[0] + info.p0.data()[0]) / 2);
				mesh.point(info.v1).data()[1] = (int)((info.p1.data()[1] + info.p0.data()[1]) / 2);
			}
			else
			{
				return;
			}
		}
		else
		{
			mesh.point(info.v1).data()[0] = remove_x;
			mesh.point(info.v1).data()[1] = remove_y;
		}

	}
	else
	{
		if (!(remain_x == 0 || remain_y == 0 || remain_x == image->height - 1 || remain_y == image->width - 1))
		{
			mesh.point(info.v1).data()[0] = (int)((info.p1.data()[0] + info.p0.data()[0]) / 2);
			mesh.point(info.v1).data()[1] = (int)((info.p1.data()[1] + info.p0.data()[1]) / 2);
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
	-- num_vertices;
	
	//调整合并后点的位置
	RegulatePosition(info);
	//调整点的颜色
	int x = mesh.point(info.v1).data()[0];
	int y = mesh.point(info.v1).data()[1];
	int remove_x = info.p0.data()[0];
	int remove_y = info.p0.data()[1];
	int remain_x = info.p1.data()[0];
	int remain_y = info.p1.data()[1];
	CvScalar remove_color = point_data[remove_x][remove_y].color;
	CvScalar remain_color = point_data[remain_x][remain_y].color;
	for (int k = 0; k < 3; ++k)
	{
		point_data[x][y].color.val[k] = (remove_color.val[k] + remain_color.val[k]) / 2;
	}

	//删除相邻的边
	do
	{
		Pair temp_pair = remove_head->pair;
		remove_head = remove_head->next;
		remove(point_pair_list.begin(), point_pair_list.end(), temp_pair);
	} while (remove_head->next);

	do
	{
		Pair temp_pair = remain_head->pair;
		remain_head = remain_head->next;
		remove(point_pair_list.begin(), point_pair_list.end(), temp_pair);
	} while (remain_head->next);
	//重新计算相连的边的代价，并加入vector
	for (auto it = mesh.voh_iter(info.v1); it; ++it)
	{
		Pair temp_pair;
		temp_pair.edge = it.handle();
		int from_x = mesh.point(mesh.from_vertex_handle(temp_pair.edge)).data()[0];
		int from_y = mesh.point(mesh.from_vertex_handle(temp_pair.edge)).data()[1];
		int to_x = mesh.point(mesh.to_vertex_handle(temp_pair.edge)).data()[0];
		int to_y = mesh.point(mesh.to_vertex_handle(temp_pair.edge)).data()[1];
		int difference = 0;
		CvScalar from_color = point_data[from_x][from_y].color;
		CvScalar to_color = point_data[to_x][to_y].color;
		for (int j = 0; j < 3; ++ j)
		{
			int temp = abs(from_color.val[j] - to_color.val[j]);
			difference += temp;
		}
		temp_pair.value = difference;
		point_pair_list.push_back(temp_pair);
	}
	make_heap(point_pair_list.begin(), point_pair_list.end(), Compare);
}
//************************************
// Method:    ReduceVertices
// FullName:  OpenmeshHelper::ReduceVertices
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: double rate 表示最后控制点占原总数的比例
// Parameter: bool visual 表示是否显示计算中的比例，方便看程序运行到什么时候
//************************************
void OpenmeshHelper::ReduceVertices(double rate, bool visual)
{
	InitPointData();
	ConnectMesh();
	CountVertices();

	InitPairList();
	LinkVertices();
	DeletePairList();

	make_heap(point_pair_list.begin(), point_pair_list.end(), Compare);

	mesh.request_edge_status();
	mesh.request_vertex_status();
	mesh.request_halfedge_status();
	mesh.request_face_status();

	while (num_vertices / num_all_vertices > rate)
		//while (point_pair_list.size() > 0)
	{
		pop_heap(point_pair_list.begin(), point_pair_list.end(), Compare);
		auto half = point_pair_list[point_pair_list.size() - 1].edge;
		point_pair_list.pop_back();

		if (mesh.is_collapse_ok(half))
		{
			CollapseEdge(half);
		}
		else
		{
			half = mesh.opposite_halfedge_handle(half);
			if (mesh.is_collapse_ok(half))
			{
				cout << "opposite" << endl;
				CollapseEdge(half);
			}
		}
		if (visual)
			cout << num_vertices / num_all_vertices << endl;
	}

	mesh.garbage_collection();

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
	cout << "Reduce Vertrices complete" << endl;
}