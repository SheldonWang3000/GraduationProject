#include "openmesh_helper.h"
#include <iostream>
using namespace std;

//************************************
// Method:    OpenmeshHelper
// FullName:  OpenmeshHelper::OpenmeshHelper
// Access:    public 
// Returns:   
// Qualifier: 
// Parameter: string input_location input�ļ���λ��
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
// ��ʼ�����ص����飬ÿ�������һ��ԭͼ���ص㣬������ɫֵ���Ƿ�Ϊ˺�ѵ㣬�Լ�openmesh��vertex����
//************************************
void OpenmeshHelper::InitPointData()
{
	point_data = new MyPoint*[image->height];
	for (int i = 0; i != image->height; ++ i)
	{
		point_data[i] = new MyPoint[image->width];
	}


	OpencvHelper helper;

	IplImage *contours = cvCreateImage(cvGetSize(image), image->depth, 1);
	helper.get_pic_contours(image, contours);
	uchar *data = (uchar *)contours->imageData;
	int step = contours->widthStep / sizeof(uchar);
	for (int i = 0; i < contours->height; ++i)
	{
		for (int j = 0; j < contours->width; ++j)
		{
			point_data[i][j].point_type = None;
			point_data[i][j].point = mesh.add_vertex(MyMesh::Point(i, j, 0));
			if (data[i * step + j] == 0)
			{
				point_data[i][j].isEdge = true;
			}
		}
	}
	ConnectMesh(false);

	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
	{
		int from_x = mesh.point(mesh.from_vertex_handle(i.handle())).data()[0];
		int from_y = mesh.point(mesh.from_vertex_handle(i.handle())).data()[1];
		int to_x = mesh.point(mesh.to_vertex_handle(i.handle())).data()[0];
		int to_y = mesh.point(mesh.to_vertex_handle(i.handle())).data()[1];
		if (point_data[from_x][from_y].isEdge && point_data[to_x][to_y].isEdge)
		{
			tear_list.push_back(i.handle());
		}
	}
	mesh.clear();
	for (int i = 0; i < image->height; ++i)
	{
		for (int j = 0; j < image->width; ++j)
		{
			point_data[i][j].point = mesh.add_vertex(MyMesh::Point(i, j, 0));
		}
	}

	cvReleaseImage(&contours);
	IplImage *area = cvCreateImage(cvGetSize(image), image->depth, 3);
	helper.get_pic_area(image, area);
	data = (uchar *)area->imageData;
	step = area->widthStep / sizeof(uchar);
	map<int, int> color_map;
	int color_idx = 0;
	for (int i = 0; i < area->height; ++ i)
	{
		for (int j = 0; j < area->width; ++ j)
		{
			int color = 0;
			CvScalar temp_color = cvGet2D(area, i, j);
			point_data[i][j].color = temp_color;
			
			
			for (int k = 0; k < 3; ++k)
			{
				color += temp_color.val[k];
				color *= 1000;
			}
			if (color_map.count(color) == 0)
			{
				color_map[color] = color_idx++;
			}
			point_data[i][j].area_idx = color_map[color];
			if (point_data[i][j].isEdge) point_data[i][j].area_idx = -1;
		}
	}
	cvReleaseImage(&area);
	IplImage *edge = cvCreateImage(cvGetSize(image), image->depth, 1);
	
	helper.get_pic_edge(image, edge);
	data = (uchar*)edge->imageData;
	step = edge->widthStep / sizeof(uchar);
	for (int i = 0; i < image->height; ++i)
	{
		for (int j = 0; j < image->width; ++j)
		{
			if (data[i * step + j] == 255)
				point_data[i][j].isEdge = true;
			else
				point_data[i][j].isEdge = false;
		}
	}
	ConnectMesh(false);
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
	}
	mesh.clear();
	for (int i = 0; i < image->height; ++i)
	{
		for (int j = 0; j < image->width; ++j)
		{
			point_data[i][j].point = mesh.add_vertex(MyMesh::Point(i, j, 0));
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
// ����������ѭ֮ǰ�Ĺ��򽫵����ӳ�����
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
// ����openmesh�ж�����������ֹͣ��
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
// ��ʼ������Լ��ϣ��������˺�ѵ����˺�ѵ�����Ӿͼ���ɫ����Ϊ����
// �˴��Ὣ���еĶ���Լ��룬Ϊ����֮����������ʱ����԰��ձ���ҵ���Ӧ�ߣ�֮����ɾ��˺�ѵ����˺�ѵ����ӵı�
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
// �����Ե�Ϊ���ĵ���ر�������Ҳ�����ҵ���i-th������ϵ�����бߣ������ӳ��������෴��������ֻ����һ��
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
// ɾ��˺�ѵ����˺�ѵ����ӵı�
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
void OpenmeshHelper::SortVertices()
{
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		int tear = 0;
		int crease = 0;
		for (auto it = mesh.voh_iter(i.handle()); it; ++it)
		{
			auto position = find(tear_list.begin(), tear_list.end(), it.handle());
			if (position != tear_list.end())
			{
				++tear;
				continue;
			}

			position = find(crease_list.begin(), crease_list.end(), it.handle());
			if (position != crease_list.end())
			{
				++crease;
			}
		}
		int x = mesh.point(i.handle()).data()[0];
		int y = mesh.point(i.handle()).data()[1];
		//TODO ������������
		if (tear == 0 && crease == 0)
		{
			point_data[x][y].point_type = Smooth;
			continue;
		}
		if (tear == 2)
		{
			point_data[x][y].point_type = Tear;
			continue;
		}
		if (crease == 2)
		{
			point_data[x][y].point_type = Crease;
			continue;
		}
		point_data[x][y].point_type = Corner;
	}
}
//************************************
// Method:    Output
// FullName:  OpenmeshHelper::Output
// Access:    public 
// Returns:   void
// Qualifier:
// Parameter: string output_location �����λ��
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
// �����ϲ�֮���λ�ã�����������ĵ��ǽ���㣬�����㲻�䣬���2���㶼���ı��ϵĵ㣬�������ͬһ����ȡ�е�
// �������ͬһ�ߣ���λ�ò��䣬���2������һ�������ıߵ㣬����λ�ü��㵽���ϣ������������ȡ�е�
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
// �ϲ��ߣ����ںϲ�֮�������ɫֵ������λ�ã���ɾ���ϲ�֮ǰ����رߣ����¼�����۲�������ر�
//************************************
void OpenmeshHelper::CollapseEdge(MyMesh::HalfedgeHandle half)
{
	auto info = OpenMesh::Decimater::CollapseInfoT<MyMesh>::CollapseInfoT(mesh, half);
	EdgeLink *remove_head = vertex_list[info.v0.idx()].head;
	EdgeLink *remain_head = vertex_list[info.v1.idx()].head;

	mesh.collapse(half);
	-- num_vertices;
	//�����ϲ�����λ��
	RegulatePosition(info);
	//���������ɫ
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

	//ɾ�����ڵı�
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
	//���¼��������ıߵĴ��ۣ�������vector
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
// Parameter: double rate ��ʾ�����Ƶ�ռԭ�����ı���
// Parameter: bool visual ��ʾ�Ƿ���ʾ�����еı��������㿴�������е�ʲôʱ��
//************************************
void OpenmeshHelper::ReduceVertices(double rate, bool visual)
{
	InitPointData();
	ConnectMesh(true);
	Output("D:/out.off");
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
		if (IsCollapseOK(half))
		{
			CollapseEdge(half);
		}
		else
		{
			half = mesh.opposite_halfedge_handle(half);
			if (IsCollapseOK(half))
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

bool OpenmeshHelper::IsCollapseOK(MyMesh::HalfedgeHandle half)
{
	return mesh.is_collapse_ok(half);

	if (!mesh.is_collapse_ok(half)) return false;
	auto info = OpenMesh::Decimater::CollapseInfoT<MyMesh>::CollapseInfoT(mesh, half); 
	int x0 = mesh.point(info.v0).data()[0];
	int y0 = mesh.point(info.v0).data()[1];
	int x1 = mesh.point(info.v1).data()[0];
	int y1 = mesh.point(info.v1).data()[1];

	mesh.point(info.v0).data()[0] = (x0 + x1) / 2;
	mesh.point(info.v0).data()[1] = (y0 + y1) / 2;
	
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

	mesh.point(info.v1).data()[0] = (x0 + x1) / 2;
	mesh.point(info.v1).data()[1] = (y0 + y1) / 2;
	
	for (auto i = mesh.vf_begin(info.v1); i != mesh.vf_end(info.v1); ++i)
	{
		auto normals = mesh.calc_face_normal(i.handle());
		if (normals.data()[2] < 0)
		{
			mesh.point(info.v1).data()[0] = x0;
			mesh.point(info.v1).data()[1] = y0;
			return false;
		}
	}

	mesh.point(info.v1).data()[0] = x1;
	mesh.point(info.v1).data()[1] = y1;
	return true;
}