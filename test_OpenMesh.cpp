#include <iostream>
#include <opencv2/opencv.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/CollapseInfoT.hh>
#include "opencv_helper.h"
using namespace std;
typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

//struct Pair
//{
//	int x, y;
//};

struct Pair
{
	OpenMesh::HalfedgeHandle edge;
	int value;
};

struct EdgeLink
{
	Pair pair;
	EdgeLink *next;
};

struct ListLink
{
	EdgeLink *head, *last;
};
void output(MyMesh mesh)
{
	try
	{
		if ( !OpenMesh::IO::write_mesh(mesh, "D:/output.off") )
		{
			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
			return;
		}
	}
	catch( std::exception& x )
	{
		std::cerr << x.what() << std::endl;
		return;
	}
}
int main()
{
	MyMesh mesh;
	MyMesh::VertexHandle vhandle[4];
	vhandle[0] = mesh.add_vertex(MyMesh::Point(0, 0, 0));
	vhandle[1] = mesh.add_vertex(MyMesh::Point(1, 0, 0));
	vhandle[2] = mesh.add_vertex(MyMesh::Point(1, 1, 0));
	vhandle[3] = mesh.add_vertex(MyMesh::Point(0, 1, 0));
	
	vector<MyMesh::VertexHandle> list;
	list.push_back(vhandle[0]);
	list.push_back(vhandle[1]);
	list.push_back(vhandle[3]);
	mesh.add_face(list);
	list.clear();
	list.push_back(vhandle[3]);
	list.push_back(vhandle[1]);
	list.push_back(vhandle[2]);
	mesh.add_face(list);


	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++ i)
	{
		if (mesh.from_vertex_handle(i.handle()) == vhandle[1] && mesh.to_vertex_handle(i.handle()) == vhandle[3])
		{
			mesh.request_edge_status();
			mesh.request_halfedge_status();
			mesh.request_face_status();
			mesh.request_vertex_status();

			if (mesh.is_collapse_ok(i.handle()))
			{
				mesh.collapse(i.handle());
				cout << "true" << endl;
			}
			else cout << "false" << endl;
			
			mesh.release_edge_status();
			mesh.release_face_status();
			mesh.release_vertex_status();
			mesh.release_halfedge_status();
			break;
		}
	}
	mesh.garbage_collection();
	output(mesh);

	//MyMesh::VertexHandle vhandle[7];
	//vhandle[0] = mesh.add_vertex(MyMesh::Point(-1, 1, 0));
	//vhandle[1] = mesh.add_vertex(MyMesh::Point(-1, 3, 0));
	//vhandle[2] = mesh.add_vertex(MyMesh::Point(0, 0, 0));
	//vhandle[3] = mesh.add_vertex(MyMesh::Point(0, 2, 0));
	//vhandle[4] = mesh.add_vertex(MyMesh::Point(0, 4, 0));
	//vhandle[5] = mesh.add_vertex(MyMesh::Point(1, 1, 0));
	//vhandle[6] = mesh.add_vertex(MyMesh::Point(1, 3, 0));

	//// Add three quad faces
	//std::vector<MyMesh::VertexHandle> face_vhandles;

	//face_vhandles.push_back(vhandle[1]);
	//face_vhandles.push_back(vhandle[0]);
	//face_vhandles.push_back(vhandle[2]);
	//face_vhandles.push_back(vhandle[3]);
	//mesh.add_face(face_vhandles);

	//face_vhandles.clear();

	//face_vhandles.push_back(vhandle[1]);
	//face_vhandles.push_back(vhandle[3]);
	//face_vhandles.push_back(vhandle[6]);
	//face_vhandles.push_back(vhandle[4]);
	//mesh.add_face(face_vhandles);

	//face_vhandles.clear();

	//face_vhandles.push_back(vhandle[3]);
	//face_vhandles.push_back(vhandle[2]);
	//face_vhandles.push_back(vhandle[5]);
	//face_vhandles.push_back(vhandle[6]);
	//mesh.add_face(face_vhandles);

	////mesh.triangulate();
	////获取顶点对，加入vector
	//vector<Pair> point_pair_list;
	//for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++ i, ++ i)
	//{
	//	int difference = 0;
	//	Pair newPair;
	//	newPair.edge = i.handle();
	//	newPair.value = difference;
	//	point_pair_list.push_back(newPair);

	//}
	//cout << "list_complete" << endl;



	////计算顶点个数，用作停止条件，当剩余顶点数为原本的30%的时候停止
	//double num_vertices = 0;
	//for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++ i)
	//{
	//	++ num_vertices;
	//}
	//cout << num_vertices << endl;
	//double num_all_vertices = num_vertices;

	//ListLink temp_list_link;
	//temp_list_link.head = nullptr;
	//temp_list_link.last = nullptr;
	//vector<ListLink> vertex_list(num_all_vertices, temp_list_link);
	//EdgeLink *temp_link;

	//for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++ i, ++ i)
	//{
	//	int vertex_idx = mesh.to_vertex_handle(i.handle()).idx();
	//	if (vertex_list[vertex_idx].head == nullptr) 
	//	{
	//		vertex_list[vertex_idx].head = new EdgeLink;
	//		vertex_list[vertex_idx].head->pair = point_pair_list[i.handle().idx() / 2];
	//		vertex_list[vertex_idx].head->next = nullptr;

	//		vertex_list[vertex_idx].last = vertex_list[vertex_idx].head;
	//	}
	//	else
	//	{
	//		temp_link = new EdgeLink;
	//		temp_link->pair = point_pair_list[i.handle().idx() / 2];
	//		temp_link->next = nullptr;
	//		vertex_list[vertex_idx].last->next = temp_link;
	//		vertex_list[vertex_idx].last = vertex_list[vertex_idx].last->next;
	//	}



	//	vertex_idx = mesh.from_vertex_handle(i.handle()).idx();
	//	if (vertex_list[vertex_idx].head == nullptr) 
	//	{
	//		vertex_list[vertex_idx].head = new EdgeLink;
	//		vertex_list[vertex_idx].head->pair = point_pair_list[i.handle().idx() / 2];
	//		vertex_list[vertex_idx].head->next = nullptr;

	//		vertex_list[vertex_idx].last = vertex_list[vertex_idx].head;
	//	}
	//	else
	//	{
	//		temp_link = new EdgeLink;
	//		temp_link->pair = point_pair_list[i.handle().idx() / 2];
	//		temp_link->next = nullptr;
	//		vertex_list[vertex_idx].last->next = temp_link;
	//		vertex_list[vertex_idx].last = vertex_list[vertex_idx].last->next;
	//	}
	//}

	//for (auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++ it, ++ it)
	//{
	//	cout << it.handle().idx() << endl;
	//	cout << mesh.from_vertex_handle(it.handle()).idx() << "->" << mesh.to_vertex_handle(it.handle()).idx() << endl;
	//}


	//int k = 0;
	//for (auto i = vertex_list.begin(); i != vertex_list.end(); ++ i)
	//{
	//	cout << k ++ << endl;
	//	EdgeLink *temp = i->head;
	//	while (temp->next)
	//	{
	//		cout << temp->pair.edge.idx() << endl;
	//		temp = temp->next;
	//	}
	//	cout << temp->pair.edge.idx() << endl;
	//	cout << endl;
	//}

	//// Now find the edge between vertex vhandle[2]
	//// and vhandle[3]
	//vector<Pair> list;
	//vector<OpenMesh::HalfedgeHandle> edge_list;
	//for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++ i)
	//{
	//	Pair temp;
	//	temp.x = mesh.from_vertex_handle(i.handle()).idx();
	//	temp.y = mesh.to_vertex_handle(i.handle()).idx();
	//	edge_list.push_back(i.handle());
	//	list.push_back(temp);
	//}
	//
	//OpenMesh::HalfedgeHandle temp;
	//for(auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) 
	//{
	//	if (mesh.from_vertex_handle(it.handle()) == vhandle[0] &&
	//		mesh.to_vertex_handle(it.handle()) == vhandle[2])
	//		temp = it.handle();

	//	if(mesh.to_vertex_handle(it.handle()) == vhandle[3] &&
	//		mesh.from_vertex_handle(it.handle()) == vhandle[2]) 
	//	{
	//		// Collapse edge
	//		mesh.request_edge_status();
	//		mesh.request_vertex_status();

	//		mesh.request_halfedge_status();

	//		auto info = OpenMesh::Decimater::CollapseInfoT<MyMesh>::CollapseInfoT(mesh, it.handle());

	//		if (mesh.is_collapse_ok(it.handle()))
	//		{
	//			mesh.collapse(it.handle());
	//			//mesh.garbage_collection();
	//			//mesh.point(info.v1).data()[1] = 1;
	//		}
	//		
	//		for (int i = 0; i < list.size(); ++ i)
	//		{
	//			//(list[i].x == 3 && list[i].y == 2) || (list[i].x == 2 && list[i].y == 3)||
	//			/*if (
	//				(list[i].x == 5 && list[i].y == 6) || (list[i].x == 6 && list[i].y == 5))
	//				continue;*/
	//			cout << list[i].x << "->" << list[i].y << endl;
	//			cout << mesh.from_vertex_handle(edge_list[i]) << "->" << mesh.to_vertex_handle(edge_list[i]) << endl;
	//			cout << mesh.point(mesh.from_vertex_handle(edge_list[i])).data()[0] << ' '
	//				<< mesh.point(mesh.from_vertex_handle(edge_list[i])).data()[1] << ' '
	//				<< mesh.point(mesh.to_vertex_handle(edge_list[i])).data()[0] << ' '
	//				<< mesh.point(mesh.to_vertex_handle(edge_list[i])).data()[1] << endl;
	//			cout << endl;
	//		}

	//		/*int kk = 0;
	//		for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++ i)
	//		{
	//			++ kk;
	//			cout << mesh.from_vertex_handle(i.handle()) << "->" << mesh.to_vertex_handle(i.handle()) << endl;
	//		}
	//		cout << "kk:" << kk << endl;*/
	//		break;
	//	}
	//}
	system("pause");
	return 0;
}