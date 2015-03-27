#include <iostream>
#include <opencv2/opencv.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/CollapseInfoT.hh>
#include "opencv_helper.h"
using namespace std;
typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

struct Pair
{
	int x, y;
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
	MyMesh::VertexHandle vhandle[7];
	vhandle[0] = mesh.add_vertex(MyMesh::Point(-1, 1, 0));
	vhandle[1] = mesh.add_vertex(MyMesh::Point(-1, 3, 0));
	vhandle[2] = mesh.add_vertex(MyMesh::Point(0, 0, 0));
	vhandle[3] = mesh.add_vertex(MyMesh::Point(0, 2, 0));
	vhandle[4] = mesh.add_vertex(MyMesh::Point(0, 4, 0));
	vhandle[5] = mesh.add_vertex(MyMesh::Point(1, 1, 0));
	vhandle[6] = mesh.add_vertex(MyMesh::Point(1, 3, 0));

	// Add three quad faces
	std::vector<MyMesh::VertexHandle> face_vhandles;

	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[3]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();

	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[6]);
	face_vhandles.push_back(vhandle[4]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();

	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[6]);
	mesh.add_face(face_vhandles);

	// Now find the edge between vertex vhandle[2]
	// and vhandle[3]
	vector<Pair> list;
	vector<OpenMesh::HalfedgeHandle> edge_list;
	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++ i)
	{
		Pair temp;
		temp.x = mesh.from_vertex_handle(i.handle()).idx();
		temp.y = mesh.to_vertex_handle(i.handle()).idx();
		edge_list.push_back(i.handle());
		list.push_back(temp);
	}
	
	OpenMesh::HalfedgeHandle temp;
	for(auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) 
	{
		if (mesh.from_vertex_handle(it.handle()) == vhandle[0] &&
			mesh.to_vertex_handle(it.handle()) == vhandle[2])
			temp = it.handle();

		if(mesh.to_vertex_handle(it.handle()) == vhandle[3] &&
			mesh.from_vertex_handle(it.handle()) == vhandle[2]) 
		{
			// Collapse edge
			mesh.request_edge_status();
			mesh.request_vertex_status();

			mesh.request_halfedge_status();

			auto info = OpenMesh::Decimater::CollapseInfoT<MyMesh>::CollapseInfoT(mesh, it.handle());

			if (mesh.is_collapse_ok(it.handle()))
			{
				mesh.collapse(it.handle());
				//mesh.garbage_collection();
				//mesh.point(info.v1).data()[1] = 1;
			}
			
			for (int i = 0; i < list.size(); ++ i)
			{
				//(list[i].x == 3 && list[i].y == 2) || (list[i].x == 2 && list[i].y == 3)||
				/*if (
					(list[i].x == 5 && list[i].y == 6) || (list[i].x == 6 && list[i].y == 5))
					continue;*/
				cout << list[i].x << "->" << list[i].y << endl;
				cout << mesh.from_vertex_handle(edge_list[i]) << "->" << mesh.to_vertex_handle(edge_list[i]) << endl;
				cout << mesh.point(mesh.from_vertex_handle(edge_list[i])).data()[0] << ' '
					<< mesh.point(mesh.from_vertex_handle(edge_list[i])).data()[1] << ' '
					<< mesh.point(mesh.to_vertex_handle(edge_list[i])).data()[0] << ' '
					<< mesh.point(mesh.to_vertex_handle(edge_list[i])).data()[1] << endl;
				cout << endl;
			}

			/*int kk = 0;
			for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++ i)
			{
				++ kk;
				cout << mesh.from_vertex_handle(i.handle()) << "->" << mesh.to_vertex_handle(i.handle()) << endl;
			}
			cout << "kk:" << kk << endl;*/
			break;
		}
	}
	output(mesh);
	system("pause");
	return 0;
}