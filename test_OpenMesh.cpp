#include <iostream>
#include <opencv2/opencv.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/CollapseInfoT.hh>
#include "Eigen/Dense"
#include "opencv_helper.h"
using namespace std;
using namespace Eigen;
typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;
void output(MyMesh mesh)
{
	try
	{
		if ( !OpenMesh::IO::write_mesh(mesh, "D:/test22.off") )
		{
			std::cout << "Cannot write mesh to file 'output.off'" << std::endl;
			return;
		}
	}
	catch( std::exception& x )
	{
		std::cerr << x.what() << std::endl;
		return;
	}
}
struct Position 
{
	double x, y;
	friend bool operator==(const Position &l, const Position &r)
	{
		return l.x == r.x && l.y == r.y;
	}
	friend bool operator<(const Position &l, const Position &r)
	{
		if (l.x != r.x)
		{
			return l.x < r.x;
		}
		return l.y < r.y;
	}
};
int main()
{
	MyMesh mesh, temp_mesh;
	//OpenMesh::IO::read_mesh(mesh, "D:/before.off");
	//map<Position, MyMesh::VertexHandle> m;
	//for (auto i = mesh.faces_begin(); i != mesh.faces_end(); ++i)
	//{
	//	vector<MyMesh::VertexHandle> list;
	//	for (auto j = mesh.fv_begin(i); j != mesh.fv_end(i); ++j)
	//	{
	//		Position temp;
	//		temp.x = mesh.point(j).data()[0];
	//		temp.y = mesh.point(j).data()[1];
	//		if (m[temp].idx() == -1)
	//		{
	//			m[temp] = temp_mesh.add_vertex(MyMesh::Point(temp.x, temp.y, 0));
	//		}
	//		list.push_back(m[temp]);
	//	}
	//	temp_mesh.add_face(list);
	//	
	//}
	//output(temp_mesh);
	// Add some vertices as in the illustration above
	MyMesh::VertexHandle vhandle[7];

	std::vector<MyMesh::VertexHandle> face_vhandles;
	vhandle[0] = mesh.add_vertex(MyMesh::Point(-1, 1, 0));
	vhandle[1] = mesh.add_vertex(MyMesh::Point(-1, 3, 0));
	vhandle[2] = mesh.add_vertex(MyMesh::Point(0, 0, 0));
	vhandle[3] = mesh.add_vertex(MyMesh::Point(0, 2, 0));
	vhandle[4] = mesh.add_vertex(MyMesh::Point(0, 4, 0));
	vhandle[5] = mesh.add_vertex(MyMesh::Point(1, 1, 0));
	vhandle[6] = mesh.add_vertex(MyMesh::Point(1, 3, 0));

	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[3]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();

	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[6]);
	mesh.add_face(face_vhandles);

	face_vhandles.clear();

	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[3]);
	mesh.add_face(face_vhandles);
	face_vhandles.clear();

	face_vhandles.push_back(vhandle[0]);
	face_vhandles.push_back(vhandle[2]);
	face_vhandles.push_back(vhandle[3]);
	mesh.add_face(face_vhandles);
	face_vhandles.clear();

	face_vhandles.push_back(vhandle[3]);
	face_vhandles.push_back(vhandle[5]);
	face_vhandles.push_back(vhandle[6]);
	mesh.add_face(face_vhandles);
	face_vhandles.clear();

	face_vhandles.push_back(vhandle[1]);
	face_vhandles.push_back(vhandle[6]);
	face_vhandles.push_back(vhandle[4]);
	mesh.add_face(face_vhandles);
	face_vhandles.clear();
	output(mesh);
	ofstream out("D:/test.off");
	out << "OFF" << endl;
	out << mesh.n_vertices() << ' ' << mesh.n_faces() << " 0" << endl;
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		out << mesh.point(i).data()[0] << " " << mesh.point(i).data()[1] << " 0" << endl;
	}
	for (auto i = mesh.faces_begin(); i != mesh.faces_end(); ++i)
	{
		out << "3 ";
		for (auto j = mesh.fv_begin(i); j != mesh.fv_end(i); ++j)
		{
			out << j.handle().idx() << ' ';
		}
		out << endl;
	}
	out.close();
	//mesh.triangulate();
	// Now find the edge between vertex vhandle[2]
	// and vhandle[3]
	/*for (auto i = mesh.faces_begin(); i != mesh.faces_end(); ++i)
	{
		cout << i.handle().idx() << endl;
	}
	cout << "-------------" << endl;
	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
	{
		cout << mesh.from_vertex_handle(i).idx() << "->" << mesh.to_vertex_handle(i) << endl;
	}
	cout << "--------------" << endl;*/
	//mesh.request_edge_status();
	//mesh.request_face_status();
	//mesh.request_halfedge_status();
	//mesh.request_vertex_status();
	//for (MyMesh::HalfedgeIter it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) 
	//{
	//	if (mesh.to_vertex_handle(it.handle()) == vhandle[3] &&
	//		mesh.from_vertex_handle(it.handle()) == vhandle[2]) 
	//	{
	//		// Collapse edge
	//		mesh.collapse(it);
	//		break;
	//	}
	//}
	//mesh.garbage_collection();
	//output(mesh);
	//for (auto i = mesh.faces_begin(); i != mesh.faces_end(); ++i)
	//{
	//	cout << i.handle().idx() << endl;
	//	cout << "++++" << endl;
	//	for (auto j = mesh.fv_begin(i); j != mesh.fv_end(i); ++j)
	//	{
	//		cout << j.handle().idx() << ' ';
	//	}
	//	cout << endl;
	//}
	//mesh.release_edge_status();
	//mesh.release_halfedge_status();
	//mesh.release_face_status();
	//mesh.release_vertex_status();
	system("pause");
	return 0;
}