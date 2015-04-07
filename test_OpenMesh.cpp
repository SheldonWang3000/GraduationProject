#include <iostream>
#include <opencv2/opencv.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/CollapseInfoT.hh>
#include "opencv_helper.h"
using namespace std;
typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

void output(MyMesh mesh)
{
	try
	{
		if ( !OpenMesh::IO::write_mesh(mesh, "D:/test.off") )
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
	list.push_back(vhandle[3]);
	list.push_back(vhandle[1]);
	
	mesh.add_face(list);
	list.clear();
	list.push_back(vhandle[3]);
	list.push_back(vhandle[2]);
	list.push_back(vhandle[1]);
	
	mesh.add_face(list);

	for (auto i = mesh.faces_begin(); i != mesh.faces_end(); ++i)
	{
		auto k = mesh.calc_face_normal(i.handle());
		if (i.handle().idx() == 1)
		{
			vector<MyMesh::HalfedgeHandle> list;
			for (auto j = mesh.fh_begin(i.handle()); j != mesh.fh_end(i.handle()); ++j)
			{
				list.push_back(j.handle());
			}
			int t = 2;
			for (auto j = mesh.fh_begin(i.handle()); j != mesh.fh_end(i.handle()); ++j)
			{
				//mesh.set_next_halfedge_handle(j.handle(), list[t--]);
			}
		}
	}

	output(mesh);

	system("pause");
	return 0;
}