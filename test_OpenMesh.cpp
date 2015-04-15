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
	MyMesh::VertexHandle vhandle[9];
	vhandle[0] = mesh.add_vertex(MyMesh::Point(0, 0, 0));
	vhandle[1] = mesh.add_vertex(MyMesh::Point(1, 0, 0));
	vhandle[2] = mesh.add_vertex(MyMesh::Point(1, 1, 0));
	vhandle[3] = mesh.add_vertex(MyMesh::Point(0, 1, 0));
	vhandle[4] = mesh.add_vertex(MyMesh::Point(-1, 1, 0));
	vhandle[5] = mesh.add_vertex(MyMesh::Point(-1, 0, 0));
	vhandle[6] = mesh.add_vertex(MyMesh::Point(-1, -1, 0));
	vhandle[7] = mesh.add_vertex(MyMesh::Point(0, -1, 0));
	vhandle[8] = mesh.add_vertex(MyMesh::Point(1, -1, 0));
	
	vector<MyMesh::VertexHandle> list;
	list.push_back(vhandle[0]);
	list.push_back(vhandle[3]);
	list.push_back(vhandle[4]);
	list.push_back(vhandle[5]);
	mesh.add_face(list);

	list.clear();
	list.push_back(vhandle[0]);
	list.push_back(vhandle[5]);
	list.push_back(vhandle[6]);
	list.push_back(vhandle[7]);
	mesh.add_face(list);

	list.clear();
	list.push_back(vhandle[0]);
	list.push_back(vhandle[7]);
	list.push_back(vhandle[8]);
	list.push_back(vhandle[1]);
	mesh.add_face(list);

	list.clear();
	list.push_back(vhandle[0]);
	list.push_back(vhandle[1]);
	list.push_back(vhandle[2]);
	list.push_back(vhandle[3]);
	mesh.add_face(list);
	mesh.triangulate();
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		if (i.handle().idx() == 0)
		{
			for (auto k = mesh.voh_begin(i.handle()); k != mesh.voh_end(i.handle()); ++k)
			{
				cout << mesh.to_vertex_handle(k.handle()).idx() << endl;
			}
		}
	}

	output(mesh);

	system("pause");
	return 0;
}