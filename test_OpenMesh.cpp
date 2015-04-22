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
enum PointType
{
	None, Crease, Tear, Corner, Smooth
};
struct MyPoint
{
	MyMesh::VertexHandle point;
	CvScalar color;
	int area_idx;
	PointType point_type;
	bool isEdge;
};
void output(MyMesh mesh)
{
	try
	{
		if ( !OpenMesh::IO::write_mesh(mesh, "D:/test.off") )
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
int main()
{
	MyMesh mesh;
	// Add some vertices as in the illustration above
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
	for (auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) {

		if (mesh.to_vertex_handle(it.handle()) == vhandle[3] &&
			mesh.from_vertex_handle(it.handle()) == vhandle[2]) {
			cout << it.handle().idx() << endl;
			// Collapse edge
			mesh.request_halfedge_status();
			mesh.request_vertex_status();
			mesh.request_face_status();
			mesh.request_edge_status();
			if (mesh.is_collapse_ok(it.handle()))
				mesh.collapse(it.handle());
			break;
		}
	}

	for (auto i = mesh.halfedges_begin(); i != mesh.halfedges_end(); ++i)
	{
		cout << i.handle().idx() << ":" << mesh.status(i).deleted() << endl;
		if (i.handle().idx() == 4)
		{
			cout << "----------" << mesh.from_vertex_handle(i) << "->" << mesh.to_vertex_handle(i) << endl;
		}
	}
	// Our mesh now looks like in the illustration above after the collapsing.
	output(mesh);

	system("pause");
	return 0;
}