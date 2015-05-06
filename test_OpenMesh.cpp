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
	MyMesh::VertexHandle vhandle[5];

	vhandle[0] = mesh.add_vertex(MyMesh::Point(-1, 1, 0));
	vhandle[1] = mesh.add_vertex(MyMesh::Point(-1, -1, 0));
	vhandle[2] = mesh.add_vertex(MyMesh::Point(1, -1, 0));
	vhandle[3] = mesh.add_vertex(MyMesh::Point(1, 1, 0));
	vhandle[4] = mesh.add_vertex(MyMesh::Point(0, 0, 0));

	
	

	// Add three quad faces
	std::vector<MyMesh::VertexHandle> list;

	list.push_back(vhandle[0]);
	list.push_back(vhandle[1]);
	list.push_back(vhandle[4]);
	mesh.add_face(list);

	//list.clear();
	//list.push_back(vhandle[1]);
	//list.push_back(vhandle[2]);
	//list.push_back(vhandle[4]);
	//mesh.add_face(list);

	//list.clear();
	//list.push_back(vhandle[2]);
	//list.push_back(vhandle[3]);
	//list.push_back(vhandle[4]);
	//mesh.add_face(list);

	//list.clear();
	//list.push_back(vhandle[3]);
	//list.push_back(vhandle[0]);
	//list.push_back(vhandle[4]);
	//mesh.add_face(list);

	/*set<MyMesh::VertexHandle> set;
	for (auto i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i)
	{
		for (auto j = mesh.vv_begin(i); j != mesh.vv_end(i); ++j)
		{
			set.insert(j);
		}
	}
	
	cout << set.size() << endl;*/
	for (auto i = mesh.faces_begin(); i != mesh.faces_end(); ++i)
	{
		for (auto j = mesh.fv_begin(i); j != mesh.fv_end(i); ++j)
		{
			cout << j.handle().idx() << endl;
		}
		for (auto j = mesh.fh_begin(i); j != mesh.fh_end(i); ++j)
		{
			cout << mesh.from_vertex_handle(j).idx() << "-->" << mesh.to_vertex_handle(j).idx() << endl;
		}
	}
	
	output(mesh);

	map<Position, MyMesh::VertexHandle> map;
	Position p;
	p.x = 1;
	p.y = 1;

	cout << map[p] << endl;
	map[p] = mesh.vertices_begin();
	cout << map[p] << endl;

	p.x = 2;
	p.y = 3;
	cout << map[p] << endl;

	system("pause");
	return 0;
}