#ifndef OPENMESH_HELPER_H
#define OPENMESH_HELPER_H
#include "opencv_helper.h"
#include <opencv2/opencv.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Decimater/CollapseInfoT.hh>
using namespace std;
typedef OpenMesh::PolyMesh_ArrayKernelT<> MyMesh;

class OpenmeshHelper
{
private:
	enum PointType
	{
		None, Crease, Tear, Corner, Smooth
	};
	struct Pair
	{
		OpenMesh::HalfedgeHandle edge;
		int value;
		friend bool operator==(const Pair &l, const Pair &r)
		{
			return l.edge.idx() == r.edge.idx();
		}
	};
	
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
	static bool Compare(Pair a, Pair b);
	vector<Pair> point_pair_list;
	vector<Pair> all_pair_list;
	vector<Pair> second_pair_list;
	vector<MyMesh::VertexHandle> control_list;
	vector<MyMesh::HalfedgeHandle> tear_list;
	vector<MyMesh::HalfedgeHandle> crease_list;
	map<Position, MyMesh::VertexHandle> control_map;
	multimap<Position, Position> tear_map;
	multimap<Position, Position> crease_map;
	MyMesh mesh;
	IplImage *image;
	MyMesh::VertexHandle **point_data;
	bool exact = false;
	double degree = 80;
	double lambda = 0.5;
	double num_vertices;
	double num_all_vertices;

	void InitPointData();
	void ReSampleColor();
	void ConnectMesh(bool isContour);
	void CountVertices();
	void InitPairList();
	void InsertCreaseEdge();
	void OptimizeColor();
	void SortVertices();
	void LoopReduce(double rate, bool visual);
	bool IsCollapseOK(MyMesh::HalfedgeHandle half);
	void RegulatePosition(OpenMesh::Decimater::CollapseInfoT<MyMesh> info, double *x, double *y);
	void CollapseEdge(MyMesh::HalfedgeHandle half);
	void OptimizePosition();
	bool IsBoundary(MyMesh::VertexHandle);
public:
	void Output(string output_location);
	void ReduceVertices(double rate, bool visual);

	OpenmeshHelper(string input_location);
	~OpenmeshHelper();
};
#endif // OPENMESH_HELPER_H
