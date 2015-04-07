#include "opencv_helper.h"
#include "openmesh_helper.h"

using namespace std;

int main()
{
	OpenmeshHelper mesh_helper("D:/lena_s.jpg");

	mesh_helper.ReduceVertices(0.25, false);
	mesh_helper.Output("D:/output.off");
	system("pause");
	return 0;
}