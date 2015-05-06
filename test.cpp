#include <cstdio>
#include <iostream>
#include "Eigen/Dense"
#include <time.h>
using namespace std;
using namespace Eigen;
int main()
{
	Eigen::Vector3f A(0, 5, 0);
	Eigen::Vector3f B(0, 3, 0);
	Eigen::Vector3f C(4, 0, 0);
	Eigen::Vector3f P(1, 0, 0);
	cout << A.cross(B) << endl;
	system("pause");
	return 0;
}