#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include <Eigen/SparseCore>
#include "Eigen/SparseQR"
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
using namespace std;
int main()
{
	Eigen::SparseMatrix < double > A(2065, 2065);

	for (int i = 0; i != 2065; ++i)
	{
		A.insert(i, i) = 1.0;
	}
	ofstream out("D:/eigen.txt");
	out << A << endl;
	system("pause");
	return 0;
}