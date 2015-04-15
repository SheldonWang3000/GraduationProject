#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
int main()
{
	MatrixXf A1 = MatrixXf::Random(3, 2);
	std::cout << "Here is the matrix A:\n" << A1 << std::endl;
	VectorXf b1 = VectorXf::Random(3);
	std::cout << "Here is the right hand side b:\n" << b1 << std::endl;
	//jacobiSvd ��ʽ:Slow (but fast for small matrices)
	std::cout << "The least-squares solution is:\n"
		<< A1.jacobiSvd(ComputeThinU | ComputeThinV).solve(b1) << std::endl;
	//colPivHouseholderQr����:fast
	std::cout << "The least-squares solution is:\n"
		<< A1.colPivHouseholderQr().solve(b1) << std::endl;
	system("pause");
}