#include <iostream>
#include <Eigen/Dense>
using Eigen::MatrixXd;
int main()
{
	MatrixXd m(1, 2);
	m(0, 0) = 3;
	m(0, 1) = 2.5;
	std::cout << m << std::endl;
	std::cout << std::endl;
	std::cout << m.transpose() << std::endl;
	system("pause");
	return 0;
}