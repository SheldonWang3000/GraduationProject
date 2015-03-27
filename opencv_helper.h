#ifndef OPENCV_HELPER_H_
#define OPENCV_HELPER_H_
#include <opencv2/opencv.hpp>
using namespace std;
class OpencvHelper
{
public:
	void pic_edge(IplImage *image, IplImage *result); 
private:
	double get_T(IplImage *image);
	double get_delta(vector<double> list, double avg);
	double get_avg(vector<double> list);
	void CheckAll(IplImage *image, double T);
	vector<double> Th;
	vector<double> Tl;
};
#endif //OPENCV_HELPER_H_