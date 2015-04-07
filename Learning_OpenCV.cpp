#include <iostream>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "opencv_helper.h"
using namespace std;
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")
IplImage *edge;
IplImage *grey;
vector<double> Th, Tl;
void rgb_to_file(uchar *data, int height, int width, int step)
{
	//int step = image->widthStep / sizeof(uchar);
	ofstream outfile("D:/edge.txt");
	int temp;
	for (int i = 0; i != height; ++ i)
	{
		for (int j = 0; j != width; ++ j)
		{	
			//CvScalar x = cvGet2D(image, i, j);
			temp = data[i * step + j];	
			outfile << temp << ' ';
		}
		outfile << endl;
	}
	outfile.close();
}
void rgb_to_file(IplImage *image)
{
	uchar *data = (uchar*) image->imageData;
	int step = image->widthStep / sizeof(uchar);
	ofstream outfile("D:/edge.txt");
	int temp;
	for (int i = 0; i != image->height; ++ i)
	{
		for (int j = 0; j != image->width; ++ j)
		{	
			//CvScalar x = cvGet2D(image, i, j);
			temp = data[i * step + j];	
			outfile << temp << ' ';
		}
		outfile << endl;
	}
	outfile.close();
}
double avgT(vector<double> list)
{
	double sum = 0;
	for (int i = 0; i != list.size(); ++ i)
	{
		sum += list[i];	
	}
	return sum / list.size();
}

void check(IplImage *image, double T)
{
	Th.clear();
	Tl.clear();
		
	double temp = 0;
	uchar *data = (uchar *)image->imageData;	
	int step = image->widthStep / sizeof(uchar);
	
	for (int i = 0; i != image->height; ++ i)
	{
		for (int j = 0; j != image->width; ++ j)
		{	
			temp = data[i * step + j];
			if (temp >= T) Th.push_back(temp);
			else Tl.push_back(temp);
		}
	}

}
double findT(IplImage *image)
{
	double T = 200;

	while (T >= 150)	
	{
		check(image, T);
		
		double uh = avgT(Th);
		double ul = avgT(Tl);
		
		T = (uh + ul) / 2;
	}

	return T;
}

double delta(vector<double> list, double avg)
{
	double temp;
	double sum = 0;
	for (int i = 0; i != list.size(); ++ i)
	{
		temp = list[i] - avg;
		temp *= temp;
		sum += temp;
	}
	return sum / list.size();
}
void grey2edge(IplImage *grey)
{
	
	edge = cvCreateImage(cvGetSize(grey), grey->depth, 1);
	cvNamedWindow("edge", CV_WINDOW_AUTOSIZE);

	double T = findT(grey);
	
	check(grey, T);
	
	double ul = avgT(Tl);
	double dl = delta(Tl, ul);

	cvCanny(grey, edge, ul - 0.3 * sqrt(dl), ul + sqrt(dl), 3);
	cvShowImage("edge", edge);

	rgb_to_file(edge);	
	cvWaitKey(0);
	cvReleaseImage(&edge);
	cvDestroyWindow("edge");
}
void BGR2Grey(IplImage *image)
{
	grey = cvCreateImage(cvGetSize(image), image->depth, 1);
	cvCvtColor(image, grey, CV_BGR2GRAY);

	cvNamedWindow("original", CV_WINDOW_AUTOSIZE);
	cvShowImage("original", image);

	cvNamedWindow("grey", CV_WINDOW_AUTOSIZE);
	cvShowImage("grey", grey);

	//rgb_to_file(image);
	//cvSaveImage("D:\\saveImage.jpg", grey);//���Ҷ�ͼ���浽����
	grey2edge(grey);

	cvWaitKey(0);
	cvReleaseImage(&image);
	cvReleaseImage(&grey);
	cvDestroyWindow("origial");
	cvDestroyWindow("gray");
}
int main()
{
	IplImage *image;
	image = cvLoadImage("D:/lena.jpg", -1);
	OpencvHelper helper;
	//IplImage *result = cvCreateImage(cvGetSize(image), image->depth, 1);
	//helper.pic_edge(image, result);
	//BGR2Grey(image);

	//rgb_to_file(result);
	IplImage *area;
	area = cvCreateImage(cvGetSize(image), image->depth, 1);
	area = cvCloneImage(image);
	cvPyrMeanShiftFiltering(image, area, 30, 20, 3);
	//2,40,3
	//rgb_to_file(area);

	cv::RNG rng = cv::theRNG();
	cv::Mat tempMat(area, 0);
	cv::Mat mask(tempMat.rows + 2, tempMat.cols + 2, CV_8UC1, cv::Scalar::all(0));
	auto tt = cv::Scalar::all(2);
	for(int i = 0; i < tempMat.rows; ++ i)    //opencvͼ��Ⱦ���Ҳ�ǻ���0������
		for(int j = 0; j < tempMat.cols; ++ j)
			if(mask.at<uchar>(i + 1, j + 1) == 0)
			{
				cv::Scalar newcolor(rng(256), rng(256), rng(256));
				floodFill(tempMat, mask, cv::Point(j, i), newcolor, 0, tt, tt);
			}
	
	

	uchar *data = (uchar*) image->imageData;
	map<int, int> color_set;
	int step = image->widthStep / sizeof(uchar);
	int temp = 0;
	int num = 0;
	for (int i = 0; i != area->height; ++ i)
	{
		for (int j = 0; j != area->width; ++ j)
		{	
			CvScalar color = cvGet2D(area, i, j);
			//temp = data[i * step + j];	
			int sum = 0;
			sum = color.val[0];
			sum *= 1000;
			sum += color.val[1];
			sum *= 1000;
			sum += color.val[2];
			
			if (color_set.count(sum) != 0)
				color_set[sum] += 1;
			else
				color_set.insert(std::map<int, int>::value_type(sum, 1));
		}
		
	}
	for (auto i = color_set.begin(); i != color_set.end(); ++ i)
	{
		if (i->second > 100)
		{
			++ num;
		}
	}
	
	cout << num << endl;
	
	//fstream out("D:/color.txt");
	//for (auto i = color_set.begin(); i != color_set.end(); ++ i)
	//{
	//	out << i->first << ' ' << i->second << endl;
	//}
	//out.close();

	cvNamedWindow("area", CV_WINDOW_AUTOSIZE);
	cvShowImage("area", area);
	cvWaitKey(0);
	cvDestroyWindow("area");
	cvDestroyWindow("result");
	cvReleaseImage(&area);
	cvReleaseImage(&image);
	return 0;
}
