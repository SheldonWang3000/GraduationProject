#include "opencv_helper.h"
#include <fstream>
using namespace std;
double OpencvHelper::get_avg(vector<double> list)
{
	double sum = 0;
	for (int i = 0; i != list.size(); ++ i)
	{
		sum += list[i];	
	}
	return sum / list.size();
}
double OpencvHelper::get_delta(vector<double> list, double avg)
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
void OpencvHelper::CheckAll(IplImage *image, double T)
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
double OpencvHelper::get_T(IplImage *image)
{
	double T = 200;
	while (T >= 150)	
	{
		CheckAll(image, T);
		double uh = get_avg(Th);
		double ul = get_avg(Tl);
		T = (uh + ul) / 2;
	}
	return T;
}
void OpencvHelper::get_pic_edge(IplImage *image, IplImage *result)
{
	IplImage *grey, *edge;
	grey = cvCreateImage(cvGetSize(image), image->depth, 1);
	cvCvtColor(image, grey, CV_BGR2GRAY);
	edge = cvCreateImage(cvGetSize(image), image->depth, 1);
	double T = get_T(grey);
	CheckAll(grey, T);
	double ul = get_avg(Tl);
	double dl = get_delta(Tl, ul);
	dl = sqrt(dl);
	cvCanny(grey, edge, ul - 0.3 * dl, ul + dl, 3);

	cvCopy(edge, result, NULL);
	//Output(result, "D:/edge");
	//cvNamedWindow("edge", CV_WINDOW_AUTOSIZE);
	//cvShowImage("edge",edge);
	//cvWaitKey(0);
	cvReleaseImage(&grey);
	cvReleaseImage(&edge);
	//cvDestroyWindow("edge");
	//
}

void OpencvHelper::get_pic_area(IplImage *image, IplImage *result)
{
	IplImage *area = cvCreateImage(cvGetSize(image), image->depth, 1);
	area = cvCloneImage(image);
	cvPyrMeanShiftFiltering(image, area, 30, 20, 3);

	// 转为灰度图  
	IplImage *pGrayImage = cvCreateImage(cvGetSize(area), IPL_DEPTH_8U, 1);
	cvCvtColor(area, pGrayImage, CV_BGR2GRAY);
	// 转为二值图  
	IplImage *pBinaryImage = cvCreateImage(cvGetSize(pGrayImage), IPL_DEPTH_8U, 1);
	cvThreshold(pGrayImage, pBinaryImage, 130, 255, CV_THRESH_BINARY);

	// 检索轮廓并返回检测到的轮廓的个数  
	CvMemStorage *pcvMStorage = cvCreateMemStorage();
	CvSeq *pcvSeq = NULL;
	cvFindContours(pBinaryImage, pcvMStorage, &pcvSeq, sizeof(CvContour), CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE, cvPoint(0, 0));

	// 画轮廓图  
	int nLevels = 5;
	// 填充成白色  
	cvRectangle(result, cvPoint(0, 0), cvPoint(result->width, result->height), CV_RGB(255, 255, 255), CV_FILLED);
	cvDrawContours(result, pcvSeq, CV_RGB(0, 0, 0), CV_RGB(0, 0, 0), 5, 2);

	cv::RNG rng = cv::theRNG();
	cv::Mat tempMat(result, 0);
	cv::Mat mask(tempMat.rows + 2, tempMat.cols + 2, CV_8UC1, cv::Scalar::all(0));
	auto tt = cv::Scalar::all(2);
	for (int i = 0; i < tempMat.rows; ++i)    //opencv图像等矩阵也是基于0索引的
		for (int j = 0; j < tempMat.cols; ++j)
			if (mask.at<uchar>(i + 1, j + 1) == 0)
			{
				cv::Scalar newcolor(rng(256), rng(256), rng(256));
				floodFill(tempMat, mask, cv::Point(j, i), newcolor, 0, tt, tt);
			}
}
void OpencvHelper::get_pic_contours(IplImage *image, IplImage *result)
{
	IplImage *area = cvCreateImage(cvGetSize(image), image->depth, 1);
	area = cvCloneImage(image);
	cvPyrMeanShiftFiltering(image, area, 30, 20, 3);

	// 转为灰度图  
	IplImage *pGrayImage = cvCreateImage(cvGetSize(area), IPL_DEPTH_8U, 1);
	cvCvtColor(area, pGrayImage, CV_BGR2GRAY);
	// 转为二值图  
	IplImage *pBinaryImage = cvCreateImage(cvGetSize(pGrayImage), IPL_DEPTH_8U, 1);
	cvThreshold(pGrayImage, pBinaryImage, 130, 255, CV_THRESH_BINARY);

	// 检索轮廓并返回检测到的轮廓的个数  
	CvMemStorage *pcvMStorage = cvCreateMemStorage();
	CvSeq *pcvSeq = NULL;
	cvFindContours(pBinaryImage, pcvMStorage, &pcvSeq, sizeof(CvContour), CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE, cvPoint(0, 0));

	// 画轮廓图  
	int nLevels = 5;
	// 填充成白色  
	cvRectangle(result, cvPoint(0, 0), cvPoint(result->width, result->height), CV_RGB(255, 255, 255), CV_FILLED);
	cvDrawContours(result, pcvSeq, CV_RGB(0, 0, 0), CV_RGB(0, 0, 0), 5, 2);
	//Output(result, "D:/contours");
}

void OpencvHelper::Output(IplImage *image, string position)
{
	string str = position + ".jpg";
	cvSaveImage(str.c_str(), image);
	uchar *data = (uchar*)image->imageData;
	int step = image->widthStep / sizeof(uchar);
	ofstream out(position + ".txt");
	for (int i = 0; i < image->height; ++i)
	{
		for (int j = 0; j < image->width; ++j)
		{
			int t = data[i * step + j];
			out << t << ' ';
		}
		out << endl;
	}
	out.close();
}
