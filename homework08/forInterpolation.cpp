// openCV Library 를 이용하여 이미지를 읽고, bilinearInterpolation 하여라.

//file I/O and resize with interpolation
#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>

using namespace std;
using namespace cv; // cv:: 생략하

void BilinearInterpolation(Mat &src, Mat &dst);
Mat bi_dst;

int main()
{
	cv::Mat src;
	string file_path;

    file_path = "/Users/junbeom/Desktop/1.jpg";
    // 원본 image 읽어오기
	src = imread(file_path, 1); // 컬러로 읽는다.

	if (src.empty()) // file path에서 image를 못찾았다면,
	{
		cout << "Cannot find an image" << std::endl;
		return -1;
	}
    // 원본 해상도 출력
    int height = src.rows, width = src.cols;
    cout << "src 해상도 : " << height << " X " << width << endl;

    // 원본 출력
	imshow("SourceImage", src);
	waitKey(0);

	// 비율 입력
	int new_height, new_width;
	cin >> new_height >> new_width;

    bi_dst = Mat(new_height, new_width, src.type(), cv::Scalar(0));

	// bilinearInterpolation
	BilinearInterpolation(src, bi_dst);


    cout << "bi-Image 해상도 : " << new_height << " X " << new_width << endl;
	// bilinearInterplolation 된 image 출력
	imshow("BIImage", bi_dst);
	waitKey(0);

	return 0;
}

void BilinearInterpolation(Mat &src, Mat &dst)
{
	///////////////////////////////////////////////
	//bilinear interpolation
    for(int i=0; i<dst.rows; i++) {
        for(int j=0; j<dst.cols; j++) {
            int px = (int)(j * src.cols / dst.cols); // 대응되는 원본 Point의 왼쪽 위 격자점 x-좌표
            int py = (int)(i * src.rows / dst.rows); // 대응되는 원본 Point의 왼쪽 위 격자점 y-좌표

            // x, y 격자 내 간격
            double fx1 = (double)(j * src.cols) / (double)dst.cols - (double)px;
            double fy1 = (double)(i * src.rows) / (double)dst.rows - (double)py;
            double fx2 = 1.0 - fx1;
            double fy2 = 1.0 - fy1;

            // fx1,fx2,fy1,fy2로 만들어진 4개 구역의 넓이
            double w1 = fx1 * fy1;
            double w2 = fx2 * fy1;
            double w3 = fx1 * fy2;
            double w4 = fx2 * fy2;

            // 격자 점(src의 pixel)의 가중치(RGB)를 이용한 bilinear interpolation
            Vec3b p1 = src.at<Vec3b>(py, px);
            Vec3b p2 = src.at<Vec3b>(py, px + 1);
            Vec3b p3 = src.at<Vec3b>(py + 1, px);
            Vec3b p4 = src.at<Vec3b>(py + 1, px + 1);

            dst.at<Vec3b>(i, j)[0] = (uchar)(w1 * (double)p4[0] + w2 * (double)p3[0] + w3 * (double)p2[0] + w4 * (double)p1[0]);
            dst.at<Vec3b>(i, j)[1] = (uchar)(w1 * (double)p4[1] + w2 * (double)p3[1] + w3 * (double)p2[1] + w4 * (double)p1[1]);
            dst.at<Vec3b>(i, j)[2] = (uchar)(w1 * (double)p4[2] + w2 * (double)p3[2] + w3 * (double)p2[2] + w4 * (double)p1[2]);
        }
    }
	///////////////////////////////////////////////
}