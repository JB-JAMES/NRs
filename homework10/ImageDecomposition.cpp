// 주어진 이미지 2개 사이의 matching point 를 추출한다.
// 데이터의 함수관계를 없애기 위해 한쪽 이미지에는 Gaussian Noise를 준다.

#include <iostream>
#include <opencv2/core.hpp>
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/features2d.hpp"

using namespace cv;
using namespace std;

Mat AddGaussianNoise(Mat &src, double average, double std);
void MatchedKeyPoint(Mat &src1, Mat &src2);

//global variable
int distThresh = 30;
vector<Point2f> Left;
vector<Point2f> Right;

int main()
{
    Mat LeftImg = imread("../2pairs_for_HW10/goodLeft.jpg");
    Mat RightImg = imread("../2pairs_for_HW10/goodRight.jpg");

    //image
//	resize(LeftImg, LeftImg, Size(), 0.25, 0.25, INTER_LINEAR);
//	resize(RightImg, RightImg, Size(), 0.25, 0.25, INTER_LINEAR);

    //Add Gaussian noise
    double average, std;

    average = 0.0;
    std = 30.0;

    RightImg = AddGaussianNoise(RightImg, average, std);

    //Extract features and matching correspondence
    MatchedKeyPoint(LeftImg, RightImg);

    // 1. Establish the feature correspondence between the two images (Output: xi, yi, xi’, yi’, i=1~N)
    cout << "N = " << Left.size() << endl;
//    for(int i=0; i<Left.size(); i++)
//        cout << i+1 << " : <" << Left[i].x << "," << Left[i].y << " >\t< " << Right[i].x << "," << Right[i].y << " >" << endl;
    for(int i=0; i<Left.size(); i++)
        cout << Left[i].x << " " << Left[i].y << " " << Right[i].x << " " << Right[i].y << endl;
    // 2. Find the parameters using the correspondence data

    return 0;
}

Mat AddGaussianNoise(Mat &src, double average, double std)
{
    Mat noise_src(src.size(), src.type());

    randn(noise_src, Scalar::all(average), Scalar::all(std));

    Mat tmp_img;
    src.convertTo(tmp_img, src.type());
    addWeighted(tmp_img, 1.0, noise_src, 1.0, 0.0, tmp_img);

    imshow("tmp", tmp_img);
    waitKey(0);

    return tmp_img;
}

void MatchedKeyPoint(Mat &src1, Mat &src2)
{
    Ptr<Feature2D> featureExtractor;
    featureExtractor = ORB::create();
    vector<KeyPoint> LeftKeypoints, RightKeypoints;
    Mat LeftDescriptors, RightDescriptors;
    Mat matchingImage;
    featureExtractor->detectAndCompute(src1, Mat(), LeftKeypoints, LeftDescriptors);
    featureExtractor->detectAndCompute(src2, Mat(), RightKeypoints, RightDescriptors);

    Ptr<DescriptorMatcher> matcher;
    vector<vector<DMatch>> knnMatches;
    vector<DMatch> matches;

    matcher = DescriptorMatcher::create("BruteForce-Hamming");
    matcher->match(LeftDescriptors, RightDescriptors, matches);
    vector<DMatch>::iterator it = matches.begin();
    for (; it != matches.end(); )
    {
        if (it->distance > distThresh)
        {
            it = matches.erase(it);
        }
        else
        {
            it++;
        }
    }

    drawMatches(src1, LeftKeypoints, src2, RightKeypoints, matches, matchingImage);

    imshow("matches", matchingImage);
    waitKey(0);

    for (int i = 0; i < matches.size(); i++)
    {
        //Left src1 point, Right src2 point
        Left.push_back(LeftKeypoints[matches[i].queryIdx].pt);
        Right.push_back(RightKeypoints[matches[i].trainIdx].pt);
    }

}