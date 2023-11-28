#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<vector>
#include<opencv2/opencv.hpp>
#include<fstream>
#include"image_point.h"
#include"camera.h"
using namespace std;
using namespace cv;
/// <summary>
/// ÄÚ¶¨Ïò
/// </summary>
class InteriorOrientation
{
public:
	void get_param_from_file(const char* file_path);
	void calculate_A_matrix(Mat_<double>& A);
	void calculate_L_matrix(Mat_<double>& L);
	void affine_interior_orientation(const char* file_path, const char* result_file_path);
public:
	Camera camera_;
	double pixel_size_;
	double width_, height_;
	vector<ImagePoint> frame_coordinate_points_;
	vector<ImagePoint> pixel_coordinate_points_;
};

