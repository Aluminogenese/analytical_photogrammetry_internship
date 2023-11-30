#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<vector>
#include<opencv2/opencv.hpp>
#include<fstream>
#include"plane_coordinates.h"
#include"camera.h"
using namespace std;
using namespace cv;
/// <summary>
/// 内定向
/// </summary>
class InteriorOrientation
{
public:
	void get_param_from_file(const char* file_path);
	void calculate_A_matrix(Mat_<double>& A);
	void calculate_L_matrix(Mat_<double>& L);
	void affine_interior_orientation(const char* file_path, const char* result_file_path);
public:
	Camera camera_;// 相机内参
	double pixel_size_;// 像素大小
	double width_, height_;// 宽度（列数），高度（行数）
	vector<PlaneCoordinates> frame_coordinate_points_;//框标坐标
	vector<PlaneCoordinates> pixel_coordinate_points_;//像素坐标
};

