#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<vector>
#include<opencv2/opencv.hpp>
#include<fstream>
#include"camera.h"
#include"plane_coordinates.h"

using namespace std;
using namespace cv;
/// <summary>
/// 相对定向
/// </summary>
class RelativeOrientation
{
public:
	void get_param_from_file(const char* image_path, Camera& interior_elements, vector<PlaneCoordinates>& image_points, double& m);
	void calculate_relative_rotation_matrix(Mat_<double>& R);
	Mat calculate_auxiliary_coordinate(int i, Mat_<double>R, Camera interior_elements, vector<PlaneCoordinates> image_points);
	void calculate_A_matrix(int i, double Bx, double N2, Mat_<double>& A, Mat_<double> right_auxiliary_coordinates);
	void calculate_L_matrix(int i, double N1, double N2, double By, Mat_<double> left_auxiliary_coordinates, Mat_<double> right_auxiliary_coordinates, Mat_<double>& L);
	void correct_relative_param(Mat_<double> X);
	bool is_tolerant(Mat_<double> X, double tolerance);
	void calculate_relative_orientation(const char* left_image_path, const char* right_image_path, const char* result_file_path);
public:
	Camera left_image_interior_elements_;// 内方位元素，x0, y0, f
	vector<PlaneCoordinates> left_image_points_;// 像点坐标
	double left_m_;// 比例尺分母

	Camera right_image_interior_elements_;
	vector<PlaneCoordinates> right_image_points_;
	double right_m_;

	Mat_<double> relative_orientation_elements_ = Mat::zeros(5, 1, CV_32F); // 相对定向元素φ, ω, κ, u, v 
};

