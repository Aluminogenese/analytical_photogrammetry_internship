#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<opencv2/opencv.hpp>
#include<vector>
#include<fstream>
#include"control_point.h"
using namespace cv;
using namespace std;
/// <summary>
/// ¾ø¶Ô¶¨Ïò
/// </summary>
class AbsoluteOrientation
{
public:
	void get_param_from_file(const char* image_path);
	void calculate_absolute_rotation_matrix(Mat_<double>& R);
	void calculate_barycentric_coordinate(vector<ControlPoint>& barycentric_coordinate, vector<ControlPoint> original_coordinate, ControlPoint& gravity_center);
	void calculate_A_matrix(int i, Mat_<double>& A, Mat_<double> model_barycentric_coordinate_rotated);
	void calculate_L_matrix(int i, Mat_<double>& A, Mat_<double> control_barycentric_coordinate_matrix, Mat_<double> model_barycentric_coordinate_rotated);
	void correct_absolute_param(Mat_<double> X);
	bool is_tolerant(Mat_<double> X, double tolerance);
	void calculate_absolute_orientation(const char* image_path, const char* result_file_path);
public:
	vector<ControlPoint> model_point_;
	vector<ControlPoint> control_point_;
	Mat_<double> absolute_orientation_elements_ = Mat::zeros(7, 1, CV_32F); //x, y, z, ¦Ë, ¦µ, ¦¸, ¦Ê
};

