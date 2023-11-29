#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<vector>
#include<opencv2/opencv.hpp>
#include<fstream>
#include"image_for_relative_orientation.h"
using namespace std;
using namespace cv;
/// <summary>
/// 相对定向
/// </summary>
class RelativeOrientation
{
public:
	void get_param_from_file(const char* image_path, ImageForRelativeOrientation& image_for_relative_orientation);
	void calculate_relative_rotation_matrix(Mat_<double>& R);
	Mat calculate_auxiliary_coordinate(int i, Mat_<double>R, ImageForRelativeOrientation image_for_relative_orientation);
	void calculate_A_matrix(int i, double Bx, double N2, Mat_<double>& A, Mat_<double> right_auxiliary_coordinates);
	void calculate_L_matrix(int i, double N1, double N2, double By, Mat_<double> left_auxiliary_coordinates, Mat_<double> right_auxiliary_coordinates, Mat_<double>& L);
	void correct_relative_param(Mat_<double> X);
	bool is_tolerant(Mat_<double> X, double tolerance);
	void calculate_relative_orientation(const char* left_image_path, const char* right_image_path, const char* result_file_path);
public:
	ImageForRelativeOrientation left_image_;
	ImageForRelativeOrientation right_image_;
	Mat_<double> relative_orientation_elements_ = Mat::zeros(5, 1, CV_32F); //fai, omega, kappa, u, v 

};

