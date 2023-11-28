#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "image_for_intersection.h"
#include<opencv2/opencv.hpp>
#include"control_point.h"
#include<fstream>
using namespace cv;
/// <summary>
/// 空间前方交会
/// </summary>
class SpaceIntersection
{
public:
	void get_param_from_file(const char* image_path, ImageForIntersection& image_for_intersection);
	void calculate_rotation_matrix(Mat_<double>& R, ImageForIntersection image_for_intersection);
	Mat coordinate_change(int i, Mat_<double>R, ImageForIntersection image_for_intersection);
	ControlPoint N1N2_intersection(int i, ImageForIntersection left_image, double N1, Mat_<double> left_auxiliary_coordinates);
	void pointfactor_space_intersection(const char* left_image_path, const char* right_image_path, const char* reslut_file_path);
public:
	ImageForIntersection left_image_;
	ImageForIntersection right_image_;
};

