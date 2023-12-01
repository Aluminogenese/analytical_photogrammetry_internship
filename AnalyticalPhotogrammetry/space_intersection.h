#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<opencv2/opencv.hpp>
#include<fstream>
#include<vector>

#include"space_coordinates.h"
#include "camera.h"
#include "plane_coordinates.h"

using namespace std;
using namespace cv;

/// <summary>
/// 空间前方交会
/// </summary>
class SpaceIntersection
{
public:
	void get_param_from_file(const char* file_path, double& m, Camera& interior_elements, Mat_<double>& exterior_elements, vector<PlaneCoordinates>& image_points);
	void calculate_rotation_matrix(Mat_<double>& R, Mat_<double> exterior_elements);
	Mat calculate_auxiliary_coordinate(int i, Mat_<double>R, Camera interior_elements, vector<PlaneCoordinates> image_points);
	SpaceCoordinates N1N2_intersection(int i, double N1, double N2, double By, Mat_<double> left_auxiliary_coordinates, Mat_<double> right_auxiliary_coordinates);
	void pointfactor_space_intersection(const char* left_image_path, const char* right_image_path, const char* reslut_file_path);
public:
	Camera left_image_interior_elements_;// 内方位元素，x0, y0, f
	Mat_<double> left_image_exterior_elements_ = Mat::zeros(6, 1, CV_32F); //外方位元素x, y, z, φ, ω, κ
	vector<PlaneCoordinates> left_image_points_;// 像点坐标
	double left_m_;// 比例尺分母

	Camera right_image_interior_elements_;
	Mat_<double> right_image_exterior_elements_ = Mat::zeros(6, 1, CV_32F);
	vector<PlaneCoordinates> right_image_points_;
	double right_m_;

};
