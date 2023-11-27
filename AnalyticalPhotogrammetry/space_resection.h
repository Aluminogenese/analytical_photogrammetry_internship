#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<vector>
#include<opencv2/opencv.hpp>
#include<fstream>
#include"control_point.h"
#include"image_point.h"
#include"camera.h"
#include"exterior_orientation_elements.h"
using namespace std;
using namespace cv;
/// <summary>
/// 空间后方交会
/// </summary>
class SpaceResection
{
public:
	void get_param_from_file(const char* camera_file_name, const char* point_file_name);
	void Initialize();
	void calculate_rotation_matrix(Mat_<double>& R);
	void calculate_A_matrix(int i, Mat_<double>R, Mat_<double>& A);
	void calculate_L_matrix(int i, Mat_<double>& L, vector<ImagePoint> approxiate_image_point);
	void correct_exterior_orientation_elements(ExteriorOrientationElements correction);
	bool if_tolerant(ExteriorOrientationElements correction, double tolerance);
	void calculate_space_resection(const char* camera_file_name, const char* point_file_name, const char* result_file_name);
public:
	vector<ImagePoint> image_point_;
	vector<ControlPoint> control_point_;
	Camera camera_;
	ExteriorOrientationElements exterior_orientation_elements_;
	double m_ = 50000.0;
};

