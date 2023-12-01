#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<vector>
#include<opencv2/opencv.hpp>
#include<fstream>
#include"space_coordinates.h"
#include"plane_coordinates.h"
#include"camera.h"
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
	void calculate_L_matrix(int i, Mat_<double>R, Mat_<double>& L);
	void correct_exterior_orientation_elements(Mat_<double> X);
	bool is_tolerant(Mat_<double> correction, double tolerance);
	void calculate_space_resection(const char* camera_file_name, const char* point_file_name, const char* result_file_name);
public:
	double m_ = 50000.0;// 比例尺分母
	Camera camera_;// 相机参数x0，y0，f
	vector<PlaneCoordinates> image_points_;// 影像坐标
	vector<SpaceCoordinates> control_points_;// 地面坐标
	Mat_<double> exterior_orientation_elements_ = Mat::zeros(6, 1, CV_32F); //外方位元素x, y, z, φ, ω, κ
};

