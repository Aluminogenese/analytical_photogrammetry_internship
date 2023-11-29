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
/// �ռ�󷽽���
/// </summary>
class SpaceResection
{
public:
	static void get_param_from_file(const char* camera_file_name, const char* point_file_name, Camera& camera, vector<ImagePoint>& image_point, vector<ControlPoint>& control_point);
	static void Initialize(double m, Camera camera, vector<ControlPoint> control_point, double& Xs, double& Ys, double& Zs, double& phi, double& omega, double& kappa);
	static void calculate_rotation_matrix(Mat_<double>& R, double phi, double omega, double kappa);
	static void calculate_A_matrix(int i, Camera camera, vector<ControlPoint> control_point, vector<ImagePoint> image_point, double Z, double phi, double omega, double kappa, Mat_<double>R, Mat_<double>& A);
	static void calculate_L_matrix(int i, Mat_<double>& L, vector<ImagePoint> image_point, vector<ImagePoint> approxiate_image_point);
	static void correct_exterior_orientation_elements(double& Xs, double& Ys, double& Zs, double& phi, double& omega, double& kappa, Mat_<double> correction);
	static bool is_tolerant(Mat_<double> correction, double tolerance);
	static void calculate_space_resection(const char* camera_file_name, const char* point_file_name, const char* result_file_name);
};

