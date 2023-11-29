// AnalyticalPhotogrammetry.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include<opencv2/opencv.hpp>
using namespace cv;
using namespace std;
//#include"space_resection.h"
//#include"space_intersection.h"
//#include"interior_orientation.h"
#include"relative_orientation.h"
int main()
{
	//SpaceResection space_resection_handler;
	//space_resection_handler.calculate_space_resection("./Data/camera.txt", "./Data/space_resection.txt", "./Result/space_resection_result.txt");

	//SpaceIntersection space_intersection_handler;
	//space_intersection_handler.pointfactor_space_intersection("./Data/space_intersection/0320.txt", "./Data/space_intersection/0319.txt", "./Result/space_intersection_result.txt");

	//InteriorOrientation interior_orientation_handler;
	//interior_orientation_handler.affine_interior_orientation("./Data/interior_orientation/0319.txt", "./Result/interior_orientation_result.txt");

	RelativeOrientation relative_orientation_handler;
	relative_orientation_handler.calculate_relative_orientation("./Data/relative_orientation/0320.txt", "./Data/relative_orientation/0319.txt", "./Result/relative_orientation_result");

}