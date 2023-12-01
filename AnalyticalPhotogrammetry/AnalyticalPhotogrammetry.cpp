#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include<opencv2/opencv.hpp>
using namespace cv;
using namespace std;
#include"space_resection.h"
#include"space_intersection.h"
#include"interior_orientation.h"
#include"relative_orientation.h"
#include"absolute_orientation.h"
int main()
{
	while (1)
	{
		string work_pattern = "";
		cout << "Choose work_patter: A(Space Resection), B(Space Intersection), C(Interior Orientation)" << endl;
		cout << "D(Relative Orientation), E(Absolute Orientation)" << endl;
		cin >> work_pattern;
		if (work_pattern == "A") {
			SpaceResection space_resection_handler;
			space_resection_handler.calculate_space_resection("./Data/space_resection/camera.txt", "./Data/space_resection/space_resection.txt", "./Result/space_resection_result.txt");
		}
		else if (work_pattern == "B") {
			SpaceIntersection space_intersection_handler;
			space_intersection_handler.pointfactor_space_intersection("./Data/space_intersection/0320.txt", "./Data/space_intersection/0319.txt", "./Result/space_intersection_result.txt");
		}
		else if (work_pattern == "C") {
			InteriorOrientation interior_orientation_handler;
			interior_orientation_handler.affine_interior_orientation("./Data/interior_orientation/0319.txt", "./Result/interior_orientation_result.txt");
		}
		else if (work_pattern == "D") {
			RelativeOrientation relative_orientation_handler;
			relative_orientation_handler.calculate_relative_orientation("./Data/relative_orientation/0320.txt", "./Data/relative_orientation/0319.txt", "./Result/relative_orientation_result.txt");
		}
		else if (work_pattern == "E") {
			AbsoluteOrientation absolute_orientation_handler;
			absolute_orientation_handler.calculate_absolute_orientation("./Data/absolute_orientation/absolute_orientation.txt", "./Result/absolute_orientation_result.txt");
		}
		else
			break;
	}
}