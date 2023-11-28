﻿// AnalyticalPhotogrammetry.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include"space_resection.h"
#include"space_intersection.h"
#include"interior_orientation.h"

int main()
{
	//SpaceResection space_resection_handler;
	//space_resection_handler.calculate_space_resection("./Data/camera.txt", "./Data/space_resection.txt", "./Result/space_resection_result.txt");

	//SpaceIntersection space_intersection_handler;
	//space_intersection_handler.pointfactor_space_intersection("./Data/space_intersection/0320.txt", "./Data/space_intersection/0319.txt", "./Result/space_intersection_result.txt");

	InteriorOrientation interior_orientation_handler;
	interior_orientation_handler.affine_interior_orientation("./Data/interior_orientation/0319.txt", "./Result/interior_orientation_result.txt");
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
