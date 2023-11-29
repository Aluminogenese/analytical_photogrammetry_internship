#pragma once
#include<vector>
#include"camera.h"
#include"image_point.h"
using namespace std;
class ImageForRelativeOrientation
{
public:
	Camera camera_;
	double m_;
	vector<ImagePoint> image_point_;
};

