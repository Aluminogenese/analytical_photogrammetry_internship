#pragma once
#include<vector>
#include"camera.h"
#include"plane_coordinates.h"
using namespace std;
class ImageForRelativeOrientation
{
public:
	Camera camera_;
	double m_;
	vector<PlaneCoordinates> image_point_;
};

