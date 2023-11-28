#pragma once
#include "camera.h"
#include<vector>
#include "exterior_orientation_elements.h"
#include "image_point.h"
using namespace std;
class ImageForIntersection
{
public:
	Camera camera_;
	double m_;
	ExteriorOrientationElements exterior_orientation_elements_;
	vector<ImagePoint> image_point_;
};

