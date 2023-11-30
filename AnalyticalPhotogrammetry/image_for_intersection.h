#pragma once
#include "camera.h"
#include<vector>
#include "exterior_orientation_elements.h"
#include "plane_coordinates.h"
using namespace std;
class ImageForIntersection
{
public:
	Camera camera_;
	double m_;
	ExteriorOrientationElements exterior_orientation_elements_;
	vector<PlaneCoordinates> image_point_;
};

