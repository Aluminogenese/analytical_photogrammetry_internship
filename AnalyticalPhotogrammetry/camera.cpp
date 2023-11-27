#include "camera.h"

Camera::Camera()
{
	x0_ = 0;
	y0_ = 0;
	f_ = 0;
}

Camera::~Camera()
{
}

Camera::Camera(double x0, double y0, double f)
{
	x0_ = x0;
	y0_ = y0;
	f_ = f;
}
