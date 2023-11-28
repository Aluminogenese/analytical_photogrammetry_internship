#include "image_point.h"

ImagePoint::ImagePoint()
{
	id_ = 0;
	x_ = 0;
	y_ = 0;
}

ImagePoint::~ImagePoint()
{
}

ImagePoint::ImagePoint(double x, double y)
{
	x_ = x;
	y_ = y;
}

ImagePoint::ImagePoint(int id, double x, double y)
{
	id_ = id;
	x_ = x;
	y_ = y;
}
