#include "control_point.h"

ControlPoint::ControlPoint()
{
	id_ = 0;
	x_ = 0;
	y_ = 0;
	z_ = 0;
}

ControlPoint::~ControlPoint()
{
}

ControlPoint::ControlPoint(double x, double y, double z)
{
	x_ = x;
	y_ = y;
	z_ = z;
}

ControlPoint::ControlPoint(int id, double x, double y, double z)
{
	id_ = id;
	x_ = x;
	y_ = y;
	z_ = z;
}
