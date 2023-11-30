#include "plane_coordinates.h"

PlaneCoordinates::PlaneCoordinates()
{
	id_ = 0;
	x_ = 0;
	y_ = 0;
}

PlaneCoordinates::~PlaneCoordinates()
{
}

PlaneCoordinates::PlaneCoordinates(double x, double y)
{
	x_ = x;
	y_ = y;
}

PlaneCoordinates::PlaneCoordinates(int id, double x, double y)
{
	id_ = id;
	x_ = x;
	y_ = y;
}
