#include "space_coordinates.h"

SpaceCoordinates::SpaceCoordinates()
{
	id_ = 0;
	x_ = 0;
	y_ = 0;
	z_ = 0;
}

SpaceCoordinates::~SpaceCoordinates()
{
}

SpaceCoordinates::SpaceCoordinates(double x, double y, double z)
{
	x_ = x;
	y_ = y;
	z_ = z;
}

SpaceCoordinates::SpaceCoordinates(int id, double x, double y, double z)
{
	id_ = id;
	x_ = x;
	y_ = y;
	z_ = z;
}
