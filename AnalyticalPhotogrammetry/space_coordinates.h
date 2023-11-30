#pragma once
/// <summary>
/// ¿Õ¼ä×ø±ê
/// </summary>
class SpaceCoordinates
{
public:
	SpaceCoordinates();
	~SpaceCoordinates();
	SpaceCoordinates(double x, double y, double z);
	SpaceCoordinates(int id, double x, double y, double z);

public:
	int id_;
	double x_;
	double y_;
	double z_;
};

