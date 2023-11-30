#pragma once
/// <summary>
/// Ïñµã×ø±ê
/// </summary>
class PlaneCoordinates
{
public:
	PlaneCoordinates();
	~PlaneCoordinates();
	PlaneCoordinates(double x, double y);
	PlaneCoordinates(int id, double x, double y);
public:
	int id_;
	double x_;
	double y_;
};

