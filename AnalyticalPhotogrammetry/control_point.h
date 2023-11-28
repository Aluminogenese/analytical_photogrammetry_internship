#pragma once
/// <summary>
/// ¿ØÖÆµã×ø±ê
/// </summary>
class ControlPoint
{
public:
	ControlPoint();
	~ControlPoint();
	ControlPoint(double x, double y, double z);
	ControlPoint(int id, double x, double y, double z);

public:
	int id_;
	double x_;
	double y_;
	double z_;
};

