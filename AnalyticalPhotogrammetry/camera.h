#pragma once

class Camera
{
public:
	Camera();
	~Camera();
	Camera(double x0, double y0, double f);
public:
	double x0_;
	double y0_;
	double f_;
};

