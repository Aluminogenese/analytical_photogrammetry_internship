#pragma once
/// <summary>
/// Ïñµã×ø±ê
/// </summary>
class ImagePoint
{
public:
	ImagePoint();
	~ImagePoint();
	ImagePoint(double x, double y);
	ImagePoint(int id, double x, double y);
public:
	int id_;
	double x_;
	double y_;
};

