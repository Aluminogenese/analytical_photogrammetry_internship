#pragma once
class ExteriorOrientationElements
{
public:
	ExteriorOrientationElements();
	~ExteriorOrientationElements();
	ExteriorOrientationElements(double Xs, double Ys, double Zs, double phi, double omega, double kappa);
public:
	double Xs_;
	double Ys_;
	double Zs_;
	double phi_;
	double omega_;
	double kappa_;
};

