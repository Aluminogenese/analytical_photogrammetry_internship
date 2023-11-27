#include "exterior_orientation_elements.h"

ExteriorOrientationElements::ExteriorOrientationElements()
{
}

ExteriorOrientationElements::~ExteriorOrientationElements()
{
}

ExteriorOrientationElements::ExteriorOrientationElements(double Xs, double Ys, double Zs, double phi, double omega, double kappa)
{
	Xs_ = Xs;
	Ys_ = Ys;
	Zs_ = Zs;
	phi_ = phi;
	omega_ = omega;
	kappa_ = kappa;
}
