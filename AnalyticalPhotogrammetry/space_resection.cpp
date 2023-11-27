#include "space_resection.h"

void SpaceResection::get_param_from_file(const char* camera_file_name, const char* point_file_name)
{
	FILE* fp_point = fopen(point_file_name, "r");
	char buffer[512] = { 0 };
	while (!feof(fp_point))
	{
		fgets(buffer, 512, fp_point);
		int point_id;
		double x=0., y=0., X=0., Y=0., Z=0.;
		ImagePoint tmp_image_point;
		ControlPoint tmp_control_point;
		sscanf(buffer, "%d %lf %lf %lf %lf %lf\n", &point_id, &x, &y, &X, &Y, &Z);
		tmp_image_point = ImagePoint(x, y);
		tmp_control_point = ControlPoint(X, Y, Z);
		image_point_.push_back(tmp_image_point);
		control_point_.push_back(tmp_control_point);
	}
	fclose(fp_point);

	FILE* fp_camera = fopen(camera_file_name, "r");
	char tmp_char[512]{};
	fgets(tmp_char, 512, fp_camera);
	sscanf(tmp_char, "%lf %lf %lf\n", &camera_.x0_, &camera_.y0_, &camera_.f_);
	camera_.f_ /= 1000.0;
	fclose(fp_camera);
}

void SpaceResection::Initialize()
{
	double Z = 0;
	int point_num = image_point_.size();
	for (vector<ControlPoint>::const_iterator it = control_point_.cbegin(); it != control_point_.cend(); it++)
	{
		exterior_orientation_elements_.Xs_ += (*it).x_;
		exterior_orientation_elements_.Ys_ += (*it).y_;
		//exterior_orientation_elements_.Zs_ += (*it).z_;
		Z += (*it).z_;
	}
	exterior_orientation_elements_.Xs_ /= point_num;
	exterior_orientation_elements_.Ys_ /= point_num;
	//exterior_orientation_elements_.Zs_ /= point_num;
	exterior_orientation_elements_.Zs_ = m_ * camera_.f_;
	exterior_orientation_elements_.phi_ = 0;
	exterior_orientation_elements_.omega_ = 0;
	exterior_orientation_elements_.kappa_ = 0;
}

void SpaceResection::calculate_rotation_matrix(Mat_<double>& R)
{
	double phi = exterior_orientation_elements_.phi_;
	double omega = exterior_orientation_elements_.omega_;
	double kappa = exterior_orientation_elements_.kappa_;

	R.at<double>(0, 0) = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
	R.at<double>(0, 1) = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
	R.at<double>(0, 2) = -sin(phi) * cos(omega);
	R.at<double>(1, 0) = cos(omega) * sin(kappa);
	R.at<double>(1, 1) = cos(omega) * cos(kappa);
	R.at<double>(1, 2) = -sin(omega);
	R.at<double>(2, 0) = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
	R.at<double>(2, 1) = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
	R.at<double>(2, 2) = cos(phi) * cos(omega);
}

void SpaceResection::calculate_A_matrix(int i, Mat_<double>R, Mat_<double>& A)
{
	double phi = exterior_orientation_elements_.phi_;
	double omega = exterior_orientation_elements_.omega_;
	double kappa = exterior_orientation_elements_.kappa_;
	double Z = R.at<double>(0, 2) * (control_point_[i].x_ - exterior_orientation_elements_.Xs_) + R.at<double>(1, 2) * (control_point_[i].y_ - exterior_orientation_elements_.Ys_) + R.at<double>(2, 2) * (control_point_[i].z_ - exterior_orientation_elements_.Zs_);

	A.at<double>(i * 2, 0) = (R.at<double>(0, 0) * camera_.f_ + R.at<double>(0, 2) * (image_point_[i].x_ - camera_.x0_)) / Z;//a11
	A.at<double>(i * 2, 1) = (R.at<double>(1, 0) * camera_.f_ + R.at<double>(1, 2) * (image_point_[i].x_ - camera_.x0_)) / Z;//a12
	A.at<double>(i * 2, 2) = (R.at<double>(2, 0) * camera_.f_ + R.at<double>(2, 2) * (image_point_[i].x_ - camera_.x0_)) / Z;//a13
	A.at<double>(i * 2, 3) = (image_point_[i].y_ - camera_.y0_) * sin(omega) - ((image_point_[i].x_ - camera_.x0_) * ((image_point_[i].x_ - camera_.x0_) * cos(kappa) - (image_point_[i].y_ - camera_.y0_) * sin(kappa)) / camera_.f_ + camera_.f_ * cos(kappa)) * cos(omega);//a14
	A.at<double>(i * 2, 4) = -camera_.f_ * sin(kappa) - (image_point_[i].x_ - camera_.x0_) * ((image_point_[i].x_ - camera_.x0_) * sin(kappa) + (image_point_[i].y_ - camera_.y0_) * cos(kappa)) / camera_.f_;//a15
	A.at<double>(i * 2, 5) = image_point_[i].y_ - camera_.y0_;//a16
	A.at<double>(i * 2 + 1, 0) = (R.at<double>(0, 1) * camera_.f_ + R.at<double>(0, 2) * (image_point_[i].y_ - camera_.y0_)) / Z;//a21
	A.at<double>(i * 2 + 1, 1) = (R.at<double>(1, 1) * camera_.f_ + R.at<double>(1, 2) * (image_point_[i].y_ - camera_.y0_)) / Z;//a22
	A.at<double>(i * 2 + 1, 2) = (R.at<double>(2, 1) * camera_.f_ + R.at<double>(2, 2) * (image_point_[i].y_ - camera_.y0_)) / Z;//a23
	A.at<double>(i * 2 + 1, 3) = -(image_point_[i].x_ - camera_.x0_) * sin(omega) - ((image_point_[i].y_ - camera_.y0_) * ((image_point_[i].x_ - camera_.x0_) * cos(kappa) - (image_point_[i].y_ - camera_.y0_) * sin(kappa)) / camera_.f_ - camera_.f_ * sin(kappa)) * cos(omega);//a24
	A.at<double>(i * 2 + 1, 4) = -camera_.f_ * cos(kappa) - ((image_point_[i].y_ - camera_.y0_) * ((image_point_[i].x_ - camera_.x0_) * sin(kappa) + (image_point_[i].y_ - camera_.y0_) * cos(kappa))) / camera_.f_;//a25
	A.at<double>(i * 2 + 1, 5) = -(image_point_[i].x_ - camera_.x0_);//a26
}

void SpaceResection::calculate_L_matrix(int i, Mat_<double>& L, vector<ImagePoint> approxiate_image_point)
{
	L.at<double>(i * 2, 0) = image_point_[i].x_ - approxiate_image_point[i].x_;
	L.at<double>(i * 2 + 1, 0) = image_point_[i].y_ - approxiate_image_point[i].y_;
}

void SpaceResection::correct_exterior_orientation_elements(ExteriorOrientationElements correction)
{
	exterior_orientation_elements_.Xs_ += correction.Xs_;
	exterior_orientation_elements_.Ys_ += correction.Ys_;
	exterior_orientation_elements_.Zs_ += correction.Zs_;
	exterior_orientation_elements_.phi_ += correction.phi_;
	exterior_orientation_elements_.omega_ += correction.omega_;
	exterior_orientation_elements_.kappa_ += correction.kappa_;
}

bool SpaceResection::if_tolerant(ExteriorOrientationElements correction, double tolerance)
{
	if (fabs(correction.phi_) < tolerance && fabs(correction.omega_) < tolerance && fabs(correction.kappa_) < tolerance)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void SpaceResection::calculate_space_resection(const char* camera_file_name, const char* point_file_name, const char* result_file_name)
{
	int iteration = 0;//迭代次数
	double tolerance = 0.01;//限差
	Mat_<double> R = Mat::zeros(3, 3, CV_32F);
	Mat_<double> correction_matrix = Mat::zeros(6, 1, CV_32F);

	ExteriorOrientationElements correction;

	get_param_from_file(camera_file_name, point_file_name);
    //Get initial parameters
	Initialize();
	int point_num = image_point_.size();

	Mat_<double> A = Mat::zeros(point_num * 2, 6, CV_32F);
	Mat_<double> L = Mat::zeros(point_num * 2, 1, CV_32F);
	vector<ImagePoint> approxiate_image_point(point_num);

	//Adjustment process with iterations
	do
	{
		//Calculate rotation matrix
		calculate_rotation_matrix(R);
		//Calculate approximation of x,y of each point
		//Form norm equation of each point
		for (int i = 0; i < point_num; i++)
		{
			//Approximation
			approxiate_image_point[i].x_ = camera_.x0_ - camera_.f_ * (R.at<double>(0, 0) * (control_point_[i].x_ - exterior_orientation_elements_.Xs_) + R.at<double>(1, 0) * (control_point_[i].y_ - exterior_orientation_elements_.Ys_) + R.at<double>(2, 0) * (control_point_[i].z_ - exterior_orientation_elements_.Zs_))
				/ (R.at<double>(0, 2) * (control_point_[i].x_ - exterior_orientation_elements_.Xs_) + R.at<double>(1, 2) * (control_point_[i].y_ - exterior_orientation_elements_.Ys_) + R.at<double>(2, 2) * (control_point_[i].z_ - exterior_orientation_elements_.Zs_));
			approxiate_image_point[i].y_ = camera_.y0_ - camera_.f_ * (R.at<double>(0, 1) * (control_point_[i].x_ - exterior_orientation_elements_.Xs_) + R.at<double>(1, 1) * (control_point_[i].y_ - exterior_orientation_elements_.Ys_) + R.at<double>(2, 1) * (control_point_[i].z_ - exterior_orientation_elements_.Zs_))
				/ (R.at<double>(0, 2) * (control_point_[i].x_ - exterior_orientation_elements_.Xs_) + R.at<double>(1, 2) * (control_point_[i].y_ - exterior_orientation_elements_.Ys_) + R.at<double>(2, 2) * (control_point_[i].z_ - exterior_orientation_elements_.Zs_));
			//Norm equation
			calculate_A_matrix(i, R, A);
			calculate_L_matrix(i, L, approxiate_image_point);
		}
		//Calculate correction according to norm equation
		//X = inv(A' * A) * A' * L
		correction_matrix = (A.t() * A).inv() * A.t() * L;
		correction = ExteriorOrientationElements(correction_matrix.at<double>(0, 0), correction_matrix.at<double>(1, 0), correction_matrix.at<double>(2, 0), correction_matrix.at<double>(3, 0),
			correction_matrix.at<double>(4, 0), correction_matrix.at<double>(5, 0));
		correct_exterior_orientation_elements(correction);
		iteration += 1;
	} while (!if_tolerant(correction, tolerance));
	Mat_<double> V = A * correction_matrix - L;
	Mat_<double> V_ = V.t() * V;
	double accuracy = sqrt(V_.at<double>(0, 0) / (point_num * 2 - 6));

	////Output final result
	cout << "--------------------------------------------" << endl;
	cout << "Space Resection Result" << endl;
	cout << "Iteration: " << iteration << endl;
	iteration = 0;
	cout << "Xs = " << fixed << setprecision(2) << exterior_orientation_elements_.Xs_ << endl;
	cout << "Ys = " << fixed << setprecision(2) << exterior_orientation_elements_.Ys_ << endl;
	cout << "Zs = " << fixed << setprecision(2) << exterior_orientation_elements_.Zs_ << endl;
	cout << "fai = " << fixed << setprecision(5) << exterior_orientation_elements_.phi_ << endl;
	cout << "omega = " << fixed << setprecision(5) << exterior_orientation_elements_.omega_ << endl;
	cout << "kappa = " << fixed << setprecision(5) << exterior_orientation_elements_.kappa_ << endl;
	cout << "Rotation Matrix:" << endl;
	cout << fixed << setprecision(5) << R << endl;
	cout << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
	cout << endl;

	ofstream outfile;
	outfile.open(result_file_name, ios::out);
	outfile << "Xs = " << fixed << setprecision(2) << exterior_orientation_elements_.Xs_ << endl;
	outfile << "Ys = " << fixed << setprecision(2) << exterior_orientation_elements_.Ys_ << endl;
	outfile << "Zs = " << fixed << setprecision(2) << exterior_orientation_elements_.Zs_ << endl;
	outfile << "fai = " << fixed << setprecision(5) << exterior_orientation_elements_.phi_ << endl;
	outfile << "omega = " << fixed << setprecision(5) << exterior_orientation_elements_.omega_ << endl;
	outfile << "kappa = " << fixed << setprecision(5) << exterior_orientation_elements_.kappa_ << endl;
	outfile << "Rotation Matrix:" << endl;
	outfile << fixed << setprecision(5) << R << endl;
	outfile << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
	outfile << endl;

}


