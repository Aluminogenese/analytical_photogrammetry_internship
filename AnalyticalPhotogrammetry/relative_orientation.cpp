#include "relative_orientation.h"

void RelativeOrientation::get_param_from_file(const char* image_path, ImageForRelativeOrientation& image_for_relative_orientation)
{
	FILE* fp = fopen(image_path, "r");
	char buffer[512] = { 0 };
	fgets(buffer, 512, fp);
	double x0, y0, f, m;
	sscanf(buffer, "%lf %lf %lf %lf\n", &f, &x0, &y0, &m);
	image_for_relative_orientation.camera_ = Camera(x0/1000.0, y0/1000.0, f/1000.0);
	image_for_relative_orientation.m_ = m;
	while (!feof(fp))
	{
		int id;
		double x, y;
		fgets(buffer, 512, fp);
		sscanf(buffer, "%d %lf %lf\n", &id, &x, &y);
		image_for_relative_orientation.image_point_.push_back(ImagePoint(id, x / 1000.0, y / 1000.0));
	}
}

void RelativeOrientation::calculate_relative_rotation_matrix(Mat_<double>& R)
{
	double phi = relative_orientation_elements_.at<double>(0, 0);
	double omega = relative_orientation_elements_.at<double>(1, 0);
	double kappa = relative_orientation_elements_.at<double>(2, 0);

	R.at<double>(0, 0) = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);//a1
	R.at<double>(0, 1) = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);//a2
	R.at<double>(0, 2) = -sin(phi) * cos(omega);//a3
	R.at<double>(1, 0) = cos(omega) * sin(kappa);//b1
	R.at<double>(1, 1) = cos(omega) * cos(kappa);//b2
	R.at<double>(1, 2) = -sin(omega);//b3
	R.at<double>(2, 0) = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);//c1
	R.at<double>(2, 1) = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);//c2
	R.at<double>(2, 2) = cos(phi) * cos(omega);//c3

}

Mat RelativeOrientation::calculate_auxiliary_coordinate(int i, Mat_<double> R, ImageForRelativeOrientation image_for_relative_orientation)
{
	Mat_<double> image_coordinate(3, 1);
	image_coordinate.at<double>(0, 0) = image_for_relative_orientation.image_point_[i].x_;
	image_coordinate.at<double>(1, 0) = image_for_relative_orientation.image_point_[i].y_;
	image_coordinate.at<double>(2, 0) = -image_for_relative_orientation.camera_.f_;
	Mat_<double> auxiliary_coordinates = R * image_coordinate;
	return auxiliary_coordinates;
}

void RelativeOrientation::calculate_A_matrix(int i, double Bx, double N2, Mat_<double>& A, Mat_<double> right_auxiliary_coordinates)
{
	double X2 = right_auxiliary_coordinates.at<double>(0, 0);
	double Y2 = right_auxiliary_coordinates.at<double>(1, 0);
	double Z2 = right_auxiliary_coordinates.at<double>(2, 0);
	A.at<double>(i, 0) = -X2 * Y2 * N2 / Z2;
	A.at<double>(i, 1) = -(Z2 + Y2 * Y2 / Z2) * N2;
	A.at<double>(i, 2) = X2 * N2;
	A.at<double>(i, 3) = Bx;
	A.at<double>(i, 4) = -Y2 * Bx / Z2;
}

void RelativeOrientation::calculate_L_matrix(int i, double N1, double N2, double By, Mat_<double> left_auxiliary_coordinates, Mat_<double> right_auxiliary_coordinates, Mat_<double>& L)
{
	double Y1 = left_auxiliary_coordinates.at<double>(1, 0);
	double Y2 = right_auxiliary_coordinates.at<double>(1, 0);

	L.at<double>(i, 0) = N1 * Y1 - N2 * Y2 - By;
}

void RelativeOrientation::correct_relative_param(Mat_<double> X)
{
	for (int i = 0; i < 5; i++)
	{
		relative_orientation_elements_.at<double>(i, 0) += X.at<double>(i, 0);
	}

}

bool RelativeOrientation::is_tolerant(Mat_<double> X, double tolerance)
{
	if (fabs(X.at<double>(0, 0)) < tolerance && fabs(X.at<double>(1, 0)) < tolerance && fabs(X.at<double>(2, 0)) < tolerance)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void RelativeOrientation::calculate_relative_orientation(const char* left_image_path, const char* right_image_path, const char* result_file_path)
{
	int iteration = 0;
	double tolerance = 3e-5;
	// ��������Ƭ���������x,y,�����ڷ�λԪ��x0,y0,f�ͱ�����m
	get_param_from_file(left_image_path, left_image_);
	get_param_from_file(right_image_path, right_image_);
	// ��������Ƭ����ת����
	Mat_<double> R_left = Mat::eye(3, 3, CV_32F);
	Mat_<double> R_right = Mat::zeros(3, 3, CV_32F);

	Mat_<double> X = Mat::zeros(5, 1, CV_32F); //dfai, domega, dkappa, du, dv 
	int point_num = left_image_.image_point_.size();
	Mat_<double> A = Mat::zeros(point_num, 5, CV_32F);
	Mat_<double> L = Mat::zeros(point_num, 1, CV_32F);
	double Bx = right_image_.image_point_[0].x_ - left_image_.image_point_[0].x_;
	double By = 0.0, Bz = 0.0;
	do
	{
		calculate_relative_rotation_matrix(R_right);
		By = Bx * relative_orientation_elements_.at<double>(3, 0);
		Bz = Bx * relative_orientation_elements_.at<double>(4, 0);
		for (int i = 0; i < point_num; i++)
		{
			// ����ռ����������ռ丨������
			Mat_<double> left_auxiliary_coordinates = calculate_auxiliary_coordinate(i, R_left, left_image_);
			Mat_<double> right_auxiliary_coordinates = calculate_auxiliary_coordinate(i, R_right, right_image_);
			// �����ͶӰϵ��
			double N1 = (Bx * right_auxiliary_coordinates.at<double>(2, 0) - Bz * right_auxiliary_coordinates.at<double>(0, 0)) / (left_auxiliary_coordinates.at<double>(0, 0) * right_auxiliary_coordinates.at<double>(2, 0) - left_auxiliary_coordinates.at<double>(2, 0) * right_auxiliary_coordinates.at<double>(0, 0));
			double N2 = (Bx * left_auxiliary_coordinates.at<double>(2, 0) - Bz * left_auxiliary_coordinates.at<double>(0, 0)) / (left_auxiliary_coordinates.at<double>(0, 0) * right_auxiliary_coordinates.at<double>(2, 0) - left_auxiliary_coordinates.at<double>(2, 0) * right_auxiliary_coordinates.at<double>(0, 0));

			calculate_A_matrix(i, Bx, N2, A, right_auxiliary_coordinates);
			calculate_L_matrix(i, N1, N2, By, left_auxiliary_coordinates, right_auxiliary_coordinates, L);
		}
		cout << "A:\n" << A << endl;
		X = (A.t() * A).inv() * A.t() * L;
		cout << "X:\n" << X << endl;
		correct_relative_param(X);
		cout << "relative_orientation_elements:\n" << relative_orientation_elements_ << endl;
		iteration += 1;

	} while (!is_tolerant(X, tolerance));

	cout << "Convergency!!!" << endl;
	cout << "Correction:" << endl;
	cout << X << endl;
	cout << "--------------------------------------------" << endl;
	cout << "Relative Orientation Result: " << endl;
	cout << "Iteration: " << iteration << endl;
	cout << "Residual:" << endl;
	Mat_<double>V = A * X - L;
	cout << V << endl;
	cout << "Five Parameters of Relative Orientation(fai, omega, kappa, u, v): " << endl;
	cout << relative_orientation_elements_.at<double>(0, 0) << " " << relative_orientation_elements_.at<double>(1, 0) << " " << relative_orientation_elements_.at<double>(2, 0) << " " << relative_orientation_elements_.at<double>(3, 0) << "  " << relative_orientation_elements_.at<double>(4, 0) << endl;
	Mat_<double>V_ = V.t() * V;
	double accuracy = sqrt(V_.at<double>(0, 0) / (point_num - 5));
	cout << "RMS Error��" << accuracy << endl;

	ofstream outfile;
	outfile.open(result_file_path, ios::out);
	outfile << "Residual:" << endl;
	outfile << V << endl;
	outfile << "Five Parameters of Relative Orientation(fai, omega, kappa, u, v): " << endl;
	outfile << relative_orientation_elements_.at<double>(0, 0) << " " << relative_orientation_elements_.at<double>(1, 0) << " " << relative_orientation_elements_.at<double>(2, 0) << " " << relative_orientation_elements_.at<double>(3, 0) << "  " << relative_orientation_elements_.at<double>(4, 0) << endl;
	outfile << "RMS Error��" << accuracy << endl;
	outfile << "Iteration: " << iteration << endl;
}
