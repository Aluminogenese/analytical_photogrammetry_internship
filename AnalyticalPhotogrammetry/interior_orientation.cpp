#include "interior_orientation.h"

void InteriorOrientation::get_param_from_file(const char* file_path)
{
	FILE* fp = fopen(file_path, "r");
	char buffer[512] = { 0 };
	// ����ע����
	for (int i = 0; i < 4; i++) {
		fgets(buffer, 512, fp);
	}
	fgets(buffer, 512, fp);
	sscanf(buffer, "%lf %lf %lf\n", &camera_.f_, &camera_.x0_, &camera_.y0_);

	fgets(buffer, 512, fp);
	sscanf(buffer, "%lf\n", &pixel_size_);

	fgets(buffer, 512, fp);
	sscanf(buffer, "%lf %lf\n", &width_, &height_);

	while (!feof(fp))
	{
		double xm, ym, xp, yp;
		fgets(buffer, 512, fp);
		sscanf(buffer, "%lf %lf %lf %lf\n", &xm, &ym, &xp, &yp);
		frame_coordinate_points_.push_back(PlaneCoordinates(xm, ym));
		pixel_coordinate_points_.push_back(PlaneCoordinates(xp, yp));
	}
	fclose(fp);
}

void InteriorOrientation::calculate_A_matrix(Mat_<double>& A)
{
	int point_num = pixel_coordinate_points_.size();
	for (int i = 0; i < point_num; i++)
	{
		A.at<double>(i * 2, 0) = 1;
		A.at<double>(i * 2, 1) = pixel_coordinate_points_[i].x_ * pixel_size_;
		A.at<double>(i * 2, 2) = -pixel_coordinate_points_[i].y_ * pixel_size_;
		A.at<double>(i * 2, 3) = 0;
		A.at<double>(i * 2, 4) = 0;
		A.at<double>(i * 2, 5) = 0;

		A.at<double>(i * 2 + 1, 0) = 0;
		A.at<double>(i * 2 + 1, 1) = 0;
		A.at<double>(i * 2 + 1, 2) = 0;
		A.at<double>(i * 2 + 1, 3) = 1;
		A.at<double>(i * 2 + 1, 4) = pixel_coordinate_points_[i].x_ * pixel_size_;
		A.at<double>(i * 2 + 1, 5) = -pixel_coordinate_points_[i].y_ * pixel_size_;
	}
}

void InteriorOrientation::calculate_L_matrix(Mat_<double>& L)
{
	int point_num = frame_coordinate_points_.size();
	for (int i = 0; i < point_num; i++)
	{
		L.at<double>(i * 2, 0) = frame_coordinate_points_[i].x_;
		L.at<double>(i * 2 + 1, 0) = frame_coordinate_points_[i].y_;
	}
}

void InteriorOrientation::affine_interior_orientation(const char* file_path, const char* result_file_path)
{
	get_param_from_file(file_path);

    double x_offset = (width_ * pixel_size_) / 2;
    double y_offset = -(height_ * pixel_size_) / 2;
    int point_num = frame_coordinate_points_.size();

    Mat_<double> A = Mat::zeros(point_num * 2, 6, CV_32F);
    Mat_<double> L = Mat::zeros(point_num * 2, 1, CV_32F);
    Mat_<double> parameters = Mat::zeros(6, 1, CV_32F);

    calculate_A_matrix(A);
    calculate_L_matrix(L);

    parameters = (A.t() * A).inv() * A.t() * L;
    Mat_<double> V = A * parameters - L;
    Mat_<double> V_ = V.t() * V;
    double accuracy = sqrt(V_.at<double>(0, 0) / (point_num * 2 - 6));

    double m0 = parameters.at<double>(0, 0);
    double m1 = parameters.at<double>(1, 0);
    double m2 = parameters.at<double>(2, 0);
    double n0 = parameters.at<double>(3, 0);
    double n1 = parameters.at<double>(4, 0);
    double n2 = parameters.at<double>(5, 0);

	vector<PlaneCoordinates> calculted_points;
	for (int i = 0; i < point_num; i++)
	{
		double x = m0 + m1 * pixel_coordinate_points_[i].x_ * pixel_size_ - m2 * pixel_coordinate_points_[i].y_ * pixel_size_;
		double y = n0 + n1 * pixel_coordinate_points_[i].x_ * pixel_size_ - n2 * pixel_coordinate_points_[i].y_ * pixel_size_;
		calculted_points.push_back(PlaneCoordinates(x, y));
	}

	cout << "--------------------------------------------" << endl;
    cout << "Interior Orientation Result" << endl;
    cout << "m0 = " << fixed << setprecision(5) << m0 << endl;
    cout << "m1 = " << fixed << setprecision(5) << m1 << endl;
    cout << "m2 = " << fixed << setprecision(5) << m2 << endl;
    cout << "n0 = " << fixed << setprecision(5) << n0 << endl;
    cout << "n1 = " << fixed << setprecision(5) << n1 << endl;
    cout << "n2 = " << fixed << setprecision(5) << n2 << endl;
    cout << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
	cout << "x\ty\tdx\tdy" << endl;
	for (int i = 0; i < calculted_points.size(); i++)
	{
		cout << fixed << setprecision(5) << calculted_points[i].x_ << "\t" << calculted_points[i].y_ << "\t" << (frame_coordinate_points_[i].x_ - calculted_points[i].x_) << " \t" << (frame_coordinate_points_[i].y_ - calculted_points[i].y_) << endl;
	}
    ofstream outfile;
    outfile.open(result_file_path, ios::out);
	outfile << "--------------------------------------------" << endl;
	outfile << "Interior Orientation Result" << endl;
    outfile << "m0 = " << fixed << setprecision(5) << m0 << endl;
    outfile << "m1 = " << fixed << setprecision(5) << m1 << endl;
    outfile << "m2 = " << fixed << setprecision(5) << m2 << endl;
    outfile << "n0 = " << fixed << setprecision(5) << n0 << endl;
    outfile << "n1 = " << fixed << setprecision(5) << n1 << endl;
    outfile << "n2 = " << fixed << setprecision(5) << n2 << endl;
    outfile << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
	outfile << "x\ty\tdx\tdy" << endl;
	for (int i = 0; i < calculted_points.size(); i++)
	{
		outfile << fixed << setprecision(5) << calculted_points[i].x_ << "\t" << calculted_points[i].y_ << "\t" << frame_coordinate_points_[i].x_ - calculted_points[i].x_ << " \t" << frame_coordinate_points_[i].y_ - calculted_points[i].y_ << endl;
	}
}

