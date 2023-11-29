#include "space_intersection.h"

void SpaceIntersection::get_param_from_file(const char* image_path, ImageForIntersection& image_for_intersection)
{
	FILE* fp = fopen(image_path, "r");
	char buffer[512] = { 0 };
	// 读入注释行
	fgets(buffer, 512, fp);
	fgets(buffer, 512, fp);

	double f = 0.0, x0 = 0.0, y0 = 0.0, m = 0.0;
	fgets(buffer, 512, fp);
	sscanf(buffer, "%lf %lf %lf %lf\n", &f, &x0, &y0, &m);
	image_for_intersection.camera_.f_ = f / 1000.0;
	image_for_intersection.camera_.x0_ = x0 / 1000.0;
	image_for_intersection.camera_.y0_ = y0 / 10000.0;
	image_for_intersection.m_ = m;
	// 读入注释行
	fgets(buffer, 512, fp);
	fgets(buffer, 512, fp);

	double phi = 0.0, omega = 0.0, kappa = 0.0, Xs = 0.0, Ys = 0.0, Zs = 0.0;
	fgets(buffer, 512, fp);
	sscanf(buffer, "%lf %lf %lf %lf %lf %lf\n", &phi, &omega, &kappa, &Xs, &Ys, &Zs);
	image_for_intersection.exterior_orientation_elements_.phi_ = phi * CV_PI / 180.0;
	image_for_intersection.exterior_orientation_elements_.omega_ = omega * CV_PI / 180.0;
	image_for_intersection.exterior_orientation_elements_.kappa_ = kappa * CV_PI / 180.0;
	image_for_intersection.exterior_orientation_elements_.Xs_ = Xs;
	image_for_intersection.exterior_orientation_elements_.Ys_ = Ys;
	image_for_intersection.exterior_orientation_elements_.Zs_ = Zs;
	// 读入注释行
	fgets(buffer, 512, fp);
	fgets(buffer, 512, fp);
	while (!feof(fp))
	{
		fgets(buffer, 512, fp);
		int id = 0;
		double x=0.0, y=0.0;
		sscanf(buffer, "%d %lf %lf", &id, &x, &y);
		image_for_intersection.image_point_.push_back(ImagePoint(id, x / 1000.0, y / 1000.0));
	}
	fclose(fp);
}

void SpaceIntersection::calculate_rotation_matrix(Mat_<double>& R, ImageForIntersection image_for_intersection)
{
	double phi = image_for_intersection.exterior_orientation_elements_.phi_;
	double omega = image_for_intersection.exterior_orientation_elements_.omega_;
	double kappa = image_for_intersection.exterior_orientation_elements_.kappa_;

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

Mat SpaceIntersection::calculate_auxiliary_coordinate(int i, Mat_<double> R, ImageForIntersection image_for_intersection)
{
	Mat_<double> image_coordinate(3, 1);
	image_coordinate.at<double>(0, 0) = image_for_intersection.image_point_[i].x_;
	image_coordinate.at<double>(1, 0) = image_for_intersection.image_point_[i].y_;
	image_coordinate.at<double>(2, 0) = -image_for_intersection.camera_.f_;
	Mat_<double> auxiliary_coordinates = R * image_coordinate;
	return auxiliary_coordinates;
}

ControlPoint SpaceIntersection::N1N2_intersection(int i, ImageForIntersection left_image, double N1, Mat_<double> left_auxiliary_coordinates)
{
	double X = left_image.exterior_orientation_elements_.Xs_ + N1 * left_auxiliary_coordinates.at<double>(0, 0);
	double Y = left_image.exterior_orientation_elements_.Ys_ + N1 * left_auxiliary_coordinates.at<double>(1, 0);
	double Z = left_image.exterior_orientation_elements_.Zs_ + N1 * left_auxiliary_coordinates.at<double>(2, 0);
	return ControlPoint(left_image.image_point_[i].id_, X, Y, Z);
}

void SpaceIntersection::pointfactor_space_intersection(const char* left_image_path, const char* right_image_path, const char* reslut_file_path)
{
	vector<ControlPoint> vec_control_point;
	// 获取左右像片已知数据x0,y0,f,Xs,Ys,Zs,φ,ω,κ
	get_param_from_file(left_image_path, left_image_);
	get_param_from_file(right_image_path, right_image_);
	// 计算基线分量
	double Bx = right_image_.exterior_orientation_elements_.Xs_ - left_image_.exterior_orientation_elements_.Xs_;
	double By = right_image_.exterior_orientation_elements_.Ys_ - left_image_.exterior_orientation_elements_.Ys_;
	double Bz = right_image_.exterior_orientation_elements_.Zs_ - left_image_.exterior_orientation_elements_.Zs_;
	// 计算每张像片的旋转矩阵
	Mat_<double> R_left = Mat::zeros(3, 3, CV_32F);
	Mat_<double> R_right = Mat::zeros(3, 3, CV_32F);
	calculate_rotation_matrix(R_left, left_image_);
	calculate_rotation_matrix(R_right, right_image_);
	int point_num = left_image_.image_point_.size();
	for (int i = 0; i < point_num; i++) 
	{
		// 由像空间坐标计算像空间辅助坐标
		Mat_<double> left_auxiliary_coordinates = calculate_auxiliary_coordinate(i, R_left, left_image_);
		Mat_<double> right_auxiliary_coordinates = calculate_auxiliary_coordinate(i, R_right, right_image_);
		// 计算点投影系数
		double N1 = (Bx * right_auxiliary_coordinates.at<double>(2, 0) - Bz * right_auxiliary_coordinates.at<double>(0, 0)) / (left_auxiliary_coordinates.at<double>(0, 0) * right_auxiliary_coordinates.at<double>(2, 0) - left_auxiliary_coordinates.at<double>(2, 0) * right_auxiliary_coordinates.at<double>(0, 0));
		double N2 = (Bx * left_auxiliary_coordinates.at<double>(2, 0) - Bz * left_auxiliary_coordinates.at<double>(0, 0)) / (left_auxiliary_coordinates.at<double>(0, 0) * right_auxiliary_coordinates.at<double>(2, 0) - left_auxiliary_coordinates.at<double>(2, 0) * right_auxiliary_coordinates.at<double>(0, 0));
		// 计算地面坐标
		ControlPoint control_point = N1N2_intersection(i, left_image_, N1, left_auxiliary_coordinates);
		vec_control_point.push_back(control_point);
	}
	////Output result
	cout << "--------------------------------------------" << endl;
	cout << "Space Intersection result of point factor N1 N2 method" << endl;
	for (int i = 0; i < vec_control_point.size(); i++)
	{
		cout << fixed << setprecision(5) << vec_control_point[i].id_ << " " << vec_control_point[i].x_ << " " << vec_control_point[i].y_ << " " << vec_control_point[i].z_ << endl;
	}
	cout << endl;
	ofstream outfile;
	outfile.open(reslut_file_path, ios::out);
	for (int i = 0; i < vec_control_point.size(); i++)
	{
		outfile << fixed << setprecision(5) << vec_control_point[i].id_ << " " << vec_control_point[i].x_ << " " << vec_control_point[i].y_ << " " << vec_control_point[i].z_ << endl;
	}
}
