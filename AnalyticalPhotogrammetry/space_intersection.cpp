#include "space_intersection.h"

void SpaceIntersection::get_param_from_file(const char* file_path, double& m, Camera& interior_elements, Mat_<double>& exterior_elements, vector<PlaneCoordinates>& image_points)
{
	FILE* fp = fopen(file_path, "r");
	char buffer[512] = { 0 };
	// 读入注释行
	fgets(buffer, 512, fp);
	fgets(buffer, 512, fp);

	double f = 0.0, x0 = 0.0, y0 = 0.0;
	fgets(buffer, 512, fp);
	sscanf(buffer, "%lf %lf %lf %lf\n", &f, &x0, &y0, &m);
	interior_elements.f_ = f / 1000.0;
	interior_elements.x0_ = x0 / 1000.0;
	interior_elements.y0_ = y0 / 1000.0;//mm→m

	// 读入注释行
	fgets(buffer, 512, fp);
	fgets(buffer, 512, fp);

	double phi = 0.0, omega = 0.0, kappa = 0.0, Xs = 0.0, Ys = 0.0, Zs = 0.0;
	fgets(buffer, 512, fp);
	sscanf(buffer, "%lf %lf %lf %lf %lf %lf\n", &phi, &omega, &kappa, &Xs, &Ys, &Zs);
	exterior_elements.at<double>(0, 0) = Xs;
	exterior_elements.at<double>(1, 0) = Ys;
	exterior_elements.at<double>(2, 0) = Zs;
	exterior_elements.at<double>(3, 0) = phi * CV_PI / 180.0;
	exterior_elements.at<double>(4, 0) = omega * CV_PI / 180.0;
	exterior_elements.at<double>(5, 0) = kappa * CV_PI / 180.0;// 度→rad
	// 读入注释行
	fgets(buffer, 512, fp);
	fgets(buffer, 512, fp);
	while (!feof(fp))
	{
		fgets(buffer, 512, fp);
		int id = 0;
		double x=0.0, y=0.0;
		sscanf(buffer, "%d %lf %lf", &id, &x, &y);
		image_points.push_back(PlaneCoordinates(id, x / 1000.0, y / 1000.0));
	}
	fclose(fp);
}

void SpaceIntersection::calculate_rotation_matrix(Mat_<double>& R, Mat_<double> exterior_elements)
{
	double phi = exterior_elements.at<double>(3, 0);
	double omega = exterior_elements.at<double>(4, 0);
	double kappa = exterior_elements.at<double>(5, 0);

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

Mat SpaceIntersection::calculate_auxiliary_coordinate(int i, Mat_<double>R, Camera interior_elements, vector<PlaneCoordinates> image_points)
{
	Mat_<double> image_coordinate(3, 1);// 像空间坐标
	image_coordinate.at<double>(0, 0) = image_points[i].x_;
	image_coordinate.at<double>(1, 0) = image_points[i].y_;
	image_coordinate.at<double>(2, 0) = -interior_elements.f_;
	Mat_<double> auxiliary_coordinates = R * image_coordinate;// 像空间辅助坐标
	return auxiliary_coordinates;
}

SpaceCoordinates SpaceIntersection::N1N2_intersection(int i, double N1, double N2, double By, Mat_<double> left_auxiliary_coordinates, Mat_<double> right_auxiliary_coordinates)
{
	double Xs1 = left_image_exterior_elements_.at<double>(0, 0);
	double Ys1 = left_image_exterior_elements_.at<double>(1, 0);
	double Zs1 = left_image_exterior_elements_.at<double>(2, 0);

	double Xs2 = right_image_exterior_elements_.at<double>(0, 0);
	double Ys2 = right_image_exterior_elements_.at<double>(1, 0);
	double Zs2 = right_image_exterior_elements_.at<double>(2, 0);


	double X1 = left_auxiliary_coordinates.at<double>(0, 0);
	double Y1 = left_auxiliary_coordinates.at<double>(1, 0);
	double Z1 = left_auxiliary_coordinates.at<double>(2, 0);

	double X2 = right_auxiliary_coordinates.at<double>(0, 0);
	double Y2 = right_auxiliary_coordinates.at<double>(1, 0);
	double Z2 = right_auxiliary_coordinates.at<double>(2, 0);

	double X = 0.5 * (Xs1 + N1 * X1 + Xs2 + N2 * X2);
	double Y = 0.5 * (Ys2 + N1 * Y1 + Ys2 + N2 * Y2);
	double Z = 0.5 * (Zs1 + N1 * Z1 + Zs2 + N2 * Z2);

	return SpaceCoordinates(left_image_points_[i].id_, X, Y, Z);
}

void SpaceIntersection::pointfactor_space_intersection(const char* left_image_path, const char* right_image_path, const char* reslut_file_path)
{
	vector<SpaceCoordinates> vec_control_points;
	// 获取左右像片已知数据x0,y0,f,Xs,Ys,Zs,φ,ω,κ
	get_param_from_file(left_image_path, left_m_, left_image_interior_elements_, left_image_exterior_elements_, left_image_points_);
	get_param_from_file(right_image_path, right_m_, right_image_interior_elements_, right_image_exterior_elements_, right_image_points_);
	// 计算基线分量
	double Bx = right_image_exterior_elements_.at<double>(0, 0) - left_image_exterior_elements_.at<double>(0, 0);
	double By = right_image_exterior_elements_.at<double>(1, 0) - left_image_exterior_elements_.at<double>(1, 0);
	double Bz = right_image_exterior_elements_.at<double>(2, 0) - left_image_exterior_elements_.at<double>(2, 0);
	// 计算每张像片的旋转矩阵
	Mat_<double> R_left = Mat::zeros(3, 3, CV_32F);
	Mat_<double> R_right = Mat::zeros(3, 3, CV_32F);
	calculate_rotation_matrix(R_left, left_image_exterior_elements_);
	calculate_rotation_matrix(R_right, right_image_exterior_elements_);

	int point_num = left_image_points_.size();
	for (int i = 0; i < point_num; i++) 
	{
		// 由像空间坐标计算像空间辅助坐标
		Mat_<double> left_auxiliary_coordinates = calculate_auxiliary_coordinate(i, R_left, left_image_interior_elements_, left_image_points_);
		Mat_<double> right_auxiliary_coordinates = calculate_auxiliary_coordinate(i, R_right, right_image_interior_elements_, right_image_points_);
		// 计算点投影系数
		double N1 = (Bx * right_auxiliary_coordinates.at<double>(2, 0) - Bz * right_auxiliary_coordinates.at<double>(0, 0)) / (left_auxiliary_coordinates.at<double>(0, 0) * right_auxiliary_coordinates.at<double>(2, 0) - left_auxiliary_coordinates.at<double>(2, 0) * right_auxiliary_coordinates.at<double>(0, 0));
		double N2 = (Bx * left_auxiliary_coordinates.at<double>(2, 0) - Bz * left_auxiliary_coordinates.at<double>(0, 0)) / (left_auxiliary_coordinates.at<double>(0, 0) * right_auxiliary_coordinates.at<double>(2, 0) - left_auxiliary_coordinates.at<double>(2, 0) * right_auxiliary_coordinates.at<double>(0, 0));
		// 计算地面坐标
		SpaceCoordinates control_point = N1N2_intersection(i, N1, N2, Bx, left_auxiliary_coordinates, right_auxiliary_coordinates);
		vec_control_points.push_back(control_point);
	}
	cout << "--------------------------------------------" << endl;
	cout << "Space Intersection result of point factor N1 N2 method" << endl;
	for (int i = 0; i < vec_control_points.size(); i++)
	{
		cout << fixed << setprecision(5) << vec_control_points[i].id_ << " " << vec_control_points[i].x_ << " " << vec_control_points[i].y_ << " " << vec_control_points[i].z_ << endl;
	}
	cout << endl;
	ofstream outfile;
	outfile.open(reslut_file_path, ios::out);
	for (int i = 0; i < vec_control_points.size(); i++)
	{
		outfile << fixed << setprecision(5) << vec_control_points[i].id_ << " " << vec_control_points[i].x_ << " " << vec_control_points[i].y_ << " " << vec_control_points[i].z_ << endl;
	}
}
