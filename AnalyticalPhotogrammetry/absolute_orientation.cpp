#include "absolute_orientation.h"

void AbsoluteOrientation::get_param_from_file(const char* image_path)
{
	FILE* fp = fopen(image_path, "r");
	char buffer[512] = { 0 };
	while (!feof(fp))
	{
		fgets(buffer, 512, fp);
		int id;
		double Xp, Yp, Zp, Xtp, Ytp, Ztp;
		sscanf(buffer, "%d %lf %lf %lf %lf %lf %lf\n", &id, &Xp, &Yp, &Zp, &Xtp, &Ytp, &Ztp);
		model_points_.push_back(SpaceCoordinates(id, Xp, Yp, Zp));
		control_points_.push_back(SpaceCoordinates(id, Xtp, Ytp, Ztp));
	}
	fclose(fp);
}

void AbsoluteOrientation::calculate_absolute_rotation_matrix(Mat_<double>& R)
{
	double phi = absolute_orientation_elements_.at<double>(4, 0);
	double omega = absolute_orientation_elements_.at<double>(5, 0);
	double kappa = absolute_orientation_elements_.at<double>(6, 0);
	
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

void AbsoluteOrientation::calculate_barycentric_coordinate(vector<SpaceCoordinates>& barycentric_coordinate, vector<SpaceCoordinates> original_coordinate, SpaceCoordinates& gravity_center)
{
	int point_num = original_coordinate.size();
	// 计算重心坐标
	double  x_center = 0.0;
	double  y_center = 0.0;
	double  z_center = 0.0;
	for (int i = 0; i < point_num; i++)
	{
		x_center += original_coordinate[i].x_;
		y_center += original_coordinate[i].y_;
		z_center += original_coordinate[i].z_;
	}
	x_center /= point_num;
	y_center /= point_num;
	z_center /= point_num;
	gravity_center.x_ = x_center;
	gravity_center.y_ = y_center;
	gravity_center.z_ = z_center;
	// 重心化
	for (int i = 0; i < point_num; i++)
	{
		double baricentric_x = original_coordinate[i].x_ - x_center;
		double baricentric_y = original_coordinate[i].y_ - y_center;
		double baricentric_z = original_coordinate[i].z_ - z_center;

		barycentric_coordinate.push_back(SpaceCoordinates(baricentric_x, baricentric_y, baricentric_z));
	}
}

void AbsoluteOrientation::calculate_A_matrix(int i, Mat_<double>& A, Mat_<double> model_barycentric_coordinate_rotated)
{	
	double lambda = absolute_orientation_elements_.at<double>(3, 0);
	double phi = absolute_orientation_elements_.at<double>(4, 0);
	double omega = absolute_orientation_elements_.at<double>(5, 0);
	double kappa = absolute_orientation_elements_.at<double>(6, 0);
	// 旋转后的重心化坐标X'_=R•X_
	double X_rotated = model_barycentric_coordinate_rotated.at<double>(0, 0);
	double Y_rotated = model_barycentric_coordinate_rotated.at<double>(1, 0);
	double Z_rotated = model_barycentric_coordinate_rotated.at<double>(2, 0);

	// 计算偏导数
	Mat_<double> partial_derivative_matrix = Mat::zeros(3, 7, CV_32F);
	partial_derivative_matrix.at<double>(0, 0) = 1;// X X
	partial_derivative_matrix.at<double>(0, 1) = 0;// X Y
	partial_derivative_matrix.at<double>(0, 2) = 0;// X Z
	partial_derivative_matrix.at<double>(0, 3) = X_rotated;// X/λ
	partial_derivative_matrix.at<double>(0, 4) = -lambda * Z_rotated; // X Φ
	partial_derivative_matrix.at<double>(0, 5) = -lambda * Y_rotated * sin(phi); // X Ω
	partial_derivative_matrix.at<double>(0, 6) = -lambda * Y_rotated * cos(phi) * cos(omega) - lambda * Z_rotated * sin(omega); // X κ
	partial_derivative_matrix.at<double>(1, 0) = 0;// Y X
	partial_derivative_matrix.at<double>(1, 1) = 1;// Y Y
	partial_derivative_matrix.at<double>(1, 2) = 0;// Y Z
	partial_derivative_matrix.at<double>(1, 3) = Y_rotated;// Y λ
	partial_derivative_matrix.at<double>(1, 4) = 0; // Y Φ
	partial_derivative_matrix.at<double>(1, 5) = lambda * X_rotated * sin(phi) - lambda * Z_rotated * cos(phi); //Y Ω
	partial_derivative_matrix.at<double>(1, 6) = lambda * X_rotated * cos(phi) * cos(omega) + lambda * Z_rotated * sin(phi) * cos(omega); // Y κ
	partial_derivative_matrix.at<double>(2, 0) = 0;// Z X
	partial_derivative_matrix.at<double>(2, 1) = 0;// Z Y
	partial_derivative_matrix.at<double>(2, 2) = 1;// Z Z
	partial_derivative_matrix.at<double>(2, 3) = Z_rotated;
	partial_derivative_matrix.at<double>(2, 4) = lambda * X_rotated; // Z Φ
	partial_derivative_matrix.at<double>(2, 5) = lambda * Y_rotated * cos(phi); // Z Ω
	partial_derivative_matrix.at<double>(2, 6) = lambda * X_rotated * sin(omega) - lambda * Y_rotated * sin(phi) * cos(omega); // Z κ
	for (int j = 0; j < 7; j++)
	{
		A.at<double>(3 * i, j) = partial_derivative_matrix.at<double>(0, j);
		A.at<double>(3 * i + 1, j) = partial_derivative_matrix.at<double>(1, j);
		A.at<double>(3 * i + 2, j) = partial_derivative_matrix.at<double>(2, j);
	}
}

void AbsoluteOrientation::calculate_L_matrix(int i, Mat_<double>& L, Mat_<double> control_barycentric_coordinate_matrix, Mat_<double> model_barycentric_coordinate_rotated)
{
	// (ΔX，ΔY，ΔZ)
	Mat_<double> delta = Mat::zeros(3, 1, CV_32F);
	delta.at<double>(0, 0) = absolute_orientation_elements_.at<double>(0, 0);
	delta.at<double>(1, 0) = absolute_orientation_elements_.at<double>(1, 0);
	delta.at<double>(2, 0) = absolute_orientation_elements_.at<double>(2, 0);

	double lambda = absolute_orientation_elements_.at<double>(3, 0);

	Mat_<double> L_i = Mat::zeros(3, 1, CV_32F);
	// L_X=Xtp_-λX'_-ΔX
	L_i = control_barycentric_coordinate_matrix - lambda * model_barycentric_coordinate_rotated - delta;

	L.at<double>(3 * i, 0) = L_i.at<double>(0, 0);
	L.at<double>(3 * i + 1, 0) = L_i.at<double>(1, 0);
	L.at<double>(3 * i + 2, 0) = L_i.at<double>(2, 0);
}

void AbsoluteOrientation::correct_absolute_param(Mat_<double> X)
{
	for (int i = 0; i < 7; i++)
	{
		absolute_orientation_elements_.at<double>(i, 0) += X.at<double>(i, 0);
	}
}

bool AbsoluteOrientation::is_tolerant(Mat_<double> X, double tolerance)
{
	if (fabs(X.at<double>(0, 0)) < tolerance && fabs(X.at<double>(1, 0)) < tolerance && fabs(X.at<double>(2, 0)) < tolerance && fabs(X.at<double>(3, 0)) < tolerance && fabs(X.at<double>(4, 0)) < tolerance && fabs(X.at<double>(5, 0)) < tolerance && fabs(X.at<double>(6, 0)) < tolerance)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void AbsoluteOrientation::calculate_absolute_orientation(const char* file_path, const char* result_file_path)
{
	int iteration = 0;
	double tolerance = 3e-5;
	// 初始化参数
	absolute_orientation_elements_.at<double>(3, 0) = 1;// λ=1
	// 读入模型点坐标和相应控制点坐标
	get_param_from_file(file_path);

	int point_num = model_points_.size();

	Mat_<double> A = Mat::zeros(3 * point_num, 7, CV_32F);
	Mat_<double> L = Mat::zeros(3 * point_num, 1, CV_32F);
	Mat_<double> X = Mat::zeros(7, 1, CV_32F);

	vector<SpaceCoordinates> model_barycentric_coordinate;//模型点重心化坐标
	vector<SpaceCoordinates> control_barycentric_coordinate;//控制点重心化坐标

	SpaceCoordinates model_gravity_center, control_gravity_center;//模型点、控制点重心坐标
	calculate_barycentric_coordinate(model_barycentric_coordinate, model_points_, model_gravity_center);
	calculate_barycentric_coordinate(control_barycentric_coordinate, control_points_, control_gravity_center);
	do
	{
		Mat_<double> R = Mat::zeros(3, 3, CV_32F);
		calculate_absolute_rotation_matrix(R);
		for (int i = 0; i < point_num; i++)
		{
			// 矩阵形式的重心化模型点坐标Xp_
			Mat_<double> model_barycentric_coordinate_matrix = Mat::zeros(3, 1, CV_32F);
			model_barycentric_coordinate_matrix.at<double>(0, 0) = model_barycentric_coordinate[i].x_;
			model_barycentric_coordinate_matrix.at<double>(1, 0) = model_barycentric_coordinate[i].y_;
			model_barycentric_coordinate_matrix.at<double>(2, 0) = model_barycentric_coordinate[i].z_;
			// 矩阵形式的重心化控制点坐标Xtp_
			Mat_<double> control_barycentric_coordinate_matrix = Mat::zeros(3, 1, CV_32F);
			control_barycentric_coordinate_matrix.at<double>(0, 0) = control_barycentric_coordinate[i].x_;
			control_barycentric_coordinate_matrix.at<double>(1, 0) = control_barycentric_coordinate[i].y_;
			control_barycentric_coordinate_matrix.at<double>(2, 0) = control_barycentric_coordinate[i].z_;
			// 旋转之后的重心化模型点坐X'p_
			Mat_<double> model_barycentric_coordinate_rotated = Mat::zeros(3, 1, CV_32F);
			// X'p_=R•Xp
			model_barycentric_coordinate_rotated = R * model_barycentric_coordinate_matrix;

			calculate_A_matrix(i, A, model_barycentric_coordinate_rotated);
			calculate_L_matrix(i, L, control_barycentric_coordinate_matrix, model_barycentric_coordinate_rotated);
		}

		X = (A.t() * A).inv() * A.t() * L;
		correct_absolute_param(X);
		iteration += 1;
	} while (!is_tolerant(X, tolerance));

	Mat_<double>V = A * X - L;

	Mat_<double>V_ = V.t() * V;
	double accuracy = sqrt(V_.at<double>(0, 0) / (3 * point_num - 7));

	cout << "--------------------------------------------" << endl;
	cout << "Absolute Orientation Result: " << endl;
	cout << "Iteration: " << iteration << endl;
	cout << "Residual:" << endl;
	cout << "Seven Parameters of Relative Orientation(ΔX ΔY ΔZ λ Φ Ω κ): " << endl;
	cout << absolute_orientation_elements_.at<double>(0, 0) + control_gravity_center.x_ << "\t" << absolute_orientation_elements_.at<double>(1, 0) + control_gravity_center.y_ << "\t" << absolute_orientation_elements_.at<double>(2, 0) + control_gravity_center.z_ << "\t" << absolute_orientation_elements_.at<double>(3, 0)
		<< "\t" << absolute_orientation_elements_.at<double>(4, 0) << "\t" << absolute_orientation_elements_.at<double>(5, 0) << "\t" << absolute_orientation_elements_.at<double>(6, 0) << endl;
	cout << "RMS Error：" << accuracy << endl;

	ofstream outfile;
	outfile.open(result_file_path, ios::out);
	outfile << "--------------------------------------------" << endl;
	outfile << "Absolute Orientation Result: " << endl;
	outfile << "Iteration: " << iteration << endl;
	outfile << "Residual:" << endl;
	outfile << V << endl;
	outfile << "Seven Parameters of Relative Orientation(ΔX ΔY ΔZ λ Φ Ω κ): " << endl;
	outfile << absolute_orientation_elements_.at<double>(0, 0) + control_gravity_center.x_ << "\t" << absolute_orientation_elements_.at<double>(1, 0) + control_gravity_center.y_ << "\t" << absolute_orientation_elements_.at<double>(2, 0) + control_gravity_center.z_ << "\t" << absolute_orientation_elements_.at<double>(3, 0)
		<< "\t" << absolute_orientation_elements_.at<double>(4, 0) << "\t" << absolute_orientation_elements_.at<double>(5, 0) << "\t" << absolute_orientation_elements_.at<double>(6, 0) << endl;

	outfile << "RMS Error：" << accuracy << endl;

	Mat_<double> R_final = Mat::zeros(3, 3, CV_32F);// 旋转矩阵
	Mat_<double> translation_quality = Mat::zeros(3, 1, CV_32F);// 平移量
	Mat_<double> control_gravity_center_matrix = Mat::zeros(3, 1, CV_32F);
	Mat_<double> result = Mat::zeros(3, 1, CV_32F);

	calculate_absolute_rotation_matrix(R_final);
	translation_quality.at<double>(0, 0) = absolute_orientation_elements_.at<double>(0, 0);
	translation_quality.at<double>(1, 0) = absolute_orientation_elements_.at<double>(1, 0);
	translation_quality.at<double>(2, 0) = absolute_orientation_elements_.at<double>(2, 0);

	control_gravity_center_matrix.at<double>(0, 0) = control_gravity_center.x_ ;
	control_gravity_center_matrix.at<double>(1, 0) = control_gravity_center.y_ ;
	control_gravity_center_matrix.at<double>(2, 0) = control_gravity_center.z_ ;

	cout << "Coordinate of ground points: " << endl;
	outfile << "Coordinate of ground points: " << endl;
	for (int i = 0; i < point_num; i++)
	{
		Mat_<double> model_point_matrix = Mat::zeros(3, 1, CV_32F);
		model_point_matrix.at<double>(0, 0) = model_barycentric_coordinate[i].x_;
		model_point_matrix.at<double>(1, 0) = model_barycentric_coordinate[i].y_;
		model_point_matrix.at<double>(2, 0) = model_barycentric_coordinate[i].z_;
		// Xtp = λRXp + Xs
		result = absolute_orientation_elements_.at<double>(3, 0) * R_final * model_point_matrix + translation_quality + control_gravity_center_matrix;
		cout << fixed << setprecision(6) << model_points_[i].id_ << "\t" << result.at<double>(0, 0) << "\t" << result.at<double>(1, 0) << "\t" << result.at<double>(2, 0) << "\t" << control_points_[i].x_ - result.at<double>(0, 0) << "\t" << control_points_[i].y_ - result.at<double>(1, 0) << "\t" << control_points_[i].z_ - result.at<double>(2, 0) << endl;
		outfile << fixed << setprecision(6) << model_points_[i].id_ << "\t" << result.at<double>(0, 0) << "\t" << result.at<double>(1, 0) << "\t" << result.at<double>(2, 0) << "\t" << control_points_[i].x_ - result.at<double>(0, 0) << "\t" << control_points_[i].y_ - result.at<double>(1, 0) << "\t" << control_points_[i].z_ - result.at<double>(2, 0) << endl;
	}
}
