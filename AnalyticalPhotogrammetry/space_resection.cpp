#include "space_resection.h"

void SpaceResection::get_param_from_file(const char* camera_file_name, const char* point_file_name, Camera& camera, vector<ImagePoint>& image_point, vector<ControlPoint>& control_point)
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
		tmp_image_point = ImagePoint(x / 1000.0, y / 1000.0);
		tmp_control_point = ControlPoint(X, Y, Z);
		image_point.push_back(tmp_image_point);
		control_point.push_back(tmp_control_point);
	}
	fclose(fp_point);

	FILE* fp_camera = fopen(camera_file_name, "r");
	char tmp_char[512]{};
	fgets(tmp_char, 512, fp_camera);
	sscanf(tmp_char, "%lf %lf %lf\n", &camera.x0_, &camera.y0_, &camera.f_);
	camera.f_ /= 1000.0;
	fclose(fp_camera);
}

void SpaceResection::Initialize(double m, Camera camera, vector<ControlPoint> control_point, double& Xs, double& Ys, double& Zs, double& phi, double& omega, double& kappa)
{
	int point_num = control_point.size();
	for (vector<ControlPoint>::const_iterator it = control_point.cbegin(); it != control_point.cend(); it++)
	{
		Xs += (*it).x_;
		Ys += (*it).y_;
	}
	Xs /= point_num;
	Ys /= point_num;
	Zs = m * camera.f_;
	phi = 0;
	omega = 0;
	kappa = 0;
}

void SpaceResection::calculate_rotation_matrix(Mat_<double>& R, double phi, double omega, double kappa)
{
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

void SpaceResection::calculate_A_matrix(int i, Camera camera, vector<ControlPoint> control_point, vector<ImagePoint> image_point, double Z, double phi, double omega, double kappa, Mat_<double>R, Mat_<double>& A)
{

	A.at<double>(i * 2, 0) = (R.at<double>(0, 0) * camera.f_ + R.at<double>(0, 2) * (image_point[i].x_ - camera.x0_)) / Z;//a11
	A.at<double>(i * 2, 1) = (R.at<double>(1, 0) * camera.f_ + R.at<double>(1, 2) * (image_point[i].x_ - camera.x0_)) / Z;//a12
	A.at<double>(i * 2, 2) = (R.at<double>(2, 0) * camera.f_ + R.at<double>(2, 2) * (image_point[i].x_ - camera.x0_)) / Z;//a13
	A.at<double>(i * 2, 3) = (image_point[i].y_ - camera.y0_) * sin(omega) - ((image_point[i].x_ - camera.x0_) * ((image_point[i].x_ - camera.x0_) * cos(kappa) - (image_point[i].y_ - camera.y0_) * sin(kappa)) / camera.f_ + camera.f_ * cos(kappa)) * cos(omega);//a14
	A.at<double>(i * 2, 4) = -camera.f_ * sin(kappa) - (image_point[i].x_ - camera.x0_) * ((image_point[i].x_ - camera.x0_) * sin(kappa) + (image_point[i].y_ - camera.y0_) * cos(kappa)) / camera.f_;//a15
	A.at<double>(i * 2, 5) = image_point[i].y_ - camera.y0_;//a16
	A.at<double>(i * 2 + 1, 0) = (R.at<double>(0, 1) * camera.f_ + R.at<double>(0, 2) * (image_point[i].y_ - camera.y0_)) / Z;//a21
	A.at<double>(i * 2 + 1, 1) = (R.at<double>(1, 1) * camera.f_ + R.at<double>(1, 2) * (image_point[i].y_ - camera.y0_)) / Z;//a22
	A.at<double>(i * 2 + 1, 2) = (R.at<double>(2, 1) * camera.f_ + R.at<double>(2, 2) * (image_point[i].y_ - camera.y0_)) / Z;//a23
	A.at<double>(i * 2 + 1, 3) = -(image_point[i].x_ - camera.x0_) * sin(omega) - ((image_point[i].y_ - camera.y0_) * ((image_point[i].x_ - camera.x0_) * cos(kappa) - (image_point[i].y_ - camera.y0_) * sin(kappa)) / camera.f_ - camera.f_ * sin(kappa)) * cos(omega);//a24
	A.at<double>(i * 2 + 1, 4) = -camera.f_ * cos(kappa) - ((image_point[i].y_ - camera.y0_) * ((image_point[i].x_ - camera.x0_) * sin(kappa) + (image_point[i].y_ - camera.y0_) * cos(kappa))) / camera.f_;//a25
	A.at<double>(i * 2 + 1, 5) = -(image_point[i].x_ - camera.x0_);//a26
}

void SpaceResection::calculate_L_matrix(int i, Mat_<double>& L, vector<ImagePoint> image_point, vector<ImagePoint> approxiate_image_point)
{
	L.at<double>(i * 2, 0) = image_point[i].x_ - approxiate_image_point[i].x_;
	L.at<double>(i * 2 + 1, 0) = image_point[i].y_ - approxiate_image_point[i].y_;
}

void SpaceResection::correct_exterior_orientation_elements(double& Xs, double& Ys, double& Zs, double& phi, double& omega, double& kappa, Mat_<double> X)
{
	Xs += X.at<double>(0, 0);
	Ys += X.at<double>(1, 0);
	Zs += X.at<double>(2, 0);
	phi += X.at<double>(3, 0);
	omega += X.at<double>(4, 0);
	kappa += X.at<double>(5, 0);
}

bool SpaceResection::if_tolerant(Mat_<double> X, double tolerance)
{
	cout << X << endl;
	if (fabs(X.at<double>(3,0)) < tolerance && fabs(X.at<double>(4, 0)) < tolerance && fabs(X.at<double>(5, 0)) < tolerance)
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
	int iteration = 0;// 迭代次数
	double tolerance = 0.01;// 限差
	Mat_<double> R = Mat::zeros(3, 3, CV_32F);// 旋转矩阵
	Mat_<double> X = Mat::zeros(6, 1, CV_32F);// 未知数ΔXs，ΔYs，ΔZs，Δφ，Δω，Δκ

	vector<ImagePoint> image_point;// 像点坐标
	vector<ControlPoint> control_point;// 控制点坐标
	Camera camera;// 相机内参x0, y0, f
	double Xs = 0, Ys = 0, Zs = 0, phi = 0, omega = 0, kappa = 0;
	//ExteriorOrientationElements exterior_orientation_elements(0,0,0,0,0,0);// 外方位元素Xs，Ys，Zs，φ，ω，κ
	double m = 50000.0;// 比例尺分母
	// 获取相机内参及像点和控制点坐标
	SpaceResection::get_param_from_file(camera_file_name, point_file_name, camera, image_point, control_point);
    // 初始化参数
	Initialize(m, camera, control_point, Xs, Ys, Zs, phi, omega, kappa);
	int point_num = image_point.size();// 像点/控制点点数
	// 初始化系数阵A与常数项系数阵L
	Mat_<double> A = Mat::zeros(point_num * 2, 6, CV_32F);
	Mat_<double> L = Mat::zeros(point_num * 2, 1, CV_32F);
	vector<ImagePoint> approxiate_image_point(point_num);// 由计算出来的外方位元素计算得到的像点坐标近似值

	//Adjustment process with iterations
	do
	{
		// 计算旋转矩阵R
		calculate_rotation_matrix(R, phi, omega, kappa);
		//计算每个像点的近似值并计算误差方程系数
		for (int i = 0; i < point_num; i++)
		{
			// 像点坐标近似值
			approxiate_image_point[i].x_ = camera.x0_ - camera.f_ * (R.at<double>(0, 0) * (control_point[i].x_ - Xs) + R.at<double>(1, 0) * (control_point[i].y_ - Ys) + R.at<double>(2, 0) * (control_point[i].z_ - Zs))
				/ (R.at<double>(0, 2) * (control_point[i].x_ - Xs) + R.at<double>(1, 2) * (control_point[i].y_ - Ys) + R.at<double>(2, 2) * (control_point[i].z_ - Zs));
			approxiate_image_point[i].y_ = camera.y0_ - camera.f_ * (R.at<double>(0, 1) * (control_point[i].x_ - Xs) + R.at<double>(1, 1) * (control_point[i].y_ - Ys) + R.at<double>(2, 1) * (control_point[i].z_ - Zs))
				/ (R.at<double>(0, 2) * (control_point[i].x_ - Xs) + R.at<double>(1, 2) * (control_point[i].y_ - Ys) + R.at<double>(2, 2) * (control_point[i].z_ - Zs));
			cout << approxiate_image_point[i].x_ << endl;
			cout << approxiate_image_point[i].y_ << endl;
			double Z = R.at<double>(0, 2) * (control_point[i].x_ - Xs) + R.at<double>(1, 2) * (control_point[i].y_ - Ys) + R.at<double>(2, 2) * (control_point[i].z_ - Zs);
			// 误差方程系数
			calculate_A_matrix(i, camera, control_point, image_point, Z, phi, omega, kappa, R, A);
			cout << A << endl;
			calculate_L_matrix(i, L, image_point, approxiate_image_point);
			cout << L << endl;
		}
		// 根据法化法方程计算外方位元素近似值改正数
		X = (A.t() * A).inv() * A.t() * L;
		cout << "X:\n" << X << endl;

		correct_exterior_orientation_elements(Xs, Ys, Zs, phi, omega, kappa, X);
		iteration += 1;
	} while (!if_tolerant(X, tolerance));
	// 计算像点坐标观测值改正数
	Mat_<double> V = A * X - L;
	Mat_<double> V_ = V.t() * V;
	double accuracy = sqrt(V_.at<double>(0, 0) / (point_num * 2 - 6));

	////Output final result
	cout << "--------------------------------------------" << endl;
	cout << "Space Resection Result" << endl;
	cout << "Iteration: " << iteration << endl;
	iteration = 0;
	cout << "Xs = " << fixed << setprecision(2) << Xs << endl;
	cout << "Ys = " << fixed << setprecision(2) << Ys << endl;
	cout << "Zs = " << fixed << setprecision(2) << Zs << endl;
	cout << "fai = " << fixed << setprecision(5) <<phi << endl;
	cout << "omega = " << fixed << setprecision(5) << omega << endl;
	cout << "kappa = " << fixed << setprecision(5) << kappa << endl;
	cout << "Rotation Matrix:" << endl;
	cout << fixed << setprecision(5) << R << endl;
	cout << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
	cout << endl;

	ofstream outfile;
	outfile.open(result_file_name, ios::out);
	outfile << "Xs = " << fixed << setprecision(2) << Xs << endl;
	outfile << "Ys = " << fixed << setprecision(2) << Ys << endl;
	outfile << "Zs = " << fixed << setprecision(2) << Zs << endl;
	outfile << "fai = " << fixed << setprecision(5) << phi << endl;
	outfile << "omega = " << fixed << setprecision(5) << omega << endl;
	outfile << "kappa = " << fixed << setprecision(5) << kappa << endl;
	outfile << "Rotation Matrix:" << endl;
	outfile << fixed << setprecision(5) << R << endl;
	outfile << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
	outfile << endl;

}


