#include "space_resection.h"

void SpaceResection::get_param_from_file(const char* camera_file_name, const char* point_file_name)
{
	// ��ȡ����������Ӧ���������
	FILE* fp_point = fopen(point_file_name, "r");
	char buffer[512] = { 0 };
	while (!feof(fp_point))
	{
		fgets(buffer, 512, fp_point);
		int point_id;
		double x=0., y=0., X=0., Y=0., Z=0.;
		PlaneCoordinates tmp_image_point;
		SpaceCoordinates tmp_control_point;
		sscanf(buffer, "%d %lf %lf %lf %lf %lf\n", &point_id, &x, &y, &X, &Y, &Z);
		tmp_image_point = PlaneCoordinates(point_id, x / 1000.0, y / 1000.0);// ��λת����mm��m
		tmp_control_point = SpaceCoordinates(point_id, X, Y, Z);
		image_points_.push_back(tmp_image_point);
		control_points_.push_back(tmp_control_point);
	}
	fclose(fp_point);
	// ��ȡ����ڲ�
	FILE* fp_camera = fopen(camera_file_name, "r");
	char tmp_char[512]{};
	fgets(tmp_char, 512, fp_camera);
	sscanf(tmp_char, "%lf %lf %lf\n", &camera_.x0_, &camera_.y0_, &camera_.f_);
	camera_.f_ /= 1000.0;// ��λת����mm��m
	fclose(fp_camera);
}

void SpaceResection::Initialize()
{
	int point_num = control_points_.size();
	for (vector<SpaceCoordinates>::const_iterator it = control_points_.cbegin(); it != control_points_.cend(); it++)
	{
		exterior_orientation_elements_.at<double>(0, 0) += (*it).x_;
		exterior_orientation_elements_.at<double>(1, 0) += (*it).y_;
	}
	exterior_orientation_elements_.at<double>(0, 0) /= point_num;// Xs
	exterior_orientation_elements_.at<double>(1, 0) /= point_num;// Ys
	exterior_orientation_elements_.at<double>(2, 0) = m_ * camera_.f_;// Zs
	exterior_orientation_elements_.at<double>(3, 0) = 0;// ��
	exterior_orientation_elements_.at<double>(4, 0) = 0;// ��
	exterior_orientation_elements_.at<double>(5, 0) = 0;// ��
}

void SpaceResection::calculate_rotation_matrix(Mat_<double>& R)
{
	double phi = exterior_orientation_elements_.at<double>(3, 0);
	double omega = exterior_orientation_elements_.at<double>(4, 0);
	double kappa = exterior_orientation_elements_.at<double>(5, 0);

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

void SpaceResection::calculate_A_matrix(int i, Mat_<double>R, Mat_<double>& A)
{
	double Xs = exterior_orientation_elements_.at<double>(0, 0);
	double Ys = exterior_orientation_elements_.at<double>(1, 0);
	double Zs = exterior_orientation_elements_.at<double>(2, 0);
	double phi = exterior_orientation_elements_.at<double>(3, 0);
	double omega = exterior_orientation_elements_.at<double>(4, 0);
	double kappa = exterior_orientation_elements_.at<double>(5, 0);

	double Z_bar = R.at<double>(0, 2) * (control_points_[i].x_ - Xs) + R.at<double>(1, 2) * (control_points_[i].y_ - Ys) + R.at<double>(2, 2) * (control_points_[i].z_ - Zs);

	A.at<double>(i * 2, 0) = (R.at<double>(0, 0) * camera_.f_ + R.at<double>(0, 2) * (image_points_[i].x_ - camera_.x0_)) / Z_bar;//a11
	A.at<double>(i * 2, 1) = (R.at<double>(1, 0) * camera_.f_ + R.at<double>(1, 2) * (image_points_[i].x_ - camera_.x0_)) / Z_bar;//a12
	A.at<double>(i * 2, 2) = (R.at<double>(2, 0) * camera_.f_ + R.at<double>(2, 2) * (image_points_[i].x_ - camera_.x0_)) / Z_bar;//a13
	A.at<double>(i * 2, 3) = (image_points_[i].y_ - camera_.y0_) * sin(omega) - ((image_points_[i].x_ - camera_.x0_) * ((image_points_[i].x_ - camera_.x0_) * cos(kappa) - (image_points_[i].y_ - camera_.y0_) * sin(kappa)) / camera_.f_ + camera_.f_ * cos(kappa)) * cos(omega);//a14
	A.at<double>(i * 2, 4) = -camera_.f_ * sin(kappa) - (image_points_[i].x_ - camera_.x0_) * ((image_points_[i].x_ - camera_.x0_) * sin(kappa) + (image_points_[i].y_ - camera_.y0_) * cos(kappa)) / camera_.f_;//a15
	A.at<double>(i * 2, 5) = image_points_[i].y_ - camera_.y0_;//a16
	A.at<double>(i * 2 + 1, 0) = (R.at<double>(0, 1) * camera_.f_ + R.at<double>(0, 2) * (image_points_[i].y_ - camera_.y0_)) / Z_bar;//a21
	A.at<double>(i * 2 + 1, 1) = (R.at<double>(1, 1) * camera_.f_ + R.at<double>(1, 2) * (image_points_[i].y_ - camera_.y0_)) / Z_bar;//a22
	A.at<double>(i * 2 + 1, 2) = (R.at<double>(2, 1) * camera_.f_ + R.at<double>(2, 2) * (image_points_[i].y_ - camera_.y0_)) / Z_bar;//a23
	A.at<double>(i * 2 + 1, 3) = -(image_points_[i].x_ - camera_.x0_) * sin(omega) - ((image_points_[i].y_ - camera_.y0_) * ((image_points_[i].x_ - camera_.x0_) * cos(kappa) - (image_points_[i].y_ - camera_.y0_) * sin(kappa)) / camera_.f_ - camera_.f_ * sin(kappa)) * cos(omega);//a24
	A.at<double>(i * 2 + 1, 4) = -camera_.f_ * cos(kappa) - ((image_points_[i].y_ - camera_.y0_) * ((image_points_[i].x_ - camera_.x0_) * sin(kappa) + (image_points_[i].y_ - camera_.y0_) * cos(kappa))) / camera_.f_;//a25
	A.at<double>(i * 2 + 1, 5) = -(image_points_[i].x_ - camera_.x0_);//a26
}

void SpaceResection::calculate_L_matrix(int i, Mat_<double>R, Mat_<double>& L)
{
	double Xs = exterior_orientation_elements_.at<double>(0, 0);
	double Ys = exterior_orientation_elements_.at<double>(1, 0);
	double Zs = exterior_orientation_elements_.at<double>(2, 0);
	// ����������ֵ
	double approxiate_image_point_x = camera_.x0_ - camera_.f_ * (R.at<double>(0, 0) * (control_points_[i].x_ - Xs) + R.at<double>(1, 0) * (control_points_[i].y_ - Ys) + R.at<double>(2, 0) * (control_points_[i].z_ - Zs))
		/ (R.at<double>(0, 2) * (control_points_[i].x_ - Xs) + R.at<double>(1, 2) * (control_points_[i].y_ - Ys) + R.at<double>(2, 2) * (control_points_[i].z_ - Zs));
	double approxiate_image_point_y = camera_.y0_ - camera_.f_ * (R.at<double>(0, 1) * (control_points_[i].x_ - Xs) + R.at<double>(1, 1) * (control_points_[i].y_ - Ys) + R.at<double>(2, 1) * (control_points_[i].z_ - Zs))
		/ (R.at<double>(0, 2) * (control_points_[i].x_ - Xs) + R.at<double>(1, 2) * (control_points_[i].y_ - Ys) + R.at<double>(2, 2) * (control_points_[i].z_ - Zs));

	L.at<double>(i * 2, 0) = image_points_[i].x_ - approxiate_image_point_x;
	L.at<double>(i * 2 + 1, 0) = image_points_[i].y_ - approxiate_image_point_y;
}

void SpaceResection::correct_exterior_orientation_elements(Mat_<double> X)
{
	for (int i = 0; i < 6; i++)
	{
		exterior_orientation_elements_.at<double>(i, 0) += X.at<double>(i, 0);
	}
}

bool SpaceResection::is_tolerant(Mat_<double> X, double tolerance)
{
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
	int iteration = 0;// ��������
	double tolerance = 3e-5;// �޲�
	Mat_<double> R = Mat::zeros(3, 3, CV_32F);// ��ת����
	Mat_<double> X = Mat::zeros(6, 1, CV_32F);// δ֪����Xs����Ys����Zs�����գ����أ�����

	// ��ȡ����ڲμ����Ϳ��Ƶ�����
	SpaceResection::get_param_from_file(camera_file_name, point_file_name);
    // ��ʼ������
	Initialize();
	int point_num = image_points_.size();// ���/���Ƶ����
	// ��ʼ��ϵ����A�볣����ϵ����L
	Mat_<double> A = Mat::zeros(point_num * 2, 6, CV_32F);
	Mat_<double> L = Mat::zeros(point_num * 2, 1, CV_32F);

	// ��������
	do
	{
		// ������ת����R
		calculate_rotation_matrix(R);
		for (int i = 0; i < point_num; i++)
		{
			// ����������ϵ���ͳ�����
			calculate_A_matrix(i, R, A);
			calculate_L_matrix(i, R, L);
		}
		// ���ݷ��������̼����ⷽλԪ�ؽ���ֵ������
		X = (A.t() * A).inv() * A.t() * L;
		// ����δ֪��
		correct_exterior_orientation_elements(X);
		iteration += 1;
	} while (!is_tolerant(X, tolerance));
	// ��������
	// �����������۲�ֵ������
	Mat_<double> V = A * X - L;
	Mat_<double> V_ = V.t() * V;
	double accuracy = sqrt(V_.at<double>(0, 0) / (point_num * 2 - 6));

	////Output final result
	cout << "--------------------------------------------" << endl;
	cout << "Space Resection Result" << endl;
	cout << "Iteration: " << iteration << endl;
	cout << "Xs = " << fixed << setprecision(2) << exterior_orientation_elements_.at<double>(0, 0) << endl;
	cout << "Ys = " << fixed << setprecision(2) << exterior_orientation_elements_.at<double>(1, 0) << endl;
	cout << "Zs = " << fixed << setprecision(2) << exterior_orientation_elements_.at<double>(2, 0) << endl;
	cout << "�� = " << fixed << setprecision(5) << exterior_orientation_elements_.at<double>(3, 0) << endl;
	cout << "�� = " << fixed << setprecision(5) << exterior_orientation_elements_.at<double>(4, 0) << endl;
	cout << "�� = " << fixed << setprecision(5) << exterior_orientation_elements_.at<double>(5, 0) << endl;
	cout << "Rotation Matrix:" << endl;
	cout << fixed << setprecision(5) << R << endl;
	cout << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
	cout << endl;

	ofstream outfile;
	outfile.open(result_file_name, ios::out);
	outfile << "--------------------------------------------" << endl;
	outfile << "Space Resection Result" << endl;
	outfile << "Iteration: " << iteration << endl;
	outfile << "Xs = " << fixed << setprecision(2) << exterior_orientation_elements_.at<double>(0, 0) << endl;
	outfile << "Ys = " << fixed << setprecision(2) << exterior_orientation_elements_.at<double>(1, 0) << endl;
	outfile << "Zs = " << fixed << setprecision(2) << exterior_orientation_elements_.at<double>(2, 0) << endl;
	outfile << "�� = " << fixed << setprecision(5) << exterior_orientation_elements_.at<double>(3, 0) << endl;
	outfile << "�� = " << fixed << setprecision(5) << exterior_orientation_elements_.at<double>(4, 0) << endl;
	outfile << "�� = " << fixed << setprecision(5) << exterior_orientation_elements_.at<double>(5, 0) << endl;
	outfile << "Rotation Matrix:" << endl;
	outfile << fixed << setprecision(5) << R << endl;
	outfile << "RMS Error: " << fixed << setprecision(5) << accuracy << endl;
	outfile << endl;
}


