#include"DynamicFormula.h"
DynamicFormula::DynamicFormula(Vector3d omega, Vector3d velocity, Matrix3d R,
	Vector3d y,double delta_t)
{
	this->w = omega;
	this->v = velocity;
	this->R = R; //R如何初始化
	this->y = y;
	this->delta_t = delta_t;
}

DynamicFormula::~DynamicFormula()
{
}
Matrix3d DynamicFormula::toDaOmegaOrY(Vector3d omega) {
	Matrix3d res;
	res(0, 0) = 0;
	res(0, 1) = -omega(2);
	res(0, 2) = omega(1);
	res(1, 0) = omega(2);
	res(1, 1) = 0;
	res(1, 2) = -omega(0);
	res(2, 0) = -omega(1);
	res(2, 1) = omega(0);
	res(2, 2) = 0;
	return res;
}
VectorXd DynamicFormula::tsfs2tf( Matrix3d Y) {
	Matrix3d Rt = R.transpose();
	Matrix3d zero = Matrix3d::Zero();
	//zero.setZero(3, 3);
	Matrix3d negRtY = zero - Rt * Y;
	MatrixXd trans(6, 6);//矩阵分块赋值
	trans.block(0, 0, 3, 3) = Rt;
	trans.block(0, 3, 3, 3) = negRtY;
	trans.block(3, 0, 3, 3) = zero;
	trans.block(3, 3, 3, 3) = Rt;
	/*
	tsfs(0) = ts(0);
	tsfs(1) = ts(1);
	tsfs(2) = ts(2);
	tsfs(3) = fs(0);
	tsfs(4) = fs(1);
	tsfs(5) = fs(2);*/
	return trans * tsfs;	
}
/*
Matrix3d DynamicFormula::computeR_() {
	Matrix3d daOmega = toDaOmegaOrY(w);
	return R * daOmega;
}*/
Vector3d DynamicFormula::computey_() {
	return R * v;
}
VectorXd DynamicFormula::computelp() {
	VectorXd wv(6);
	wv.block(0, 0, 3, 1) = w;
	wv.block(3, 0, 3, 1) = v;
	return K * wv;
}
VectorXd DynamicFormula::computelp_() {
	Matrix3d Y = toDaOmegaOrY(y);
	VectorXd tf = tsfs2tf(Y);
	Vector3d l = vec62Vec31(lp);
	Vector3d p = vec62Vec32(lp);
	Vector3d a = l.cross(w) + p.cross(v);
	Vector3d b = p.cross(w);
	VectorXd ab(6);
	ab.block(0, 0, 3, 1) = a;
	ab.block(3, 0, 3, 1) = b;
	return ab + tf;
}
Matrix3d DynamicFormula::computeNextR() {
	Vector3d Rw = R *w;// 
	//double length = 1;sqrt(Rw(0)*Rw(0) + Rw(1)*Rw(1) + Rw(2)*Rw(2))
	Quaterniond q_w(0, Rw(0)* delta_t, Rw(1)*delta_t, Rw(2)* delta_t);
	q_w = q_w * q;
	q_w.w() = q_w.w()*0.5;
	q_w.x() = q_w.x()*0.5;
	q_w.y() = q_w.y()*0.5;
	q_w.z() = q_w.z()*0.5;
	//cout << "q_w:" << q_w.coeffs() << endl;//计算经过旋转后的orientation 
	q.w() = q_w.w() + q.w();
	q.x() = q_w.x() + q.x();
	q.y() = q_w.y() + q.y();
	q.z() = q_w.z() + q.z();
	q.normalize();
	return q.toRotationMatrix();
}
float* DynamicFormula::GetRotAndTransData() {
	static float data[16];//!!!!
	data[0] = R(0,0);
	data[1] = R(1,0);
	data[2] = R(2,0);
	data[3] = 0;

	data[4] = R(0,1);
	data[5] = R(1,1);
	data[6] = R(2,1);
	data[7] = 0;

	data[8] = R(0,2);
	data[9] = R(1,2);
	data[10] = R(2,2);
	data[11] = 0;

	data[12] = 0;
	data[13] = 0;
	data[14] = 0;
	data[15] = 1;
	return data;
}
Vector3d DynamicFormula::computeNexty( Vector3d y_) {
	Vector3d temp_deltay = delta_t  * y_;
//	cout << "temp_deltay" << temp_deltay << endl;
//	cout << "(" << y(0) << "," << y(1) << "," << y(2) << ")" << endl;
	return y + temp_deltay;
}
VectorXd DynamicFormula::computeNextlp() {
	/*
	Vector3d l = lp.block(0, 0, 3, 1);
	Vector3d p = lp.block(3, 0, 3, 1);
	Vector3d l_ = lp_.block(0, 0, 3, 1);
	Vector3d p_ = lp_.block(3, 0, 3, 1);
	Quaterniond q_l(0, l_(0)*delta_t, l_(1) * delta_t, l_(2) * delta_t);
	q_l = q_l * l;
	q_l.w() = q_l.w()*0.5;
	q_l.x() = q_l.x()*0.5;
	q_l.y() = q_l.y()*0.5;
	q_l.z() = q_l.z()*0.5;
	cout << "q_l:" << q_l.coeffs() << endl;//计算经过旋转后的orientation 
	q.w() = q_l.w() + l.w();
	q.x() = q_l.x() + l.x();
	q.y() = q_l.y() + q.y();
	q.z() = q_l.z() + q.z();
	q.normalize();
	return q.toRotationMatrix();*/
	return lp +  lp_* delta_t ;
}
VectorXd DynamicFormula::computeNextwv() {
	MatrixXd Kinv = K.inverse();
	//cout << "Kinv" << Kinv << endl;
	VectorXd res=Kinv * lp;
	//w = res.block(0, 0, 3, 1);
	//v = res.block(3, 0, 3, 1);
	return res;
}
Vector3d DynamicFormula::vec62Vec31(VectorXd wv) {//也可以用于wv
	return wv.block(0, 0, 3, 1);
}
Vector3d DynamicFormula::vec62Vec32(VectorXd wv) {
	return wv.block(3, 0, 3, 1);
}

void DynamicFormula::nextTime() {
	lp_ = computelp_();
//	cout << "oldlp" << lp << endl;
//	cout << "nextlp_:" << lp_<< endl;
	lp=computeNextlp();//lp
//	cout << "lp:" << lp << endl;
	VectorXd tempwv= 0.5*computeNextwv();
	//cout << "tempwv" << tempwv << endl;
//	cout << "w:"<<w(0)<<" " << w(1) << " " << w(2) << endl;
//	cout << "v:"<<v(0) << " " << v(1) << " " << v(2) << endl;
	Vector3d y_ = computey_();
//	cout << "y_" << y_ << endl;
	y=computeNexty(y_);//y
//	cout << "y:" << y << endl;
	R=computeNextR();//R
	cout << "R:" << R << endl;
	w = tempwv.block(0, 0, 3, 1);//w
	v = tempwv.block(3, 0, 3, 1);//v
//	cout << "w:" << w(0) << " " << w(1) << " " << w(2) << endl;
//	cout << "v:" << v(0) << " " << v(1) << " " << v(2) << endl;
}
void DynamicFormula::set_tsfs(Vector3d ts,Vector3d fs) {

	//this->tsfs
}