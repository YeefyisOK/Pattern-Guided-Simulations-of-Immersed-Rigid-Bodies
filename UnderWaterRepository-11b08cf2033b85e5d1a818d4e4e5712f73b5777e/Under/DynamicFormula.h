#pragma once
#include <Eigen/Dense>
#include <Eigen/Cholesky>  
#include <Eigen/LU>  
#include <Eigen/QR>  
#include <Eigen/SVD>  
//#include"Quaternion.h"
#include <iostream>
using namespace std;
using namespace Eigen;
class DynamicFormula
{
public:
	DynamicFormula(Vector3d omega, Vector3d velocity, Matrix3d R,
		Vector3d y, double delta_t);
	~DynamicFormula();
	void setK(MatrixXd newK) {
		K = newK;
	}
	void setR(Matrix3d newR) {
		R = newR;
	}
	Matrix3d toDaOmegaOrY(Vector3d omega);
	VectorXd tsfs2tf(Matrix3d Y);
	//Matrix3d computeR_();
	Vector3d computey_();
	VectorXd computelp();
	VectorXd computelp_();
	Matrix3d computeNextR();
	Vector3d computeNexty(Vector3d y_);
	VectorXd computeNextlp();
	VectorXd computeNextwv();
	Vector3d vec62Vec31(VectorXd ab);
	Vector3d vec62Vec32(VectorXd ab);
	void nextTime();
	void set_tsfs(Vector3d ts,Vector3d fs);
	float* GetRotAndTransData();
//private:
	//初始时刻的 6+1个量
	Vector3d w;
	Vector3d v;
	Matrix3d R;//R初始化一个正交矩阵
	Quaterniond q;//R转化四元数
	Vector3d y;
	VectorXd tsfs;
	VectorXd lp;
	VectorXd lp_;

	MatrixXd K;
	double delta_t;
	//平移旋转需要的两个参数
	Quaterniond delta_q;
	double theta;

};