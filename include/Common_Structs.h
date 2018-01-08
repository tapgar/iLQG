/*
 * Common_Structs.h
 *
 *  Created on: Jun 21, 2017
 *      Author: tapgar
 */

#ifndef COMMON_STRUCTS_H_
#define COMMON_STRUCTS_H_

#include <Eigen/Dense>

#define nX nQ + nQd
#define LOGGING true
#define PAUSE_VIS false

using namespace Eigen;

struct Traj_Pt {
	Matrix<double, nX, 1, DontAlign> x;
	Matrix<double, nU, 1, DontAlign> u;
	Traj_Pt()
	{
		x = Matrix<double, nX, 1, DontAlign>::Zero();
		u = Matrix<double, nU, 1, DontAlign>::Zero();
	}
};

static const double m_dDeltaSimTime_s = 0.001;
static const double m_dControlTime_s = 0.01;
static const int N = 100;
static const int K = 30;

static constexpr double m_deps = 1e-6;          // finite-difference epsilon


struct Gains {
	Matrix<double, nU, nX> fb;
	Matrix<double, nU, 1> ff;
	Gains()
	{
		fb = Matrix<double, nU, nX>::Zero();
		ff = Matrix<double, nU, 1>::Zero();
	}
};

struct Cost {
	double l;
	Matrix<double, nX, 1> lx;
	Matrix<double, nU, 1> lu;
	Matrix<double, nX, nX> lxx;
	Matrix<double, nU, nU> luu;
	Matrix<double, nU, nX> lux;
};

struct Derivs {
	Cost C;
	Matrix<double, nX, nX> A;
	Matrix<double, nX, nU> B;
};

typedef Matrix<double, nX, 1> StateMatrix;
typedef Matrix<double, nU, 1> ActionMatrix;
typedef Matrix<double, nU, nU> ActionCovMatrix;
typedef Matrix<double, nL, nL> ControllableMatrix;


#endif /* COMMON_STRUCTS_H_ */
