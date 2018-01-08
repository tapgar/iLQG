/*
 * main.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#include "Hopper.h"
//#include "PI2.h"
#include "iLQG.h"
#include <iostream>

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
	bool save_vid = false;
	if (argc > 1)
	{
		save_vid = true;
	}
	Hopper* hopper = new Hopper(save_vid);

	cout << "here" << endl;
	iLQG* ilqg = new iLQG(hopper);

	cout << "here" << endl;
//	PI2* pi2 = new PI2(hopper);

	vector<ActionMatrix> u;
	vector<StateMatrix> x;
	Matrix<double, nU, nX> fb = Matrix<double, nU, nX>::Zero();

	cout << "here" << endl;
	ilqg->GetTrajectory(&x, &u);

	int c = 0;
	while (c++ < 1)
	{
#if LOGGING
		hopper->MPC_Start(c);
#endif
		cout << "here" << endl;
		printf("MPC Iteration: %d\n", c);
		ilqg->RunMPC(hopper->getModel(), hopper->getData(), u);

		cout << "here" << endl;
		ilqg->GetTrajectory(&x, &u);
		cout << "here" << endl;
		ilqg->GetController(&fb);

		cout << "here" << endl;
		//cout << fb << endl;

		hopper->DrawTargetTrajectory(x);
		hopper->Step(x[1], u[0], fb);
		printf("\n\n\n\n");
#if LOGGING
		hopper->MPC_End();
#endif
		//usleep(1E6);
	}
//	int c = 0;
//
//	while (c++ < 20000)
//	{
//		dip->Step();
//	}
}


