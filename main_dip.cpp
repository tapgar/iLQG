/*
 * main.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#include "DIP.h"
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
	DIP* dip = new DIP(save_vid);

	iLQG* ilqg = new iLQG(dip);

	vector<ActionMatrix> u;
	vector<StateMatrix> x;
	Matrix<double, nU, nX> fb;

	ilqg->GetTrajectory(&x, &u);

	int c = 0;
	while (c++ < 500)
	{
		printf("MPC Iteration: %d\n", c);
		ilqg->RunMPC(dip->getModel(), dip->getData(), u);

		ilqg->GetTrajectory(&x, &u);
		ilqg->GetController(&fb);

		cout << fb << endl;

		dip->DrawTargetTrajectory(x);
		dip->Step(x[1], u[0], fb);
		printf("\n\n\n\n");

		//usleep(1E6);
	}
//	int c = 0;
//
//	while (c++ < 20000)
//	{
//		dip->Step();
//	}
}


