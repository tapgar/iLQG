/*
 * main.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#include "SLIP.h"
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
	SLIP* slip = new SLIP(save_vid);

	iLQG* ilqg = new iLQG(slip);

	vector<ActionMatrix> u;
	vector<StateMatrix> x;
	Matrix<double, nU, nX> fb;

	ilqg->GetTrajectory(&x, &u);

	int c = 0;
	while (c++ < 1000)
	{
#if LOGGING
		slip->MPC_Start(c);
#endif
		printf("MPC Iteration: %d\n", c);
		ilqg->RunMPC(slip->getModel(), slip->getData(), u);

		ilqg->GetTrajectory(&x, &u);
		ilqg->GetController(&fb);

		//cout << fb << endl;

		slip->DrawTargetTrajectory(x);
		slip->Step(x[1], u[0], fb);
		printf("\n\n\n\n");
#if LOGGING
		slip->MPC_End();
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


