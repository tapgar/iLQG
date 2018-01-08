/*
 * main.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: tapgar
 */

#include "SLIP.h"
#include <iostream>

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
	bool save_vid = false;
	if (argc != 7)
	{
		printf("o fuck 7 args please\n");
		return -1;
	}

	double vtarg = atof(argv[1]);
	double zdes = atof(argv[2]);
	double k = atof(argv[3]);
	double C = atof(argv[4]);
	double T = atof(argv[5]);
	double Kpt = atof(argv[6]);
	double Kdt = atof(argv[7]);

	SLIP* slip = new SLIP(save_vid);

	while (true)
	{
		slip->RunRaibert(vtarg, zdes, k, C, T, Kpt, Kdt);
	}
}


