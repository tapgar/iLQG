/*
 * Environment.cpp
 *
 *  Created on: Jul 5, 2017
 *      Author: tapgar
 */

#include "Environment.h"

using namespace Eigen;
using namespace std;

void Environment::Step(Matrix<double, nX, 1> targ_x, Matrix<double, nU, 1> u, Matrix<double, nU, nX> K) {

	int nSteps = int(m_dControlTime_s/m_dDeltaSimTime_s);
	int iter = 0;

	Matrix<double, nX, 1> cur_x = Matrix<double, nX, 1>::Zero();
	Matrix<double, nU, 1> ufb = Matrix<double, nU, 1>::Zero();

//	double Fx, Fy, Fz;
////	m_Vis->GetDisturbanceForce(&Fx, &Fy, &Fz);
////
//
//	Fx = 0.0;
//
//	if (nMPCDir < 10)
//		Fx = 50;
//
//	mj_Data->xfrc_applied[6] = Fx;
//	mj_Data->xfrc_applied[7] = Fy;
//	mj_Data->xfrc_applied[8] = Fz;

	while (iter++ < nSteps)
	{
		mj_step1(mj_Model, mj_Data);
		GetStateFromMujoco(mj_Data, &cur_x);

		for (int i = 0; i < nX; i++)
			fileFollow << targ_x(i,0) << ",";
		for (int i = 0; i < nX; i++)
			fileFollow << cur_x(i,0) << ",";
		for (int i = 0; i < nU; i++)
			fileFollow << u(i,0) << ",";

		ufb = u + K*(cur_x - targ_x);
		for (int i = 0; i < nU; i++)
			fileFollow << ufb(i,0) << ",";
		fileFollow << endl;
		UpdateControl(mj_Data, ufb);
		//printf("%f\t%f\n", mj_Data->xfrc_applied[6], mj_Data->xfrc_applied[8]);
//		for (int i = 0; i < nU; i++)
//		{
//			mj_Data->ctrl[i] = ufb(i,0);
//		}

		printf("FB Error\n");
		for (int i = 0; i < nX; i++)
			printf("%f\t", cur_x(i,0) - targ_x(i,0));
		printf("\n");

		mj_step2(mj_Model, mj_Data);
	}

	GetStateFromMujoco(mj_Data, &cur_x);

	for (int i = 0; i < nX; i++)
		fileFollow << targ_x(i,0) << ",";
	for (int i = 0; i < nX; i++)
		fileFollow << cur_x(i,0) << ",";
	for (int i = 0; i < nU; i++)
		fileFollow << u(i,0) << ",";
	for (int i = 0; i < nU; i++)
		fileFollow << ufb(i,0) << ",";
	fileFollow << endl;


	m_Vis->Draw(mj_Data);
}


#if LOGGING
void Environment::MPC_Start(int mpc_iter)
{
	nMPCDir = mpc_iter;

	char temp_buff[500];

	sprintf(temp_buff, "out/MPC%d", nMPCDir);

	mkdir(temp_buff, 0700);


	UpdateFN(temp_buff, fileTrajx_name);
	fileTrajx.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileTraju_name);
	fileTraju.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileTrajc_name);
	fileTrajc.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileEvalTrajx_name);
	fileEvalTrajx.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileEvalTraju_name);
	fileEvalTraju.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileEvalTrajc_name);
	fileEvalTrajc.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileHP_name);
	fileHP.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileA_name);
	fileA.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileB_name);
	fileB.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, filelx_name);
	filelx.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, filelxx_name);
	filelxx.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, filelu_name);
	filelu.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, filelux_name);
	filelux.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileluu_name);
	fileluu.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileQx_name);
	fileQx.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileQxx_name);
	fileQxx.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileQu_name);
	fileQu.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileQux_name);
	fileQux.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileQuu_name);
	fileQuu.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileVx_name);
	fileVx.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileVxx_name);
	fileVxx.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, filefb_name);
	filefb.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileff_name);
	fileff.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileEnvVars_name);
	fileEnvVars.open (temp_buff, std::ofstream::out);

	UpdateFN(temp_buff, fileFollow_name);
	fileFollow.open (temp_buff, std::ofstream::out);

	filepi2xs.clear();
	filepi2us.clear();
	filepi2cs.clear();
	filepi2Ss.clear();
	filepi2Ps.clear();

	for (int k = 0; k <= K; k++)
	{
		UpdateFN(temp_buff, filepi2x_name, k);
		shared_ptr<ofstream> tempFile(new std::ofstream);
		tempFile->open(temp_buff);
		filepi2xs.push_back(tempFile);

		UpdateFN(temp_buff, filepi2u_name, k);
		shared_ptr<ofstream> tempFile2(new std::ofstream);
		tempFile2->open(temp_buff);
		filepi2us.push_back(tempFile2);

		UpdateFN(temp_buff, filepi2c_name, k);
		shared_ptr<ofstream> tempFile3(new std::ofstream);
		tempFile3->open(temp_buff);
		filepi2cs.push_back(tempFile3);

		if (k == K)
			break;

		UpdateFN(temp_buff, filepi2S_name, k);
		shared_ptr<ofstream> tempFile4(new std::ofstream);
		tempFile4->open(temp_buff);
		filepi2Ss.push_back(tempFile4);

		UpdateFN(temp_buff, filepi2P_name, k);
		shared_ptr<ofstream> tempFile5(new std::ofstream);
		tempFile5->open(temp_buff);
		filepi2Ps.push_back(tempFile5);
	}
}

void Environment::UpdateFN(char* buff, std::string orig_fn)
{
	sprintf(buff, "out/MPC%d/%s", nMPCDir, orig_fn.c_str());
}

void Environment::UpdateFN(char* buff, std::string orig_fn, int iter)
{
	sprintf(buff, "out/MPC%d/%s_%d.csv", nMPCDir, orig_fn.c_str(), iter);
}

void Environment::MPC_End()
{
	fileTrajx.close();
	fileTraju.close();
	fileTrajc.close();
	fileEvalTrajx.close();
	fileEvalTraju.close();
	fileEvalTrajc.close();
	fileHP.close();
	fileA.close();
	fileB.close();
	filelx.close();
	filelxx.close();
	filelu.close();
	filelux.close();
	fileluu.close();
	fileQx.close();
	fileQxx.close();
	fileQu.close();
	fileQux.close();
	fileQuu.close();
	fileVx.close();
	fileVxx.close();
	fileff.close();
	filefb.close();
	fileEnvVars.close();
	fileFollow.close();

	for (int k = 0; k <= K; k++)
	{
		(*filepi2xs[k]).close();
		(*filepi2us[k]).close();
		(*filepi2cs[k]).close();
		if (k == K)
			break;
		(*filepi2Ss[k]).close();
		(*filepi2Ps[k]).close();
	}
}

void Environment::MPC_DirExist()
{

}

void Environment::LogTrajectory(Matrix<double, nX, 1> x, Matrix<double, nU, 1> u, double c)
{
	for (int i = 0; i < nX; i++)
		fileTrajx << x(i,0) << ",";
	for (int i = 0 ; i < nU; i++)
		fileTraju << u(i,0) << ",";
	fileTrajc << c << ",";

	fileTrajx << "\n";
	fileTraju << "\n";
	fileTrajc << "\n";
}

void Environment::LogEvalTrajectory(Matrix<double, nX, 1> x, Matrix<double, nU, 1> u, double c, int numUpdates)
{
	for (int i = 0; i < nX; i++)
		fileEvalTrajx << x(i,0) << ",";
	for (int i = 0 ; i < nU; i++)
		fileEvalTraju << u(i,0) << ",";
	fileEvalTrajc << c << "," << numUpdates;

	fileEvalTrajx << "\n";
	fileEvalTraju << "\n";
	fileEvalTrajc << "\n";

}

void Environment::LogHyperParams(double lambda, double alpha)
{
	fileHP << lambda << "," << alpha << "\n";
}


void Environment::LogPI2Trajectory(Matrix<double, nX, 1> x, Matrix<double, nU, 1> u, double c, double S, double P, int iter)
{
	for (int i = 0; i < nX; i++)
		(*filepi2xs[iter]) << x(i,0) << ",";
	for (int i = 0 ; i < nU; i++)
		(*filepi2us[iter]) << u(i,0) << ",";
	(*filepi2cs[iter]) << c << ",";

	(*filepi2xs[iter]) << "\n";
	(*filepi2us[iter]) << "\n";
	(*filepi2cs[iter]) << "\n";

	if (iter >= K)
		return;
	(*filepi2Ss[iter]) << S << ",\n";
	(*filepi2Ps[iter]) << P << ",\n";

}


void Environment::LogForwardPass(Matrix<double, nX, nX> A, Matrix<double, nX, nU> B, Matrix<double, nX, 1> lx, Matrix<double, nX, nX> lxx,
		Matrix<double, nU, nX> lux, Matrix<double, nU, nU> luu, Matrix<double, nU, 1> lu)
{
	for (int i = 0; i < nX; i++)
	{
		filelx << lx(i,0) << ",";
		for (int j = 0; j < nX; j++)
		{
			fileA << A(i,j) << ",";
			filelxx << lxx(i,j) << ",";
		}

		for (int j = 0; j < nU; j++)
			fileB << B(i,j) << ",";

	}

	for (int i = 0; i < nU; i++)
	{
		filelu << lu(i,0) << ",";
		for (int j = 0; j < nX; j++)
			filelux << lux(i,j) << ",";

		for (int j = 0; j < nU; j++)
			fileluu << luu(i,j) << ",";
	}

	filelx << "\n";
	filelxx << "\n";
	fileA << "\n";
	filelu << "\n";
	fileB << "\n";
	filelux << "\n";
	fileluu << "\n";
}

void Environment::LogBackwardPass(Matrix<double, nX, 1> Vx, Matrix<double, nX, nX> Vxx, Matrix<double, nX, 1> Qx, Matrix<double, nX, nX> Qxx,
				Matrix<double, nU, nX> Qux, Matrix<double, nU, nU> Quu, Matrix<double, nU, 1> Qu, Matrix<double, nU, nX> fb, Matrix<double, nU, 1> ff, int back_iter)
{
	fileVx << back_iter << ",";
	fileVxx << back_iter << ",";
	fileQx << back_iter << ",";
	fileQxx << back_iter << ",";
	fileQux << back_iter << ",";
	fileQu << back_iter << ",";
	fileQuu << back_iter << ",";
	filefb << back_iter << ",";
	fileff << back_iter << ",";

	for (int i = 0; i < nX; i++)
	{
		fileVx << Vx(i,0) << ",";
		fileQx << Qx(i,0) << ",";
		for (int j = 0; j < nX; j++)
		{
			fileVxx << Vxx(i,j) << ",";
			fileQxx << Qxx(i,j) << ",";
		}

	}

	for (int i = 0; i < nU; i++)
	{
		fileff << ff(i,0) << ",";
		fileQu << Qu(i,0) << ",";
		for (int j = 0; j < nX; j++)
		{
			filefb << fb(i,j) << ",";
			fileQux << Qux(i,j) << ",";
		}

		for (int j = 0; j < nU; j++)
		{
			fileQuu << Quu(i,j) << ",";
		}
	}

	fileVx << "\n";
	fileQx << "\n";
	fileVxx << "\n";
	fileQxx << "\n";
	fileff << "\n";
	fileQu << "\n";
	filefb << "\n";
	fileQux << "\n";
	fileQuu << "\n";
}

#endif
