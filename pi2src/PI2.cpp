/*
 * PI2.cpp
 *
 *  Created on: Jul 13, 2017
 *      Author: tapgar
 */



#include "PI2.h"

using namespace Eigen;
using namespace std;

PI2::PI2(Environment* pEnv)
{
	m_pEnv = pEnv;

	ActionMatrix null_ctrl = ActionMatrix::Zero();
	StateMatrix null_state = StateMatrix::Zero();

	CovMat = Matrix<double, nP, nP>::Zero();
	R = CovMat;
	Rinv = R;

	PI_Derivs null_PI_derivs;
	null_PI_derivs.B = Matrix<double, nL, nU>::Zero();
	null_PI_derivs.Phi = Matrix<double, nL, 1>::Zero();

	ControllableMatrix null_ctrl_mat = ControllableMatrix::Zero();
	Path null_path;
	null_path.cntrl = Matrix<double, nP, 1>::Zero();
	for (int i = 0; i < N; i++)
	{
		null_path.state.push_back(null_state);
		null_path.cost.push_back(0.0);
		null_path.S.push_back(0.0);
		null_path.P.push_back(0.0);

		null_path.Derivs.push_back(null_PI_derivs);

		null_path.H.push_back(null_ctrl_mat);
		null_path.Hinv.push_back(null_ctrl_mat);

		unew.push_back(null_ctrl);
		xnew.push_back(null_state);
		cnew.push_back(0.0);
	}
	xnew.push_back(null_state);
	cnew.push_back(0.0);
	null_path.state.push_back(null_state); // one more for final state
	null_path.cost.push_back(0.0);

	for (int i = 0; i < K; i++)
		m_Paths.push_back(null_path);

	//FD matrix
	Matrix<double, N + 1, N> A = Matrix<double, N + 1, N>::Zero();
	for (int i = 0; i < N + 1; i++)
	{
		if (i < N)
			A(i,i) = 1.0/m_dControlTime_s;
		if (i > 0)
			A(i,i-1) = -0.7/m_dControlTime_s;
	}
	A = A * 0.05;

	Matrix<double, N, N> tempR = A.transpose()*A;

	Matrix<double, nU, 1> Rweight;
	Rweight << 0.0001, 0.0001, 0.0001;

	R = Matrix<double, nP, nP>::Zero();
	for (int i = 0; i < nU; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				R(i*N + j, i*N + k) = tempR(j,k)*Rweight(i,0);
			}
		}
	}
	Rinv = R.inverse();

	lambda = 1/h;

	noise_factor = 0.5;

	CovMat = noise_factor*Rinv;

	//cout << CovMat << endl;

	for (int n = 0; n < N; n++)
	{
		BasisMatrix tempG = BasisMatrix::Zero();
		for (int i = 0; i < nU; i++)
			tempG(i, i*N + n) = 1.0;
		G.push_back(tempG);
	}

}

void PI2::RunMPC(const mjModel* m, const mjData* dmain, vector<ActionMatrix> init_cntrl)
{

	int iter = 0;
	int Num_iter = 10;

	unew = init_cntrl;
	while (iter++ < Num_iter)
	{
		Rollouts(m, dmain, unew);
		UpdateNominal();
		NominalRollout(m, dmain);

#if LOGGING
		for (int n = 0; n < N; n++)
		{
			for (int k = 0; k < K; k++)
			{
				ActionMatrix u;
				for (int j = 0; j < nU; j++)
					u(j,0) = m_Paths[k].cntrl(N*j + n,0);
				m_pEnv->LogPI2Trajectory(m_Paths[k].state[n],u,m_Paths[k].cost[n],m_Paths[k].S[n],m_Paths[k].P[n],k);
			}
			m_pEnv->LogPI2Trajectory(xnew[n],unew[n],cnew[n],0.0,0.0,K);
		}
#endif
	}



}

void PI2::Rollouts(const mjModel* m, const mjData* dmain, vector<ActionMatrix> init_cntrl)
{

	struct timeval start;
	gettimeofday(&start, NULL);

	long useconds = start.tv_sec*1E6 + start.tv_usec;

	Matrix<double, nP, 1> null_ctrl = Matrix<double, nP, 1>::Zero();
	EigenMultivariateNormal<double> normX_solver(null_ctrl, CovMat, false, useconds);

	for (int k = 0; k < K; k++)
	{
		m_Paths[k].cntrl = normX_solver.samples(1);
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < nU; j++)
			{
				m_Paths[k].cntrl(j*N + i,0) += init_cntrl[i](j,0);

				if (m_Paths[k].cntrl(j*N + i,0) < 0.95*m->actuator_ctrlrange[j*2])
					m_Paths[k].cntrl(j*N + i,0) = 0.95*m->actuator_ctrlrange[j*2];
				if (m_Paths[k].cntrl(j*N + i,0) > 0.95*m->actuator_ctrlrange[j*2 +1])
					m_Paths[k].cntrl(j*N + i,0) = 0.95*m->actuator_ctrlrange[j*2 +1];
			}
		}
		Rollout(m, dmain, k);
	}
}

void PI2::NominalRollout(const mjModel* m, const mjData* dmain)
{
	mjData* d = mj_makeData(m);
	mj_copyData(d, m, dmain);

	m_pEnv->GetStateFromMujoco(d, &xnew[0]);

	int nSimSteps = int(m_dControlTime_s/m_dDeltaSimTime_s);

//	printf("Updated Controls\n");

	for (int i = 0; i < N; i++)
	{
		m_pEnv->UpdateControl(d, unew[i]);

//		cout  << unew[i].transpose() << endl;

		//update simulation
		for (int j = 0; j < nSimSteps; j++)
			mj_step(m, d);
//
//		for (int j = 0; j < nU; j++)
//			printf("%f\t", d->ctrl[j]);
//		printf("\n");

		m_pEnv->Draw(d);

		//Calculate state cost
		m_pEnv->GetStateCost(m, d, xnew[i], &cnew[i]);

		//copy mujoco state
		m_pEnv->GetStateFromMujoco(d, &xnew[i+1]);
	}

	m_pEnv->GetStateCost(m, d, xnew[N], &cnew[N], true);

	mj_deleteData(d);
}

void PI2::Rollout(const mjModel* m, const mjData* dmain, int path_idx)
{

	mjData* d = mj_makeData(m);
	mj_copyData(d, m, dmain);

	m_pEnv->GetStateFromMujoco(d, &m_Paths[path_idx].state[0]);

	int nSimSteps = int(m_dControlTime_s/m_dDeltaSimTime_s);

	printf("Rolling out %d\n", path_idx);

	for (int i = 0; i < N; i++)
	{
//		cout << "Noisy Control" << endl;
//		cout << m_Paths[path_idx].cntrl[i].transpose() << endl;

		ActionMatrix u;
		for (int j = 0; j < nU; j++)
		{
			u(j,0) = m_Paths[path_idx].cntrl(N*j + i,0);
		}

		m_pEnv->UpdateControl(d, u);

		UpdateDerivatives(m, d, path_idx, i);


//		cout << "Updated?" << endl;
//		for (int j = 0; j < nU; j++)
//			printf("%f\t", d->ctrl[j]);
//		printf("\n");

		//update simulation
		for (int j = 0; j < nSimSteps; j++)
			mj_step(m, d);

		//m_pEnv->Draw(d);

//		ActionMatrix duL = lambda*m_Paths[path_idx].Derivs[i].B.inverse()*m_Paths[path_idx].H[i]*m_Paths[path_idx].Derivs[i].Phi;
//
//		for (int j = 0; j < nU; j++)
//		{
//			if (fabs(duL(j,0)) > 20.0)
//			{
//				cout << "ul:\n" << duL << endl;
//				cout << "B:\n" << m_Paths[path_idx].Derivs[i].B << endl;
//				cout << "Binv:\n" << m_Paths[path_idx].Derivs[i].B.inverse() << endl;
//				cout << "Phi:\n" << m_Paths[path_idx].Derivs[i].Phi << endl;
//				cout << "H:\n" << m_Paths[path_idx].H[i] << endl;
//
//
//				usleep(10E6);
//			}
//		}



		//Calculate state cost
		m_pEnv->GetStateCost(m, d, m_Paths[path_idx].state[i], &m_Paths[path_idx].cost[i]);

		//copy mujoco state
		m_pEnv->GetStateFromMujoco(d, &m_Paths[path_idx].state[i+1]);
	}

	m_pEnv->GetStateCost(m, d, m_Paths[path_idx].state[N], &m_Paths[path_idx].cost[N], true);

	mj_deleteData(d);
}

void PI2::UpdateDerivatives(const mjModel* m, const mjData* dmain, int path_idx, int time_idx)
{
	bool bPhiWarning = false;
	m_pEnv->JacobianControl(m, dmain, &m_Paths[path_idx].Derivs[time_idx].B);
	//m_Paths[path_idx].Derivs[time_idx].B *= m_dControlTime_s;

	m_Paths[path_idx].H[time_idx] = m_Paths[path_idx].Derivs[time_idx].B*G[time_idx]*Rinv*G[time_idx].transpose()*m_Paths[path_idx].Derivs[time_idx].B.transpose();

//	cout << " B matrix: \n" << m_Paths[path_idx].Derivs[time_idx].B << endl;
//	cout << " H matrix: \n" << m_Paths[path_idx].H[time_idx] << endl;

	m_Paths[path_idx].Hinv[time_idx] = m_Paths[path_idx].H[time_idx].inverse();

	for (int i = 0; i < nL; i++)
	{
		Matrix<double, nL, nU> B = Matrix<double, nL, nU>::Zero();
		m_pEnv->JacobianControl(m, dmain, &B, i);
		//B *= m_dControlTime_s;

		Matrix<double, nL, nL> dHdxi = m_Paths[path_idx].Hinv[time_idx]*((B*G[time_idx]*Rinv*G[time_idx].transpose()*B.transpose() - m_Paths[path_idx].H[time_idx])/m_deps);

		m_Paths[path_idx].Derivs[time_idx].Phi(i,0) = 0.0;

		for (int j = 0; j < nL; j++)
			m_Paths[path_idx].Derivs[time_idx].Phi(i,0) += dHdxi(j,j);

		if (fabs(m_Paths[path_idx].Derivs[time_idx].Phi(i,0)) > 10.0)
		{
			bPhiWarning = true;
			cout << "\nPhi Warning at index: " << i << endl;

			cout << "Original B:\n";
			cout << m_Paths[path_idx].Derivs[time_idx].B << endl;

			cout << "Perturbed B:\n";
			cout << B << endl;

			cout << "A fluke?\n";
			m_pEnv->JacobianControl(m, dmain, &B, i);
			cout << B << endl;

			cout << "dHdxi: " << endl;
			cout << dHdxi << endl;

			cout << "new H: \n";
			cout << B*G[time_idx]*Rinv*G[time_idx].transpose()*B.transpose() << endl;

			cout << "prev H: \n";
			cout << m_Paths[path_idx].H[time_idx] << endl;

			cout << "H inv\n";
			cout << m_Paths[path_idx].Hinv[time_idx] << endl;

			printf("Old vel:\n");
			for (int i = 0; i < nQd; i++)
				printf("%f\t", dmain->qvel[i]);
			printf("\n");
		}
	}
	m_Paths[path_idx].Derivs[time_idx].Phi *= 0.5;

	if (bPhiWarning)
	{
		cout << "Phi: \n";
		cout << m_Paths[path_idx].Derivs[time_idx].Phi << endl;
		//usleep(10E6);
		printf("\n\n\n");
	}

}

void PI2::UpdateNominal()
{
	int window_size = 2;
	double dTotal = 0.0;
	Matrix<double, nP, 1> uLnorm = Matrix<double, nP, 1>::Zero();
	MatrixXd CovMatnorm = Matrix<double, nP, nP>::Zero();
	Matrix<double, nP, 1> theta = Matrix<double, nP, 1>::Zero();

	vector<int> samples(N);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < nU; j++)
		{
			theta(j*N + i,0) = unew[i](j,0);
		}
	}

	for (int i = 0; i < N; i++)
	{
		double minS = 1E10;
		double maxS = -1E10;
		for (int k = 0; k < K; k++)
		{
			//calculate S for path k at time i
			CalculateS(k, i);

//			printf("S(%d,%d):\t%f\n", k, i, m_Paths[k].S[i]);

			if (m_Paths[k].S[i] < minS)
				minS = m_Paths[k].S[i];
			if (m_Paths[k].S[i] > maxS)
				maxS = m_Paths[k].S[i];
		}

//		printf("MinS: %f\tMaxS: %f\n", minS, maxS);

		double p_total = 0.0;
		for (int k = 0; k < K; k++)
		{
			CalculateP(k, i, minS, maxS);
			p_total += m_Paths[k].P[i];
		}

		unew[i] = ActionMatrix::Zero();
		Matrix<double, nP, 1> uLtot = Matrix<double, nP, 1>::Zero();
		MatrixXd CovMattot = Matrix<double, nP, nP>::Zero();
		for (int k = 0; k < K; k++)
		{
			m_Paths[k].P[i] /= p_total;

			Matrix<double, nP, 1> uL = Matrix<double, nP, 1>::Zero();

			//uL = m_Paths[k].cntrl[i] - lambda*m_Paths[k].Derivs[i].B.inverse()*m_Paths[k].H[i]*m_Paths[k].Derivs[i].Phi;
			uL = Rinv*G[i].transpose()*m_Paths[k].Derivs[i].B.transpose()*m_Paths[k].Hinv[i]*(m_Paths[k].Derivs[i].B*G[i]*m_Paths[k].cntrl - lambda*m_Paths[k].H[i]*m_Paths[k].Derivs[i].Phi);
//			cout << "P:\t" << m_Paths[k].P[i] << endl;
			//cout << "U:\t" << uL.transpose() << endl;
			uLtot += m_Paths[k].P[i]*uL;
			MatrixXd temp = (m_Paths[k].cntrl - theta);
			CovMattot = CovMattot +  m_Paths[k].P[i]*temp*temp.transpose();
		}


		for (int k = i - window_size; k < i + window_size; k++)
		{
			if (k < 0 || k >= N)
				continue;

			samples[k]++;
			for (int j = 0; j < nU; j++)
				unew[k](j,0) += uLtot(j*N + i,0);
		}

		uLnorm += double(N-i)*uLtot;
		CovMatnorm += double(N-i)*CovMattot;

		dTotal += double(N-i);
	}

	for (int i = 0; i < N; i++)
	{
		unew[i] /= double(samples[i]);
	}

//	for (int i = 0; i < N; i++)
//	{
//		for (int j = 0; j < nU; j++)
//			unew[i](j,0) = uLnorm(j*N + i,0)/dTotal;
//	}

	CovMat = CovMat*0.15 + CovMatnorm*0.85/dTotal;

	//CovarianceMatrixAdaptation();
}

void PI2::CalculateS(int path_idx, int time_idx)
{
	double state_cost_sum = accumulate(m_Paths[path_idx].cost.begin()+time_idx, m_Paths[path_idx].cost.end(), 0.0);
//	printf("State cost sum: %f\n", state_cost_sum);

	m_Paths[path_idx].S[time_idx] = state_cost_sum;

	Matrix<double, 1, 1> res;
	for (int j = time_idx; j < N; j++)
	{

		res = m_dControlTime_s*0.25*(m_Paths[path_idx].cntrl.transpose()*G[j].transpose()*m_Paths[path_idx].Derivs[j].B.transpose()*m_Paths[path_idx].Hinv[j]*m_Paths[path_idx].Derivs[j].B*G[j]*m_Paths[path_idx].cntrl);
//		printf("Res: %f\n", res(0,0));
		m_Paths[path_idx].S[time_idx] += res(0,0);

//		cout << "H: \n" << m_Paths[path_idx].H[j] << endl;
//		cout << "B:\n" << m_Paths[path_idx].Derivs[j].B << endl;

		double detH = m_Paths[path_idx].H[j].determinant();
//		printf("det H: %f\n", detH);
		if (detH < 0.0)
		{
			printf("Alert!!!\n");
			usleep(10E6);
		}
		m_Paths[path_idx].S[time_idx] += 0.5*lambda*log(detH);
//		printf("log det h: %f\n", 0.5*lambda*log(detH));
	}
}

void PI2::CalculateP(int path_idx, int time_idx, double minS, double maxS)
{
	m_Paths[path_idx].P[time_idx] = exp(-h*(m_Paths[path_idx].S[time_idx] - minS)/(maxS-minS));
}

void PI2::CovarianceMatrixAdaptation()
{
	return;
//	double dTotal = 0.0;
//	Matrix<double, nU, nU> tempCov = CovMat;
//
//	CovMat = Matrix<double, nU, nU>::Zero();
//
//	for (int i = 0; i < N; i++)
//	{
//		Matrix<double, nU, nU> CovMat_i = Matrix<double, nU, nU>::Zero();
//		for (int k = 0; k < K; k++)
//		{
//			CovMat_i += m_Paths[k].P[i]*((m_Paths[k].cntrl[i] - unew[i])*(m_Paths[k].cntrl[i] - unew[i]).transpose());
//		}
//		CovMat += CovMat_i;//double(N-i)*CovMat_i;
//		dTotal += 1.0;//double(N-i);
//	}
//
//	CovMat /= dTotal;
//
//	CovMat = CovMat*0.75 + tempCov*0.25; //keep exploration high
////	Rinv = CovMat*(1/lambda);
////	R = Rinv.inverse();
//
//	cout << "New Covariance:\n" << CovMat << endl;
//	cout << "New Rinv:\n" << Rinv << endl;
//	cout << "New R:\n" << R << endl;
}
