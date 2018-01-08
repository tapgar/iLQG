#include "iLQG.h"

iLQG::iLQG(Environment* env) {

	m_pEnv = env;

	cout << "here" << endl;
	for (int i = 0; i < N+1; i++)
	{
		Derivs tempDeriv;
		m_Derivs.push_back(tempDeriv);
	}

	cout << "here" << endl;
	Gains null_gains;
	Matrix<double, nX, 1> null_x = Matrix<double, nX, 1>::Zero();

	cout << "here" << endl;
    for (int i = 0; i < N; i++)
	{
		m_Controller.push_back(null_gains);
		xnew.push_back(null_x);
		cnew.push_back(0.0);
	}

	cout << "here" << endl;
    m_pEnv->GetInitU(&unew);

	cout << "here" << endl;
    for (int i = 0; i < N; i++)
    	cout << unew[i].transpose() << endl;

	cout << "here" << endl;
    m_LastController = null_gains;
	xnew.push_back(null_x);
	cnew.push_back(0.0);
	u_init = unew;
	x_init = xnew;

	cout << "here" << endl;
}

iLQG::~iLQG() {
	// TODO Auto-generated destructor stub
}

void iLQG::RunMPC(const mjModel* m, const mjData* dmain, vector<ActionMatrix> u)
{
	static bool bFirstIter = true;
	double prev_cost = 0.00;
	int numUpdates = 0;

	vector<StateMatrix> x_res = x_init;
	vector<ActionMatrix> u_res = u_init;
	vector<double> c_res = cnew;

	if (bFirstIter)
	{
		//1) do initial rollout to get first trajectory
		Rollout(m , dmain, x_init, u_init, &x_res, &u_res, &c_res, m_Controller, 1.0, true);
		u = u_res;
		prev_cost = accumulate(c_res.begin(), c_res.end(), 0.0);
	}

	vector<StateMatrix> x = x_res;
	if (!bFirstIter)
	{
		rotate(u.begin(), u.begin()+1, u.end());
		u[N-1] = u_init[N-1];

		Rollout(m , dmain, x, u, &x_res, &u_res, &c_res, m_Controller, 0.0, true);
		prev_cost = accumulate(c_res.begin(), c_res.end(), 0.0);
		u = u_res;
		x = x_res;
	}

	unew = u_res;
	xnew = x_res;
	cnew = c_res;

#if LOGGING
	for (int i = 0; i < N; i++)
	{
		m_pEnv->LogTrajectory(x_res[i], u_res[i], c_res[i]);
	}
	m_pEnv->LogTrajectory(x_res[N], u_init[N-1], c_res[N]);
#endif

	bFirstIter = false;

	int Num_Iters = 10;

	static bool bLastIncreased = false;

	static double lambda = 1.0;
	static double dlambda = 1.0;

	if (bLastIncreased)
		DecreaseLambda(&lambda, &dlambda);

	double zmin = 0;

	int iter = 0;


	while (iter++ < Num_Iters)
	{
		//perform backward pass
		bool bBackwardPassDone = false;
		int nBackPasses = 0;
		double dV1 = 0.0;
		double dV2 = 0.0;

		bLastIncreased = true;

		while (!bBackwardPassDone)
		{
//			bool bDiverged = false;
//			for (int lam = -3; lam <= 3; lam++)
//			{
//				nBackPasses++;
//				BackwardPass(m, u, pow(10.0, double(lam)), &dV1, &dV2, nBackPasses);
//			}

			bool bDiverged = BackwardPass(m, u, lambda, &dV1, &dV2, nBackPasses);

			if (bDiverged)
			{
				if (IncreaseLambda(&lambda, &dlambda))
					break;
			}
			else
				bBackwardPassDone = true;
		}
//		break;

//		double dSum = 0.0;
//		for (int i = 0; i < N; i++)
//		{
//			double dMax = -100000;
//			for (int k = 0; k < nU; k++)
//			{
//				double dVal = abs(m_Controller[i].ff(k,0))/(abs(u[i](k,0)) + 1.0);
//				if (dVal > dMax)
//					dMax = dVal;
//			}
//			dSum += dMax;
//		}

//		if (dSum/double(N) < 1e-4 && lambda < 1e-5)
//		{
//			DecreaseLambda(&lambda, &dlambda);
//			printf("Small gradient termination\n");
//			break;
//		}


		bool bForwardPassDone = false;
		double dCostSum = 0.0;
		//perform forward pass (w/ line search)
		if (bBackwardPassDone)
		{
//			for (double i = -3; i < 3; i += 3.0/44)
//			{
//				double alpha = pow(10.0, -i);
//
//				Rollout(m, dmain, x, u, &x_res, &u_res, &c_res, m_Controller, alpha, false);
//
//#if LOGGING
//	for (int j = 0; j < N; j++)
//	{
//		m_pEnv->LogEvalTrajectory(x_res[j], u_res[j], c_res[j], numUpdates);
//	}
//	m_pEnv->LogEvalTrajectory(x_res[N], u_init[N-1], c_res[N], numUpdates);
//#endif
//
//				dCostSum = accumulate(c_res.begin(), c_res.end(), 0.0);
//
//				double dV = -alpha*(dV1 + alpha*dV2);
//
//				printf("COST (prev/next/expected change): %f\t%f\t%f\n", prev_cost, dCostSum, dV);
//			}
			for (double i = 0; i < 3; i += 3.0/11)
			{
				double alpha = pow(10.0, -i);

				Rollout(m, dmain, x, u, &x_res, &u_res, &c_res, m_Controller, alpha, false);

				#if LOGGING
					for (int j = 0; j < N; j++)
					{
						m_pEnv->LogEvalTrajectory(x_res[j], u_res[j], c_res[j], numUpdates);
					}
					m_pEnv->LogEvalTrajectory(x_res[N], u_init[N-1], c_res[N], numUpdates);
				#endif

				dCostSum = accumulate(c_res.begin(), c_res.end(), 0.0);

				double dV = -alpha*(dV1 + alpha*dV2);

				printf("COST (prev/next/expected change/lambda): %f\t%f\t%f\t%f\n", prev_cost, dCostSum, dV, lambda);

				double z = (prev_cost - dCostSum)/dV;
				if (dV < 0)
				{
					printf("Warning!!! negative expected cost\n");
					break;
				}
				else if (z > zmin)
				{
#if LOGGING
		m_pEnv->LogHyperParams(lambda, alpha);
		numUpdates++;
#endif
					//perform another rollout while calculating derivatives
					Rollout(m, dmain, x, u, &x_res, &u_res, &c_res, m_Controller, alpha, true);
					bForwardPassDone = true;
					dCostSum = accumulate(c_res.begin(), c_res.end(), 0.0);

					break;
				}
//				else if (i == 0 && dV < 1e-5)
//				{
//					printf("Small expected cost decrease\n");
//					break;
//				}
			}
		}

		if (bForwardPassDone)
		{
			//iterate and copy new u
			iter++;
			u = u_res;
			x = x_res;
			unew = u_res;
			xnew = x_res;
			cnew = c_res;

#if LOGGING
			for (int i = 0; i < N; i++)
			{
				m_pEnv->LogTrajectory(x_res[i], u_res[i], c_res[i]);
			}
			m_pEnv->LogTrajectory(x_res[N], u_init[N-1], c_res[N]);
#endif

			if (prev_cost - dCostSum < 1e-7)
			{
				printf("Cost update < 1e-7\n");
				break;
			}

			m_LastController = m_Controller[0]; //save in case you need a good fb policy

			prev_cost = dCostSum;

			DecreaseLambda(&lambda, &dlambda);
			bLastIncreased = false;
		}
		else
		{
			//Update lambda and rerun backward pass
			if (IncreaseLambda(&lambda, &dlambda))
			{
				printf("lambda exceeded max allowed value\n");
				break;
			}
		}

	}
}

/*
 * Run rollout for defined number of iterations following open tape u and controller.
 * dmain is the starting state of the robot
 * 		a local copy of dmain is made initially for advancing the simulation forward
 * 		that local copy may be passed to the jacobian func but will be unchanged
 */
void iLQG::Rollout(const mjModel* m, const mjData* dmain, vector<StateMatrix> x, vector<ActionMatrix> u, vector<StateMatrix>* x_res, vector<ActionMatrix>* u_res, vector<double>* c_res, vector<Gains> cntrl, double alpha, bool bCalcDerivs)
{
	mjData* d = mj_makeData(m);
	mj_copyData(d, m, dmain);

	m_pEnv->GetStateFromMujoco(d, &((*x_res)[0]));

	int nSimSteps = int(m_dControlTime_s/m_dDeltaSimTime_s);

	Matrix<double, nX, 1> cur_x = Matrix<double, nX, 1>::Zero();

	for (int i = 0; i < N; i++)
	{
		//Extract state from data
		m_pEnv->GetStateFromMujoco(d, &cur_x);

#if PAUSE_VIS

		printf("FB Term\n");
		cout << (cntrl[i].fb*(cur_x - x[i])).transpose() << endl;
		printf("FF Term\n");
		cout << (alpha*cntrl[i].ff).transpose() << endl;
		printf("Previous U\n");
		cout << u[i].transpose() << endl;

#endif

		//update u following update rule
		(*u_res)[i] = u[i] + cntrl[i].fb*(cur_x - x[i]) + alpha*cntrl[i].ff;
//		(*u_res)[i] = u[i] + alpha*cntrl[i].ff;

		if (alpha < 1e-10)
			(*u_res)[i] = u[i];

		Matrix<double, nU, 2> bounds = Matrix<double, nU, 2>::Zero();
		m_pEnv->GetBounds(&bounds);

		for (int j = 0; j < nU; j++)
		{
			(*u_res)[i](j,0) = min(max((*u_res)[i](j,0), bounds(j,0)), bounds(j,1));
		}

		m_pEnv->UpdateControl(d, (*u_res)[i]);

		//calculate jacobians
		if (bCalcDerivs)
		{
			for (int k = 0; k < nU; k++)
				d->userdata[k] = (*u_res)[i](k,0);
			UpdateDerivatives(m, d, i);
#if LOGGING
			m_pEnv->LogEnvVars(m, d);
#endif
		}

		//update simulation
		for (int j = 0; j < nSimSteps; j++)
		{
			mj_step1(m, d);

//			//Extract state from data
//			m_pEnv->GetStateFromMujoco(d, &cur_x);
//			//update u following update rule
//			ActionMatrix new_u = u[i]  + cntrl[i].fb*(cur_x - x[i]) + alpha*cntrl[i].ff;
//
//			m_pEnv->GetBounds(&bounds);
//
//			for (int k = 0; k < nU; k++)
//			{
//				new_u(k,0) = min(max(new_u(k,0), bounds(k,0)), bounds(k,1));
//			}
//
//			m_pEnv->UpdateControl(d, new_u);

			mj_step2(m, d);
		}

		//copy mujoco state
		m_pEnv->GetStateFromMujoco(d, &((*x_res)[i+1]));

		//Calculate state action cost
		m_pEnv->GetStateActionCost(m, d, (*x_res)[i], (*u_res)[i], &((*c_res)[i]));

		if (bCalcDerivs)
		{
//			m_pEnv->Draw(d);
		}
	}

	if (bCalcDerivs)
	{
		m_pEnv->UpdateControl(d, (*u_res)[N-1]);
		for (int k = 0; k < nU; k++)
			d->userdata[k] = (*u_res)[N-1](k,0);
		UpdateDerivatives(m, d, N);
#if LOGGING
			m_pEnv->LogEnvVars(m, d);
#endif

	}
	//Calculate state action cost
	m_pEnv->GetStateActionCost(m, d, (*x_res)[N], (*u_res)[N-1], &((*c_res)[N]), true);

	mj_deleteData(d);
}

bool iLQG::BackwardPass(const mjModel* m, vector<ActionMatrix> u, double lambda, double *dV1, double *dV2, int back_iter )
{

	Vx = m_Derivs[N].C.lx;
	Vxx = m_Derivs[N].C.lxx;

	for (int i = N-1; i >= 0; i--)
	{
		printf("Backpass iter: %d\n", i);
		Qu = m_Derivs[i].C.lu + m_Derivs[i].B.transpose()*Vx;
		Qx = m_Derivs[i].C.lx + m_Derivs[i].A.transpose()*Vx;
		Qxx = m_Derivs[i].C.lxx + m_Derivs[i].A.transpose()*Vxx*m_Derivs[i].A;
		Qux = m_Derivs[i].C.lux + m_Derivs[i].B.transpose()*Vxx*m_Derivs[i].A;
		Quu = m_Derivs[i].C.luu + m_Derivs[i].B.transpose()*Vxx*m_Derivs[i].B;


		//regularization of Vxx
		Matrix<double, nX, nX> VxxF = Vxx;
		for (int j = 0; j < nX; j++)
			VxxF(j,j) += 0.0;//lambda;

		Matrix<double, nU, nX> QuxReg = m_Derivs[i].C.lux + m_Derivs[i].B.transpose()*VxxF*m_Derivs[i].A;
		Matrix<double, nU, nU> QuuReg = m_Derivs[i].C.luu + m_Derivs[i].B.transpose()*VxxF*m_Derivs[i].B;

		for (int j = 0; j < nU; j++)
			QuuReg(j,j) += lambda;

		bool bUseboxQP = true;

		if (bUseboxQP)
		{
			vector<bool> clamped(nU);
			for (int j = 0; j < nU; j++)
				clamped[j] = false;

			Matrix<double, nU, 1> kff = Matrix<double, nU, 1>::Zero();
			MatrixXd Hfree;
			Matrix<double, nU, 2> bounds = Matrix<double, nU, 2>::Zero();
			m_pEnv->GetBounds(&bounds);
			for (int j = 0; j < nU; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					bounds(j,k) -= u[i](j,0);
				}
			}

//			for (int k = 0; k < nU; k++)
//			{
//				printf("%f\t", xnew[i](m_pEnv->getMotorPosIdx(k),0));
//			}
//			printf("\n");

			bool bDiverged = BoxQP(QuuReg, Qu, m_Controller[min(N-1, i+1)].ff, bounds, &kff, &Hfree, &clamped);

			if (bDiverged)
			{
				printf("Wellll fuck not sure what to do now...\n");
				return true;
			}

			m_Controller[i].ff = kff;

			m_Controller[i].fb = Matrix<double, nU, nX>::Zero();
			bool bAllClamped = true;
			int numFree = 0;
			for (int j = 0; j < nU; j++)
			{
				if (!clamped[j])
				{
					bAllClamped = false;
					numFree++;
				}
			}

			if (!bAllClamped)
			{
				MatrixXd QuxFree;
				QuxFree.resize(numFree, nX);
				int cj = 0;
				for (int j = 0; j < nU; j++)
				{
					if (clamped[j])
						continue;
					for (int k = 0; k < nX; k++)
						QuxFree(cj,k) = QuxReg(j,k);
					cj++;
				}

				//MatrixXd fbFree = -Hfree.inverse()*((Hfree.transpose()).inverse()*QuxFree);
				MatrixXd fbFree = -Hfree*QuxFree;

				cj = 0;
				for (int j = 0; j < nU; j++)
				{
					if (clamped[j])
						continue;
					for (int k = 0; k < nX; k++)
						m_Controller[i].fb(j,k) = fbFree(cj,k);
					cj++;
				}
			}

		}
		else
		{

			Matrix<double, nU, nU> L = QuuReg.inverse();

			m_Controller[i].fb = -L*QuxReg;
			m_Controller[i].ff = -L*Qu;
		}

		Vx = Qx + m_Controller[i].fb.transpose()*Quu*m_Controller[i].ff + m_Controller[i].fb.transpose()*Qu + Qux.transpose()*m_Controller[i].ff;
		Vxx = Qxx + m_Controller[i].fb.transpose()*Quu*m_Controller[i].fb + m_Controller[i].fb.transpose()*Qux + Qux.transpose()*m_Controller[i].fb;
		Vxx = 0.5*(Vxx + Vxx.transpose());

#if LOGGING
		m_pEnv->LogBackwardPass(Vx, Vxx, Qx, Qxx, Qux, Quu, Qu, m_Controller[i].fb, m_Controller[i].ff, back_iter);
#endif

		Matrix<double, 1, 1> res;

		res = m_Controller[i].ff.transpose()*Qu;
		(*dV1) += res(0,0);
		res = 0.5*(m_Controller[i].ff.transpose()*Quu*m_Controller[i].ff);
		(*dV2) += res(0,0);
	}

	return false;
}

void iLQG::UpdateDerivatives(const mjModel* m, const mjData* dmain, int iter)
{
	Matrix<double, nX, nX> A = Matrix<double, nX, nX>::Zero();
	Matrix<double, nX, nU> B = Matrix<double, nX, nU>::Zero();
	Cost C;
	m_pEnv->Jacobian(m, dmain, &A, &B, &C, iter == N);
	m_pEnv->IntegrateJacobians(dmain, &A, &B);

#if LOGGING
	m_pEnv->LogForwardPass(A, B, C.lx, C.lxx, C.lux, C.luu, C.lu);
#endif

	m_Derivs[iter].A = A;
	m_Derivs[iter].B = B;
	m_Derivs[iter].C = C;
}

bool iLQG::IncreaseLambda(double* lambda, double *dlambda)
{
	(*dlambda) *= lambda_factor;

	if ((*dlambda) < lambda_factor)
		(*dlambda) = lambda_factor;

	(*lambda) *= (*dlambda);
	if ((*lambda) < lambda_min)
		(*lambda) = lambda_min;
	if ((*lambda) > lambda_max)
		return true;

	return false;
}

bool iLQG::DecreaseLambda(double* lambda, double *dlambda)
{
	//decrease lambda
	(*dlambda) /= lambda_factor;
	if ((*dlambda) > 1/lambda_factor)
		(*dlambda) = 1/lambda_factor;
	(*lambda) *= (*dlambda);
	if ((*lambda) < lambda_min)
	{
		(*lambda) = 0;
		return true; //not sure if this will ever be needed
	}
	return false;
}

void iLQG::ClampAction(Matrix<double, nU, 2> bounds, Matrix<double, nU, 1>* out)
{
	for (int i = 0; i < nU; i++)
	{
		if ((*out)(i,0) > bounds(i,1))
			(*out)(i,0) = bounds(i,1);
		if ((*out)(i,0) < bounds(i,0))
			(*out)(i,0) = bounds(i,0);
	}
}

void iLQG::IsClamped(Matrix<double, nU, 1> kff, Matrix<double, nU, 1> grad, Matrix<double, nU, 2> bounds, vector<bool>* isClamped)
{
	for (int i = 0; i < nU; i++)
	{
		if (fabs(kff(i,0) - bounds(i,1)) < 1e-6 && grad(i,0) < 0.0)
			(*isClamped)[i] = true;
		else if (fabs(kff(i,0) - bounds(i,0)) < 1e-6 && grad(i,0) > 0.0)
			(*isClamped)[i] = true;
		else
			(*isClamped)[i] = false;
	}
}

bool iLQG::BoxQP(Matrix<double, nU, nU> H, Matrix<double, nU, 1> g, Matrix<double, nU, 1> next_kff, Matrix<double, nU, 2> bounds, Matrix<double, nU, 1>* kff, MatrixXd* Hfree, vector<bool>* clamped)
{

	int maxIter = 100; //max number of iterations
	double minGrad = 1e-8; //min norm of non-fixed gradient
	double minRelImprove = 1e-8;  //minimum relative improvement
	double stepDec = 0.6;  //factor for decreasing the step size
	double minStep = 1e-22; //minimal stepsize for linesearch
	double Armijo = 0.1; //Armijo parameter

	(*kff) = next_kff;
	ClampAction(bounds, kff);

	Matrix<double, 1, 1> value = (*kff).transpose()*g + 0.5*(*kff).transpose()*H*(*kff);
	Matrix<double, 1, 1> oldValue = value;

	int iter = 0;
	while (iter++ < maxIter)
	{
		if (iter > 1 && (oldValue(0,0) - value(0,0)) < minRelImprove*fabs(oldValue(0,0)))
		{
			//minimum relative improvement
			//printf("Complete: Minimum relative improvement\n");
			return false;
		}
		oldValue = value;

		//get gradient
		Matrix<double, nU, 1> grad = g + H*(*kff);

		//find clamped dimensions
		vector<bool> oldClamped = (*clamped);
		IsClamped((*kff), grad, bounds, clamped);

		bool allClamped = true;
		bool factorize = false;
		int numFreeDims = 0;
		for (int j = 0; j < nU; j++)
		{
			if (!(*clamped)[j])
			{
				numFreeDims++;
				allClamped = false;
			}

			if ((*clamped)[j] != oldClamped[j])
				factorize = true;
		}

//		printf("Iter: %d\n",iter);
//		cout << "Kff orig: \n" << (*kff).transpose() << endl;
//		for (int i = 0; i < nU; i++)
//			(*clamped)[i] ? printf("true,") : printf("false,");
//		printf("\n");

		if (allClamped)
		{
			printf("Complete: All outputs are clamped\n");
			return false;
		}

		(*Hfree).resize(numFreeDims, numFreeDims);

		MatrixXd subH = (*Hfree);
		if (factorize || iter == 1)
		{
			int ci = 0;
			for (int i = 0; i < nU; i++)
			{
				if ((*clamped)[i])
					continue;

				int cj = 0;
				for (int j = 0; j < nU; j++)
				{
					if ((*clamped)[j])
						continue;
					subH(ci, cj) = H(i,j);
					cj++;
				}
				ci++;
			}

			(*Hfree) = subH.inverse();
			//cout << "Hfree: " << (*Hfree) << endl;
			//need to check conditioning here... and return true
		}

		//check gradient norm of free dims
		MatrixXd gradfree;
		gradfree.resize(numFreeDims,1);
		int ci = 0;
		for (int i = 0; i < nU; i++)
		{
			if ((*clamped)[i])
				continue;

			gradfree(ci,0) = grad(i,0);
			ci++;
		}

		double gradNorm = gradfree.norm();
		if (gradNorm < minGrad)
		{
			//printf("Reached minimum gradient\n");
			return false;
		}


		//get search direction
		Matrix<double, nU, 1> kff_zeroed = (*kff);
		for (int i = 0; i < nU; i++)
		{
			if (!(*clamped)[i])
				kff_zeroed(i,0) = 0.0;
		}
		Matrix<double, nU, 1> gradClamped = g + H*kff_zeroed;
		MatrixXd kff_free;
		kff_free.resize(numFreeDims, 1);
		MatrixXd grad_free = kff_free;

		ci = 0;
		for (int i = 0; i < nU; i++)
		{
			if ((*clamped)[i])
				continue;
			kff_free(ci,0) = (*kff)(i,0);
			grad_free(ci,0) = gradClamped(i,0);
			ci++;
		}

		Matrix<double, nU, 1> search_dir = Matrix<double, nU, 1>::Zero();

		//MatrixXd free_search_dir = -(*Hfree).inverse()*(((*Hfree).transpose()).inverse()*grad_free) - kff_free;
		MatrixXd free_search_dir = -(*Hfree)*grad_free - kff_free;
		ci = 0;
		for (int i = 0; i < nU; i++)
		{
			if ((*clamped)[i])
				continue;
			search_dir(i,0) = free_search_dir(ci,0);
			ci++;
		}

		//sanity check for descent direction
		Matrix<double, 1, 1> sdotg = search_dir.transpose()*grad;
		if (sdotg(0,0) > 0.0)
		{
			printf("SdotG: %f\n", sdotg(0,0));
			cout << "Kff orig: \n" << (*kff) << endl;
			for (int i = 0; i < nU; i++)
				(*clamped)[i] ? printf("true,") : printf("false,");
			printf("\n");
			cout << "Gradient: \n" << grad << endl;
			cout << "Search Direction: \n" << search_dir << endl;
			cout << "Kff free: \n" << kff_free << endl;
			cout << "H free: \n" << (*Hfree) << endl;
			printf("Oh shit this isn't supposed to happen!!!!\n");
			return true;
		}

		//linesearch
		double step = 1.0;
		Matrix<double, nU, 1> kff_c = (*kff) + step*search_dir;
		ClampAction(bounds, &kff_c);
		Matrix<double, 1, 1> value_c = kff_c.transpose()*g + 0.5*kff_c.transpose()*H*kff_c;
		while ((value_c(0,0) - oldValue(0,0))/(step*sdotg(0,0)) < Armijo)
		{
			step *= stepDec;
			kff_c = (*kff) + step*search_dir;
			ClampAction(bounds, &kff_c);
			value_c = kff_c.transpose()*g + 0.5*kff_c.transpose()*H*kff_c;
			if (step < minStep)
			{
				//printf("Min step size reached\n");
				break;
			}
		}

		(*kff) = kff_c;
		value = value_c;
	}
	return false;
}

