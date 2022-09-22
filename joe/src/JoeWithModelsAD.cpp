#include "../include/JoeWithModelsAD.h"
#include "../../common/include/myMem.h"
#include "adolc/adolc_sparse.h"

//*************** Initialize and run adjoint solver ************
void JoeWithModels_AD::runAdjoint(int rest)
//***************************************************************
{
    if(mpi_rank==0)
    {
        cout<<"======================="<<endl;
        cout<<"ADJOINT SOLVER DEPLOYED"<<endl;
        cout<<"======================="<<endl;
    }

    // read mesh or restart file
    if(rest) initializeFromRestartFile(getStringParam("RESTART"));

	// Make sure ghost cells are populated
	updateCvDataG1G2(rho, REPLACE_DATA);
    updateCvDataG1G2(rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(rhoE, REPLACE_DATA);

    for(int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
       	updateCvDataG1G2(scalarTranspEqVector[iScal].phi, REPLACE_DATA);

    updateCvDataG1G2(psi_rho, REPLACE_DATA);
    updateCvDataG1G2(psi_rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(psi_rhoE, REPLACE_DATA);

    initializeAdjoint();

    //initialize_turb_adjoint();
    //initialize_comb_adjoint();

    // Exercise some caution in this. If you are using any variable other than
    // the five conserved variables, make sure you work out the dependencies
    calcFunctionalDerivative();

    string tIntName = getStringParam("TIME_INTEGRATION");
 
    if      (tIntName == "FORWARD_EULER")          solveADSystemExplicitEuler();
    else if (tIntName == "RK")                     solveADSystemExplicitRK();
    else if (tIntName == "BACKWARD_EULER")         solveADSystemImplicitEuler();
    else if (tIntName == "BACKWARD_EULER_COUPLED") solveADSystemImplicitEuler_Coupled();
    else
    {
       if (mpi_rank == 0)
       {
         cerr << "ERROR: wrong time integration scheme specified !" << endl;
         cerr << "available integration schemes are: FORWARD_EULER, RK, BACKWARD_EULER, BDF2" << endl;
       }
    }
}

//*************** Derivative of objective functional ************
void JoeWithModels_AD::calcFunctionalDerivative()
//***************************************************************
{
   int tag;
   double dependent_var;
   double functional_total;
   int nScal = scalarTranspEqVector.size();


  //START Gradient of cost function wrt flow
   tag = mpi_rank;
   trace_on(tag, 1);

    initialize_adoubles();


  // Independent variables

     for (int icv = 0; icv < ncv_ggff; icv++)
        { rho_AD[icv]      <<= rho[icv] ;
          rhou_AD[icv][0]  <<= rhou[icv][0] ;
          rhou_AD[icv][1]  <<= rhou[icv][1] ;
          rhou_AD[icv][2]  <<= rhou[icv][2] ;
          rhoE_AD[icv]     <<= rhoE[icv] ;
	     for(int iScal=0; iScal<nScal; iScal++)
		scalarTranspEqVector_AD[iScal].phi[icv] <<= scalarTranspEqVector[iScal].phi[icv];
        }

    //cout<<"NOT CALLING INITIAL HOOKS"<<endl;
    initialHookScalarRansTurbModel_AD();
    initialHookScalarRansCombModel_AD(1);

    
    //TBD_AD: For this to be complete, calcRansturbvisc should also be added
    calcStateVariables_AD(rho_AD, rhou_AD, rhoE_AD);

    calcMaterialProperties_AD(rho_AD, rhou_AD, rhoE_AD);

    setBC_AD(rho_AD, rhou_AD, rhoE_AD);
    
   // Functional  ---> To be defined by the user

      calcFunctional_AD(rho_AD,rhou_AD,rhoE_AD) ;

   // Dependent variable

      functional >>= dependent_var;

      destroy_adoubles();

      trace_off();

      if(mpi_rank==0) {
      cout<<"==============================================="<<endl;
      cout<<"Tape stats for Functional derivative evaluation"<<endl;
      cout<<"==============================================="<<endl;
      print_tapestats(tag);
      }

      // Calculate dJ/dU using reverse mode

      double *dummy_vector = new double[1];
      dummy_vector[0] = 1.0;

      reverse(tag, 1,(5+nScal)*ncv_ggff, 0, dummy_vector, functional_der);

      MPI_Allreduce(&dependent_var, &functional_total, 1, MPI_DOUBLE, MPI_SUM, mpi_comm);

	  if(mpi_rank==0){
		  printf("Cost Function J = %.16e\n", functional_total);
		  FILE *fpj;
		  if ((fpj=fopen("objective.dat","wt"))==NULL){
			  printf("Error: cannot open file objective.dat \n");
			  throw(-1);
		  }
		  fprintf(fpj, "%.12le", functional_total);
		  fclose(fpj);
	  }
      removeTape(tag, 1);

      delete  [] dummy_vector ;

      /*
  int nScal = scalarTranspEqVector.size();

     for (int i=0; i<(5+nScal)*ncv ; i++)
	functional_der[i]=0.;

     for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;

      if (getParam(param, zone->getName()))
      {
        string name = zone->getName();
         if (name == "outlet")
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            assert( icv0 >= 0 );
              double nVec[3];
              double area = normVec3d(nVec, fa_normal[ifa]);
	      functional_der[(5+nScal)*(icv0+1)-1] = fa_normal[ifa][0];
       }
      }
     }
    }
    */

}

//*************** Derivative of objective functional ************
void JoeWithModels_AD::calcResidualDerivative(double (*A)[5][5], double ***AScal, int flagImplicit)
//***************************************************************
{

  int nScal = scalarTranspEqVector.size();

 // Flow residual R
  REALA *rhs_rho_AD       = new REALA[ncv_gg];
  REALA (*rhs_rhou_AD)[3] = new REALA[ncv_gg][3];
  REALA *rhs_rhoE_AD      = new REALA[ncv_gg];

  REALAS **rhs_rhoScal_AD = NULL;
  if (nScal > 0) getMem2D(&rhs_rhoScal_AD, 0, nScal-1, 0, ncv_gg-1, (char*)("rhs_rhoScal_AD"));


   int tag;
   double dependent_var;

  //START Gradient of Residual wrt flow
   tag = mpi_rank;
   trace_on(tag, 1);

    initialize_adoubles();

  // Independent variables

     for (int icv = 0; icv < ncv_gg; icv++)
        { rho_AD[icv]      <<= rho[icv] ;
          rhou_AD[icv][0]  <<= rhou[icv][0] ;
          rhou_AD[icv][1]  <<= rhou[icv][1] ;
          rhou_AD[icv][2]  <<= rhou[icv][2] ;
          rhoE_AD[icv]     <<= rhoE[icv] ;
	     for(int iScal=0; iScal<nScal; iScal++)
		scalarTranspEqVector_AD[iScal].phi[icv] <<= scalarTranspEqVector[iScal].phi[icv];
        }

    initialHookScalarRansTurbModel_AD();
    initialHookScalarRansCombModel_AD(0);

   // Jacobian

      calcResidual_AD(rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD, rhs_rhoScal_AD, rho_AD, rhou_AD, rhoE_AD, A, AScal, flagImplicit);

    //Dependent variables
    
    for (int icv = 0; icv < ncv_gg; icv++)
        { rhs_rho_AD[icv]      >>= dependent_var ;
	//if(fabs(dependent_var)>1.e-8 && icv<ncv_g) cout<<0<<" "<<mpi_rank<<" "<<x_cv[icv][0]<<" "<<x_cv[icv][2]<<" "<<dependent_var<<endl;
          rhs_rhou_AD[icv][0]  >>= dependent_var ;
	//if(fabs(dependent_var)>1.e-8 && icv<ncv_g) cout<<1<<" "<<mpi_rank<<" "<<x_cv[icv][0]<<" "<<x_cv[icv][2]<<" "<<dependent_var<<endl;
          rhs_rhou_AD[icv][1]  >>= dependent_var ;
	//if(fabs(dependent_var)>1.e-8 && icv<ncv_g) cout<<2<<" "<<mpi_rank<<" "<<x_cv[icv][0]<<" "<<x_cv[icv][2]<<" "<<dependent_var<<endl;
          rhs_rhou_AD[icv][2]  >>= dependent_var ;
	//if(fabs(dependent_var)>1.e-8 && icv<ncv_g) cout<<3<<" "<<mpi_rank<<" "<<x_cv[icv][0]<<" "<<x_cv[icv][2]<<" "<<dependent_var<<endl;
          rhs_rhoE_AD[icv]     >>= dependent_var ;
	//if(fabs(dependent_var)>1.e-8 && icv<ncv_g) cout<<4<<" "<<mpi_rank<<" "<<x_cv[icv][0]<<" "<<x_cv[icv][2]<<" "<<dependent_var<<endl;
	     for(int iScal=0; iScal<nScal; iScal++) {
		rhs_rhoScal_AD[iScal][icv] >>= dependent_var;
	//if(fabs(dependent_var)>1.e-8 && icv<ncv_g) cout<<5+iScal<<" "<<mpi_rank<<" "<<x_cv[icv][0]<<" "<<x_cv[icv][2]<<" "<<dependent_var<<endl;
		}
        }
	

      //destroy_adoubles();

      trace_off();

      if(mpi_rank==0) {
      cout<<"==============================================="<<endl;
      cout<<"Tape stats for Residual derivative evaluation"<<endl;
      cout<<"==============================================="<<endl;
      print_tapestats(tag);
      }


    // ---------------------------------------------------------------------------------
    // calculate residual
    // ---------------------------------------------------------------------------------

    double *myResidual = new double[5+nScal];
    double *Residual   = new double[5+nScal];

    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }

    for (int icv = 0; icv < ncv; icv++) {
        myResidual[0] += fabs(rhs_rho_AD[icv].value());
        myResidual[1] += fabs(rhs_rhou_AD[icv][0].value());
        myResidual[2] += fabs(rhs_rhou_AD[icv][1].value());
        myResidual[3] += fabs(rhs_rhou_AD[icv][2].value());
        myResidual[4] += fabs(rhs_rhoE_AD[icv].value());
        for (int iScal = 0; iScal < nScal; iScal++)
        myResidual[5+iScal] += fabs(rhs_rhoScal_AD[iScal][icv].value());
    }

    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

    showResidue(Residual,-1000000);


      delete [] rhs_rhoScal_AD ;
      delete [] rhs_rhoE_AD ;
      delete [] rhs_rhou_AD ;
      delete [] rhs_rho_AD ;

}


//*************** Calculate derivatives of residual, BCs etc  ************
void JoeWithModels_AD::calcResidual_AD(REALA *rhs_rho_AD, 
 REALA (*rhs_rhou_AD)[3],REALA *rhs_rhoE_AD, REALAS **rhs_rhoScal_AD, REALQ *rho_AD, 
 REALQ (*rhou_AD)[3], REALQ *rhoE_AD, double (*A)[5][5], double ***AScal, int flagImplicit)
//************************************************************************
{

    calcStateVariables_AD(rho_AD, rhou_AD, rhoE_AD);

    calcMaterialProperties_AD(rho_AD, rhou_AD, rhoE_AD);
	
    setBC_AD(rho_AD, rhou_AD, rhoE_AD);

    if(mu_ref>0.0 || sndOrder==true)
    calcCv2Grad_AD(grad_u, vel, limiterNavierS, sos, epsilonSDWLS);

    if(mu_ref>0.0)
    calcRansTurbViscMuet_AD(rho_AD, rhou_AD);

    calcRhs_AD(rhs_rho_AD, rhs_rhou_AD, rhs_rhoE_AD, rhs_rhoScal_AD, rho_AD, rhou_AD, rhoE_AD, A, AScal, flagImplicit);

}



//*************** Derivative of objective functional ************
void JoeWithModels_AD::calcResidualDerivative_Coupled(double ***A, int flagImplicit)
//***************************************************************
{

  int nScal = scalarTranspEqVector.size();

  REALA **rhs_AD = NULL;

   //getMem2D(&rhs_AD, 0, ncv_gg-1,     0, 5+nScal-1,               "JoeWithModels_AD::calcResidualDerivativeCoupled -> rhs_AD", true);
   getMem2D(&rhs_AD, 0, ncv_gg-1,     0, 5+nScal-1,               (char*)("rhs_AD"));

   int tag;
   double dependent_var;

  //START Gradient of Residual wrt flow
   tag = mpi_rank;
   trace_on(tag, 1);

    initialize_adoubles();

  // Independent variables

     for (int icv = 0; icv < ncv_gg; icv++)
        { rho_AD[icv]      <<= rho[icv] ;
          rhou_AD[icv][0]  <<= rhou[icv][0] ;
          rhou_AD[icv][1]  <<= rhou[icv][1] ;
          rhou_AD[icv][2]  <<= rhou[icv][2] ;
          rhoE_AD[icv]     <<= rhoE[icv] ;
	     for(int iScal=0; iScal<nScal; iScal++)
		scalarTranspEqVector_AD[iScal].phi[icv] <<= scalarTranspEqVector[iScal].phi[icv];
        }

    initialHookScalarRansTurbModel_AD();
    initialHookScalarRansCombModel_AD(0);

   // Jacobian

      calcResidual_AD_Coupled(rhs_AD, rho_AD, rhou_AD, rhoE_AD, A, flagImplicit);

    /*for (int icv = 0; icv < ncv_gg; icv++)
       for (int i = 0; i < 5+nScal; i++) 
         rhs_AD[icv][i]     = rhs_AD[icv][i]/cv_volume[icv];
     */

    //Dependent variables
    for (int icv = 0; icv < ncv_gg; icv++)
       for (int i = 0; i < 5+nScal; i++) {
         rhs_AD[icv][i]      >>= dependent_var ;
            //    if(fabs(dependent_var)>1.e-8 && icv<ncv_g) cout<<i<<" "<<icv<<" "<<x_cv[icv][0]<<" "<<x_cv[icv][1]<<" "<<x_cv[icv][2]<<" "<<dependent_var<<endl;
       }

      //destroy_adoubles();
      trace_off();

    // ---------------------------------------------------------------------------------
    // calculate residual
    // ---------------------------------------------------------------------------------
 
    double *myResidual = new double[5+nScal];
    double *Residual   = new double[5+nScal];

    for (int i = 0; i < 5+nScal; i++)
    {
      myResidual[i] = 0.0;
      Residual[i] = 0.0;
    }

    for (int icv = 0; icv < ncv; icv++)
      for (int i = 0; i < 5+nScal; i++)
        myResidual[i] += fabs(rhs_AD[icv][i].value());
       // myResidual[i] += fabs(rhs_AD[icv][i].value()*cv_volume[icv]);

    MPI_Allreduce(myResidual, Residual, 5+nScal, MPI_DOUBLE, MPI_SUM, mpi_comm);

    showResidue(Residual,-1000000);

        



      if(mpi_rank==0) {
      cout<<"==============================================="<<endl;
      cout<<"Tape stats for Residual derivative evaluation"<<endl;
      cout<<"==============================================="<<endl;
      print_tapestats(tag);
      }

      delete [] rhs_AD ;

}


//*************** Calculate derivatives of residual, BCs etc  ************
void JoeWithModels_AD::calcResidual_AD_Coupled(REALA **rhs_AD, REALQ *rho_AD, 
 REALQ (*rhou_AD)[3], REALQ *rhoE_AD, double ***A, int flagImplicit)
//************************************************************************
{

    calcStateVariables_AD(rho_AD, rhou_AD, rhoE_AD);

    calcMaterialProperties_AD(rho_AD, rhou_AD, rhoE_AD);
	
    setBC_AD(rho_AD, rhou_AD, rhoE_AD);

    if(mu_ref>0.0 || sndOrder==true)
    calcCv2Grad_AD(grad_u, vel, limiterNavierS, sos, epsilonSDWLS);

    if(mu_ref>0.0)
    calcRansTurbViscMuet_AD(rho_AD, rhou_AD);

    calcRhs_AD_Coupled(rhs_AD, rho_AD, rhou_AD, rhoE_AD, A, flagImplicit);

}


//*************** Explicit Euler solution of AD system  ************
void JoeWithModels_AD::solveADSystemExplicitEuler()
//******************************************************************
{
    //setDebugCv();

    int nScal = scalarTranspEqVector.size();
    int nAdjvars = 5+nScal;


    double (*dummyA)[5][5] = NULL;
    double ***dummyAScal = NULL;

    double *adj_residue   = new double[nAdjvars];
    double *my_adj_residue = new double[nAdjvars];


    calcResidualDerivative(dummyA,dummyAScal,0);

    //dumpDebugState("1");

   // Initialize Scalar adjoint variables

    initialize_turb_adjoint();
    initialize_comb_adjoint();

    int tag=mpi_rank;

	 for (int icv = 0; icv < ncv_gg; icv++)
    {
       int icv_i = nAdjvars*icv; 
       adj_vars[icv_i]   = psi_rho[icv];
       adj_vars[icv_i+1] = psi_rhou[icv][0];
       adj_vars[icv_i+2] = psi_rhou[icv][1];
       adj_vars[icv_i+3] = psi_rhou[icv][2];
       adj_vars[icv_i+4] = psi_rhoE[icv];
           for(int iScal=0; iScal<nScal; iScal++) {
                adj_vars[icv_i+5+iScal] = scalarTranspEqVector_psi[iScal].phi[icv];
	   }
    } 


  int done = 0;
  int iter = 0;

  if (nsteps == 0)    done = 1;

  while (done != 1)
  {
    iter++;
    double dtMin = calcDt(cfl);

   // V^T dR/dU
    fos_reverse(tag, nAdjvars * ncv_gg, nAdjvars * ncv_gg, adj_vars, residual_der);

    // ---------------------------------------------------------------
    // advance adjoint vector
    // ---------------------------------------------------------------

    double tmp, dadj_vars[nAdjvars] ;

    for (int i = 0; i < nAdjvars; i++){
           adj_residue[i]=0.;
        my_adj_residue[i]=0.;
      }   


    for (int icv = 0; icv < ncv; icv++)
    {
      tmp = local_dt[icv]/cv_volume[icv];
      
      int ii = nAdjvars*icv ;

      for (int i = 0; i < 5; i++){
         dadj_vars[i] = tmp*(residual_der[ii]+functional_der[ii]) ;
         my_adj_residue[i] += fabs(dadj_vars[i]);
         adj_vars[ii] += dadj_vars[i];
	 ii++;
      }
	for (int iScal = 0; iScal < nScal; iScal++) {
         dadj_vars[5+iScal] = tmp*(residual_der[ii]+functional_der[ii]);
         my_adj_residue[5+iScal] += fabs(dadj_vars[5+iScal]);
         adj_vars[ii] += dadj_vars[5+iScal];
         ii++ ;
         }

    }

    MPI_Allreduce(my_adj_residue, adj_residue, nAdjvars, MPI_DOUBLE, MPI_SUM, mpi_comm);


     for (int icv = 0; icv < ncv; icv++)
    {
       int ii = nAdjvars*icv;
       psi_rho[icv]    =adj_vars[ii];
       psi_rhou[icv][0]=adj_vars[ii+1];
       psi_rhou[icv][1]=adj_vars[ii+2];
       psi_rhou[icv][2]=adj_vars[ii+3];
       psi_rhoE[icv]   =adj_vars[ii+4];
        for(int iScal=0; iScal<nScal; iScal++)
                scalarTranspEqVector_psi[iScal].phi[icv]= adj_vars[ii+5+iScal] ;

    }

    updateCvDataG1G2(psi_rho, REPLACE_DATA);
    updateCvDataG1G2(psi_rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(psi_rhoE, REPLACE_DATA);

    for (int icv = 0; icv < ncv_gg; icv++)
    {
       int ii = nAdjvars*icv;
       adj_vars[ii]   = psi_rho[icv];
       adj_vars[ii+1] = psi_rhou[icv][0];
       adj_vars[ii+2] = psi_rhou[icv][1];
       adj_vars[ii+3] = psi_rhou[icv][2];
       adj_vars[ii+4] = psi_rhoE[icv];
    }


        for(int iScal=0; iScal<nScal; iScal++) {
              double *psi_phi = scalarTranspEqVector_psi[iScal].phi;
              updateCvDataG1G2(psi_phi, REPLACE_DATA);

                for (int icv = 0; icv < ncv_gg; icv++){
                        int ii = nAdjvars*icv+5+iScal;
                        adj_vars[ii] = psi_phi[icv];
               		scalarTranspEqVector_psi[iScal].phi[icv]=psi_phi[icv];

                }

        }




     if(iter%check_interval==0) {

              showResidue(adj_residue,iter);

              calcFunctionalGradient(adj_vars);

      }

      writeData(iter);

      if ((write_restart > 0) && (iter % write_restart == 0))
      writeRestart(iter);

    if ((nsteps != -1) && (iter >= nsteps))
      done = 1;
  }

     writeRestart();
     finalHook_AD();



}

//*************** Explicit Euler solution of AD system  ************
void JoeWithModels_AD::solveADSystemExplicitRK()
//******************************************************************
{
    //setDebugCv();

    int nScal = scalarTranspEqVector.size();
    int nAdjvars = 5+nScal;
    int nStage = 4;

    double (*dummyA)[5][5] = NULL;
    double ***dummyAScal = NULL;

    double *adj_residue   = new double[nAdjvars];
    double *my_adj_residue = new double[nAdjvars];

  double *psi_rho0 = new double[ncv];
  double (*psi_rhou0)[3] = new double[ncv][3];
  double *psi_rhoE0 = new double[ncv];
  double **psi_rhoScal0 = NULL;
  if (nScal > 0) getMem2D(&psi_rhoScal0, 0, nScal-1, 0, ncv-1, (char*)("psi_rhoScal0"));


    calcResidualDerivative(dummyA,dummyAScal,0);


   // Initialize Scalar adjoint variables

    initialize_turb_adjoint();
    initialize_comb_adjoint();

    int tag=mpi_rank;

	 for (int icv = 0; icv < ncv_gg; icv++)
    {
       int icv_i = nAdjvars*icv; 
       adj_vars[icv_i]   = psi_rho[icv] ;
       adj_vars[icv_i+1] = psi_rhou[icv][0];
       adj_vars[icv_i+2] = psi_rhou[icv][1];
       adj_vars[icv_i+3] = psi_rhou[icv][2];
       adj_vars[icv_i+4] = psi_rhoE[icv];
           for(int iScal=0; iScal<nScal; iScal++)
                adj_vars[icv_i+5+iScal] = scalarTranspEqVector_psi[iScal].phi[icv];
    } 


  int done = 0;
  int iter = 0;

  if (nsteps == 0)    done = 1;

  while (done != 1)
  {
    iter++;
    double dtMin = calcDt(cfl);

    for (int icv = 0; icv < ncv; icv++)
    {
       psi_rho0[icv]     = psi_rho[icv] ;
       psi_rhou0[icv][0] = psi_rhou[icv][0];
       psi_rhou0[icv][1] = psi_rhou[icv][1];
       psi_rhou0[icv][2] = psi_rhou[icv][2];
       psi_rhoE0[icv]    = psi_rhoE[icv];
           for(int iScal=0; iScal<nScal; iScal++)
                psi_rhoScal0[iScal][icv] = scalarTranspEqVector_psi[iScal].phi[icv];
    }

    for (int i = 0; i < nAdjvars; i++){
           adj_residue[i]=0.;
        my_adj_residue[i]=0.;
      }   

    for (int iStage = 1; iStage <= nStage ; iStage++) {

   // V^T dR/dU
    fos_reverse(tag, nAdjvars * ncv_gg, nAdjvars * ncv_gg, adj_vars, residual_der);

    // ---------------------------------------------------------------
    // advance adjoint vector
    // ---------------------------------------------------------------

    double tmp;

     for (int icv = 0; icv < ncv; icv++)
    {
       tmp = local_dt[icv]/(cv_volume[icv]*(nStage-iStage+1.));
       int ii = nAdjvars*icv;
       psi_rho[icv]    =psi_rho0[icv]    +tmp*(residual_der[ii]  +functional_der[ii]);
       psi_rhou[icv][0]=psi_rhou0[icv][0]+tmp*(residual_der[ii+1]+functional_der[ii+1]);
       psi_rhou[icv][1]=psi_rhou0[icv][1]+tmp*(residual_der[ii+2]+functional_der[ii+2]);
       psi_rhou[icv][2]=psi_rhou0[icv][2]+tmp*(residual_der[ii+3]+functional_der[ii+3]);
       psi_rhoE[icv]   =psi_rhoE0[icv]   +tmp*(residual_der[ii+4]+functional_der[ii+4]);

	if(iStage==nStage) {
        my_adj_residue[0] += tmp*fabs(residual_der[ii]    +functional_der[ii]);
        my_adj_residue[1] += tmp*fabs(residual_der[ii+1]  +functional_der[ii+1]);
        my_adj_residue[2] += tmp*fabs(residual_der[ii+2]  +functional_der[ii+2]);
        my_adj_residue[3] += tmp*fabs(residual_der[ii+3]  +functional_der[ii+3]);
        my_adj_residue[4] += tmp*fabs(residual_der[ii+4]  +functional_der[ii+4]);
	}
        for(int iScal=0; iScal<nScal; iScal++) {
                scalarTranspEqVector_psi[iScal].phi[icv]= psi_rhoScal0[iScal][icv]+tmp*(residual_der[ii+5+iScal]+functional_der[ii+5+iScal]) ;
                if(iStage==nStage) my_adj_residue[5+iScal] += tmp*fabs(residual_der[ii+5+iScal]+functional_der[ii+5+iScal]);
        }
    }

    // POPULATE GHOSTS AND TRANSFER

    updateCvDataG1G2(psi_rho, REPLACE_DATA);
    updateCvDataG1G2(psi_rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(psi_rhoE, REPLACE_DATA);
    for (int icv = 0; icv < ncv_gg; icv++)
    {
       int ii = nAdjvars*icv;
       adj_vars[ii]   = psi_rho[icv];
       adj_vars[ii+1] = psi_rhou[icv][0];
       adj_vars[ii+2] = psi_rhou[icv][1];
       adj_vars[ii+3] = psi_rhou[icv][2];
       adj_vars[ii+4] = psi_rhoE[icv];
    }
        for(int iScal=0; iScal<nScal; iScal++) {
              double *psi_phi = scalarTranspEqVector_psi[iScal].phi;
              updateCvDataG1G2(psi_phi, REPLACE_DATA);
                for (int icv = 0; icv < ncv_gg; icv++){
                        int ii = nAdjvars*icv+5+iScal;
                        adj_vars[ii] = psi_phi[icv];
               		scalarTranspEqVector_psi[iScal].phi[icv]=psi_phi[icv];
                }
        }
    }

    MPI_Allreduce(my_adj_residue, adj_residue, nAdjvars, MPI_DOUBLE, MPI_SUM, mpi_comm);

     if(iter%check_interval==0) {

              showResidue(adj_residue,iter);

              calcFunctionalGradient(adj_vars);

      }

      writeData(iter);

      if ((write_restart > 0) && (iter % write_restart == 0))
      writeRestart(iter);

    if ((nsteps != -1) && (iter >= nsteps))
      done = 1;
  }

     writeRestart();
     finalHook_AD();

    delete [] psi_rhoScal0 ;
    delete [] psi_rhoE0 ;
    delete [] psi_rhou0 ;
    delete [] psi_rho0 ;


}


//********** Implicit Euler solution of adjoint system  ************
void JoeWithModels_AD::solveADSystemImplicitEuler()
//******************************************************************
{
      int nScal = scalarTranspEqVector.size();
      int nAdjvars = 5+nScal;

    double (*A)[5][5]       = new double[nbocv_s][5][5];
    double (*rhs_adj)[5]    = new double[ncv_gg][5];
    double (*dadj_vars)[5]  = new double[ncv_gg][5];

  double ***AScal      = NULL;  if (nScal > 0) getMem3D(&AScal,      0, nScal-1, 0, 5, 0, nbocv_s-1, (char*)("AScal"));
  double **dadjScal    = NULL;  if (nScal > 0) getMem2D(&dadjScal,   0, nScal-1, 0, ncv_gg-1, (char*)("dadjScal"));
  double **rhs_adjScal = NULL;  if (nScal > 0) getMem2D(&rhs_adjScal,0, nScal-1, 0, ncv_gg-1, (char*)("rhs_adjScal"));

    double *adj_residue   = new double[nAdjvars];
    double *my_adj_residue = new double[nAdjvars];
  //-----------------------------------
  // some parameters (This should be replaced later with adjoint specific parameters)
  //------------------------------------

  double underRelax = getDoubleParam("UNDER_RELAXATION", "0.3");

  if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
  {
    ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
    if (mpi_rank == 0)
      cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
              " to parameter map!" << endl;
  }
  int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
  double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
  double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");
    
  if (!checkParam("CFL_RAMP"))
  {
    ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
    if (mpi_rank == 0)
      cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
  } 

  int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");
    
  double resid_energ_th = getDoubleParam("RESID_ENERG_TH", 1.0e-20);

   // Initialize Jacobian matrices to zero
    for (int noc=0; noc<nbocv_s; noc++)   
      for (int i=0; i<5; i++)
      for (int j=0; j<5; j++)
        A[noc][i][j] = 0.0;

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      for (int i = 0; i <= 5; i++)
        for (int noc = 0; noc < nbocv_s; noc++)
          AScal[iScal][i][noc] = 0.0;
    }

    calcResidualDerivative(A,AScal,1);

   // Initialize Scalar adjoint variables

   initialize_turb_adjoint();
   initialize_comb_adjoint();

    for (int icv = 0; icv < ncv_gg; icv++)
    {  
       int icv_i = nAdjvars*icv; 
       adj_vars[icv_i]   = psi_rho[icv] ;
       adj_vars[icv_i+1] = psi_rhou[icv][0];
       adj_vars[icv_i+2] = psi_rhou[icv][1];
       adj_vars[icv_i+3] = psi_rhou[icv][2];
       adj_vars[icv_i+4] = psi_rhoE[icv];
	   for(int iScal=0; iScal<nScal; iScal++)
                adj_vars[icv_i+5+iScal] = scalarTranspEqVector_psi[iScal].phi[icv];
    }

  int done = 0;
  int iter = 0;

   if (initial_flowfield_output == "YES")
   writeData(0);

  if (nsteps == 0)    done = 1;

  while (done != 1)
  {
    iter++;
    if ((iter >= startIncCFL) && (iter%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
    double dt_min = calcDt(cfl);

   // V^T dR/dU
    int tag=mpi_rank;
    fos_reverse(tag, nAdjvars * ncv_gg, nAdjvars * ncv_gg, adj_vars, residual_der);

    // ---------------------------------------------------------------
    // advance adjoint vector
    // ---------------------------------------------------------------

    /*
    // ******* TBD_AD: JUST FOR EPISTEMIC UQ PROBLEM ***************
    if((iter%10)==0 && mpi_rank==0) cout<<" JUNKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK "<<endl;
	double a         = 1.173958365311995e+03;
        double a2        = a*a;
        double To        = 423.4159779614325;
        double To2       = To*To;
        double To3       = To2*To;
        double To4       = To3*To;
        double base_fact = exp(a/To);
	double q_release = 3e5;
	double alpha     = 1000.;



    for (int icv = 0; icv < ncv; icv++)
    {
        double s_by_f     = alpha*(1.-scalarTranspEqVector[0].phi[icv])*rho[icv]*0.0625*cv_volume[icv];
        double df_by_dT   = base_fact*(-a/To2);
	       df_by_dT  += base_fact*(temp[icv].value()-To)*(a2/To4+2*a/To3);
	adouble gam1 = gamma[icv]-1.;
	adouble dT_by_drho = gam1/rho[icv]/RoM[icv]*(-scalarTranspEqVector[0].phi[icv]*q_release/rho[icv]
			   +0.5*(rhou[icv][0]*rhou[icv][0]+rhou[icv][1]*rhou[icv][1]+rhou[icv][2]*rhou[icv][2])/rho[icv]/rho[icv])-temp[icv]/rho[icv];
	adouble dT_by_drhou = -gam1/rho[icv]/RoM[icv]*rhou[icv][0]/rho[icv];
	adouble dT_by_drhov = -gam1/rho[icv]/RoM[icv]*rhou[icv][1]/rho[icv];
	adouble dT_by_drhow = -gam1/rho[icv]/RoM[icv]*rhou[icv][2]/rho[icv];
	adouble dT_by_drhoE =  gam1/rho[icv]/RoM[icv];
	adouble dT_by_drhoLam =  gam1/rho[icv]/RoM[icv]/rho[icv]*q_release;

         int ii = nAdjvars*icv;
	double mult_const = adj_vars[ii+5]*s_by_f*df_by_dT;

	
         residual_der[ii]   += dT_by_drho.value()*mult_const;
         residual_der[ii+1] += dT_by_drhou.value()*mult_const;
         residual_der[ii+2] += dT_by_drhov.value()*mult_const;
         residual_der[ii+3] += dT_by_drhow.value()*mult_const;
         residual_der[ii+4] += dT_by_drhoE.value()*mult_const;
         residual_der[ii+5] += dT_by_drhoLam.value()*mult_const;
	
    }
    // ******* TBD_AD: JUST FOR EPISTEMIC UQ PROBLEM ***************
    */

    for (int icv = 0; icv < ncv; icv++)
    {

         int ii = nAdjvars*icv;

      for (int i = 0; i < 5; i++){
         rhs_adj[icv][i] = residual_der[ii]+functional_der[ii];
         ii++ ;
      }
         for (int iScal = 0; iScal < nScal; iScal++) {
         rhs_adjScal[iScal][icv] = residual_der[ii]+functional_der[ii];
         ii++ ;
         }

      double tmp = cv_volume[icv]/local_dt[icv];

         for (int i=0; i<5; i++) {
          A[nbocv_i[icv]][i][i] += tmp;    // Add diagonal part ( vol/dt + A )
          rhs_adj[icv][i] *= underRelax;
         }


         for (int iScal = 0; iScal < nScal; iScal++) {
          	AScal[iScal][5][nbocv_i[icv]] += tmp;
          	rhs_adjScal[iScal][icv] *= underRelax;
         }


    }


      // ---------------------------------------------------------------------------------
      // solve lin system 
      // ---------------------------------------------------------------------------------



       for (int iScal = 0; iScal < nScal; iScal++)  {  // prepare rhs and A
   		 solveLinSysScalar(dadjScal[iScal], AScal[iScal][5], rhs_adjScal[iScal],
                        scalarTranspEqVector[iScal].phiZero,
                        scalarTranspEqVector[iScal].phiZeroRel,
                        scalarTranspEqVector[iScal].phiMaxiter,
                        scalarTranspEqVector[iScal].getName());

    		for (int icv = 0; icv < ncv; icv++)
    		{

        		int noc_f = nbocv_i[icv];
        		int noc_l = nbocv_i[icv + 1] - 1;

        		// move the other implicit terms to the RHS
        		for (int noc = noc_f; noc <= noc_l; noc++) {
                 		for (int i=0; i<5; i++)
                        		rhs_adj[icv][i] -= AScal[iScal][i][noc]*dadjScal[iScal][nbocv_v[noc]];

        		}

    		}
    }

      solveCoupledLinSysNS(dadj_vars, A, rhs_adj, zeroAbsLS, zeroRelLS, maxIterLS);  // solve linear system

      // ---------------------------------------------------------------------------------
      // update vec
      // ---------------------------------------------------------------------------------


      for (int i = 0; i < nAdjvars; i++){
	   adj_residue[i]=0.;
	my_adj_residue[i]=0.;
      }


    for (int icv = 0; icv < ncv; icv++)
    {
         int ii = nAdjvars*icv ;
      for (int i = 0; i < 5; i++){
         adj_vars[ii] += dadj_vars[icv][i];
         my_adj_residue[i] += fabs(dadj_vars[icv][i]);
	 ii++;
      }
         for (int iScal = 0; iScal < nScal; iScal++) {
         adj_vars[ii] += dadjScal[iScal][icv];
         my_adj_residue[5+iScal] += fabs(dadjScal[iScal][icv]);
         ii++ ;
         }

      double tmp = cv_volume[icv]/local_dt[icv];

         for (int i=0; i<5; i++)
          A[nbocv_i[icv]][i][i] -= tmp;    // Remove diagonal part ( vol/dt + A )

         for (int iScal = 0; iScal < nScal; iScal++) 
          	AScal[iScal][5][nbocv_i[icv]] -= tmp;
    }

    MPI_Allreduce(my_adj_residue, adj_residue, nAdjvars, MPI_DOUBLE, MPI_SUM, mpi_comm);

   

    for (int icv = 0; icv < ncv; icv++)
    {  
       int ii = nAdjvars*icv;
       psi_rho[icv]    =adj_vars[ii];
       psi_rhou[icv][0]=adj_vars[ii+1];
       psi_rhou[icv][1]=adj_vars[ii+2];
       psi_rhou[icv][2]=adj_vars[ii+3];
       psi_rhoE[icv]   =adj_vars[ii+4];
	for(int iScal=0; iScal<nScal; iScal++)
                scalarTranspEqVector_psi[iScal].phi[icv]= adj_vars[ii+5+iScal] ;

    }

    updateCvDataG1G2(psi_rho, REPLACE_DATA);
    updateCvDataG1G2(psi_rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(psi_rhoE, REPLACE_DATA);

    for (int icv = 0; icv < ncv_gg; icv++)
    {  
       int ii = nAdjvars*icv;
       adj_vars[ii]   = psi_rho[icv];
       adj_vars[ii+1] = psi_rhou[icv][0];
       adj_vars[ii+2] = psi_rhou[icv][1];
       adj_vars[ii+3] = psi_rhou[icv][2];
       adj_vars[ii+4] = psi_rhoE[icv];
    }

	
	for(int iScal=0; iScal<nScal; iScal++) {
	      double *psi_phi = scalarTranspEqVector_psi[iScal].phi;
              updateCvDataG1G2(psi_phi, REPLACE_DATA);

    		for (int icv = 0; icv < ncv_gg; icv++){
       			int ii = nAdjvars*icv+5+iScal;
       			adj_vars[ii] = psi_phi[icv];
		
                }
        }
	

     if(iter%(check_interval)==0) { 
              if ((mpi_rank == 0) && (iter%(check_interval*10) == 0))
        cout << "\ndone step: "<< iter << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

              showResidue(adj_residue,iter);
     	if(iter%(check_interval*1)==0) calcFunctionalGradient(adj_vars);
      }

      copy_turb_adjoint();
      copy_comb_adjoint();

      writeData(iter);

      if ((write_restart > 0) && (iter % write_restart == 0))
      writeRestart(iter);

    if ((nsteps != -1) && (iter >= nsteps))
      done = 1;
  }

     calcFunctionalGradient(adj_vars);
     writeRestart();
     finalHook_AD();

    delete [] A ;
    delete [] rhs_adj;
    delete [] dadj_vars;
}

//********** Implicit Euler solution of adjoint system  ************
void JoeWithModels_AD::solveADSystemImplicitEuler_Coupled()
//******************************************************************
{
      int nScal = scalarTranspEqVector.size();
      for (int iScal = 0; iScal < nScal; iScal++)
    	 scalarTranspEqVector[iScal].coupling = "COUPLED";

      int nAdjvars = 5+nScal;

    //double ***A;         getMem3D(&A,        0, nbocv_s-1, 0, nAdjvars-1, 0, nAdjvars-1, "JoeWithModels_AD::solveADSystemImplicitEuler_Coupled -> A",    true);
  //double **rhs_adj;      getMem2D(&rhs_adj,  0, ncv_gg-1,  0, nAdjvars-1, "JoeWithModels_AD::solveADSystemImplicitEuler_Coupled -> rhs_adj",    true);
  //double **dadj_vars;    getMem2D(&dadj_vars,0, ncv_gg-1,  0, nAdjvars-1, "JoeWithModels_AD::solveADSystemImplicitEuler_Coupled -> dadj_vars",  true);
    double ***A;         getMem3D(&A,        0, nbocv_s-1, 0, nAdjvars-1, 0, nAdjvars-1, (char*)("A"));
    double **rhs_adj;      getMem2D(&rhs_adj,  0, ncv_gg-1,  0, nAdjvars-1, (char*)("rhs_adj"));
    double **dadj_vars;    getMem2D(&dadj_vars,0, ncv_gg-1,  0, nAdjvars-1, (char*)("dadj_vars"));

    double *adj_residue   = new double[nAdjvars];
    double *my_adj_residue = new double[nAdjvars];
  //-----------------------------------
  // some parameters (This should be replaced later with adjoint specific parameters)
  //------------------------------------

  double underRelax = getDoubleParam("UNDER_RELAXATION", "0.3");

  if (!checkParam("LINEAR_SOLVER_NS_TRESHOLDS"))
  {
    ParamMap::add("LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2");    // add default values
    if (mpi_rank == 0)
      cout << "WARNING: added keyword \"LINEAR_SOLVER_NS_TRESHOLDS  MAX_ITER=30  ABS_RESID=1.0e-8  REL_RESID=1.0e-2\"" <<
              " to parameter map!" << endl;
  }
  int maxIterLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getInt("MAX_ITER");
  double zeroAbsLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("ABS_RESID");
  double zeroRelLS = getParam("LINEAR_SOLVER_NS_TRESHOLDS")->getDouble("REL_RESID");
    
  if (!checkParam("CFL_RAMP"))
  {
    ParamMap::add("CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0");    // add default: no increase of CFL number!
    if (mpi_rank == 0)
      cout << "WARNING: added \"CFL_RAMP AFTER_ITER=100  INTERVAL_ITER=10  FACTOR_CFL=1.0  MAX_CFL=100.0\" to parameter map" << endl;
  } 

  int startIncCFL = getParam("CFL_RAMP")->getInt("AFTER_ITER");
  int intervalIncCFL = getParam("CFL_RAMP")->getInt("INTERVAL_ITER");
  double incCFL = getParam("CFL_RAMP")->getDouble("FACTOR_CFL");
  double maxCFL = getParam("CFL_RAMP")->getDouble("MAX_CFL");
    
  double resid_energ_th = getDoubleParam("RESID_ENERG_TH", 1.0e-20);

   // Initialize Jacobian matrices to zero
    for (int noc=0; noc<nbocv_s; noc++)   
      for (int i=0; i<nAdjvars; i++)
        for (int j=0; j<nAdjvars; j++)
           A[noc][i][j] = 0.0;

    calcResidualDerivative_Coupled(A,1);









   // Initialize Scalar adjoint variables

   initialize_turb_adjoint();
   initialize_comb_adjoint();

    for (int icv = 0; icv < ncv_gg; icv++)
    {  
       int icv_i = nAdjvars*icv; 
       adj_vars[icv_i]   = psi_rho[icv] ;
       adj_vars[icv_i+1] = psi_rhou[icv][0];
       adj_vars[icv_i+2] = psi_rhou[icv][1];
       adj_vars[icv_i+3] = psi_rhou[icv][2];
       adj_vars[icv_i+4] = psi_rhoE[icv];
	   for(int iScal=0; iScal<nScal; iScal++)
                adj_vars[icv_i+5+iScal] = scalarTranspEqVector_psi[iScal].phi[icv];
    }

  int done = 0;
  int iter = 0;

   if (initial_flowfield_output == "YES")
   writeData(0);

  if (nsteps == 0)    done = 1;

  while (done != 1)
  {
    iter++;
    if ((iter >= startIncCFL) && (iter%intervalIncCFL == 0) && (cfl < maxCFL))      cfl *= incCFL;
    double dt_min = calcDt(cfl);

   // V^T dR/dU
    int tag=mpi_rank;
    fos_reverse(tag, nAdjvars * ncv_gg, nAdjvars * ncv_gg, adj_vars, residual_der);

    // ---------------------------------------------------------------
    // advance adjoint vector
    // ---------------------------------------------------------------


    for (int icv = 0; icv < ncv; icv++)
    {

         int ii = nAdjvars*icv;

      for (int i = 0; i < 5; i++){
         rhs_adj[icv][i] = residual_der[ii]+functional_der[ii];
         ii++ ;
      }
         for (int iScal = 0; iScal < nScal; iScal++) {
         rhs_adj[icv][5+iScal] = residual_der[ii]+functional_der[ii];
         ii++ ;
         }

      double tmp = cv_volume[icv]/local_dt[icv];

         for (int i=0; i<nAdjvars; i++) {
          A[nbocv_i[icv]][i][i] += tmp;    // Add diagonal part ( vol/dt + A )
          rhs_adj[icv][i] *= underRelax;
         }
    }

      // solve linear system

       solveCoupledLinSysNSCoupled(dadj_vars, A, rhs_adj, zeroAbsLS, zeroRelLS, maxIterLS, nScal); 

      // ---------------------------------------------------------------------------------
      // update vec
      // ---------------------------------------------------------------------------------


      for (int i = 0; i < nAdjvars; i++){
	   adj_residue[i]=0.;
	my_adj_residue[i]=0.;
      }


    for (int icv = 0; icv < ncv; icv++)
    {
         int ii = nAdjvars*icv ;
      for (int i = 0; i < nAdjvars; i++){
         // DEBUG dadj_vars[icv][i] = rhs_adj[icv][i]/underRelax*local_dt[icv]/cv_volume[icv];
         //if(i<5) adj_vars[ii] += dadj_vars[icv][i]; //DEBUG
         adj_vars[ii] += dadj_vars[icv][i];
         my_adj_residue[i] += fabs(dadj_vars[icv][i]);
	 ii++;
      }

      double tmp = cv_volume[icv]/local_dt[icv];

         for (int i=0; i<nAdjvars; i++)
          A[nbocv_i[icv]][i][i] -= tmp;    // Remove diagonal part ( vol/dt + A )

    }

    MPI_Allreduce(my_adj_residue, adj_residue, nAdjvars, MPI_DOUBLE, MPI_SUM, mpi_comm);

    for (int icv = 0; icv < ncv; icv++)
    {  
       int ii = nAdjvars*icv;
       psi_rho[icv]    =adj_vars[ii];
       psi_rhou[icv][0]=adj_vars[ii+1];
       psi_rhou[icv][1]=adj_vars[ii+2];
       psi_rhou[icv][2]=adj_vars[ii+3];
       psi_rhoE[icv]   =adj_vars[ii+4];
	for(int iScal=0; iScal<nScal; iScal++)
                scalarTranspEqVector_psi[iScal].phi[icv]= adj_vars[ii+5+iScal] ;

    }

    updateCvDataG1G2(psi_rho, REPLACE_DATA);
    updateCvDataG1G2(psi_rhou, REPLACE_ROTATE_DATA);
    updateCvDataG1G2(psi_rhoE, REPLACE_DATA);

    for (int icv = 0; icv < ncv_gg; icv++)
    {  
       int ii = nAdjvars*icv;
       adj_vars[ii]   = psi_rho[icv];
       adj_vars[ii+1] = psi_rhou[icv][0];
       adj_vars[ii+2] = psi_rhou[icv][1];
       adj_vars[ii+3] = psi_rhou[icv][2];
       adj_vars[ii+4] = psi_rhoE[icv];
    }

	
	for(int iScal=0; iScal<nScal; iScal++) {
	      double *psi_phi = scalarTranspEqVector_psi[iScal].phi;
              updateCvDataG1G2(psi_phi, REPLACE_DATA);

    		for (int icv = 0; icv < ncv_gg; icv++){
       			int ii = nAdjvars*icv+5+iScal;
       			adj_vars[ii] = psi_phi[icv];
		
                }
        }
	

     if(iter%(check_interval)==0) { 
              if ((mpi_rank == 0) && (iter%(check_interval*10) == 0))
        cout << "\ndone step: "<< iter << ", cfl: " << cfl << ", min. dt: " << dt_min << " time: " << time << endl;

              showResidue(adj_residue,iter);
     	if(iter%(check_interval*1)==0) calcFunctionalGradient(adj_vars);
      }

      copy_turb_adjoint();
      copy_comb_adjoint();

      writeData(iter);

      if ((write_restart > 0) && (iter % write_restart == 0))
      writeRestart(iter);

    if ((nsteps != -1) && (iter >= nsteps))
      done = 1;
  }

     calcFunctionalGradient(adj_vars);
     writeRestart();
     finalHook_AD();

    delete [] A ;
    delete [] rhs_adj;
    delete [] dadj_vars;
}

//*************** Print Stats in tape ************
void JoeWithModels_AD::print_tapestats(int tag)
//************************************************
{
   int unsigned long count[100];
   tapestats(tag, count);
   printf("Number of independents         = %6d\n", count[0]);
   printf("Number of dependents           = %6d\n", count[1]);
   printf("Max number of live active vars = %6d\n", count[2]);
   printf("Size of value stack            = %6d\n", count[3]);
   printf("Buffer size                    = %6d\n", count[4]);
   cout<<"==============================================="<<endl;
}


void JoeWithModels_AD::showResidue(double *rhsResid,int step)
{
  int nScal = scalarTranspEqVector.size();
  temporalHook_AD(step);
  // residual label at every 10 output steps
  if ((mpi_rank == 0) && ((step%(check_interval*10) == 0) || (step == 1)))
  {
    printf("            psi_rho      psi_rhou-X     psi_rhou-Y    psi_rhou-Z   psi_rhoE      ");
    for (int iScal = 0; iScal < nScal; iScal++)
      printf("%12s", scalarTranspEqVector[iScal].getName());
    cout << endl;
  }

  if (((step%(check_interval) == 0) && (step >= 0)))
  {
	  if (checkParam("WRITE_TURB_ADJOINTS")){
		  write_turb_adjoints_AD();
	  }
  }
  // residual value at each output step
  if (mpi_rank == 0)
  {
    printf("RESID: %6d %12.4e %12.4e %12.4e %12.4e %12.4e", step, rhsResid[0], rhsResid[1], rhsResid[2], rhsResid[3], rhsResid[4]);
    for (int iScal = 0; iScal < nScal; iScal++)
      printf("%12.4e", rhsResid[5+iScal]);
    cout << endl;
  }
}




void JoeWithModels_AD::calcRhs_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, REALAS **rhs_rhoScal, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double (*A)[5][5], double ***AScal, int flagImplicit)
{
  // set RHS to zero
  for (int icv = 0; icv < ncv_gg; icv++)
  {
    rhs_rho[icv] = 0.0;
    for (int i = 0; i < 3; i++)
      rhs_rhou[icv][i] = 0.0;
    rhs_rhoE[icv] = 0.0;
  }

  // set scalars RHS to zero
  
  for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
    for (int icv = 0; icv < ncv_gg; icv++)
      rhs_rhoScal[iScal][icv] = 0.0;
  


  // =======================================================================================
  // NAVIER-STOKES
  // =======================================================================================

  // compute Euler Flux for NS and scalars
  calcEulerFlux_AD(rhs_rho, rhs_rhou, rhs_rhoE, rhs_rhoScal, rho, rhou, rhoE, A, AScal, flagImplicit);
  

  // compute viscous Flux for NS
  if (mu_ref > 0.0)
    calcViscousFluxNS_AD(rhs_rho, rhs_rhou, rhs_rhoE, rho, rhou, rhoE, A, flagImplicit);


  // add source terms to RHS of Navier-Stokes equations
  sourceHook_AD(rhs_rho, rhs_rhou, rhs_rhoE, A);
  sourceHookRansTurb_AD(rhs_rho, rhs_rhou, rhs_rhoE, A);
  sourceHookRansComb_AD(rhs_rho, rhs_rhou, rhs_rhoE, A);

  
  // =======================================================================================
  // SCALARS
  // =======================================================================================
  for (int iScal = 0; iScal < scalarTranspEqVector.size(); iScal++)
  {

    // compute viscous Flux for scalars and add source terms to RHS of scalar equations
    if (AScal == NULL)
    {
      if (mu_ref > 0.0)
        calcViscousFluxScalar_new_AD(rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal], scalarTranspEqVector_AD[iScal], flagImplicit);

      sourceHookScalarRansTurb_new_AD(rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
      sourceHookScalarRansComb_new_AD(rhs_rhoScal[iScal], NULL, scalarTranspEqVector[iScal].getName(), flagImplicit);
    }
    else
    {
      if (mu_ref > 0.0)
        calcViscousFluxScalar_new_AD(rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal], scalarTranspEqVector_AD[iScal], flagImplicit);

      sourceHookScalarRansTurb_new_AD(rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
      sourceHookScalarRansComb_new_AD(rhs_rhoScal[iScal], AScal[iScal][5], scalarTranspEqVector[iScal].getName(), flagImplicit);
    }
  }

}

void JoeWithModels_AD::calcEulerFlux_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, REALAS **rhs_rhoScal, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double (*A)[5][5], double ***AScal, int flagImplicit)
{
  int nScal = scalarTranspEqVector.size();

  double (*Apl)[5] = NULL;
  double (*Ami)[5] = NULL;

  if (flagImplicit)
  {
    Apl = new double[5][5];
    Ami = new double[5][5];

    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
        Apl[i][j] = Ami[i][j] = 0.0;
  }

  // Implicit matrix for scalars: the definition is not easy, so ask an expert (e.g., R. Pecnik)
  double (*AplScal)[6] = NULL;
  double (*AmiScal)[6] = NULL;

  if (flagImplicit)
  {
    if (nScal > 0) AplScal = new double[nScal][6];
    if (nScal > 0) AmiScal = new double[nScal][6];

    for (int iScal = 0; iScal < nScal; iScal++)
      for (int i = 0; i <= 5; i++)
        AplScal[iScal][i] = AmiScal[iScal][i] = 0.0;
  }

  REALQ Frho, Frhou[3], FrhoE;
  REALQS *FrhoScal     = NULL;         if (nScal > 0) FrhoScal     = new REALQS[nScal];
  REALQS *Scalar0      = NULL;         if (nScal > 0) Scalar0      = new REALQS[nScal];            // cell face if second order
  REALQS *Scalar1      = NULL;         if (nScal > 0) Scalar1      = new REALQS[nScal];
  REALQS *ScalCV0      = NULL;         if (nScal > 0) ScalCV0      = new REALQS[nScal];            // cell face if second order
  REALQS *ScalCV1      = NULL;         if (nScal > 0) ScalCV1      = new REALQS[nScal];
  double *ScalConvTerm = NULL;         if (nScal > 0) ScalConvTerm = new double[nScal];            // 0 if convective term no considered, otherwise 1


  // count how many cells switched back to first order due to extrapolated negative density or pressure/temperature at the faces
  int CountReducedOrder = 0;
  int myCountReducedOrder = 0;

  // save the index of kine if defined and save convTerm for speed
  int kine_Index = getScalarTransportIndex("kine");
  for (int iScal = 0; iScal < nScal; iScal++)
    ScalConvTerm[iScal] = double(scalarTranspEqVector[iScal].convTerm);

  // =============================================================================================
  // compute gradients, with boundary values
  // =============================================================================================
  if (sndOrder == true)
  {
    calcCv2Grad_AD(grad_rho, rho, limiterNavierS, rho, epsilonSDWLS);
#ifdef temp_reconstruction
    calcCv2Grad_AD(grad_temp, temp, limiterNavierS, temp, epsilonSDWLS);
#else
    calcCv2Grad_AD(grad_p, press, limiterNavierS, press, epsilonSDWLS);
#endif

    
    for (int iScal = 0; iScal < nScal; iScal++)
    {
      if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
      {
        // For the scalar reconstruction, Grad(rho*Phi) is required
        // Gradients of rho*Phi are saved in grad_phi temporarily
        // Gradients of rho*Phi are also limited like rho with alpha_rho
        // Boundary face values rho_bfa*Phi_fa are saved in Phi_fa temporarily
        REALQS *rhoPhi = new REALQS[ncv_ggff];
        REALQS *phi = scalarTranspEqVector_AD[iScal].phi;
        REALQS (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;

        // Compute rho*Phi
        for (int icv = 0; icv < ncv_ggff; icv++)
          rhoPhi[icv] = rho[icv] * phi[icv];

        // Compute gradients of rho*Phi and limit based on rho*Phi
        calcCv2Grad_AD(grad_phi, rhoPhi, limiterNavierS, rhoPhi, epsilonSDWLS);

        delete [] rhoPhi;
      }
      else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
      {
        REALQS *phi = scalarTranspEqVector_AD[iScal].phi;
        REALQS (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;
        calcCv2Grad_AD(grad_phi, phi, limiterNavierS, phi, epsilonSDWLS);
      }

      else
      {
        cerr << "### JoeWithModels::calcEulerFlux_new => Wrong reconstruction type for scalars! ###" << endl;
        throw(-1);
      }
    }

  }

//===================================================================
//====================================================================
//====================================================================
//	ROUTINES FROM jj_AD.cpp
//====================================================================
//====================================================================
//====================================================================


  // ===============================================================================================
  // cycle through internal faces, assembling flux to both sides
  // ===============================================================================================
  //for (int ifa = nfa_b; ifa < nfa; ifa++)
  for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
  { 
   if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
   {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    assert( icv0 >= 0 && icv0 < ncv_gg);
    assert( icv1 >= 0 && icv1 < ncv_gg);

    int noc00, noc01, noc11, noc10;
    if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
      getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    // face unit normal and area...
    REALX nVec[3] = {0.0, 0.0, 0.0};
    REALX area = normVec3d(nVec, fa_normal[ifa]);

    // .............................................................................................
    // reconstruction of variables at faces: rho, u, T or P, scalars
    // .............................................................................................
    REALQ rho0 = rho[icv0];
    REALQ u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
    REALQ p0 = press[icv0];
    REALQ T0 = temp[icv0];
    REALQ h0 = enthalpy[icv0];
    REALQ gam0 = gamma[icv0];
    REALQ R0 = RoM[icv0];
    REALQS kineCV0 = 0.0;          // cell center
    REALQS kineFA0 = 0.0;          // cell face if second order

    REALQ rho1 = rho[icv1];
    REALQ u1[3] = {vel[icv1][0], vel[icv1][1], vel[icv1][2]};
    REALQ p1 = press[icv1];
    REALQ T1 = temp[icv1];
    REALQ h1 = enthalpy[icv1];
    REALQ gam1 = gamma[icv1];
    REALQ R1 = RoM[icv1];
    REALQS kineCV1 = 0.0;
    REALQS kineFA1 = 0.0;

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      REALQS *phi = scalarTranspEqVector_AD[iScal].phi;
      ScalCV0[iScal] = Scalar0[iScal] = phi[icv0];
      ScalCV1[iScal] = Scalar1[iScal] = phi[icv1];
    }

    if (sndOrder == true)
    {
      REALX r0[3] = {0.0, 0.0, 0.0}, r1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);

      // ----------------------------------------
      // left side
      // ----------------------------------------
      rho0 += vecDotVec3d_AD(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
      T0 += vecDotVec3d_AD(r0, grad_temp[icv0]);
      if ((T0 <= 0.0) || (rho0 <= 0.0))
      {
        T0 = temp[icv0];
        rho0 = rho[icv0];
        myCountReducedOrder++;
      }
#else
      p0 += vecDotVec3d_AD(r0, grad_p[icv0]);
      if ((p0 <= 0.0) || (rho0 <= 0.0))
      {
        p0 = press[icv0];
        rho0 = rho[icv0];
        myCountReducedOrder++;
      }
#endif
      else
      {
        for (int i = 0; i < 3; i++)
          u0[i] += vecDotVec3d_AD(r0, grad_u[icv0][i]);
        for (int iScal = 0; iScal < nScal; iScal++)
        {
          REALQS (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;

          if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
            Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d_AD(r0, grad_phi[icv0])) / rho0;
          else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
            Scalar0[iScal] += vecDotVec3d_AD(r0, grad_phi[icv0]);
        }
      }

      // ----------------------------------------
      // right side
      // ----------------------------------------
      rho1 += vecDotVec3d_AD(r1, grad_rho[icv1]);
#ifdef temp_reconstruction
      T1 += vecDotVec3d_AD(r1, grad_temp[icv1]);
      if ((T1 <= 0.0) || (rho1 <= 0.0))
      {
        T1 = temp[icv1];
        rho1 = rho[icv1];
        myCountReducedOrder++;
      }
#else
      p1 += vecDotVec3d_AD(r1, grad_p[icv1]);
      if ((p1 <= 0.0) || (rho1 <= 0.0))
      {
        p1 = press[icv1];
        rho1 = rho[icv1];
        myCountReducedOrder++;
      }
#endif
      else
      {
        for (int i = 0; i < 3; i++)
          u1[i] += vecDotVec3d_AD(r1, grad_u[icv1][i]);
        for (int iScal = 0; iScal < nScal; iScal++)
        {
          REALQS (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;

          if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
            Scalar1[iScal] = (rho[icv1] * Scalar1[iScal] + vecDotVec3d_AD(r1, grad_phi[icv1])) / rho1;    // CHANGE FOR SCALAR AD

          else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
            Scalar1[iScal] += vecDotVec3d_AD(r1, grad_phi[icv1]);
        }
      }


      // .............................................................................................
      // calculation of other variables at faces: p/T, h, R, gam
      // .............................................................................................
#ifdef temp_reconstruction
      calcThermoProp_T_AD(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
      calcThermoProp_T_AD(p1, h1, R1, gam1, rho1, T1, Scalar1, nScal);
#else
      calcThermoProp_p_AD(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
      calcThermoProp_p_AD(T1, h1, R1, gam1, rho1, p1, Scalar1, nScal);
#endif

    }


    if (kine_Index > -1)   // save kine if defined
    {
      kineCV0 = ScalCV0[kine_Index];         // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA0 = Scalar0[kine_Index];         // cell face left, if second order
      kineCV1 = ScalCV1[kine_Index];         // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA1 = Scalar1[kine_Index];         // cell face right,
    }

    // .............................................................................................
    // calculation of Euler Flux explicit using HLLC
    // .............................................................................................
    calcEulerFlux_HLLC_AD(Frho, Frhou, FrhoE, FrhoScal,
                       rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
                       rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
                       area, nVec, nScal, 0.0);

    // icv0 is always valid...
    rhs_rho[icv0] -= Frho;
    for (int i = 0; i < 3; i++)
      rhs_rhou[icv0][i] -= Frhou[i];
    rhs_rhoE[icv0] -= FrhoE;
    for (int iScal = 0; iScal < nScal; iScal++)
      rhs_rhoScal[iScal][icv0] -= ScalConvTerm[iScal] * FrhoScal[iScal];

    // icv1 can be ghost...  but still valid
      rhs_rho[icv1] += Frho;
      for (int i = 0; i < 3; i++)
        rhs_rhou[icv1][i] += Frhou[i];
      rhs_rhoE[icv1] += FrhoE;
      for (int iScal = 0; iScal < nScal; iScal++)
        rhs_rhoScal[iScal][icv1] += ScalConvTerm[iScal] * FrhoScal[iScal];

    // .............................................................................................
    // calculate implicit matrix using HLLC
    // .............................................................................................
    if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
    {
      calcEulerFluxMatrices_HLLC_AD(Apl, Ami, AplScal, AmiScal,
                rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
                rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, kineCV1,
                area, nVec, nScal, 0.0);

      if (icv0 < ncv)  // if icv0 is internal...
      {
      for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
      {
        A[noc00][i][j] += Apl[j][i];
        A[noc01][i][j] -= Apl[j][i];
      }
      }

      if (icv1 < ncv)  // if icv1 is internal...
      {
        for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
        {
          A[noc11][i][j] -= Ami[j][i];
          A[noc10][i][j] += Ami[j][i];
        }
      }

      for (int iScal = 0; iScal < nScal; iScal++)
        for (int i = 0; i <= 5; i++)
        {
          if (icv0 < ncv)
          {
          AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
          AScal[iScal][i][noc01] -= ScalConvTerm[iScal] * AplScal[iScal][i];
          }

          if (icv1 < ncv)
          {
            AScal[iScal][i][noc11] -= ScalConvTerm[iScal] * AmiScal[iScal][i];
            AScal[iScal][i][noc10] += ScalConvTerm[iScal] * AmiScal[iScal][i];
          }
        }
    }

   }
  }

  // ===============================================================================================
  // cycle through boundary faces, assembling flux
  // ===============================================================================================
  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;

      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION OR WALL BOUNDARY, ATTENTION: ONLY FIRST ORDER!!!
        // .............................................................................................
        if (param->getString() == "SYMMETRY" || (param->getString() == "WALL"))
        {
	   for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];

            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];

            REALX nVec[3] = {0.0, 0.0, 0.0};
            REALX area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
             Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv1];

            REALQS kineFA = 0.0;
            if (kine_Index > -1)
              kineFA = scalarTranspEqVector_AD[kine_Index].phi[icv1];

            calcEulerFlux_HLLC_AD(Frho, Frhou, FrhoE, FrhoScal,
                rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                area, nVec, nScal, 0.0);

            rhs_rho[icv0] -= Frho;
            for (int i = 0; i < 3; i++)
              rhs_rhou[icv0][i] -= Frhou[i];
            rhs_rhoE[icv0] -= FrhoE;
            for (int iScal = 0; iScal < nScal; iScal++)
              rhs_rhoScal[iScal][icv0] -= ScalConvTerm[iScal] * FrhoScal[iScal];

            if (flagImplicit && icv0 < ncv)
            {
              calcEulerFluxMatrices_HLLC_AD(Apl, NULL, AplScal, AmiScal,
                rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                area, nVec, nScal, 0.0);

              int noc00 = nbocv_i[icv0]; // icv0's diagonal
              for (int i = 0; i < 5; i++)
              for (int j = 0; j < 5; j++)
                A[noc00][i][j] += Apl[j][i];

              for (int iScal = 0; iScal < nScal; iScal++)
                for (int i = 0; i <= 5; i++)
                  AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
            }
          }
        }
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), NEUMANN, ...)
        // .............................................................................................
        else
        {
	   for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];

            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];

            REALX nVec[3] = {0.0, 0.0, 0.0};
            REALX area = normVec3d(nVec, fa_normal[ifa]);

            REALQ rho0 = rho[icv0];
            REALQ u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
            REALQ p0 = press[icv0];
            REALQ T0 = temp[icv0];
            REALQ h0 = enthalpy[icv0];
            REALQ gam0 = gamma[icv0];
            REALQ R0 = RoM[icv0];
            REALQS kineCV0 = 0.0;           // cell center
            REALQS kineFA0 = 0.0;           // cell face

            REALQS kineFA1 = 0.0;

            for (int iScal = 0; iScal < nScal; iScal++)
            {
	      ScalCV0[iScal] = Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv0];
              Scalar1[iScal] = scalarTranspEqVector_AD[iScal].phi[icv1];
            }

            if (sndOrder == true)
            {
              REALX r0[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);

              // left side
              rho0 += vecDotVec3d_AD(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
              T0 += vecDotVec3d_AD(r0, grad_temp[icv0]);
              if ((T0 <= 0.0) || (rho0 <= 0.0))
              {
                T0 = temp[icv0];
                rho0 = rho[icv0];
                myCountReducedOrder++;
              }
#else
              p0 += vecDotVec3d_AD(r0, grad_p[icv0]);
              if ((p0 <= 0.0) || (rho0 <= 0.0))
              {
                p0 = press[icv0];
                rho0 = rho[icv0];
                myCountReducedOrder++;
              }
#endif
              else
              {
                for (int i = 0; i < 3; i++)
                  u0[i] += vecDotVec3d_AD(r0, grad_u[icv0][i]);
                for (int iScal = 0; iScal < nScal; iScal++)
                {
                  REALQS (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;
                  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
                    Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d_AD(r0, grad_phi[icv0])) / rho0;    // CHANGE FOR SCALAR AD
                  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
                    Scalar0[iScal] += vecDotVec3d_AD(r0, grad_phi[icv0]);
                }
              }

              // calculation of other variables at faces: p/T, h, R, gam
#ifdef temp_reconstruction
              calcThermoProp_T_AD(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
#else
              calcThermoProp_p_AD(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
#endif
            }

            if (kine_Index > -1)   // save kine if defined
            {
              kineCV0 = ScalCV0[kine_Index];          // cell center
              kineFA0 = Scalar0[kine_Index];                                 // cell face
              kineFA1 = Scalar1[kine_Index];
            }

            calcEulerFlux_HLLC_AD(Frho, Frhou, FrhoE, FrhoScal,
                               rho0,      u0,        p0,          T0,         h0,             R0,        gam0,        Scalar0, kineFA0,
                               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar1, kineFA1,
                               area, nVec, nScal, 0.0);

            rhs_rho[icv0] -= Frho;
            for (int i = 0; i < 3; i++)
              rhs_rhou[icv0][i] -= Frhou[i];
            rhs_rhoE[icv0] -= FrhoE;
            for (int iScal = 0; iScal < nScal; iScal++)
              rhs_rhoScal[iScal][icv0] -= ScalConvTerm[iScal] * FrhoScal[iScal];

            if (flagImplicit && icv0 < ncv)
            {
              calcEulerFluxMatrices_HLLC_AD(Apl, NULL, AplScal, AmiScal,
		        rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, kineCV0,
                        rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar1, kineFA1,
                        area, nVec, nScal, 0.0);

              int noc00 = nbocv_i[icv0]; // icv0's diagonal
              for (int i = 0; i < 5; i++)
              for (int j = 0; j < 5; j++)
                A[noc00][i][j] += Apl[j][i];

              for (int iScal = 0; iScal < nScal; iScal++)
                for (int i = 0; i <= 5; i++)
                  AScal[iScal][i][noc00] += ScalConvTerm[iScal] * AplScal[iScal][i];
            }
          }
        }
      }
    }
  }

  // output the number of times switched back to first order at faces
  MPI_Allreduce(&myCountReducedOrder, &CountReducedOrder, 1, MPI_INT, MPI_SUM, mpi_comm);
  if ((CountReducedOrder > 0) && (mpi_rank == 0))
    cout << "Switched back to first order at " << CountReducedOrder << " face(s)" << endl;

  if (Apl != NULL)  delete [] Apl;
  if (Ami != NULL)  delete [] Ami;

  if (AplScal != NULL) delete [] AplScal;
  if (AmiScal != NULL) delete [] AmiScal;

  if (nScal > 0) delete [] FrhoScal;
  if (nScal > 0) delete [] Scalar0;
  if (nScal > 0) delete [] Scalar1;
  if (nScal > 0) delete [] ScalCV0;
  if (nScal > 0) delete [] ScalCV1;
  if (nScal > 0) delete [] ScalConvTerm;
}

void JoeWithModels_AD::calcViscousFluxNS_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double (*A)[5][5], int flagImplicit)
{
  double (*A0)[5];
  double (*A1)[5];

  if (flagImplicit)
  {
    A0 = new double[5][5];
    A1 = new double[5][5];

    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
        A0[i][j] = A1[i][j] = 0.0;
  }
  else
    A0 = A1 = NULL;

  REALQ Frhou[3] = {0.0, 0.0, 0.0}, FrhoE = 0.0;

  // save the index of kine if defined
  int kine_index = -1;
      kine_index = getScalarTransportIndex("kine");


  // ====================================================================
  //        compute gradients, with boundary values
  // ====================================================================

  calcCv2Grad_AD(grad_enthalpy, enthalpy, limiterNavierS, enthalpy, epsilonSDWLS);

  // ====================================================================
  // cycle through internal faces, assembling flux and matrix
  // ====================================================================
  //for (int ifa = nfa_b; ifa < nfa; ifa++)
  for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
  {
   if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
   {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];

    int noc00, noc01, noc11, noc10;
    if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
      getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    // face unit normal and area...
    REALX nVec[3] = {0.0, 0.0, 0.0};
    REALX area = normVec3d(nVec, fa_normal[ifa]);
    REALX sVec[3] = {0.0, 0.0, 0.0};
    vecMinVec3d(sVec, x_cv[icv1], x_cv[icv0]);
    REALX smag = normVec3d(sVec);

    double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
    vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
    vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
    double w0 = sqrt(vecDotVec3d(dx0, dx0));
    double w1 = sqrt(vecDotVec3d(dx1, dx1));
    double ws = w0 + w1;
    w0 /= ws;
    w1 /= ws;
    REALQ uAux_fa[3] = { w1*vel[icv0][0]+ w0*vel[icv1][0],
                         w1*vel[icv0][1]+ w0*vel[icv1][1],
                         w1*vel[icv0][2]+ w0*vel[icv1][2]};

    // kinetic energy, if defined
    REALQS kine0 = 0.0;
    REALQS kine1 = 0.0;
    REALQS kine_fa = 0.0;
    if (kine_index > -1)
    {
      kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];
      kine1 = scalarTranspEqVector_AD[kine_index].phi[icv1];
      kine_fa = w1*kine0 + w0*kine1;
    }

    // calculate viscous flux
    addViscFlux_AD(Frhou, FrhoE, A0, A1,
              rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], kine0,
              rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], kine1,
              mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa,
              area, nVec, smag, sVec);

    if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
    {
      for (int i=0; i<5; i++)
      for (int j=0; j<5; j++)
      {
        A[noc00][i][j] -= A0[j][i];
        A[noc01][i][j] += A0[j][i];
      }

      if (icv1 < ncv)  // if icv1 is internal...
      {
        for (int i=0; i<5; i++)
        for (int j=0; j<5; j++)
        {
          A[noc11][i][j] += A1[j][i];
          A[noc10][i][j] -= A1[j][i];
        }
      }
    }

    // icv0 is always valid...
    for (int i = 0; i < 3; i++)
      rhs_rhou[icv0][i] -= Frhou[i];
    rhs_rhoE[icv0] -= FrhoE;

    // icv1 can be ghost... But this may be required for the Adjoint
    if (icv1 < ncv_gg)
    {
      for (int i = 0; i < 3; i++)
        rhs_rhou[icv1][i] += Frhou[i];
      rhs_rhoE[icv1] += FrhoE;
    }
   }
  }


  // ====================================================================
  // cycle through boundary faces, assembling flux and matrix
  // ====================================================================
  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;
      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION
        // .............................................................................................
        if (param->getString() == "SYMMETRY" || (param->getString() == "NOTHING"))
        {
	       if (kine_index > -1)
                    {
          		for (int index = 0; index < zone->faVec.size(); ++index) {
          			int ifa = zone->faVec[index];
                                int icv0 = cvofa[ifa][0];
                                int icv1 = cvofa[ifa][1];
                                assert( icv0 >= 0 );
                                double nVec[3] = {0.0, 0.0, 0.0};
                                double area = normVec3d(nVec, fa_normal[ifa]);
      
                                adouble *phi = scalarTranspEqVector_AD[kine_index].phi;
                                adouble kine_fa  = phi[icv1];
          
                                adouble tmp = 1.0/3.0*(rho[icv0] + rho[icv1])*kine_fa;
     
                               for (int i = 0; i < 3; i++)
                                   rhs_rhou[icv0][i] -= tmp*fa_normal[ifa][i];
                               }
                    }
        }
        // .............................................................................................
        // WALL BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "WALL")
        {
          for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];

            REALX nVec[3] = {0.0, 0.0, 0.0};
            REALX area = normVec3d(nVec, fa_normal[ifa]);
            REALX sVec[3] = {0.0, 0.0, 0.0};
            vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
            REALX smag = fabs(vecDotVec3d(sVec, nVec));  // project sVec to wall face normal

            REALQS kine0 = 0.0;
            REALQS kine1 = 0.0;
            REALQS kine_fa = 0.0;
            if (kine_index > -1)
              kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];

            // calculate viscous flux
            addViscFlux_AD(Frhou, FrhoE, A0, NULL, 
                      rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  kine0,
                      rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1],  kine1,
                      mul_fa[ifa], 0.0, lamOcp_fa[ifa], kine_fa, vel[icv1],
                      area, nVec, smag, nVec);  /* <- use nVec here instead of sVec, to avoid inaccurate correction*/

            if (flagImplicit && icv0<ncv)
            {
              int noc00 = nbocv_i[icv0]; // icv0's diagonal

              for (int i=0; i<5; i++)
                for (int j=0; j<5; j++)
                  A[noc00][i][j] -= A0[j][i];
            }

            for (int i = 0; i < 3; i++)
              rhs_rhou[icv0][i] -= Frhou[i];
            rhs_rhoE[icv0] -= FrhoE;
          }
        }
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS
        // .............................................................................................
        else
        {
          for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];

            REALX nVec[3] = {0.0, 0.0, 0.0};
            REALX area = normVec3d(nVec, fa_normal[ifa]);
            REALX sVec[3] = {0.0, 0.0, 0.0};
            vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
            REALX smag = normVec3d(sVec);

            REALQS kine0 = 0.0;
            REALQS kine1 = 0.0;
            REALQS kine_fa = 0.0;
            if (kine_index > -1)
            {
              kine0 = scalarTranspEqVector_AD[kine_index].phi[icv0];
              kine1 = scalarTranspEqVector_AD[kine_index].phi[icv1];
              kine_fa = kine1;
            }

            // calculate viscous flux
            addViscFlux_AD(Frhou, FrhoE, A0, NULL,
                rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0],  kine0,
                rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1],  kine1,
                mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, vel[icv1],
                area, nVec, smag, sVec);

            if (flagImplicit && icv0<ncv)
            {
              int noc00 = nbocv_i[icv0]; // icv0's diagonal

              for (int i=0; i<5; i++)
                for (int j=0; j<5; j++)
                  A[noc00][i][j] -= A0[j][i];
            }

            for (int i = 0; i < 3; i++)
              rhs_rhou[icv0][i] -= Frhou[i];
            rhs_rhoE[icv0] -= FrhoE;
          }
        }
      }
    }
  }

  if (A0  != NULL)  delete [] A0;
  if (A1  != NULL)  delete [] A1;
}

void JoeWithModels_AD::setBC_AD(REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE)
{
  static int first = 1;
  int bc_err = 0;

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;

      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // HOOK BOUNDARY CONDITION
        // .............................................................................................
        if (param->getString() == "HOOK")
        {
          if ((first) && (mpi_rank == 0))
            cout << "Applying HOOK                to zone: "<< zone->getName() << endl;

          boundaryHook_AD(temp, vel, press, &(*zone));

          setScalarBC_AD(&(*zone));
          ComputeBCProperties_T_AD(&(*zone));
        }
        // .............................................................................................
        // CBC BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "CBC")
        {
          double u_bc[3], T_bc, p_bc;

          for (int i=0; i<3; i++)
            u_bc[i] = param->getDouble(i+2);
          T_bc = param->getDouble(5);
          p_bc = param->getDouble(6);

          if ((first)&&(mpi_rank == 0))
            cout << "Applying CBC                 to zone: "<< zone->getName() <<"\t u_bc: "<< u_bc[0]<< " "<< u_bc[1]<< " "
                 << u_bc[2] <<" T_bc: " << T_bc << " p_bc: " << p_bc << endl;

	  for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            //assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

            REALX nVec[3];
            REALX area = normVec3d(nVec, fa_normal[ifa]);

            if (vecDotVec3d(u_bc, nVec) > 0.0)   // outlet
            {
              REALQ velMagN = vecDotVec3d_AD(vel[icv0], nVec);
              REALQ mach = fabs(velMagN)/sos[icv0];

	      temp[icv1] = temp[icv0];
              for (int i=0; i<3; i++)
                vel[icv1][i] = vel[icv0][i];

              if (mach >= 1.0)  press[icv1] = press[icv0];
              else              press[icv1] = p_bc;
            }
            else      // inlet
            {
              temp[icv1] = T_bc;
              for (int i=0; i<3; i++)
                vel[icv1][i] = u_bc[i];
              press[icv1] = p_bc;
            }
          }

          setScalarBC_AD(&(*zone));
          ComputeBCProperties_T_AD(&(*zone));
       }
        // .............................................................................................
        // CBC SUBSONIC INLET BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "CBC_SUBSONIC_INLET")
        {
          double angleU[3], Ttot, htot, ptot;

          for (int i=0; i<3; i++)
            angleU[i] = param->getDouble(i+2);
          Ttot = param->getDouble(5);
          ptot = param->getDouble(6);

          if ((first)&&(mpi_rank == 0))
            cout << "Applying CBC_SUBSONIC        to zone: "<< zone->getName() <<"\t angleU: "<<angleU[0]<<" "<<angleU[1]<<" "<<angleU[2]
                 << " Ttot: " << Ttot << " Ptot: "<< ptot << endl;

	  for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            //assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

            REALQ u[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
            REALQ wPow2 = vecDotVec3d_AD(u, u);            // velocity squared
            REALQ velMag = sqrt(wPow2);                   // apply angle to extrapolated velocity
            for (int i=0; i<3; i++)
              vel[icv1][i] = angleU[i]*velMag;
            temp[icv1] = Ttot;                              // save temporary total temperature (to compute then total enthalpy)
          }

          setScalarBC_AD(&(*zone));
          ComputeBCProperties_T_AD(&(*zone));                  // total enthalpy from total temperature

	  for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];
            int icv1 = cvofa[ifa][1];
            REALQ wPow2 = vecDotVec3d_AD(vel[icv1], vel[icv1]);         // velocity squared
            enthalpy[icv1] -= 0.5* wPow2;                                   // static enthalpy
          }
          
          ComputeBCProperties_H_AD(&(*zone));                  // static temperature and thermo properties from static enthalpy
          
          // Assumes isentropic relations to determine static pressure (= constant cp)
          // At first approximation ok, but could be improved; should for now be considered in defining p_bc
	  for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];
            int icv1 = cvofa[ifa][1];
            press[icv1] = ptot*pow(temp[icv1]/Ttot, gamma[icv1]/(gamma[icv1]-1.0));
          }
        }
        // .............................................................................................
        // CBC SUBSONIC OUTLET BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "CBC_SUBSONIC_OUTLET")
        {
          double p_bc = param->getDouble(2);

          if ((first)&&(mpi_rank == 0))
            cout << "Applying CBC_SUBSONIC_OUTLET to zone: "<< zone->getName() << "\t pOut: "<< p_bc << endl;

	    for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            //assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

            // Assumes that the temperature is constant across outlet boundary (extrapolation) while p_bc is prescribed
             press[icv1] = p_bc;
            for (int i=0; i<3; i++)
              vel[icv1][i] = vel[icv0][i];

            temp[icv1] = temp[icv0];
          }

          setScalarBC_AD(&(*zone));
          ComputeBCProperties_T_AD(&(*zone));
        }
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "SYMMETRY")
        {
          if ((first)&&(mpi_rank == 0))
            cout << "Applying SYMMETRY            to zone: "<< zone->getName() << endl;

	    for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            //assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

            REALX nVec[3];
            REALX area = normVec3d(nVec, fa_normal[ifa]);

            // flip u, APPROXIMATION ---> take velocity at the cell center
            REALQ u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
            REALQ un = vecDotVec3d_AD(nVec, u0);
            for (int i = 0; i < 3; i++)
              vel[icv1][i] = u0[i] - 1.0*un*nVec[i];
            //assert(fabs(vecDotVec3d(vel_bfa[ifa], nVec)) < 1.0e-10);

            temp[icv1]  = temp[icv0];
            press[icv1] = press[icv0];
          }

          setScalarBC_AD(&(*zone));
          ComputeBCProperties_T_AD(&(*zone));
        }
        // .............................................................................................
        // NEUMANN BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "NEUMANN")
        {
          if ((first)&&(mpi_rank == 0))
            cout << "Applying NEUMANN             to zone: "<< zone->getName() << endl;

	    for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
	    //assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

            for (int i = 0; i < 3; i++)
              vel[icv1][i] = vel[icv0][i];

            temp[icv1] = temp[icv0];
            press[icv1] = press[icv0];
          }

          setScalarBC_AD(&(*zone));
          ComputeBCProperties_T_AD(&(*zone));
        }
        // .............................................................................................
        // WALL BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "WALL")
        {
          int i=0;
          double T_bc = 0.0;
          if ((i = param->findString("TEMP")) != 0)
            T_bc = param->getDouble(i+1);

          if ((first)&&(mpi_rank == 0))
          {
            if (T_bc > 0.0)    cout << "Applying WALL isothermal     to zone: "<< zone->getName() << "\t Temp_BC: " << T_bc << endl;
            else               cout << "Applying WALL adiabatic      to zone: "<< zone->getName() << endl;
          }

	    for (int index = 0; index < zone->faVec.size(); ++index) {
            int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
	    int icv1 = cvofa[ifa][1];
            //assert((icv0<ncv)&&(icv1>=ncv_gg)&&(icv1<ncv_ggf));

            for (int i = 0; i < 3; i++)
              vel[icv1][i] = 0.0;

            if (T_bc > 0.0)   temp[icv1] = T_bc;           // wall temperature
            else              temp[icv1] = temp[icv0];     // adiabatic wall

            press[icv1] = press[icv0];                      // Assumes zero pressure gradient at the wall
          }

          setScalarBC_AD(&(*zone));
          ComputeBCProperties_T_AD(&(*zone));
        }
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS
        // .............................................................................................
        else
        {
          if (mpi_rank == 0)
            cerr << "Error: cannot recognize bc for zone: " << zone->getName() << " bc: " << param->getString() << endl;
          bc_err = 1;
        }
      }
      else
      {
        if (mpi_rank == 0)
          cerr << "Error: no bc set for: "<< zone->getName() << endl;
        bc_err = 1;
      }
    }

   // update density at boundary using EOS
  for (int ifa = 0; ifa < nfa_b; ifa++)
  {
    int icv1 = cvofa[ifa][1];
    rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
    sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
  }
  for (int ifa = nfa; ifa < nfa_b2; ifa++)
  {
    int icv1 = cvofa[ifa][1];
    rho[icv1] = press[icv1]/(RoM[icv1]*temp[icv1]);
    sos[icv1] = sqrt(gamma[icv1]*press[icv1]/RoM[icv1]);
  }


  if (bc_err != 0)
    throw(-1);

  first = 0;
}


void JoeWithModels_AD::calcRhs_AD_Coupled(REALA **rhs, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double ***A, int flagImplicit)
{
    int nScal = scalarTranspEqVector.size();
    int nAdjvars = 5+nScal;

  // set RHS to zero
  for (int icv = 0; icv < ncv_gg; icv++)
    for (int i = 0; i < nAdjvars; i++)
         rhs[icv][i] = 0.0;

  // compute Euler Flux for NS and scalars
  calcFluxCoupled_AD(rhs, rho, rhou, rhoE, A, nScal, flagImplicit);

  // add source terms to RHS of Navier-Stokes equations
  sourceHookCoupled_AD(rhs, A, flagImplicit);
  sourceHookRansTurbCoupled_AD(rhs, A, flagImplicit);
  sourceHookRansCombCoupled_AD(rhs, A, flagImplicit);

}


//===================================================================
//====================================================================
//====================================================================
//                 ROUTINES FROM J.cpp
//====================================================================
//====================================================================
//====================================================================


void JoeWithModels_AD::calcFluxCoupled_AD(REALA **rhs, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double ***A, int nScal, int flagImplicit)
{
  double **Apl = NULL;
  double **Ami = NULL;
  double **A0  = NULL;
  double **A1  = NULL;
  
  if (flagImplicit)
  {
    //getMem2D(&Apl, 0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> Apl", true);
    //getMem2D(&Ami, 0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> Ami", true);
    getMem2D(&Apl, 0, 5+nScal-1, 0, 5+nScal-1, (char*)("Apl"));
    getMem2D(&Ami, 0, 5+nScal-1, 0, 5+nScal-1, (char*)("Ami"));
    for (int i = 0; i < 5+nScal; i++)
      for (int j = 0; j < 5+nScal; j++)
      {
        Apl[i][j] = 0.0;
        Ami[i][j] = 0.0;
      }
    if (mu_ref > 0.0)
    {
      //getMem2D(&A0,  0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> A0 ", true);
      //getMem2D(&A1,  0, 5+nScal-1, 0, 5+nScal-1, "JoeWithModels::calcFluxCoupled -> A1 ", true);
      getMem2D(&A0,  0, 5+nScal-1, 0, 5+nScal-1, (char*)("A0"));
      getMem2D(&A1,  0, 5+nScal-1, 0, 5+nScal-1, (char*)("A1"));
      for (int i = 0; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
        {
          A0[i][j] = 0.0;
          A1[i][j] = 0.0;
        }
    }
  }

  adouble *EulerFlux   = new adouble[5+nScal];
  adouble *ViscousFlux = new adouble[5+nScal];
  
  for (int i = 0; i < 5+nScal; i++)
  {
    EulerFlux[i]   = 0.0;
    ViscousFlux[i] = 0.0;
  }
  
  adouble *Scalar0        = NULL;         if (nScal > 0) Scalar0       = new adouble[nScal];            // cell face if 2nd order, cell center if 1st order
  adouble *Scalar1        = NULL;         if (nScal > 0) Scalar1       = new adouble[nScal];            // cell face if 2nd order, cell center if 1st order
  adouble *ScalCV0        = NULL;         if (nScal > 0) ScalCV0       = new adouble[nScal];            // cell center
  adouble *ScalCV1        = NULL;         if (nScal > 0) ScalCV1       = new adouble[nScal];            // cell center
  adouble (*gradScal0)[3] = NULL;         if (nScal > 0) gradScal0     = new adouble[nScal][3];         // gradient of scalars
  adouble (*gradScal1)[3] = NULL;         if (nScal > 0) gradScal1     = new adouble[nScal][3];         // gradient of scalars
  adouble *dpress_dscal0  = NULL;         if (nScal > 0) dpress_dscal0 = new adouble[nScal];            // derivative of pressure with respect to scalars
  adouble *dpress_dscal1  = NULL;         if (nScal > 0) dpress_dscal1 = new adouble[nScal];            // derivative of pressure with respect to scalars
  adouble *diffScal       = NULL;         if (nScal > 0) diffScal      = new adouble[nScal];            // diffusivity of scalars
  double *ConvTerm       = NULL;          if (nScal > 0) ConvTerm      = new double[nScal];             // 0 if convective term not considered, otherwise 1
  double *DiffTerm       = NULL;          if (nScal > 0) DiffTerm      = new double[nScal];             // 0 if diffusive  term not considered, otherwise 1
  
  // count how many cells switched back to first order due to extrapolated negative density or pressure/temperature at the faces
  int CountReducedOrder = 0;
  int myCountReducedOrder = 0;

  // save the index of kine if defined and save convTerm for speed
  int kine_Index = getScalarTransportIndex("kine");
    
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    ConvTerm[iScal] = (double)scalarTranspEqVector[iScal].convTerm;
    DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
    dpress_dscal0[iScal] = 0.0;
    dpress_dscal1[iScal] = 0.0;
  }
  
  // =============================================================================================
  // compute gradients, with boundary values
  // =============================================================================================
  if (sndOrder == true)
  {
    calcCv2Grad_AD(grad_rho, rho, limiterNavierS, rho, epsilonSDWLS);                    // density gradients
#ifdef temp_reconstruction
    calcCv2Grad_AD(grad_temp, temp, limiterNavierS, temp, epsilonSDWLS);                   // temperature gradients
#else
    calcCv2Grad_AD(grad_p, press, limiterNavierS, press, epsilonSDWLS);                    // pressure gradients
#endif
  }
  
  if (mu_ref > 0.0)
    calcCv2Grad_AD(grad_enthalpy, enthalpy, limiterNavierS, enthalpy, epsilonSDWLS);       // enthalpy gradients
  
  for (int iScal = 0; iScal < nScal; iScal++)
  {
    string scalName(scalarTranspEqVector[iScal].getName());

    if ((scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE") && (sndOrder == true))
    {
      // For the scalar reconstruction, Grad(rho * phi) is required
      adouble *phi              = scalarTranspEqVector_AD[iScal].phi;
      adouble *rhoPhi           = scalarTranspEqVector_AD[iScal].rhophi;
      adouble (*grad_rhophi)[3] = scalarTranspEqVector_AD[iScal].grad_rhophi;

      // Compute rho * phi
      for (int icv = 0; icv < ncv_ggff; icv++)
        rhoPhi[icv] = rho[icv] * phi[icv];

      // Compute gradients of rho*Phi and limit based on rho*Phi
      calcCv2Grad_AD(grad_rhophi, rhoPhi, limiterNavierS, rhoPhi, epsilonSDWLS);     // scalars gradients
    }

    if ( ((scalarTranspEqVector[iScal].reconstruction == "STANDARD") && (sndOrder == true)) || ((mu_ref > 0.0) /*&& (scalarTranspEqVector[iScal].diffTerm == 1.0) */) )
    {
      adouble *phi           = scalarTranspEqVector_AD[iScal].phi;
      adouble (*grad_phi)[3] = scalarTranspEqVector_AD[iScal].grad_phi;
      calcCv2Grad_AD(grad_phi, phi, limiterNavierS, phi, epsilonSDWLS);                 // scalars gradients
    }
    
    if ((mu_ref > 0.0) /*&& (scalarTranspEqVector[iScal].diffTerm == 1.0)*/)
    {
      diffusivityHookScalarRansTurb_AD(scalName);
      diffusivityHookScalarRansComb_AD(scalName);
    }
  }
  
  if (flagImplicit)
  {
    pressureDerivativeHookScalarRansTurb_AD();      // compute pressure derivative dP/dScal
    pressureDerivativeHookScalarRansComb_AD();
  }


  // ===============================================================================================
  // cycle through internal faces, assembling flux to both sides
  // ===============================================================================================
  //for (int ifa = nfa_b; ifa < nfa; ifa++)
  for (int ifa = nfa_b; ifa < nfa_b2gg; ifa++) // Ensure we sweep through internal faces of ghosts too
  { 
   if(ifa<nfa || (ifa >= nfa_b2 && ifa < nfa_b2gg))
   {
    int icv0 = cvofa[ifa][0];
    int icv1 = cvofa[ifa][1];
    assert( icv0 >= 0 );
    assert( icv1 >= 0 );

    int noc00, noc01, noc11, noc10;
    if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
      getImplDependencyIndex(noc00, noc01, noc11, noc10, icv0, icv1);

    // face unit normal and area...
    double nVec[3] = {0.0, 0.0, 0.0};
    double area = normVec3d(nVec, fa_normal[ifa]);

    // .............................................................................................
    // reconstruction of variables at faces: rho, u, T or P, scalars
    // .............................................................................................
    adouble rho0 = rho[icv0];
    adouble u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
    adouble p0 = press[icv0];
    adouble T0 = temp[icv0];
    adouble h0 = enthalpy[icv0];
    adouble gam0 = gamma[icv0];
    adouble R0 = RoM[icv0];
    adouble kineCV0 = 0.0;          // cell center
    adouble kineFA0 = 0.0;          // cell face if second order

    adouble rho1 = rho[icv1];
    adouble u1[3] = {vel[icv1][0], vel[icv1][1], vel[icv1][2]};
    adouble p1 = press[icv1];
    adouble T1 = temp[icv1];
    adouble h1 = enthalpy[icv1];
    adouble gam1 = gamma[icv1];
    adouble R1 = RoM[icv1];
    adouble kineCV1 = 0.0;
    adouble kineFA1 = 0.0;
    
    adouble kine_fa = 0.0;          // interpolated kinetic energy at face

    for (int iScal = 0; iScal < nScal; iScal++)
    {
      adouble *phi = scalarTranspEqVector_AD[iScal].phi;
      ScalCV0[iScal] = Scalar0[iScal] = phi[icv0];
      ScalCV1[iScal] = Scalar1[iScal] = phi[icv1];
    }

    if (sndOrder == true)
    {
      double r0[3] = {0.0, 0.0, 0.0};
      double r1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(r1, x_fa[ifa], x_cv[icv1]);

      // ----------------------------------------
      // left side
      // ----------------------------------------
      rho0 += vecDotVec3d_AD(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
      T0 += vecDotVec3d_AD(r0, grad_temp[icv0]);
      if ((T0 <= 0.0) || (rho0 <= 0.0))
      {
        T0 = temp[icv0];
        rho0 = rho[icv0];
        myCountReducedOrder++;
      }
#else
      p0 += vecDotVec3d_AD(r0, grad_p[icv0]);
      if ((p0 <= 0.0) || (rho0 <= 0.0))
      {
        p0 = press[icv0];
        rho0 = rho[icv0];
        myCountReducedOrder++;
      }
#endif
      else
      {
        for (int i = 0; i < 3; i++)
          u0[i] += vecDotVec3d_AD(r0, grad_u[icv0][i]);
        for (int iScal = 0; iScal < nScal; iScal++)
        {
          if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
            Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_rhophi[icv0])) / rho0;
          else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
            Scalar0[iScal] += vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_phi[icv0]);
        }
      }

      // ----------------------------------------
      // right side
      // ----------------------------------------
      rho1 += vecDotVec3d_AD(r1, grad_rho[icv1]);
#ifdef temp_reconstruction
      T1 += vecDotVec3d_AD(r1, grad_temp[icv1]);
      if ((T1 <= 0.0) || (rho1 <= 0.0))
      {
        T1 = temp[icv1];
        rho1 = rho[icv1];
        myCountReducedOrder++;
      }
#else
      p1 += vecDotVec3d_AD(r1, grad_p[icv1]);
      if ((p1 <= 0.0) || (rho1 <= 0.0))
      {
        p1 = press[icv1];
        rho1 = rho[icv1];
        myCountReducedOrder++;
      }
#endif
      else
      {
        for (int i = 0; i < 3; i++)
          u1[i] += vecDotVec3d_AD(r1, grad_u[icv1][i]);
        for (int iScal = 0; iScal < nScal; iScal++)
        {
          if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
            Scalar1[iScal] = (rho[icv1] * Scalar1[iScal] + vecDotVec3d_AD(r1, scalarTranspEqVector_AD[iScal].grad_rhophi[icv1])) / rho1;
          else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
            Scalar1[iScal] += vecDotVec3d_AD(r1, scalarTranspEqVector_AD[iScal].grad_phi[icv1]);
        }
      }
     
      
      // .............................................................................................
      // calculation of other variables at faces: p/T, h, R, gam
      // .............................................................................................
#ifdef temp_reconstruction
      calcThermoProp_T_AD(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
      calcThermoProp_T_AD(p1, h1, R1, gam1, rho1, T1, Scalar1, nScal);
#else
      calcThermoProp_p_AD(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);
      calcThermoProp_p_AD(T1, h1, R1, gam1, rho1, p1, Scalar1, nScal);
#endif

    }


    if (kine_Index > -1)   // save kine if defined
    {
      kineCV0 = ScalCV0[kine_Index];                     // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA0 = Scalar0[kine_Index];                     // cell face left, if second order
      kineCV1 = ScalCV1[kine_Index];                     // cell center left for implicit side (Jacobi is computed with cell center)
      kineFA1 = Scalar1[kine_Index];                     // cell face right,
      
      double dx0[3] = {0.0, 0.0, 0.0};
      double dx1[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));
      kine_fa  = (w1 * kineCV0 + w0 * kineCV1) / (w0 + w1);  // cell face interpolated
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    // .............................................................................................
    // calculate Euler Flux explicit using HLLC
    // .............................................................................................
    calcEulerFluxCoupled_HLLC_AD(EulerFlux,
             rho0, u0, p0, T0, h0, R0, gam0, Scalar0, kineFA0,
             rho1, u1, p1, T1, h1, R1, gam1, Scalar1, kineFA1,
             area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

    // icv0 is always valid...
    for (int i = 0; i < 5+nScal; i++)
      rhs[icv0][i] -= EulerFlux[i];
    
    // icv1 can be ghost... but still valid
      for (int i = 0; i < 5+nScal; i++)
        rhs[icv1][i] += EulerFlux[i];


    // .............................................................................................
    // calculate Euler implicit matrix using HLLC
    // .............................................................................................
    if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
    {
      for (int iScal = 0; iScal < nScal; iScal++)
        if (scalarTranspEqVector_AD[iScal].dpress_dphi != NULL)  
        {
          dpress_dscal0[iScal] = scalarTranspEqVector_AD[iScal].dpress_dphi[icv0];
          dpress_dscal1[iScal] = scalarTranspEqVector_AD[iScal].dpress_dphi[icv1];
        }

      calcEulerFluxMatricesCoupled_HLLC_AD(Apl, Ami,
               rho[icv0], vel[icv0], press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
               rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], ScalCV1, dpress_dscal1, kineCV1,
               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
//      calcEulerFluxMatricesCoupled_HLLC(Apl, Ami,
//               rho0, u0, p0, T0, h0, R0, gam0, Scalar0, dpress_dscal0, kineFA0,
//               rho1, u1, p1, T1, h1, R1, gam1, Scalar1, dpress_dscal1, kineFA1,
//               area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

      for (int i = 0; i < 5+nScal; i++)
      for (int j = 0; j < 5+nScal; j++)
      {
        A[noc00][i][j] += Apl[j][i];
        A[noc01][i][j] -= Apl[j][i];
      }

      if (icv1 < ncv)  // if icv1 is internal...
        for (int i = 0; i < 5+nScal; i++)
        for (int j = 0; j < 5+nScal; j++)
        {
          A[noc11][i][j] -= Ami[j][i];
          A[noc10][i][j] += Ami[j][i];
        }
    }

    // .............................................................................................
    // calculate viscous Flux explicit and implicit matrix
    // .............................................................................................
    if (mu_ref > 0.0)
    {
      double sVec[3] = {0.0, 0.0, 0.0};
      vecMinVec3d(sVec, x_cv[icv1], x_cv[icv0]);
      double smag = normVec3d(sVec);
      double alpha = vecDotVec3d(nVec, sVec);
      assert((alpha > 0.0) && (alpha < 1.0+1.0E-12));             // alpha should now contain s dot n...
      adouble uAux_fa[3] = { 0.5 * (vel[icv0][0] + vel[icv1][0]),
                            0.5 * (vel[icv0][1] + vel[icv1][1]),
                            0.5 * (vel[icv0][2] + vel[icv1][2])};
      
      for (int iScal = 0; iScal < nScal; iScal++)
      {
        for (int i = 0; i < 3; i++)
        {
          gradScal0[iScal][i] = scalarTranspEqVector_AD[iScal].grad_phi[icv0][i];
          gradScal1[iScal][i] = scalarTranspEqVector_AD[iScal].grad_phi[icv1][i];
        }
        diffScal[iScal] = scalarTranspEqVector_AD[iScal].diff[ifa];
      }      
      
      calcViscousFluxCoupled_AD(ViscousFlux, A0, A1,
                 rho[icv0], vel[icv0], grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                 rho[icv1], vel[icv1], grad_u[icv1], enthalpy[icv1], grad_enthalpy[icv1], temp[icv1], RoM[icv1], gamma[icv1], ScalCV1, gradScal1, dpress_dscal1, kineCV1,
                 mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kine_fa, uAux_fa, diffScal, DiffTerm, 
                 area, nVec, smag, sVec, alpha, nScal);
  
      // icv0 is always valid...
      for (int i = 0; i < 5+nScal; i++)
        rhs[icv0][i] -= ViscousFlux[i];
      
      // icv1 can be ghost... but still valid
        for (int i = 0; i < 5+nScal; i++)
          rhs[icv1][i] += ViscousFlux[i]; 
      
      if (flagImplicit && icv0 < ncv && icv1 < ncv_g)
      {
        for (int i = 0; i < 5+nScal; i++) 
        for (int j = 0; j < 5+nScal; j++)
        {
          A[noc00][i][j] += A0[j][i];
          A[noc01][i][j] -= A0[j][i];
        }
  
        if (icv1 < ncv)  // if icv1 is internal...
          for (int i = 0; i < 5+nScal; i++) 
          for (int j = 0; j < 5+nScal; j++)
          {
            A[noc11][i][j] -= A1[j][i];
            A[noc10][i][j] += A1[j][i];
          }
      }
    }
   }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ===============================================================================================
  // cycle through boundary faces, assembling flux
  // ===============================================================================================
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;

      if (getParam(param, zone->getName()))
      {
        // .............................................................................................
        // SYMMETRY BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
        // .............................................................................................
        if (param->getString() == "SYMMETRY")
        {
          for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];          // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
              Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv0];

            adouble kineFA = 0.0;
            if (kine_Index > -1)
              kineFA = scalarTranspEqVector_AD[kine_Index].phi[icv0];

            // Euler flux
            calcEulerFluxCoupled_HLLC_AD(EulerFlux,
                rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
//                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                     area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];

            if (flagImplicit && icv0<ncv)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector_AD[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector_AD[iScal].dpress_dphi[icv0];

              calcEulerFluxMatricesCoupled_HLLC_AD(Apl, NULL,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
//                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
                       area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");          // more "correct" but apparently less stable...

              for (int i = 0; i < 5+nScal; i++)
              for (int j = 0; j < 5+nScal; j++)
                A[noc00][i][j] += Apl[j][i];
            }

            // Viscous flux: only 2/3 viscosity times trace of strain rate tensor and 2/3 * rho * kine
            // Viscous flux
            if (mu_ref > 0.0)
            {
              adouble tmp = 2.0 / 3.0 * ((mul_fa[ifa] + mut_fa[ifa]) * (grad_u[icv0][0][0] + grad_u[icv0][1][1] + grad_u[icv0][2][2]) + rho[icv1] * kineFA);
              rhs[icv0][1] -= area * tmp * nVec[0];
              rhs[icv0][2] -= area * tmp * nVec[1];
              rhs[icv0][3] -= area * tmp * nVec[2];
              
            if (flagImplicit && icv0<ncv)
              {
                // No implicit term considered here!
              }
            }
            
          }
        }
        
        // .............................................................................................
        // NEUMANN BOUNDARY CONDITION, ATTENTION: ONLY FIRST ORDER!!!
        // ASSUMES THAT IF NSE HAS NEUMANN BC, ALL SCALARS HAVE NEUMANN BC!!!
        // .............................................................................................
        else if (param->getString() == "NEUMANN")
        {
          for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];          // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
              Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv1];

            adouble kineFA = 0.0;
            if (kine_Index > -1)
              kineFA = scalarTranspEqVector_AD[kine_Index].phi[icv1];

            // Euler flux
            calcEulerFluxCoupled_HLLC_AD(EulerFlux,
                     rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                     rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];

            if (flagImplicit && icv0<ncv)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector_AD[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector_AD[iScal].dpress_dphi[icv0];

              calcEulerFluxMatricesCoupled_HLLC_AD(Apl, Ami,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

              for (int i = 0; i < 5+nScal; i++)
              for (int j = 0; j < 5+nScal; j++)
                A[noc00][i][j] += Apl[j][i]; // + Ami[i][j];   // !!! Apl + Ami for Neumann boundary condition: Qr = Ql -> dF/dQ = dF/dQr + dF/dQl but more unstable
            }

            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
              double alpha = 1.0;
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector_AD[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = 0.0;
              }
              
              calcViscousFluxCoupled_AD(ViscousFlux, A0, A1,
                         rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kineFA, vel[icv1], diffScal, DiffTerm, 
                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */

              for (int i = 0; i < 5+nScal; i++)
                rhs[icv0][i] -= ViscousFlux[i];
              
            if (flagImplicit && icv0<ncv)
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += A0[j][i]; // + A1[i][j];   // !!!! Apl + Ami for Neumann boundary condition: Qr = Ql -> dF/dQ = dF/dQr + dF/dQl
              }
            }            
          }
        }
        
        // .............................................................................................
        // WALL BOUNDARY CONDITION
        // .............................................................................................
        else if (param->getString() == "WALL")
        {
          // Set DiffTerm flag for scalars to 0 if Neumann BC for viscous flux
          for (int iScal = 0; iScal < nScal; iScal++)
          {
            double dummy;
            string scalName(scalarTranspEqVector[iScal].getName());
            if (!(scalarZoneIsHook(zone->getName(), scalName) || scalarZoneIsDirichlet(dummy, zone->getName(), scalName)))
              DiffTerm[iScal] = 0.0;
          }

          for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0];            // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              ScalCV0[iScal]  = scalarTranspEqVector_AD[iScal].phi[icv0];
              Scalar0[iScal]  = scalarTranspEqVector_AD[iScal].phi[icv1];
            }
            
            adouble kineCV0 = 0.0;
            adouble kineFA = 0.0;
            if (kine_Index > -1)
              kineCV0 = scalarTranspEqVector_AD[kine_Index].phi[icv0];
                        
            // Euler flux
            calcEulerFluxCoupled_HLLC_AD(EulerFlux,
                     rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                     rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, kineFA,
                     area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
//                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];

            if (flagImplicit && icv0<ncv)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector_AD[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector_AD[iScal].dpress_dphi[icv0];

              calcEulerFluxMatricesCoupled_HLLC_AD(Apl, NULL,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
                       rho[icv1], vel[icv1], press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar0, dpress_dscal0, kineFA,
                       area, nVec, nScal, ConvTerm, 0.0, "ONLY_PRESSURE");
//                       area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

              for (int i = 0; i < 5+nScal; i++)
              for (int j = 0; j < 5+nScal; j++)
                A[noc00][i][j] += Apl[j][i];
            }
            
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
              double alpha = 1.0;
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector_AD[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = scalarTranspEqVector_AD[iScal].grad_phi[icv0][i];
              }
              
              // Neumann BC for scalars is enforced by setting DiffTerm to 0
              calcViscousFluxCoupled_AD(ViscousFlux, A0, NULL,
                         rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0],    gamma[icv0], ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                         rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1],    gamma[icv1], Scalar0, gradScal0, dpress_dscal0, kineFA,
                         mul_fa[ifa], 0.0, lamOcp_fa[ifa], kineFA, vel[icv1], diffScal, DiffTerm, 
                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */

              for (int i = 0; i < 5+nScal; i++)
                rhs[icv0][i] -= ViscousFlux[i];
              
            if (flagImplicit && icv0<ncv)
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += A0[j][i];
              }
            }
          }
          
          // Set back DiffTerm flag for scalars to original setting
          for (int iScal = 0; iScal < nScal; iScal++)
            DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
        }
        
        // .............................................................................................
        // OTHER BOUNDARY CONDITIONS (HOOK, DIRICHLET(CBC), ...)
        // .............................................................................................
        else
        {
          // Set DiffTerm flag for scalars to 0 if Neumann BC for viscous flux
          for (int iScal = 0; iScal < nScal; iScal++)
          {
            double dummy;
            string scalname(scalarTranspEqVector[iScal].getName());
            if (!(scalarZoneIsHook(zone->getName(), scalname) || scalarZoneIsDirichlet(dummy, zone->getName(), scalname)))
              DiffTerm[iScal] = 0.0;
          }

          for (int index = 0; index < zone->faVec.size(); ++index) {
          int ifa = zone->faVec[index];
            int icv0 = cvofa[ifa][0];
            int icv1 = cvofa[ifa][1];
            assert( icv0 >= 0 );
            int noc00 = nbocv_i[icv0]; // icv0's diagonal

            double nVec[3] = {0.0, 0.0, 0.0};
            double area = normVec3d(nVec, fa_normal[ifa]);

            adouble rho0 = rho[icv0];
            adouble u0[3] = {vel[icv0][0], vel[icv0][1], vel[icv0][2]};
            adouble p0 = press[icv0];
            adouble T0 = temp[icv0];
            adouble h0 = enthalpy[icv0];
            adouble gam0 = gamma[icv0];
            adouble R0 = RoM[icv0];
            adouble kineCV0 = 0.0;           // cell center
            adouble kineFA0 = 0.0;           // cell face

            adouble kineFA1 = 0.0;

            for (int iScal = 0; iScal < nScal; iScal++)
            {
              ScalCV0[iScal] = Scalar0[iScal] = scalarTranspEqVector_AD[iScal].phi[icv0];
              Scalar1[iScal] = scalarTranspEqVector_AD[iScal].phi[icv1];
            }

            if (sndOrder == true)
            {
              double r0[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(r0, x_fa[ifa], x_cv[icv0]);

              // left side
              rho0 += vecDotVec3d_AD(r0, grad_rho[icv0]);
#ifdef temp_reconstruction
              T0 += vecDotVec3d_AD(r0, grad_temp[icv0]);
              if ((T0 <= 0.0) || (rho0 <= 0.0))
              {
                T0 = temp[icv0];
                rho0 = rho[icv0];
                myCountReducedOrder++;
              }
#else
              p0 += vecDotVec3d_AD(r0, grad_p[icv0]);
              if ((p0 <= 0.0) || (rho0 <= 0.0))
              {
                p0 = press[icv0];
                rho0 = rho[icv0];
                myCountReducedOrder++;
              }
#endif
              else
              {
                for (int i = 0; i < 3; i++)
                  u0[i] += vecDotVec3d_AD(r0, grad_u[icv0][i]);
                for (int iScal = 0; iScal < nScal; iScal++)
                {
                  if (scalarTranspEqVector[iScal].reconstruction == "CONSERVATIVE")
                    Scalar0[iScal] = (rho[icv0] * Scalar0[iScal] + vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_rhophi[icv0])) / rho0;
                  else if (scalarTranspEqVector[iScal].reconstruction == "STANDARD")
                    Scalar0[iScal] += vecDotVec3d_AD(r0, scalarTranspEqVector_AD[iScal].grad_phi[icv0]);
                }
              }

              // calculation of other variables at faces: p/T, h, R, gam
#ifdef temp_reconstruction
              calcThermoProp_T_AD(p0, h0, R0, gam0, rho0, T0, Scalar0, nScal);
#else
              calcThermoProp_p_AD(T0, h0, R0, gam0, rho0, p0, Scalar0, nScal);   // TBD WHAT IS THIS FOR
#endif
            }

            if (kine_Index > -1)   // save kine if defined
            {
              kineCV0 = ScalCV0[kine_Index];          // cell center
              kineFA0 = Scalar0[kine_Index];          // cell face
              kineFA1 = Scalar1[kine_Index];
            }

            calcEulerFluxCoupled_HLLC_AD(EulerFlux,
                     rho0,         u0,           p0,         T0,         h0,         R0,           gam0,         Scalar0, kineFA0,
                     rho[icv1],    vel[icv1], press[icv1],  temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1],  Scalar1, kineFA1,
                     area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

            for (int i = 0; i < 5+nScal; i++)
              rhs[icv0][i] -= EulerFlux[i];
            
            if (flagImplicit && icv0<ncv)
            {
              for (int iScal = 0; iScal < nScal; iScal++)
                if (scalarTranspEqVector_AD[iScal].dpress_dphi != NULL)  
                  dpress_dscal0[iScal] = scalarTranspEqVector_AD[iScal].dpress_dphi[icv0];

              // HACK:: Q_bfa might depend on Q0 so that the Jacobi matrix Apl should be adapted, not yet implemented
              calcEulerFluxMatricesCoupled_HLLC_AD(Apl, NULL,
                        rho[icv0],    vel[icv0],    press[icv0], temp[icv0], enthalpy[icv0], RoM[icv0], gamma[icv0], ScalCV0, dpress_dscal0, kineCV0,
                        rho[icv1],    vel[icv1],    press[icv1], temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1], Scalar1, dpress_dscal0, kineFA1,
                        area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");
//              calcEulerFluxMatricesCoupled_HLLC_AD(Apl, NULL,
//                        rho0,         u0,           p0,          T0,         h0,         R0,           gam0,         Scalar0, dpress_dscal0, kineFA0,
//                        rho[icv1], vel[icv1], press[icv1],  temp[icv1], enthalpy[icv1], RoM[icv1], gamma[icv1],      Scalar1, dpress_dscal0, kineFA1,
//                        area, nVec, nScal, ConvTerm, 0.0, "ALL_TERMS");

              for (int i = 0; i < 5+nScal; i++)
              for (int j = 0; j < 5+nScal; j++)
                A[noc00][i][j] += Apl[j][i];
            }

            
            // Viscous flux
            if (mu_ref > 0.0)
            {
              double sVec[3] = {0.0, 0.0, 0.0};
              vecMinVec3d(sVec, x_fa[ifa], x_cv[icv0]);
//              double smag = fabs(vecDotVec3d(sVec, nVec));                // project sVec to wall face normal
//              double alpha = 1.0;
              double smag = normVec3d(sVec);
              double alpha = vecDotVec3d(nVec, sVec);
              
              for (int iScal = 0; iScal < nScal; iScal++)
              {
                diffScal[iScal] = scalarTranspEqVector_AD[iScal].diff[ifa];
                for (int i = 0; i < 3; i++)
                  gradScal0[iScal][i] = scalarTranspEqVector_AD[iScal].grad_phi[icv0][i];
              }
              
              // Neumann BC for scalars is enforced by setting DiffTerm to 0
              calcViscousFluxCoupled_AD(ViscousFlux, A0, NULL,
                         rho[icv0],    vel[icv0],    grad_u[icv0], enthalpy[icv0], grad_enthalpy[icv0], temp[icv0], RoM[icv0], gamma[icv0], ScalCV0, gradScal0, dpress_dscal0, kineCV0,
                         rho[icv1],    vel[icv1],    grad_u[icv0], enthalpy[icv1], grad_enthalpy[icv0], temp[icv1], RoM[icv1], gamma[icv1], Scalar0, gradScal0, dpress_dscal0, kineFA1,
                         mul_fa[ifa], mut_fa[ifa], lamOcp_fa[ifa], kineFA1, vel[icv1], diffScal, DiffTerm, 
//                         area, nVec, smag, nVec, alpha, nScal);  /* <- use nVec here instead of sVec, to avoid inaccurate correction */
                         area, nVec, smag, sVec, alpha, nScal);  

              for (int i = 0; i < 5+nScal; i++)
                rhs[icv0][i] -= ViscousFlux[i];
              
              if (flagImplicit && icv0<ncv)
              {
                for (int i = 0; i < 5+nScal; i++)
                for (int j = 0; j < 5+nScal; j++)
                  A[noc00][i][j] += A0[j][i];
              }
            }
          }
          
          // Set back DiffTerm flag for scalars to original setting
          for (int iScal = 0; iScal < nScal; iScal++)
            DiffTerm[iScal] = (double)scalarTranspEqVector[iScal].diffTerm;
        }
      }
    }
  }

  // output the number of times switched back to first order at faces
  MPI_Allreduce(&myCountReducedOrder, &CountReducedOrder, 1, MPI_INTEGER, MPI_SUM, mpi_comm);
  if ((CountReducedOrder > 0) && (mpi_rank == 0))
    cout << "Switched back to first order at " << CountReducedOrder << " face(s)" << endl;

  if (Apl != NULL) {freeMem2D(Apl, 0, 5+nScal-1, 0, 5+nScal-1);   Apl = NULL;}
  if (Ami != NULL) {freeMem2D(Ami, 0, 5+nScal-1, 0, 5+nScal-1);   Ami = NULL;}
  if (A0  != NULL) {freeMem2D(A0,  0, 5+nScal-1, 0, 5+nScal-1);   A0  = NULL;}
  if (A1  != NULL) {freeMem2D(A1,  0, 5+nScal-1, 0, 5+nScal-1);   A1  = NULL;}
  
  if (EulerFlux     != NULL) {delete [] EulerFlux;        EulerFlux     = NULL;}
  if (ViscousFlux   != NULL) {delete [] ViscousFlux;      ViscousFlux   = NULL;}
  
  if (Scalar0       != NULL) {delete [] Scalar0;          Scalar0       = NULL;}
  if (Scalar1       != NULL) {delete [] Scalar1;          Scalar1       = NULL;}
  if (ScalCV0       != NULL) {delete [] ScalCV0;          ScalCV0       = NULL;}
  if (ScalCV1       != NULL) {delete [] ScalCV1;          ScalCV1       = NULL;}
  if (gradScal0     != NULL) {delete [] gradScal0;        gradScal0     = NULL;}
  if (gradScal1     != NULL) {delete [] gradScal1;        gradScal1     = NULL;}
  if (dpress_dscal0 != NULL) {delete [] dpress_dscal0;    dpress_dscal0 = NULL;}
  if (dpress_dscal1 != NULL) {delete [] dpress_dscal1;    dpress_dscal1 = NULL;}
  if (diffScal      != NULL) {delete [] diffScal;         diffScal      = NULL;}
  if (ConvTerm      != NULL) {delete [] ConvTerm;         ConvTerm      = NULL;}
  if (DiffTerm      != NULL) {delete [] DiffTerm;         DiffTerm      = NULL;}
}
