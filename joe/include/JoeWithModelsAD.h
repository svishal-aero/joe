#pragma once

#include "../../common/include/MiscUtils.h"
#include "../../common/include/tc_vec3d.h"

#include "UgpWithCvCompFlowAD.h"
#include "JoeWithModels.h"
#include "adolc/adolc.h"

class JoeWithModels_AD: virtual public JoeWithModels, virtual public UgpWithCvCompFlow_AD
{
  private:

  int icv_debug;

  void setDebugCv() {

    icv_debug = -1;
     double x_debug[3] = {46.2617,0.5,0.376642};

    for (int icv = 0; icv < ncv_ggff; ++icv) {
      double d2 = 0.0;
      FOR_I3 {
        double dx = x_cv[icv][i] - x_debug[i];
        d2 += dx*dx;
      }
      if (sqrt(d2) < 1.0E-4) {
        cout << "found it on rank: " << mpi_rank << " " << icv<< " "<<x_cv[icv][0] << " " << x_cv[icv][1] << " " << x_cv[icv][2];
        if (icv < ncv)
          cout << " internal cv" << endl;
        else if (icv < ncv_g)
          cout << " ghost level 1 cv" << endl;
        else if (icv < ncv_gg)
          cout << " ghost level 2 cv" << endl;
        else if (icv < ncv_ggf)
          cout << " fake level 1 cv" << endl;
        else
          cout << " fake level 2 cv" << endl;
        icv_debug = icv;
      }
    }
    MPI_Pause((char*)("did we find one?"));

  }

   void dumpDebugState(char * message) {

    int nScal =  scalarTranspEqVector.size();
    if (icv_debug >= 0) {
      cout << "On processor "<<mpi_rank<<" "<<"debug state: " << rho_AD[icv_debug] << " " << rhou_AD[icv_debug][0] << " " << rhou_AD[icv_debug][1] << " " << rhoE_AD[icv_debug] << endl;
	for(int iScal=0; iScal<nScal; iScal++)
            cout<<"Scalar "<<iScal<<" "<<scalarTranspEqVector_AD[iScal].phi[icv_debug]<<endl;
      }
    MPI_Pause(message);
  }


 public:
  /**
   * constructor, pass ParamMap
   */
  JoeWithModels_AD(ParamMap &p) : JoeWithModels(p), UgpWithCvCompFlow_AD(p) { init();}

  /**
   * constructor, pass name of joe's input file
   */
  JoeWithModels_AD(char *name) : JoeWithModels(name), UgpWithCvCompFlow_AD(name) { init();}

  virtual ~JoeWithModels_AD() {init();}


  virtual void init()
  {
    if (mpi_rank == 0)
      cout << "JoeWithModels_AD()"<< endl;
  }


public:

  void runAdjoint(int rest);

  void calcFunctionalDerivative();

  virtual void initialHook_AD() {/* empty */};
  virtual void temporalHook_AD() {/* empty */};
  virtual void temporalHook_AD(int step_val) {/* empty */};

  virtual void calcFunctional_AD(REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE) {/* empty */};

  void calcResidualDerivative(double (*A)[5][5], double ***AScal, int flagImplicit);

  void calcResidual_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3], REALA *rhs_rhoE_AD, REALAS **rhs_rhoScal_AD,
                            REALQ *rho_AD, REALQ (*rhou_AD)[3], REALQ *rhoE_AD, double (*A)[5][5], double ***AScal, int flagImplicit);

  virtual void boundaryHook_AD(REALQ *T_fa, REALQ (*vel_fa)[3], REALQ *p_fa, FaZone *zone) {/* empty */} 

  virtual void sourceHook_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, double (*A)[5][5]) {
	  if(checkParam("AXISYMMETRIC")){
		  if((step%100==0) && mpi_rank == 0)
			  printf("Using axisymmetric\n");
		  axisymmetric_source_AD(rhs_rho, rhs_rhou, rhs_rhoE, A);
	  }
  }
  
  virtual void axisymmetric_source_AD(adouble *rhs_rho, adouble (*rhs_rhou)[3], adouble *rhs_rhoE, double (*A)[5][5])  { 
	  for(int icv = 0; icv < ncv; icv++){
		  adouble _rho, _u, _v, _w, _p, _mu;
		  double _fac;
		  _rho = rho_AD[icv];
		  _u = rhou_AD[icv][0]/_rho;
		  _v = rhou_AD[icv][1]/_rho;
		  _w = rhou_AD[icv][2]/_rho;
		  _p =  (GAMMA-1.0)*(rhoE_AD[icv]- 0.5*(rhou_AD[icv][0]*rhou_AD[icv][0]+rhou_AD[icv][1]*rhou_AD[icv][1]
											 +rhou_AD[icv][2]*rhou_AD[icv][2])/rho_AD[icv]);
		  
		  _fac = 1.0/max(x_cv[icv][1], 1e-10)*cv_volume[icv];
		  
		  // inviscid contribution
		  rhs_rho[icv] -= _rho*_v*_fac;
		  rhs_rhou[icv][0] -= _rho*_u*_v*_fac;
		  rhs_rhou[icv][1] -= _rho*_v*_v*_fac;
		  rhs_rhou[icv][2] -= 0.0;
		  rhs_rhoE[icv] -= (rhoE_AD[icv] + _p)*_v*_fac;


		  int noc00 = nbocv_i[icv];
		  A[noc00][1][1] += _v.value()*_fac;
		  A[noc00][2][2] += 2.0*_v.value()*_fac;
		  A[noc00][4][4] += (GAMMA)*_v.value()*_fac;
		  
		  // viscous contribution
		  adouble muLam   = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);
		  adouble muTurb   = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
		  _mu = muLam + muTurb;
		  
		  adouble _tauxy, _tmp, _tauyy;
		  _tauxy = _mu*(grad_u[icv][0][1] + grad_u[icv][1][0]);
		  _tmp = 2.0/3.0*_mu*(grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2]);
		  _tauyy = 2.0*_mu*grad_u[icv][1][1] - _tmp - _v*_mu/max(x_cv[icv][1], 1e-10)*2.0/3.0;
		  
		  //double lamovercp = InterpolateAtCellCenterFromFaceValues(lamOcp_fa, icv);
		  adouble keff = muLam/0.72 + muTurb/PrTurb;
		  adouble qr = -keff*GAMMA*RoM[icv]/(GAMMA-1)*grad_temp[icv][1];
		  
		  rhs_rho[icv] += 0.0;
		  rhs_rhou[icv][0] += _tauxy*_fac;
		  rhs_rhou[icv][1] += _tauyy*_fac;
		  rhs_rhou[icv][2] += 0.0;
		  rhs_rhoE[icv] += (-qr + _tauxy*_u + _tauyy*_v)*_fac;
	  }
  }

  virtual void sourceHookRansTurb_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, double (*A)[5][5]) {/* empty */}
  virtual void sourceHookRansComb_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, double (*A)[5][5]) {/* empty */}
  //virtual void sourceHookCoupled_AD(REALA **rhs,          double ***A,  int flagImplicit) {/* empty */}
  //virtual void sourceHookRansTurbCoupled_AD(REALA **rhs,          double ***A,  int flagImplicit) {/* empty */}
  //virtual void sourceHookRansCombCoupled_AD(REALA **rhs,          double ***A,  int flagImplicit) {/* empty */}

  void calcEulerFlux_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, REALAS **rhs_rhoScal, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double (*A)[5][5], double ***AScal, int flagImplicit);

  void calcViscousFluxNS_AD(REALA *rhs_rho, REALA (*rhs_rhou)[3], REALA *rhs_rhoE, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double (*A)[5][5], int flagImplicit);

  void setBC_AD(REALQ *rho_AD, REALQ (*rhou_AD)[3], REALQ *rhoE_AD);

  void calcRhs_AD(REALA *rhs_rho_AD, REALA (*rhs_rhou_AD)[3], REALA *rhs_rhoE_AD, REALAS **rhs_rhoScal_AD,
                            REALQ *rho_AD, REALQ (*rhou_AD)[3], REALQ *rhoE_AD, double (*A)[5][5], double ***AScal, int flagImplicit);

  void solveADSystemExplicitEuler();

  void solveADSystemExplicitRK();

  void solveADSystemImplicitEuler();

  void solveADSystem_Sparse();
  void solveADSystem_Sparse1();
	void solveADSystem_IO();


  virtual void calcFunctionalGradient(double *adj_vars) {/* empty */};

  void print_tapestats(int tag);

  void showResidue(double *rhsResid, int step);

  void solveADSystemImplicitEuler_Coupled();

  void calcResidualDerivative_Coupled(double ***A, int flagImplicit);

  void calcResidual_AD_Coupled(REALA **rhs_AD, REALQ *rho_AD, REALQ (*rhou_AD)[3], REALQ *rhoE_AD, double ***A, int flagImplicit);

  void calcRhs_AD_Coupled(REALA **rhs, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double ***A, int flagImplicit);

  void calcFluxCoupled_AD(REALA **rhs, REALQ *rho, REALQ (*rhou)[3], REALQ *rhoE, double ***A, int nScal, int flagImplicit);

  virtual void finalHook_AD() {/*empty */}
};
