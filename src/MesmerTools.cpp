//MesmerTools.cpp
#include "MesmerTools.h"


namespace mesmer
{

  // translation contribution for the partition function of two molecules
  double translationalContribution(const double m1, const double m2, const double beta){
    // Translational contribution
    // 2.0593e19 = conversion factor,  1e-6*(((cm-1 -> j)/(h2*na)))^3/2
    // double tp_C = 2.0593e19 * pow(2. * M_PI ,1.5);

    return (tp_C * pow(m1 * m2 / ((m1 + m2) * beta), 1.5));
  }

  double canonicalPartitionFunction(const vector<double>& DOS, const vector<double>& Ene, const double beta){

    double CanPrtnFn(0.0) ;
	
	// add the function to calculate the partition function and the derivatives and output the thermodynamic data to the .test file
	// Q denotes the partition function, which is the same as CanPrtnFn in the original code
	// Q1 denotes the first derivatives
	// Q2 denotes the second derivatives
	// units: cal, mol, K, i.e. Cp and S for cal/mol/K, H for cal/mol
    for (size_t i(0), j(DOS.size()-1); i < DOS.size(); i++, j--) {
      if (DOS[j] > 0.0)
	  {
        CanPrtnFn += exp( log(DOS[j]) - beta*Ene[j] ) ;
	  }
    }

	return CanPrtnFn;

  }

  double canonicalMeanEnergy(const vector<double>& DOS, const vector<double>& Ene, const double beta){

    double meanEnergy(0.0), CanPrtnFn(0.0) ;
    for (size_t i(0), j(DOS.size()-1); i < DOS.size(); i++, j--) {
      if (DOS[j] > 0.0) {
        double tmp  = exp( log(DOS[j]) - beta*Ene[j] ) ;
        CanPrtnFn  += tmp ;
        meanEnergy += Ene[j]*tmp ;
      }
    }
    return meanEnergy/CanPrtnFn ;

  }

  // add the function to calculate the partition function and the derivatives
  // Q denotes the partition function, which is the same as CanPrtnFn in the original code
  // Q1 denotes the first derivatives
  // Q2 denotes the second derivatives
  // units: cal, mol, K, i.e. Cp and S for cal/mol/K, H for cal/mol
  bool canonicalTestPartitionFunction(const vector<double>& DOS, const vector<double>& Ene, const double beta, double* prtnFn)
  {
	  double Q = 0.0;
	  double Q1 = 0.0;
	  double Q2 = 0.0;

	  for (size_t i(0), j(DOS.size()-1); i < DOS.size(); i++, j--) {
		  if (DOS[j] > 0.0)
		  {
			  Q += DOS[j]*exp(-beta*Ene[j]);
			  Q1 += -Ene[j]/boltzmann_RCpK*DOS[j]*exp(-beta*Ene[j]);
			  Q2 += pow(Ene[j]/boltzmann_RCpK,2)*DOS[j]*exp(-beta*Ene[j]);
		  }
	  }

	  prtnFn[0] = Q;
	  prtnFn[1] = Q1/Q;
	  prtnFn[2] = Q2/Q - pow(Q1/Q, 2);

	  return true;
  }

  // Function to calculate the thermodynamic data and output to the .test file
  // prtnFn[0] is the partition function z1*z2*...*zj*...*zn
  // prtnFn[1] denotes for sum(z'[j]/z[j])
  // prtnFn[2] denotes for sum((z'[j]/z[j])')=sum(z''[j]/z[j]-(z'[j]/z[j])^2)
  // z'[j] is dz/d(1/T)
  void thermodynamicCalc(const double* prtnFn,const double beta, double MW)
  {
	  double temp = 1.0/boltzmann_RCpK/beta;
	  double S, Cp, HmH0;

	  S = idealGasC/Calorie_in_Joule*(2.5+1.5*log(2*M_PI*MW/1000)-4*log(AvogadroC)-3*log(PlancksConstant_in_JouleSecond)-log(atm_in_pascal)+2.5*log(idealGasC*temp)+log(prtnFn[0])-prtnFn[1]/temp);
	  Cp = idealGasC/Calorie_in_Joule*(prtnFn[2]/temp/temp+2.5);
	  HmH0 = idealGasC/Calorie_in_Joule*(-prtnFn[1]+temp*2.5);
	  ctest << "temperature\tQ, S, Cp and H(T)-H(0):\t" << temp << "\t" << prtnFn[0] << "\t" << S << "\t" << Cp << "\t" << HmH0 << endl;
  }


  //
  // Inserts leading zeros to cellDOS and cellEne vector to accounts for the graining integrity.
  //
  void shiftCells(int MaximumCell, int cellOffset, const vector<double>& cellDOS, const vector<double>& cellEne, std::vector<double>& shiftedCellDOS, std::vector<double>& shiftedCellEne){
    for(int i = 0; i < cellOffset; ++i){
      shiftedCellDOS.push_back(0.0);
      shiftedCellEne.push_back(0.0);
    }
    for(int i = cellOffset, j = 0; i < MaximumCell; ++i, ++j){
      shiftedCellDOS.push_back(cellDOS[j]);
      shiftedCellEne.push_back(cellEne[j]);
    }
  }

  //
  // Calculate the average grain energy and then number of states per grain.
  //
  void calcGrainAverages(const int MaximumGrain, const int GrainSize, const vector<double>& shiftedCellDOS, const vector<double>& shiftedCellEne, vector<double>& grainDOS, vector<double>& grainEne)
  {
    grainEne.clear() ;
    grainDOS.clear() ;
    grainEne.resize(MaximumGrain, 0.) ;
    grainDOS.resize(MaximumGrain, 0.) ;

    // Check that there are enough cells.
    if (GrainSize < 1) {
      throw (std::runtime_error("The number of Cells is insufficient to produce requested number of Grains.")); 
    }

    int idx1 = 0 ;
    int idx2 = 0 ;
    for (int i = 0 ; i < MaximumGrain ; ++i ) {

      int idx3(idx1);

      // Calculate the number of states in a grain.
      double gNOS = 0.0 ;
      for (int j = 0 ; j < GrainSize ; ++j, ++idx1 ){
        gNOS += shiftedCellDOS[idx1] ;
      }

      // Calculate average energy of the grain if it contains sum states.
       if ( gNOS > 0.0 ){
        double gSE = 0.0 ; // grain sum of state energy
        for (int j = 0 ; j < GrainSize ; ++j, ++idx3 ){
          gSE += shiftedCellEne[idx3] * shiftedCellDOS[idx3] ;
        }
        grainDOS[idx2] = gNOS ;
        grainEne[idx2] = gSE/gNOS ;
        idx2++ ;
      }

    }

    // Issue warning if number of grains produced is less that requested.

    if ( idx2 != MaximumGrain ) {
      cinfo << "Number of grains produced is not equal to that requested" << endl
        << "Number of grains requested: " << MaximumGrain << endl
        << "Number of grains produced : " << idx2 << endl;
    }
  }

}
