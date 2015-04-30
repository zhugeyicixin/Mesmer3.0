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
	double temp = 0.0;
	double Q, Q1, Q2, S, Cp, HmH0;
	Q = 0.0;
	Q1 = 0.0;
	Q2 = 0.0;
	temp = 1.0/boltzmann_RCpK/beta;
	
    for (size_t i(0), j(DOS.size()-1); i < DOS.size(); i++, j--) {
      if (DOS[j] > 0.0)
	  {
        CanPrtnFn += exp( log(DOS[j]) - beta*Ene[j] ) ;
		Q += DOS[j]*exp(-beta*Ene[j]);
		Q1 += -Ene[j]/boltzmann_RCpK*DOS[j]*exp(-beta*Ene[j]);
		Q2 += pow(Ene[j]/boltzmann_RCpK,2)*DOS[j]*exp(-beta*Ene[j]);
	  }
    }
	S = 8.314/4.184*log(Q)-8.314/4.184*Q1/Q/temp+1.5*8.314/4.184+8.314/4.184*(11.16970828-log(1.013e5)+2.5*log(8.314*temp));
	Cp = 8.314/4.184/temp/temp*(Q2/Q-(Q1*Q1/Q/Q))+8.314/4.184*2.5;
	HmH0 = -8.314/4.184*Q1/Q+8.314/4.184*temp*2.5;
	ctest << "temperature\tQ, S, Cp and HmH0:\t" << temp << "\t" << Q << "\t" << S << "\t" << Cp << "\t" << HmH0 << endl;
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
