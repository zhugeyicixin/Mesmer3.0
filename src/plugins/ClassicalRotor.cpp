#include "../DensityOfStates.h"
#include "../MolecularComponents.h"

using namespace std;
namespace mesmer
{
  class ClassicalRotor : public DensityOfStatesCalculator
  {
  public:

    // Function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, size_t MaximumCell);

    // Function to calculate contribution to canonical partition function.
    virtual double canPrtnFnCntrb(gDensityOfStates* gdos, double beta) ;

	// Function to calculate contribution to canonical partition function and the derivatives.
	virtual bool canTestPrtnFnCntrb(gDensityOfStates* gdos, double beta, double* prtnFn) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) ;

    // Constructor which registers with the list of DensityOfStatesCalculators in TopPlugin
    // This class calculates a complete DOS: it is not an extra class. 
    ClassicalRotor(const char* id) : m_id(id){ Register(); }

    virtual ~ClassicalRotor() {}
    virtual const char* getID()  { return m_id; }

    // Included only in a subset of DensityOfStatesCalculators.
    // Otherwise the baseclass function returns false.
    virtual bool includesRotations(){return true;}

    virtual ClassicalRotor* Clone() { return new ClassicalRotor(*this); }
  private:
    const char* m_id;
  } ;

  //************************************************************
  //Global instance, defining its id (usually the only instance)
  ClassicalRotor theClassicalRotor("ClassicalRotors");
  //************************************************************


  // Provide a function to define particular counts of the DOS of a molecule.
  bool ClassicalRotor::countCellDOS(gDensityOfStates* pDOS, size_t MaximumCell)
  {
    vector<double> cellEne;
    getCellEnergies(MaximumCell, cellEne);
    vector<double> cellDOS(MaximumCell, 0.0) ;

    //
    // Initialize density of states array using calculated rotational
    // density of state from inverse Laplace transform of rotors.
    //
    vector<double> rotConst;
    RotationalTop rotorType = pDOS->get_rotConsts(rotConst);
    double sym = pDOS->get_Sym();
    double qele = pDOS->getSpinMultiplicity();
    double cnt = 0.;

    switch (rotorType){
      case NONLINEAR: //3-D symmetric/asymmetric/spherical top
        cnt = qele * sqrt(4./(rotConst[0] * rotConst[1] * rotConst[2]))/sym ;
        for (size_t i(0) ; i < MaximumCell ; ++i )
          cellDOS[i] = cnt*sqrt(cellEne[i]) ;
        break;
      case LINEAR: //2-D linear
        cnt = qele / (rotConst[0] * sym);
        for (size_t i(0) ; i < MaximumCell ; ++i )
          cellDOS[i] = cnt ;
        break;
      default: // Assume atom.
        cellDOS[0] = qele  ;
        break;
    }

    // Electronic excited states.
    vector<double> eleExc;
    pDOS->getEleExcitation(eleExc);
	vector<double> tmpCellDOS(cellDOS);
    for (size_t j(0) ; j < eleExc.size() ; ++j){
      size_t nr = int(eleExc[j]) ;
      if (nr < MaximumCell) {
        for (size_t i(0) ; i < MaximumCell - nr ; i++ ) {
          tmpCellDOS[i + nr] = tmpCellDOS[i + nr] + cellDOS[i] ;
        }
      }
    }

    pDOS->setCellDensityOfStates(tmpCellDOS) ;

    return true;
  }

  // Calculate contribution to canonical partition function.
  double ClassicalRotor::canPrtnFnCntrb(gDensityOfStates* gdos, double beta) {

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);
    double sym = gdos->get_Sym();

    double qtot(1.0) ; 
    qtot *= double(gdos->getSpinMultiplicity());

	switch(rotorType){
      case NONLINEAR://3-D symmetric/asymmetric/spherical top
        qtot *= (sqrt(M_PI/(rotConst[0] * rotConst[1] * rotConst[2]))*(pow(beta,-1.5))/sym) ;
        break;
      case LINEAR://2-D linear
        qtot /= (rotConst[0]*sym*beta) ;
        break;
      default:
        break; // Assume atom.
    }

    // Electronic excited states.
    vector<double> eleExc;
    gdos->getEleExcitation(eleExc);
    for (size_t j(0) ; j < eleExc.size() ; ++j){
	  qtot += qtot*exp(-beta*eleExc[j]) ;
    }

	return qtot ;
  }
  
  // Function to calculate contribution to canonical partition function and the derivatives.
  // prtnFn[0] is the partition function z1*z2*...*zj*...*zn
  // prtnFn[1] denotes for sum(z'[j]/z[j])
  // prtnFn[2] denotes for sum((z'[j]/z[j])')=sum(z''[j]/z[j]-(z'[j]/z[j])^2)
  // z'[j] is dz/d(1/T)
  bool ClassicalRotor::canTestPrtnFnCntrb(gDensityOfStates* gdos, double beta, double* prtnFn)
  {
	  prtnFn[0] = 1.0;
	  prtnFn[1] = 0.0;
	  prtnFn[2] = 0.0;

	  vector<double> rotConst;
	  RotationalTop rotorType = gdos->get_rotConsts(rotConst);
	  double sym = gdos->get_Sym();

	  double temp = 1.0/boltzmann_RCpK/beta;

	  prtnFn[0] *= double(gdos->getSpinMultiplicity());

	  switch(rotorType){
	  case NONLINEAR://3-D symmetric/asymmetric/spherical top
		  prtnFn[0] *= (sqrt(M_PI/(rotConst[0] * rotConst[1] * rotConst[2]))*(pow(beta,-1.5))/sym) ;
		  prtnFn[1] += -1.5*temp;
		  prtnFn[2] += 1.5*temp*temp;
		  break;
	  case LINEAR://2-D linear
		  prtnFn[0] /= (rotConst[0]*sym*beta) ;
		  prtnFn[1] += -temp;
		  prtnFn[2] += temp*temp;
		  break;
	  default:
		  break; // Assume atom.
	  }

	  // Electronic excited states.
	  // this code could be accelerated because of the recalculated items if needed
	  // the current state is clearer for the physical meaning
	  vector<double> eleExc;
	  gdos->getEleExcitation(eleExc);
	  for (size_t j(0) ; j < eleExc.size() ; ++j){
		  prtnFn[0] += prtnFn[0]*exp(-beta*eleExc[j]) ;
		  prtnFn[1] += -(eleExc[j]/boltzmann_RCpK) * exp(-beta*eleExc[j]) / (1+exp(-beta*eleExc[j]));
		  prtnFn[2] += pow(eleExc[j]/boltzmann_RCpK, 2) * exp(-beta*eleExc[j]) / (1+exp(-beta*eleExc[j]));
	  }

	  return true;
  }

  // Function to return the number of degrees of freedom associated with this count.
  unsigned int ClassicalRotor::NoDegOfFreedom(gDensityOfStates* gdos) {

    vector<double> rotConst;
    RotationalTop rotorType = gdos->get_rotConsts(rotConst);

    unsigned int nDOF(0) ;
    switch(rotorType){
      case NONLINEAR:
        nDOF = 3 ;
        break;
      case LINEAR:
        nDOF = 2 ;
        break;
      default:
        // Assume atom.
        break; 
    }

    return nDOF ;
  }

}//namespace
