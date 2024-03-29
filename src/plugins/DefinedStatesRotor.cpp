
//-------------------------------------------------------------------------------------------
//
// DefinedStatesRotor.cpp
//
// This file contains the implementation of the methods for importing and testing the 
// density of states calculated or determine by other means. An example where this class
// might be deployed are the rotary-electronic states of OH which are not well approximated
// by treating these modes as decoupled, the coupling being described by one of Hund's
// coupling cases.
//
//-------------------------------------------------------------------------------------------

#include "../DensityOfStates.h"
#include "../MolecularComponents.h"
#include "../Molecule.h"
#include "../Constants.h"

using namespace std;
using namespace Constants;

namespace mesmer
{
  class DefinedStatesRotor : public DensityOfStatesCalculator
  {
  public:
    //Read data from XML. 
    virtual bool ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC=NULL);

    // Function to define particular counts of the DOS of a molecule.
    virtual bool countCellDOS(gDensityOfStates* mol, size_t MaximumCell);

    // Function to calculate contribution to canonical partition function.
    virtual double canPrtnFnCntrb(gDensityOfStates* gdos, double beta) ;

	// Function to calculate contribution to canonical partition function and the derivatives.
	virtual bool canTestPrtnFnCntrb(gDensityOfStates* gdos, double beta, double* prtnFn) ;

    // Function to return the number of degrees of freedom associated with this count.
    virtual unsigned int NoDegOfFreedom(gDensityOfStates* gdos) ;

    // Constructor which registers with the list of DensityOfStatesCalculators in the base class
    DefinedStatesRotor(const char* id) : m_id(id), 
      m_energyLevels(), 
      m_degeneracies() 
    { Register(); }

    virtual ~DefinedStatesRotor() {}
    virtual const char* getID()  { return m_id; }
    virtual bool includesRotations(){return true;}
    virtual DefinedStatesRotor* Clone() { return new DefinedStatesRotor(*this); }

  private:
    const char* m_id;
    vector<double> m_energyLevels ;   // Energies of the states and their associated degeneracies.
    vector<int>    m_degeneracies ;

  } ;

  //-------------------------------------------------------------
  //Global instance, defining its id
  DefinedStatesRotor theDefinedStatesRotor("DefinedStatesRotors");
  //-------------------------------------------------------------

  //Read data from XML and store in this instance.
  bool DefinedStatesRotor::ReadParameters(gDensityOfStates* gdos, PersistPtr ppDOSC)
  {
    // Read in defined states data.

    PersistPtr pp = ppDOSC->XmlMoveTo("me:States") ;

    if (pp) {

      const char* p = pp->XmlReadValue("units", optional);
      string units = p ? p : "cm-1";

      while(pp = pp->XmlMoveTo("me:State"))
      {
        double energy = pp->XmlReadDouble("energy", optional);
        energy = getConvertedEnergy(units, energy);
        m_energyLevels.push_back(energy) ;
        int degeneracy = pp->XmlReadInteger("degeneracy", optional);
        m_degeneracies.push_back(degeneracy) ;
      }

    } else {

      // Default : free rotor.

      cinfo << "No rotory-electronic states defined for " << gdos->getHost()->getName() << "." << endl;

      m_energyLevels.push_back(0.0) ;
      m_degeneracies.push_back(1) ;
    }

    return true;
  }

  //
  // Calculate quantum mechanical 1D rotor densities of states of a free 
  // rotor and convolve them with the main density of states.
  //
  bool DefinedStatesRotor::countCellDOS(gDensityOfStates* pDOS, size_t MaximumCell)
  {

    vector<double> cellDOS(MaximumCell, 0.0) ;

    // Shift eigenvalues by the zero point energy and form initial denisity of states array.

    double zeroPointEnergy(m_energyLevels[0]) ;
    for (size_t i(0) ; i < m_energyLevels.size() ; i++ ) {
      size_t nr = int(m_energyLevels[i] - zeroPointEnergy) ;
      if (nr < MaximumCell) {
        cellDOS[nr] = double(m_degeneracies[i]) ;
      }
    }

    // Initialise denisty of states.   

    pDOS->setCellDensityOfStates(cellDOS) ;

    return true;

  }

  //
  // Provide a function to calculate contribution to canonical partition function.
  //
  double DefinedStatesRotor::canPrtnFnCntrb(gDensityOfStates* gdos, double beta)
  {
    double Qrot(0.0) ;

    double zeroPointEnergy(m_energyLevels[0]) ; 
    for (size_t i(0) ; i < m_energyLevels.size() ; i++ ) {
      Qrot += double(m_degeneracies[i])*exp(-beta*(m_energyLevels[i] - zeroPointEnergy)) ;
    }

    return Qrot ;
  }

  // Function to calculate contribution to canonical partition function and the derivatives.
  // prtnFn[0] is the partition function z1*z2*...*zj*...*zn
  // prtnFn[1] denotes for sum(z'[j]/z[j])
  // prtnFn[2] denotes for sum((z'[j]/z[j])')=sum(z''[j]/z[j]-(z'[j]/z[j])^2)
  // z'[j] is dz/d(1/T)
  bool DefinedStatesRotor::canTestPrtnFnCntrb(gDensityOfStates* gdos, double beta, double* prtnFn)
  {
	  double Q = 0.0;
	  double Q1 = 0.0;
	  double Q2 = 0.0;

	  double zeroPointEnergy(m_energyLevels[0]) ; 
	  for (size_t i(0) ; i < m_energyLevels.size() ; i++ ) {
		  Q += double(m_degeneracies[i])*exp(-beta*(m_energyLevels[i] - zeroPointEnergy)) ;
		  Q1 += -(m_energyLevels[i]-zeroPointEnergy)/boltzmann_RCpK * double(m_degeneracies[i])*exp(-beta*(m_energyLevels[i] - zeroPointEnergy));
		  Q2 += pow((m_energyLevels[i]-zeroPointEnergy)/boltzmann_RCpK, 2) * double(m_degeneracies[i])*exp(-beta*(m_energyLevels[i] - zeroPointEnergy));
	  }

	  prtnFn[0] = Q;
	  prtnFn[1] = Q1/Q;
	  prtnFn[2] = Q2/Q - pow(Q1/Q, 2);

	  return true;
  }

  // Function to return the number of degrees of freedom associated with this count.
  unsigned int DefinedStatesRotor::NoDegOfFreedom(gDensityOfStates* gdos) {

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
