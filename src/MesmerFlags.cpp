#include "MesmerFlags.h"

namespace mesmer{
  MesmerFlags::MesmerFlags() : testDOSEnabled(false),
    testRateConstantEnabled(false),
    microRateEnabled(false),
    grainBoltzmannEnabled(false),
    grainDOSEnabled(false),
    grainTSsosEnabled(false),
    cyclePrintGrainDOS(false),
    cellDOSEnabled(false),
    cyclePrintCellDOS(false),
    reactionOCSEnabled(false),
    kfEGrainsEnabled(false),
    kbEGrainsEnabled(false),
    TunnellingCoeffEnabled(false),
    CrossingCoeffEnabled(false),
    cellFluxEnabled(false),
    grainFluxEnabled(false),
    rateCoefficientsOnly(false),
    useTheSameCellNumber(false),
    grainedProfileEnabled(false),
    speciesProfileEnabled(false),
    viewEvents(false),
    shortestTimeOfInterest(1.0e100),
    maxEvolutionTime(0.),
    printEigenValuesNum(0),
	printEigenVectors(false),
    printReactionOperatorNum(0),
    allowSmallerDEDown(false),
    print_TabbedMatrices(true),
    showCollisionOperator(0),
    useDOSweightedDT(false),
    punchSymbols(),
    punchNumbers(),
    punchFileName(),
	overwriteXmlAnalysis(false),
	autoSetMaxEne(false),
	//default poulation threshold for automatic calculation of energy above the top barrier
	popThreshold(1.0e-5),
    bForceMacroDetailedBalance(false)
  {}
}//namespace

