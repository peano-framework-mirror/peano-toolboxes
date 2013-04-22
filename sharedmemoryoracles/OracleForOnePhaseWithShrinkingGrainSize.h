// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _SHARED_MEMORY_ORACLE_FOR_ONE_PHASE_WITH_SHRINKING_GRAIN_SIZE_H_
#define _SHARED_MEMORY_ORACLE_FOR_ONE_PHASE_WITH_SHRINKING_GRAIN_SIZE_H_


#include "tarch/logging/Log.h"
#include "peano/datatraversal/autotuning/OracleForOnePhase.h"
#include "tarch/timing/Measurement.h"


#include <map>


namespace sharedmemoryoracles {
  class OracleForOnePhaseWithShrinkingGrainSize;
}


/**
 * Oracle With Shrinking Grain Size
 *
 * This is a very simple oracle that runs through three different states:
 *
 * - First, it tries to find out a reasonable time per unknown in serial mode.
 * - Second, it sets the optimal grain size to half of the biggest grain size
 *   found so far, and halves this grain size iteratively as long as the
 *   runtime improves.
 * - If the grain size halfing leads to an increase in runtime, the grain size
 *   is doubled and this is the optimal grain size returned on each request.
 *
 * @author Tobias Weinzierl
 */
class sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize: public peano::datatraversal::autotuning::OracleForOnePhase {
  private:
    static tarch::logging::Log                           _log;

    const bool                                           _pipelineDescendProcessing;
    const bool                                           _pipelineAscendProcessing;
    const peano::datatraversal::autotuning::MethodTrace  _methodTrace;

    bool                                                 _oracleIsSearching;
    int                                                  _biggestProblemSize;
    int                                                  _currentGrainSize;
    tarch::timing::Measurement                           _currentMeasurement;
    double                                               _previousMeasuredTime;
    double                                               _lastProblemSize;
  public:
    /**
     * Oracle Constructor
     */
    OracleForOnePhaseWithShrinkingGrainSize( bool pipelineDescendProcessing = false, bool pipelineAscendProcessing = false, const peano::datatraversal::autotuning::MethodTrace& methodTrace = peano::datatraversal::autotuning::NumberOfDifferentMethodsCalling);

    virtual ~OracleForOnePhaseWithShrinkingGrainSize();

    virtual std::pair<int,bool> parallelise(int problemSize);
    virtual void parallelSectionHasTerminated(double elapsedCalendarTime);
    virtual void plotStatistics() const;

    virtual void informAboutElapsedTimeOfLastTraversal(double elapsedTime);

    /**
     * For this oracle type, the adapter number is completely irrelevant.
     */
    virtual peano::datatraversal::autotuning::OracleForOnePhase* createNewOracle(int adapterNumber, const peano::datatraversal::autotuning::MethodTrace& methodTrace) const;
};


#endif
