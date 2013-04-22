#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"
#include "peano/utils/Globals.h"
#include "tarch/Assertions.h"

#include <cstdlib>
#include <limits>


tarch::logging::Log  sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::_log( "sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize" );


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::OracleForOnePhaseWithShrinkingGrainSize(
  bool                                                  pipelineDescendProcessing,
  bool                                                  pipelineAscendProcessing,
  const peano::datatraversal::autotuning::MethodTrace&  methodTrace
):
  _pipelineDescendProcessing(pipelineDescendProcessing),
  _pipelineAscendProcessing(pipelineAscendProcessing),
  _methodTrace(methodTrace),
  _oracleIsSearching(true),
  _biggestProblemSize(0),
  _currentGrainSize(std::numeric_limits<int>::max()),
  _currentMeasurement(1.0e-2), // magic constant
  _previousMeasuredTime(-1.0),
  _lastProblemSize(0.0) {
}


std::pair<int,bool> sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelise(int problemSize) {
  if (_biggestProblemSize<problemSize) _biggestProblemSize = problemSize;

  _lastProblemSize = problemSize;

  if (problemSize <= _currentGrainSize) {
    return std::pair<int,bool>(0,_oracleIsSearching);
  }
  else {
    return std::pair<int,bool>(_currentGrainSize,_oracleIsSearching);
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::parallelSectionHasTerminated(double elapsedCalendarTime) {
  assertion(_oracleIsSearching);
  assertion(_lastProblemSize!=0.0);

  _currentMeasurement.setValue(elapsedCalendarTime/_lastProblemSize);

  if (_currentMeasurement.isAccurateValue()) {
    _currentMeasurement.increaseAccuracy(2.0);

    // first phase has finished
    if (_biggestProblemSize < _currentGrainSize) {
      assertion(_previousMeasuredTime==-1.0);
      _currentGrainSize = _biggestProblemSize / 2 + 1;
    }
    else if (
      _previousMeasuredTime > _currentMeasurement.getValue() &&
      _currentGrainSize >= 2
    ) {
      _currentGrainSize /= 2;
    }
    else if (
      _previousMeasuredTime > _currentMeasurement.getValue()
    ) {
      _oracleIsSearching  = false;
    }
    else {
      _oracleIsSearching  = false;
      _currentGrainSize  *= 2;
    }

   _previousMeasuredTime = _currentMeasurement.getValue();
    _currentMeasurement.erase();
  }
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::plotStatistics() const {
  if (_biggestProblemSize < _currentGrainSize) {
    logInfo(
      "plotStatistics()",
      "(method-trace=" << toString(_methodTrace) <<
      ",still determining serial runtime" <<
      ")"
    );
  }
  else if (_oracleIsSearching) {
    logInfo(
      "plotStatistics()",
      "(method-trace=" << toString(_methodTrace) <<
      ",grain-size=" << _currentGrainSize <<
      ",biggest-problem-size=" << _biggestProblemSize <<
      ",t[prev]=" << _previousMeasuredTime <<
      ",t[current]=" << _currentMeasurement.getValue() <<
      ")"
    );
  }
  else if (_biggestProblemSize==_currentGrainSize) {
    logInfo(
      "plotStatistics()",
      "(method-trace=" << toString(_methodTrace) <<
      ",does not scale, oracle is not searching anymore" <<
      ")"
    );
  }
  else {
    logInfo(
      "plotStatistics()",
      "(method-trace=" << toString(_methodTrace) <<
      ",grain-size=" << _currentGrainSize <<
      ",biggest-problem-size=" << _biggestProblemSize <<
      ",t=" << _previousMeasuredTime <<
      ",oracle is not searching anymore" <<
      ")"
    );
  }
}


sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::~OracleForOnePhaseWithShrinkingGrainSize() {
}


peano::datatraversal::autotuning::OracleForOnePhase* sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::createNewOracle(int adapterNumber, const peano::datatraversal::autotuning::MethodTrace& methodTrace) const {
  return new OracleForOnePhaseWithShrinkingGrainSize(_pipelineDescendProcessing, _pipelineAscendProcessing,methodTrace);
}


void sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize::informAboutElapsedTimeOfLastTraversal(double elapsedTime) {
}
