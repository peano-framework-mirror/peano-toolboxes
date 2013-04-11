#include "mpibalancing/OracleForOnePhasePerfectBalancing.h"

#include "tarch/Assertions.h"
#include "tarch/la/ScalarOperations.h"
#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "peano/parallel/loadbalancing/Oracle.h"

using namespace peano::parallel::loadbalancing;


tarch::logging::Log mpibalancing::OracleForOnePhasePerfectBalancing::_log( "mpibalancing::OracleForOnePhasePerfectBalancing" );


bool mpibalancing::OracleForOnePhasePerfectBalancing::_forkHasFailed = false;



mpibalancing::OracleForOnePhasePerfectBalancing::OracleForOnePhasePerfectBalancing(bool joinsAllowed, OracleState& oracleState):
  _joinsAllowed(joinsAllowed),
  _oracleState(oracleState)
{
    _commandFromMaster = Continue;
}


mpibalancing::OracleForOnePhasePerfectBalancing::OracleForOnePhasePerfectBalancing(bool joinsAllowed):
  _joinsAllowed(joinsAllowed),
  _oracleState(*(new OracleState))
{
    _commandFromMaster = Continue;
}


mpibalancing::OracleForOnePhasePerfectBalancing::~OracleForOnePhasePerfectBalancing() {
}


void mpibalancing::OracleForOnePhasePerfectBalancing::receivedStartCommand(const LoadBalancingFlag& commandFromMaster ) {
    logInfo("receivedStartCommand(LoadBalancingFlag)", "got command from Master: " << commandFromMaster);
  _commandFromMaster = commandFromMaster;
}


peano::parallel::loadbalancing::LoadBalancingFlag mpibalancing::OracleForOnePhasePerfectBalancing::getCommandForWorker( int workerRank, bool forkIsAllowed, bool joinIsAllowed ) {
  logTraceInWith4Arguments( "getCommandForWorker(int,bool)", workerRank, forkIsAllowed, joinIsAllowed, _joinsAllowed );
  
  unsigned int totalNumberOfWorkers = peano::parallel::loadbalancing::Oracle::getInstance().getNumberOfWorkers();
  std::vector<WorkerInfo>& _workerSet = _oracleState.workerSet;

  LoadBalancingFlag result = Continue;

    // TODO: what about parts which belong to the master? how do we get this information?
   logInfo("getCommandForWorker(int,bool)", "number of Workers:" << totalNumberOfWorkers << " | _workerSet.size(): " << _workerSet.size());

    if (totalNumberOfWorkers < _workerSet.size()) {
        _workerSet.clear();
    } else if (totalNumberOfWorkers == _workerSet.size()) { // wait until all workers have replied
        // 0 = min, 1=avg, 2=max
        double previous_total_all_stats[3];
        double previous_local_all_stats[3];

        double previous_total_below_stats[3];
        double previous_local_below_stats[3];
        int totalNumberOfWorkersTotalBelow = 0;
        int totalNumberOfWorkersLocalBelow = 0;

        double previous_total_above_stats[3];
        double previous_local_above_stats[3];
        int totalNumberOfWorkersTotalAbove = 0;
        int totalNumberOfWorkersLocalAbove = 0;
        int maxLevel = 0;

        previous_total_all_stats[0] = 0.0;
        previous_total_all_stats[1] = 0.0;
        previous_total_all_stats[2] = 0.0;
 
        previous_local_all_stats[0] = 0.0;
        previous_local_all_stats[1] = 0.0;
        previous_local_all_stats[2] = 0.0;
 
        previous_total_below_stats[0] = -1.0;
        previous_total_below_stats[1] = 0.0;
        previous_total_below_stats[2] = -1.0;
 
        previous_local_below_stats[0] = -1.0;
        previous_local_below_stats[1] = 0.0;
        previous_local_below_stats[2] = -1.0;

        previous_total_above_stats[0] = -1.0;
        previous_total_above_stats[1] = 0.0;
        previous_total_above_stats[2] = -1.0;
 
        previous_local_above_stats[0] = -1.0;
        previous_local_above_stats[1] = 0.0;
        previous_local_above_stats[2] = -1.0;

        double master_localWorkload = 0.0;
        //double master_totalWorkload = 0.0;
        double master_localWorkload_additional = 0.0;
 
        std::vector<WorkerInfo>::iterator it = _workerSet.begin();
        if (it != _workerSet.end()) {
            master_localWorkload = it->getLocalWorkloadOfParent();
            //master_totalWorkload = it->getTotalWorkloadOfParent();
        }

        /*// initialize stats with first process 
        it = _workerSet.begin();
        if (it != _workerSet.end()) {
            previous_total_all_stats[0] = it->getTotalWorkload();
            previous_total_all_stats[1] = it->getTotalWorkload();
            previous_total_all_stats[2] = it->getTotalWorkload();
 
            previous_local_all_stats[0] = it->getLocalWorkload();
            previous_local_all_stats[1] = it->getLocalWorkload();
            previous_local_all_stats[2] = it->getLocalWorkload();

            it++;
        }*/
        
        // initialize stats with master
        previous_total_all_stats[0] = 0.0;
        previous_total_all_stats[1] = 0.0;
        previous_total_all_stats[2] = 0.0;

        previous_local_all_stats[0] = master_localWorkload;
        previous_local_all_stats[1] = master_localWorkload;
        previous_local_all_stats[2] = master_localWorkload;

        totalNumberOfWorkers++;
  
        // compute stats over all processes first
        it = _workerSet.begin();
        while (it != _workerSet.end()) {
            previous_total_all_stats[0] = fmin(previous_total_all_stats[0], it->getTotalWorkload());
            previous_total_all_stats[1] += it->getTotalWorkload() ;
            previous_total_all_stats[2] = fmax(previous_total_all_stats[2], it->getTotalWorkload());
 
            previous_local_all_stats[0] = fmin(previous_local_all_stats[0], it->getLocalWorkload());
            previous_local_all_stats[1] += it->getLocalWorkload();
            previous_local_all_stats[2] = fmax(previous_local_all_stats[2], it->getLocalWorkload());

            if (it->isMarkedForRemoval()) {
                master_localWorkload_additional += it->getLocalWorkload();
            }

            if (it->getMaxLevel() > maxLevel) maxLevel = it->getMaxLevel();

            it++;
        }
 
        //const int MIN_STABILITY = 2*maxLevel+1;
        //int min_stability = MIN_STABILITY;
 
        const int MIN_STABILITY = 1;
        int min_stability = MIN_STABILITY;

        previous_total_all_stats[1] = previous_total_all_stats[1] / totalNumberOfWorkers;
        previous_local_all_stats[1] = previous_local_all_stats[1] / totalNumberOfWorkers;
 
        // now compute stats for all processes below or above main average
        if (master_localWorkload >= previous_local_all_stats[1]) {
            previous_local_above_stats[0] = master_localWorkload;
            previous_local_above_stats[1] = master_localWorkload;
            previous_local_above_stats[2] = master_localWorkload;
            totalNumberOfWorkersLocalAbove++;
        } else {
            previous_local_below_stats[0] = master_localWorkload;
            previous_local_below_stats[1] = master_localWorkload;
            previous_local_below_stats[2] = master_localWorkload;
            totalNumberOfWorkersLocalBelow++;
        }
 
        std::vector<WorkerInfo>::iterator it_found = _workerSet.end();
        it = _workerSet.begin();
        while (it != _workerSet.end()) {
            if (it->getTotalWorkload() >= previous_total_all_stats[1]) {
                // ---- total above
                
                if (previous_total_above_stats[0] == -1)  previous_total_above_stats[0] = it->getTotalWorkload();
                else previous_total_above_stats[0] = fmin(previous_total_above_stats[0], it->getTotalWorkload());
 
                previous_total_above_stats[1] +=  it->getTotalWorkload();

                if (previous_total_above_stats[2] == -1)  previous_total_above_stats[2] = it->getTotalWorkload();
                else previous_total_above_stats[2] = fmax(previous_total_above_stats[2], it->getTotalWorkload());
                
                totalNumberOfWorkersTotalAbove++;
            } else {
                // ---- total below
                if (previous_total_below_stats[0] == -1)  previous_total_below_stats[0] = it->getTotalWorkload();
                else previous_total_below_stats[0] = fmin(previous_total_below_stats[0], it->getTotalWorkload());
 
                previous_total_below_stats[1] +=  it->getTotalWorkload();

                if (previous_total_below_stats[2] == -1)  previous_total_below_stats[2] = it->getTotalWorkload();
                else previous_total_below_stats[2] = fmax(previous_total_below_stats[2], it->getTotalWorkload());

                totalNumberOfWorkersTotalBelow++;
            }
 
            if (it->getLocalWorkload() >= previous_local_all_stats[1]) {
                // ---- local above
                
                if (previous_local_above_stats[0] == -1)  previous_local_above_stats[0] = it->getLocalWorkload();
                else previous_local_above_stats[0] = fmin(previous_local_above_stats[0], it->getLocalWorkload());
 
                previous_local_above_stats[1] +=  it->getLocalWorkload();

                if (previous_local_above_stats[2] == -1)  previous_local_above_stats[2] = it->getLocalWorkload();
                else previous_local_above_stats[2] = fmax(previous_local_above_stats[2], it->getLocalWorkload());
                
                totalNumberOfWorkersLocalAbove++;
            } else {
                // ---- local below
                if (previous_local_below_stats[0] == -1)  previous_local_below_stats[0] = it->getLocalWorkload();
                else previous_local_below_stats[0] = fmin(previous_local_below_stats[0], it->getLocalWorkload());
 
                previous_local_below_stats[1] +=  it->getLocalWorkload();

                if (previous_local_below_stats[2] == -1)  previous_local_below_stats[2] = it->getLocalWorkload();
                else previous_local_below_stats[2] = fmax(previous_local_below_stats[2], it->getLocalWorkload());

                totalNumberOfWorkersLocalBelow++;
            }
 
            if (it->getStability() < min_stability) {
                min_stability = it->getStability();
            }
            if (it->getRank() == workerRank) it_found = it;
            
            it++;
        }
  
        previous_total_above_stats[1] = previous_total_above_stats[1] / totalNumberOfWorkersTotalAbove;
        previous_total_below_stats[1] = previous_total_below_stats[1] / totalNumberOfWorkersTotalBelow;
        previous_local_above_stats[1] = previous_local_above_stats[1] / totalNumberOfWorkersLocalAbove;
        previous_local_below_stats[1] = previous_local_below_stats[1] / totalNumberOfWorkersLocalBelow;


    if (it_found != _workerSet.end())  {
        if (it_found->getLocalWorkload() == previous_total_all_stats[2]) {
            result = ForkGreedy;
        }
    } else {
        // INCONSISTENCY
        if (it_found == _workerSet.end()) {
            result = Continue;
            _workerSet.clear();
        }

    }


        // sanity check for loadbalancing decision
        if (it_found != _workerSet.end()) {
            //if (!_oracleState.isGridBalanced) {
            //    result = Continue;
            //    it_found->resetStability();
            //}

            switch (result) {
                case Join:
                    if (joinIsAllowed && _joinsAllowed) {
                        it_found->markForRemoval();
                    } else {
                        result = Continue;
                    }
                    break;
                case ForkGreedy:
                    if (!forkIsAllowed || _forkHasFailed) {
                        logInfo("getCommandForWorker(int,bool)", "i want to fork sooooo hard, but i can't!!!!!");
                        result = Continue;
                    } else {
                        //it_found->markForRemoval();
                    }
                    break;
                default:
                    break;
            }
        }
    } else {
        // worker not known
        result = Continue;
    }
 
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      result = ForkAllChildrenAndBecomeAdministrativeRank;
//        assertion( !_forkHasFailed );
//        if (!forkIsAllowed || _forkHasFailed) {
//                result = Continue;
//        } else {
//            std::vector<WorkerInfo>::iterator it = _workerSet.begin();
//            if (it->getLocalWorkload() == it->getTotalWorkload() /*&& _oracleState.isGridBalanced */) { // replacement for removed hasWorkerAWorker
//                result = ForkGreedy;
//            } else {
//                result = Continue;
//            }
//        }
     }


    // @todo Muss ich selber nochmal ueberlegen
//    if (_forkHasFailed) {
//        // api change
//        _forkHasFailed = false;
//    }


  logTraceOutWith1Argument( "getCommandForWorker(int,bool)", toString(result) );
  return result;
}


void mpibalancing::OracleForOnePhasePerfectBalancing::receivedTerminateCommand(
      int     workerRank,
      double  waitedTime,
      double  workerNumberOfInnerVertices,
      double  workerNumberOfBoundaryVertices,
      double  workerNumberOfOuterVertices,
      double  workerNumberOfInnerCells,
      double  workerNumberOfOuterCells,
      int     workerMaxLevel,
      double  workerLocalWorkload,
      double  workerTotalWorkload,
      int     currentLevel,
      double  parentCellLocalWorkload,
      double  parentCellTotalWorkload,
      const tarch::la::Vector<DIMENSIONS,double>& boundingBoxOffset,
      const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize)
{
    std::vector<WorkerInfo>& _workerSet = _oracleState.workerSet;

    // find worker in current set of worker information
    std::vector<WorkerInfo>::iterator it = _workerSet.begin();
    std::vector<WorkerInfo>::iterator it_found = _workerSet.end();
    bool found = false;
    while (!found && (it != _workerSet.end())) {
        if (workerRank == it->getRank() ) {
            found = true;
            it_found = it;
        } else {
            it++;
        }
    }

    if (it_found != _workerSet.end()) {
        // already inside our set of worker information
        // => update workerCells
        if (it_found->isMarkedForRemoval()) {
            // this is the last result
           _workerSet.erase(it_found);
        } else {
            logInfo("receivedTerminateCommand(existing): ", "load vector: " << workerLocalWorkload << " | " << workerTotalWorkload << " | " << parentCellLocalWorkload << " | " << parentCellTotalWorkload);
            it_found->update(workerLocalWorkload, workerTotalWorkload, parentCellLocalWorkload, parentCellTotalWorkload, workerMaxLevel);
        }
    } else {
        // worker is new
        // => add to set of worker information
        logInfo("receivedTerminateCommand(new): ", "load vector: " << workerLocalWorkload << " | " << workerTotalWorkload << " | " << parentCellLocalWorkload << " | " << parentCellTotalWorkload);
        WorkerInfo temp(workerRank, workerLocalWorkload, workerTotalWorkload, parentCellLocalWorkload, parentCellTotalWorkload, workerMaxLevel);
        _workerSet.push_back(temp);
    }
}

void mpibalancing::OracleForOnePhasePerfectBalancing::plotStatistics() {
}


mpibalancing::OracleForOnePhasePerfectBalancing::OracleForOnePhase* mpibalancing::OracleForOnePhasePerfectBalancing::createNewOracle(int adapterNumber) const {
  return new OracleForOnePhasePerfectBalancing(_joinsAllowed,_oracleState);
}


void mpibalancing::OracleForOnePhasePerfectBalancing::forkFailed() {
  _forkHasFailed = true;
}

int mpibalancing::OracleForOnePhasePerfectBalancing::getCoarsestRegularInnerAndOuterGridLevel() const {
    return 3;
}

/* --------------------- */
mpibalancing::WorkerInfo::WorkerInfo() {
    this->_rank = 0;
    this->_stability = 0;
    this->_localWorkload = 0.0;
    this->_totalWorkload = 0.0;
    this->_localWorkloadOfParent = 0.0;
    this->_totalWorkloadOfParent = 0.0;
    this->_removal = false;
    this->_maxLevel = 0;
}

mpibalancing::WorkerInfo::WorkerInfo(int rank) {
    WorkerInfo();
    std::cout << "removal: " << this->_removal << std::endl;
    this->_rank = rank;
}

mpibalancing::WorkerInfo::WorkerInfo(int rank, double localWorkload, double totalWorkload, double localWorkloadOfParent, double totalWorkloadOfParent, int maxLevel) {
    WorkerInfo();
    this->_rank = rank;
    update(localWorkload, totalWorkload, localWorkloadOfParent, totalWorkloadOfParent, maxLevel);
}

int mpibalancing::WorkerInfo::getRank() {
    return this->_rank;
}

double mpibalancing::WorkerInfo::getLocalWorkload() {
    return this->_localWorkload;
}

double mpibalancing::WorkerInfo::getTotalWorkload() {
    return this->_totalWorkload;
}

double mpibalancing::WorkerInfo::getLocalWorkloadOfParent() {
    return this->_localWorkloadOfParent;
}

double mpibalancing::WorkerInfo::getTotalWorkloadOfParent() {
    return this->_totalWorkloadOfParent;
}


void mpibalancing::WorkerInfo::update(double localWorkload, double totalWorkload, double localWorkloadOfParent, double totalWorkloadOfParent, int maxLevel) {
    if ( (localWorkload != _localWorkload) ||
         (totalWorkload != _totalWorkload) ) {
        //this->resetStability();
    } else {
        //this->_stability++;
    }
    this->_stability = this->_stability + 1;
    this->_localWorkload = localWorkload;
    this->_totalWorkload = totalWorkload;
    this->_localWorkloadOfParent = localWorkloadOfParent;
    this->_totalWorkloadOfParent = totalWorkloadOfParent;
    this->_maxLevel = maxLevel;
}

bool mpibalancing::WorkerInfo::isIdle() {
    return this->_localWorkload < 1.0;
}

int mpibalancing::WorkerInfo::getStability() {
    return this->_stability;
}

void mpibalancing::WorkerInfo::resetStability() {
    this->_stability = 0;
}

bool mpibalancing::WorkerInfo::isMarkedForRemoval() {
    return this->_removal;
}

void mpibalancing::WorkerInfo::markForRemoval() {
    this->_removal = true;
}

int mpibalancing::WorkerInfo::getMaxLevel() {
    return this->_maxLevel;
}


