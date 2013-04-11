// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _ORACLE_FOR_ONE_PHASE_PERFECT_BALANCING_H_
#define _ORACLE_FOR_ONE_PHASE_PERFECT_BALANCING_H_

#include "peano/parallel/loadbalancing/OracleForOnePhase.h"
#include "tarch/logging/Log.h"


#include <set>
#include <vector>
#include <algorithm>


namespace mpibalancing {
  class OracleForOnePhasePerfectBalancing;
  class WorkerInfo;
     struct _OracleState;
      typedef struct _OracleState OracleState;
}


class mpibalancing::WorkerInfo {
    public:
        WorkerInfo();
        WorkerInfo(int rank);
        WorkerInfo(int rank, double localWorkload, double totalWorkload, 
                             double localWorkloadOfParent, double totalWorkloadOfParent,
                             int maxLevel);
        int getRank();
        double getLocalWorkload();
        double getTotalWorkload();
        double getLocalWorkloadOfParent();
        double getTotalWorkloadOfParent();
        void update(double localWorkload, double totalWorkload, double localWorkloadOfParent, double totalWorkloadOfParent, int maxLevel);
 
        bool isIdle();
        int getStability();
        void resetStability();

        bool isMarkedForRemoval();
        void markForRemoval();

        int getMaxLevel();
    private:
        int _rank;
        int _stability;
        int _maxLevel;
        double _localWorkload;
        double _totalWorkload;
        double _localWorkloadOfParent;
        double _totalWorkloadOfParent;
        bool _removal;
};

struct mpibalancing::_OracleState {
    std::vector<mpibalancing::WorkerInfo> workerSet;
};



/**
 * @author Roland Wittmann
 */
class mpibalancing::OracleForOnePhasePerfectBalancing: public peano::parallel::loadbalancing::OracleForOnePhase {
  public:
    /**
     * Constructor
     *
     * !!! OracleState
     *
     * The oracle state is shared by all instances of this oracle no matter
     * which adapter uses them. Please c
     *
     * @param joinsAllowed  Set true if you wanna enable joins
     * @param oracleState   Global oracle state shared by all oracles of this type.
     */
    OracleForOnePhasePerfectBalancing(bool joinsAllowed, OracleState& oracleState);

    /**
     * Constructor
     *
     * This is a brief fix that basically creates one state object on the heap and
     * passes it to this oracle instance (and implicitly to all others). If you don't
     * want to change oracles throughout your computation, this constructor is fine.
     * Please note however that it is not 100% clean code, as it induces a memory
     * leak.
     */
    OracleForOnePhasePerfectBalancing(bool joinsAllowed);


    virtual ~OracleForOnePhasePerfectBalancing();

    virtual void receivedStartCommand(const peano::parallel::loadbalancing::LoadBalancingFlag& commandFromMaster );

    /**
     * This operation is not const, as it might update some local stuff.
     */
    virtual peano::parallel::loadbalancing::LoadBalancingFlag getCommandForWorker( int workerRank, bool forkIsAllowed, bool joinIsAllowed );

    /**
     * Information about termination call
     *
     * Is called whenever the master receives the acknowledgement from the
     * worker that the latter has finished its local traversal. Provides
     * statistics that you might want to bookkeep for load balancing. All t
     * the workload doubles stem from the cells' workload. See
     * Cell::setNodeWorkload() to adopt these figures to your own workload
     * model.
     *
     * @see receivedStartCommand()
     * @see getCommandForWorker()
     * @see Oracle
     *
     * @param workerRank   Rank of the worker that has just reported that it
     *                     finished its traversal.
     * @param waitedTime   Time (user time) that the master had to wait until
     *                     the worker delivered its finish message. This time
     *                     is zero, if the message was already in the MPI queue
     *                     when the master checked for the worker.
     * @param workerNumberOfInnerVertices  Number of inner vertices handled by
     *                     this worker. If you require the total number, you
     *                     have to feed your oracle manually within
     *                     endIteration(). Here, you have access to the state
     *                     object holding the total numbers.
     * @param workerNumberOfBoundaryVertices Number of boundary vertices
     *                     handled by this worker.
     * @param workerNumberOfOuterVertices Number of outer vertices handled by
     *                     this worker.
     * @param workerNumberOfInnerCells Number of inner cells handled by
     *                     this worker.
     * @param workerNumberOfOuterCells Number of outer cells handled by this
     *                     worker.
     * @param workerMaxLevel Maximum level handled by the worker. If you
     *                     compare it to current level, you have information
     *                     about the height of the worker tree.
     * @param workerLocalWorkload The workload handled by the worker. This is
     *                     the worker's local workload, i.e. if the worker has
     *                     forked itself again, workload of these children is
     *                     not contained within this figure.
     * @param workerTotalWorkload The workload represented by the worker. This
     *                     number is equal to workerLocalWorkload if the worker
     *                     has not forked again, i.e. if it does not have any
     *                     children. Otherwise, total workload comprises both
     *                     the worker's workload and those of its children.
     * @param currentLevel Current level, i.e. level of the root cell of the
     *                     worker partition.
     * @param parentCellLocalWorkload Local workload of the next coarser cell
     *                     on this rank. This number does not comprise workload
     *                     of any worker, i.e. it does in particular neither
     *                     comprise workerLocalWorkload nor
     *                     workerTotalWorkload.
     * @param parentCellTotalWorkload Total workload of the next coarser cell
     *                     comprising the ranks workload of the whole tree
     *                     induced by this cell plus the workloads of all
     *                     forked subtrees.
     * @param boundingBoxOffset Bounding box of this worker subtree.
     * @param boundingBoxSize   Bounding box of this worker subtree.
     */
    virtual void receivedTerminateCommand(
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
      const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize
    );

    /**
     * Plot something to info log device.
     */
    virtual void plotStatistics();

    /**
     * Clone this oracle. This operation is used by the singleton whenever a
     * piece of code asks for parallelisation that never asked before.
     *
     * @param adapterNumber Number of your adapter. Have a closer look to your
     *        repository's state if you want to find out which adapters are
     *        mapped to which state. You can even use the toString() operation
     *        there to map this parameter to a string.
     */
    virtual OracleForOnePhase* createNewOracle(int adapterNumber) const;

    virtual void forkFailed();
 
    virtual int getCoarsestRegularInnerAndOuterGridLevel() const;

  private:
    /**
     * Logging device
     */
    static tarch::logging::Log  _log;

    /**
     * If a fork failed, all the oracles should stop to ask for further forks.
     * Wouldn't make sense and just slow down the application.
     */
    static bool                 _forkHasFailed;

    /**
     * Global flag set at construction time.
     */
    bool                        _joinsAllowed;

    OracleState& _oracleState;
    
    peano::parallel::loadbalancing::LoadBalancingFlag _commandFromMaster;
};



#endif // _ORACLE_FOR_ONE_PHASE_BINARY_PARTITIONING_H_
