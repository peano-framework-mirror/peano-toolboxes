// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _MPIBALANCING_FAIR_NODE_POOL_STRATEGY_H_
#define _MPIBALANCING_FAIR_NODE_POOL_STRATEGY_H_


#ifdef Parallel
#include <mpi.h>
#endif

#include "tarch/parallel/NodePoolStrategy.h"
#include "tarch/logging/Log.h"

#include <list>


namespace mpibalancing {
  class FairNodePoolStrategy;
}


/**
 * Fair Node Pool Strategy
 *
 * Also a very simple node pool strategy. Different to the default strategy
 * answering node requests on a fcfs basis, this class collects multiple
 * requests and does some bookkeeping how many classes have already got
 * additional workers. It than answers those nodes first that got the smallest
 * number of additional workers so far.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.6 $
 */
class mpibalancing::FairNodePoolStrategy: public tarch::parallel::NodePoolStrategy {
  protected:
    /**
     * Copy from FCFS but enriched by counter how many rank have already
     * requested an update.
     *
     * @author Tobias Weinzierl
     * @version $Revision: 1.4 $
     */
    class NodePoolListEntry {
      public:
        /**
         * Represents the state of the worker, i.e. whether it is idle or busy.
         */
        enum State {
          IDLE,
          WORKING
        };

      private:
        /**
         * Holds the rank of the process represented by this object.
         */
        int         _rank;

        double      _bookedWorkers;

        /**
         * Holds the state of the process.
         */
        State       _state;

        /**
         * Machine name
         */
        std::string _name;

      public:
        /**
         * Construct one entry. By default this entry corresponds to an idle worker.
         */
        NodePoolListEntry( int rank, const std::string& name );

        virtual ~NodePoolListEntry();

        /**
         * Activates the node. Precondition: Node is idle. Thus, the local min level
         * is overwritten by the argument level and the state is set to working.
         */
        void activate();

        /**
         * The local rank is set to 0 and the state is switched to idle.
         */
        void deActivate();

        /**
         * @return Rank of process.
         */
        int getRank() const;

        /**
         * @return Name of the node the process is running on.
         */
        std::string getNodeName() const;

        /**
         * @return Is the node idle?
         */
        bool isIdle() const;

        /**
         * An element is smaller if and only if it is idle and the subsequent node
         * than is not idle.
         *
         * @return Object is smaller
         */
        bool operator<( const NodePoolListEntry& than ) const;

        /**
         * Two entries are equal if and only if their rank equals.
         */
        bool operator==( const NodePoolListEntry& than ) const;

        /**
         * Create string representation.
         */
        void toString(std::ostream& out) const;

        /**
         * Return string representation.
         */
        std::string toString() const;

        /**
         * We could increment the worker counter, but we actually add the level.
         * This way, ranks that already deployed very fine grid levels are not
         * as important as workers trying to fork rather coarse areas.
         */
        void addNewWorker();

        /**
         * Halves all the entries
         */
        void reduceNumberOfBookedWorkers();

        int getNumberOfBookedWorkers() const;
    };

    typedef std::list<NodePoolListEntry>   NodeContainer;

    /**
     * Logging Device
     */
    static tarch::logging::Log _log;

    /**
     * Tag on which the node pool works
     */
    int _tag;

    /**
     * The ist the list of active nodes. Every entry corresponds to one node.
     * If the entry is set, the node is working already and the server is not
     * allowed to deploy another job on this node. If the entry isn't set, there
     * is a job request message in the queue and the server is allowed to send
     * a job. Therefore, in the beginning, all the entries are set. For the very
     * first entry, corresponding to the server node, the invariant holds, that
     * this entry is set always.
     */
    NodeContainer _nodes;

    void logQueue( const RequestQueue& queue ) const;

    int getWorkersOfNode( int rank ) const;

    void updateNodeWeights();
  public:
  /**
   * Constructor
   *
   * Construct all the attributes.
   */
    FairNodePoolStrategy();
    virtual ~FairNodePoolStrategy();

    virtual void setNodePoolTag(int tag);
    virtual tarch::parallel::messages::WorkerRequestMessage extractElementFromRequestQueue(RequestQueue& queue);
    virtual void fillWorkerRequestQueue(RequestQueue& queue);
    virtual void addNode(const tarch::parallel::messages::RegisterAtNodePoolMessage& node );
    virtual void removeNode( int rank );
    virtual int getNumberOfIdleNodes() const;
    virtual void setNodeIdle( int rank );
    virtual int reserveNode(int forMaster);
    virtual bool isRegisteredNode(int rank) const;
    virtual bool isIdleNode(int rank) const;
    virtual void clearRegisteredNodes();
    virtual int getNumberOfRegisteredNodes() const;
    virtual std::string toString() const;
    virtual bool hasIdleNode(int forMaster) const;
    virtual int removeNextIdleNode();
};

#endif
