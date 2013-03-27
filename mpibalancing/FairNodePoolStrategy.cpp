#include "mpibalancing/FairNodePoolStrategy.h"


#include <sstream>
#include <limits>


tarch::logging::Log mpibalancing::FairNodePoolStrategy::_log( "mpibalancing::FairNodePoolStrategy" );


mpibalancing::FairNodePoolStrategy::FairNodePoolStrategy():
  NodePoolStrategy(),
  _tag(-1),
  _nodes() {
}


mpibalancing::FairNodePoolStrategy::~FairNodePoolStrategy() {
}


void mpibalancing::FairNodePoolStrategy::fillWorkerRequestQueue(RequestQueue& queue) {
  #ifdef Parallel
  assertion( _tag >= 0 );
  while ( tarch::parallel::messages::WorkerRequestMessage::isMessageInQueue(_tag, true) ) {
    tarch::parallel::messages::WorkerRequestMessage message;
    message.receive(MPI_ANY_SOURCE,_tag, true);
    queue.push_back( message );
  }
  #endif
}


void mpibalancing::FairNodePoolStrategy::logQueue( const RequestQueue& queue ) const {
  #ifdef Parallel
  if (queue.empty()) {
  _log.debug( "logQueue()", "queue is empty" );
  }
  else {
    std::ostringstream msg;
    msg << "queue: ";

    for (RequestQueue::const_iterator p = queue.begin(); p != queue.end(); p++ ) {
      msg << p->getSenderRank() << ",";
    }
    _log.debug( "logQueue()", msg.str() );
  }
  #endif
}


tarch::parallel::messages::WorkerRequestMessage mpibalancing::FairNodePoolStrategy::extractElementFromRequestQueue(RequestQueue& queue) {
  assertion( !queue.empty() );

  RequestQueue::iterator pResultInQueue;
  int                    workersOfRankCurrentlyAnswered = std::numeric_limits<int>::max();

  #ifdef Parallel
  for (RequestQueue::iterator p = queue.begin(); p != queue.end(); p++) {
    if (getWorkersOfNode(p->getSenderRank()) < workersOfRankCurrentlyAnswered) {
      workersOfRankCurrentlyAnswered = getWorkersOfNode(p->getSenderRank());
      pResultInQueue                 = p;
    }
  }

#endif
  tarch::parallel::messages::WorkerRequestMessage result = *pResultInQueue;
  queue.erase(pResultInQueue);

  return result;
}


int mpibalancing::FairNodePoolStrategy::getWorkersOfNode( int rank ) const {
  int result = -1;

  for (NodeContainer::const_iterator p=_nodes.begin(); p!=_nodes.end(); p++) {
    if (p->getRank()==rank) {
      result = p->getNumberOfBookedWorkers();
    }
  }

  assertion( result>=0 );
  return result;
}


void mpibalancing::FairNodePoolStrategy::addNode(const tarch::parallel::messages::RegisterAtNodePoolMessage& node) {
  #ifdef Parallel
  assertion( !isRegisteredNode(node.getSenderRank()) );

  logTraceInWith1Argument( "addNode(...)", node.getSenderRank() );
  NodePoolListEntry newEntry(
    node.getSenderRank(),
    tarch::parallel::StringTools::convert(node.getNodeName())
  );
  _nodes.push_back( newEntry ) ;
  _nodes.sort();
  logTraceOutWith1Argument( "addNode(...)", newEntry.toString() );
  #endif
}


void mpibalancing::FairNodePoolStrategy::removeNode( int rank ) {
  assertion( isRegisteredNode(rank) );

  for (
    NodeContainer::iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      #ifdef Debug
      _log.debug( "removeNode(int)", "remove entry " + p->toString() );
      #endif
      _nodes.erase(p);
      _nodes.sort();
      return;
    }
  }
}


bool mpibalancing::FairNodePoolStrategy::hasIdleNode(int forMaster) const {
  return !_nodes.empty() &&
         _nodes.front().isIdle();
}


int mpibalancing::FairNodePoolStrategy::removeNextIdleNode() {
  assertion1( !_nodes.empty(), _nodes.size() );
  assertion3( FairNodePoolStrategy::hasIdleNode(-1), _nodes.size(), _nodes.front().toString(), _nodes.begin()->isIdle() );

  int result = _nodes.front().getRank();
  _nodes.pop_front();
  return result;
}


int mpibalancing::FairNodePoolStrategy::getNumberOfIdleNodes() const {
  int result = 0;
  NodeContainer::const_iterator p = _nodes.begin();
  while (p != _nodes.end()&& p->isIdle() ) {
  p++;
  result++;
  }
  return result;
}


void mpibalancing::FairNodePoolStrategy::setNodeIdle( int rank ) {
  for (
    NodeContainer::iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      p->deActivate();
    }
  }

  _nodes.sort();
}


bool mpibalancing::FairNodePoolStrategy::isRegisteredNode(int rank) const {
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank ) {
      return true;
    }
  }
  return false;
}


bool mpibalancing::FairNodePoolStrategy::isIdleNode(int rank) const {
  assertion1( isRegisteredNode(rank), rank );
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    if ( p->getRank() == rank && p->isIdle() ) {
      return true;
    }
  }
  return false;
}


void mpibalancing::FairNodePoolStrategy::clearRegisteredNodes() {
  _nodes.clear();
}


int mpibalancing::FairNodePoolStrategy::getNumberOfRegisteredNodes() const {
  return static_cast<int>( _nodes.size() );
}


std::string mpibalancing::FairNodePoolStrategy::toString() const {
  std::ostringstream result;
  for (
    NodeContainer::const_iterator p = _nodes.begin();
    p != _nodes.end();
    p++
  ) {
    p->toString(result);
  }
  return result.str();
}


int mpibalancing::FairNodePoolStrategy::reserveNode(int forMaster) {
  assertion1(hasIdleNode(forMaster),forMaster);

  NodePoolListEntry result = _nodes.front();
  _nodes.pop_front();

  logDebug( "getFreeNode(int)", "found free node: " << result.toString() );

  result.activate();
  _nodes.push_back(result);
  _nodes.sort();

  for (NodeContainer::iterator p=_nodes.begin(); p!=_nodes.end(); p++ ) {
    if (p->getRank()==forMaster) {
      p->addNewWorker();
    }
  }

  updateNodeWeights();

  return result.getRank();
}


void mpibalancing::FairNodePoolStrategy::updateNodeWeights() {
  const double reductionOfWeights = 1.0 / (getNumberOfRegisteredNodes() - getNumberOfIdleNodes() + 1.0);
  for (NodeContainer::iterator p=_nodes.begin(); p!=_nodes.end(); p++ ) {
    p->reduceNumberOfBookedWorkers(reductionOfWeights);
  }
}


void mpibalancing::FairNodePoolStrategy::setNodePoolTag(int tag) {
  _tag = tag;
}


mpibalancing::FairNodePoolStrategy::NodePoolListEntry::NodePoolListEntry(int rank, const std::string& name):
  _rank(rank),
  _bookedWorkers(0),
  _state(WORKING),
  _name(name) {
}


mpibalancing::FairNodePoolStrategy::NodePoolListEntry::~NodePoolListEntry() {
}


std::string mpibalancing::FairNodePoolStrategy::NodePoolListEntry::getNodeName() const {
  return _name;
}


std::string mpibalancing::FairNodePoolStrategy::NodePoolListEntry::toString() const {
  std::ostringstream out;
  toString(out);
  return out.str();
}


void mpibalancing::FairNodePoolStrategy::NodePoolListEntry::toString(std::ostream& out) const {
  out << "(rank:" << _rank;
  switch (_state) {
    case IDLE:    out << ",state:idle";     break;
    case WORKING: out << ",state:working";  break;
  }
  out << ",booked-workers:" << _bookedWorkers;
  out << ",name:" << _name << ")";
}


bool mpibalancing::FairNodePoolStrategy::NodePoolListEntry::operator==( const NodePoolListEntry& than ) const {
  return _rank==than._rank;
}


void mpibalancing::FairNodePoolStrategy::NodePoolListEntry::activate() {
  assertionEquals1( _state, IDLE, toString() );
  _state = WORKING;
}


void mpibalancing::FairNodePoolStrategy::NodePoolListEntry::deActivate() {
  _state         = IDLE;
  _bookedWorkers = 0;
}


int mpibalancing::FairNodePoolStrategy::NodePoolListEntry::getRank() const {
  return _rank;
}


bool mpibalancing::FairNodePoolStrategy::NodePoolListEntry::isIdle() const {
  return _state == IDLE;
}


bool mpibalancing::FairNodePoolStrategy::NodePoolListEntry::operator<( const mpibalancing::FairNodePoolStrategy::NodePoolListEntry& than ) const {
  return isIdle() && !than.isIdle();
}


void mpibalancing::FairNodePoolStrategy::NodePoolListEntry::addNewWorker() {
  _bookedWorkers+=1.0;
}


int mpibalancing::FairNodePoolStrategy::NodePoolListEntry::getNumberOfBookedWorkers() const {
  return _bookedWorkers;
}


void mpibalancing::FairNodePoolStrategy::NodePoolListEntry::reduceNumberOfBookedWorkers(double value) {
  assertion( value>=0.0 );
  _bookedWorkers -= value;
  if (_bookedWorkers<0.0) {
    _bookedWorkers = 0.0;
  }
}
