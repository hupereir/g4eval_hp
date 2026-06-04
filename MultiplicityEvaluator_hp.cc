#include "MultiplicityEvaluator_hp.h"

#include <ffarawobjects/MicromegasRawHitContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <ffarawobjects/Gl1RawHit.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

//_____________________________________________________________________
namespace
{

  //! range adaptor to be able to use range-based for loop
  template<class T> class range_adaptor
  {
    public:
    range_adaptor( const T& range ):m_range(range){}
    const typename T::first_type& begin() {return m_range.first;}
    const typename T::second_type& end() {return m_range.second;}
    private:
    T m_range;
  };
}

//_____________________________________________________________________
MultiplicityEvaluator_hp::MultiplicityEvaluator_hp( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int MultiplicityEvaluator_hp::Init(PHCompositeNode* topNode )
{
  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "MultiplicityEvaluator_hp::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "MultiplicityEvaluator_hp::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  // add container to output tree
  auto newNode = new PHIODataNode<PHObject>( new Container, "MultiplicityEvaluator_hp::Container","PHObject");

  // overwrite split level for easier offline browsing
  newNode->SplitLevel(99);
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MultiplicityEvaluator_hp::process_event(PHCompositeNode* topNode)
{

  // current event multiplicity
  MultiplicityStruct current_mult_struct;

  // get BCO from gl1
  auto gl1rawhit = findNode::getClass<Gl1RawHit>(topNode,"GL1RAWHIT");
  if( gl1rawhit )
  {
    current_mult_struct._gtm_bco = gl1rawhit->get_bco()&0xFFFFFFFFFF;
  } else {
    std::cout << "MultiplicityEvaluator_hp::process_event - GL1RAWHIT not found" << std::endl;
  }

  // raw hits multiplicity
  auto rawhitcontainer = findNode::getClass<MicromegasRawHitContainer>(topNode, "MICROMEGASRAWHIT" );
  if( rawhitcontainer )
  {
    current_mult_struct._rawhits =  rawhitcontainer->get_nhits();
  } else {
    std::cout << "MultiplicityEvaluator_hp::process_event - MICROMEGASRAWHIT not found" << std::endl;
  }

  // calibrated hits multiplicity
  auto hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if( hitsetcontainer )
  {
    // loop over all TPOT hitsets
    for( const auto& [hitsetkey, hitset]:range_adaptor( hitsetcontainer->getHitSets(TrkrDefs::micromegasId ) ) )
    {
      current_mult_struct._hits += hitset->size();
    }
  } else {
    std::cout << "MultiplicityEvaluator_hp::process_event - TRKR_HITSET not found" << std::endl;
  }

  // clusters
  auto cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if( cluster_map )
  {
    // loop over TPOT hitset keys
    for( const auto& hitsetkey:cluster_map->getHitSetKeys(TrkrDefs::micromegasId) )
    {
      const auto& range = cluster_map->getClusters(hitsetkey);
      current_mult_struct._clusters += std::distance( range.first, range.second );
    }
  } else {
    std::cout << "MultiplicityEvaluator_hp::process_event - TRKR_CLUSTER not found" << std::endl;
  }

  // update container
  if( m_container )
  {
    m_container->set_previous_multiplicity( m_container->current_multiplicity() );
    m_container->set_current_multiplicity( current_mult_struct );
  } else {
    std::cout << "MultiplicityEvaluator_hp::process_event - m_container not found" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________
int MultiplicityEvaluator_hp::End(PHCompositeNode*)
{ return Fun4AllReturnCodes::EVENT_OK; }
