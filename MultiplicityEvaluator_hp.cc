#include "MultiplicityEvaluator_hp.h"

#include <ffarawobjects/MicromegasRawHit.h>
#include <ffarawobjects/MicromegasRawHitContainer.h>
#include <ffarawobjects/Gl1Packet.h>
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
  MultiplicityStruct current_mult;
  MultiplicityStruct::Array det_mult;

  // get BCO from gl1
  auto gl1rawhit = findNode::getClass<Gl1Packet>(topNode,"GL1RAWHIT");
  if( gl1rawhit )
  {

    current_mult._gtm_bco = gl1rawhit->getBCO()&0xFFFFFFFFFF;
    for( auto&& mult:det_mult ) { mult._gtm_bco = gl1rawhit->getBCO()&0xFFFFFFFFFF; }

  } else {

    std::cout << "MultiplicityEvaluator_hp::process_event - GL1RAWHIT not found" << std::endl;

  }

  auto get_det_id = []( int layer, int tile ) { return (layer-55)*MicromegasDefs::m_ntiles + tile; };

  // raw hits multiplicity
  auto rawhitcontainer = findNode::getClass<MicromegasRawHitContainer>(topNode, "MICROMEGASRAWHIT" );
  if( rawhitcontainer )
  {

    current_mult._rawhits =  rawhitcontainer->get_nhits();
    for( unsigned int i = 0; i < rawhitcontainer->get_nhits(); ++i )
    {
      // get hit, fee, hitsetkey, tile and layer
      const auto& rawhit = rawhitcontainer->get_hit(i);
      const int fee = m_mapping.get_new_fee_id(rawhit->get_fee());
      const TrkrDefs::hitsetkey hitsetkey = m_mapping.get_hitsetkey(fee);

      const int layer = int(TrkrDefs::getLayer(hitsetkey));
      const int tile = int(MicromegasDefs::getTileId(hitsetkey));
      const int detid = get_det_id( layer, tile );

      // increment count
      ++det_mult[detid]._rawhits;

      // add samples
      const auto sample = rawhit->get_sample_begin();
      current_mult._rawhit_samples.emplace_back( sample );
      det_mult[detid]._rawhit_samples.emplace_back( sample );

      if( sample < 1024 )
      {
        ++current_mult._rawhits_truncated;
        ++det_mult[detid]._rawhits_truncated;
      }
    }

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
      // increment total hit multiplicity
      current_mult._hits += hitset->size();

      // increment per detector hit multiplicity
      const int layer = int(TrkrDefs::getLayer(hitsetkey));
      const int tile = int(MicromegasDefs::getTileId(hitsetkey));
      const int detid = get_det_id( layer, tile );
      det_mult[detid]._hits += hitset->size();

      // loop over hits
      for( const auto& [hitkey,hit]:range_adaptor( hitset->getHits() ) )
      {
        const auto sample = MicromegasDefs::getSample(hitkey);
        current_mult._hit_samples.emplace_back( sample );
        det_mult[detid]._hit_samples.emplace_back( sample );
      }

    }

  } else {

    // std::cout << "MultiplicityEvaluator_hp::process_event - TRKR_HITSET not found" << std::endl;

  }

  // clusters
  auto cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if( cluster_map )
  {

    // loop over TPOT hitset keys
    for( const auto& hitsetkey:cluster_map->getHitSetKeys(TrkrDefs::micromegasId) )
    {
      const auto& range = cluster_map->getClusters(hitsetkey);
      const auto& clusters = std::distance( range.first, range.second );
      current_mult._clusters += clusters;

      // increment per detector hit multiplicity
      const int layer = int(TrkrDefs::getLayer(hitsetkey));
      const int tile = int(MicromegasDefs::getTileId(hitsetkey));
      const int detid = get_det_id( layer, tile );
      det_mult[detid]._clusters += clusters;
    }

  } else {

    // std::cout << "MultiplicityEvaluator_hp::process_event - TRKR_CLUSTER not found" << std::endl;

  }

  // update container
  auto container = findNode::getClass<Container>(topNode, "MultiplicityEvaluator_hp::Container");
  if( container )
  {
    container->set_previous_multiplicity( container->current_multiplicity() );
    container->set_current_multiplicity( current_mult );

    container->set_previous_det_multiplicity( container->current_det_multiplicity() );
    container->set_current_det_multiplicity( det_mult );

  } else {
    std::cout << "MultiplicityEvaluator_hp::process_event - m_container not found" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________________
int MultiplicityEvaluator_hp::End(PHCompositeNode*)
{ return Fun4AllReturnCodes::EVENT_OK; }
