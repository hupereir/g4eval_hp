#include "ClearTrackSeeds_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>

//_____________________________________________________________________
ClearTrackSeeds_hp::ClearTrackSeeds_hp( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int ClearTrackSeeds_hp::Init(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int ClearTrackSeeds_hp::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int ClearTrackSeeds_hp::process_event(PHCompositeNode* topNode)
{
  // load track container from topnode
  const auto track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if( !track_map )
  {
    std::cout << "ClearTrackSeeds_hp::process_event - unable to load track map" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // keep track of seed pointers
  std::set<TrackSeed*> seeds;

  // loop over tracks
  for( const auto& [key,track]:*track_map )
  {
    if( track )
    {
      seeds.insert( track->get_tpc_seed() );
      seeds.insert( track->get_silicon_seed() );
    }
  }

  // delete seeds
  if( Verbosity() )
  { std::cout << "ClearTrackSeeds_hp::process_event - deleting " << seeds.size() << " seeds" << std::endl; }

  for (auto* s : seeds)
  { delete s; }


  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int ClearTrackSeeds_hp::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }
