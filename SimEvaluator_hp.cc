#include "SimEvaluator_hp.h"

#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/SvtxTrack.h>

#include <Geant4/G4SystemOfUnits.hh>
#include <HepMC/GenEvent.h>
#include <TDatabasePDG.h>

#include <algorithm>
#include <bitset>
#include <iostream>
#include <numeric>

//_____________________________________________________________________
namespace
{

  // print vectors
  template<class T>
  void print_vector( const std::string& name, const std::vector<T> values )
  {
    std::cout << "  " << name << "[" << values.size() << "] = {" << std::endl << "    ";
    int current = 0;
    int total = 0;
    for( const auto& value:values )
    {
      std::cout << value;
      if( total < int( values.size()-1 ) ) std::cout << ", ";
      ++current;
      ++total;
      if( current == 10 ) { std::cout << std::endl << "    "; current = 0; }
    }
    std::cout << " };" << std::endl;
  };

  //! convenient class to use range for loop on a pair of iterators, as returned by most sphenix container classes
  template<class T>
  class range_adaptor
  {
    public:
    range_adaptor( const std::pair<T,T>& range ) : m_range( range ) {};
    const T& begin() const { return m_range.first; }
    const T& end() const { return m_range.second; }

    private:
    const std::pair<T,T> m_range;
  };

  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  template<class T> inline T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  //! pt
  template<class T> inline T get_pt( T px, T py ) { return std::sqrt( square(px) + square(py) ); }

  //! p
  template<class T> inline T get_p( T px, T py, T pz ) { return std::sqrt( square(px) + square(py) + square(pz) ); }

  //! eta
  template<class T> inline T get_eta( T p, T pz ) { return std::log( (p+pz)/(p-pz) )/2; }

  //! true if particle is primary
  inline bool is_primary( PHG4Particle* particle )
  { return particle->get_parent_id() == 0; }

  //_____________________________________________________________________
  SimEvaluator_hp::EventStruct create_event(PHHepMCGenEventMap* eventmap)
  {
    SimEvaluator_hp::EventStruct eventStruct;
    if( eventmap )
    {
      // keep track of the main underlying event, associated to key 0
      PHHepMCGenEvent* mainevent = nullptr;

      // check event keys
      for( const auto& pair:eventmap->get_map() )
      {
        std::cout << "SimEvaluator_hp::create_event - adding event with key: " << pair.first << std::endl;
        ++eventStruct._nevt;
        if( pair.first < 0 ) ++eventStruct._nevt_bg;
        else ++eventStruct._nevt_active;

        if( pair.first == 0 ) mainevent = pair.second;
      }

      // store centrality and reaction plane angle
      if( mainevent )
      {
        auto hi = mainevent->getEvent()->heavy_ion();
        if( hi )
        {
          eventStruct._bimp = hi->impact_parameter();
          eventStruct._rplane = hi->event_plane_angle();
          std::cout << "SimEvaluator_hp::create_event - impact parameter: " << eventStruct._bimp << std::endl;
        } else {
          std::cout << "SimEvaluator_hp::create_event - unable to load heavy ion data" << std::endl;
        }
      } else {
        std::cout << "SimEvaluator_hp::create_event - no event found for key 0" << std::endl;
      }
    }
    return eventStruct;
  }

  //! create track struct from struct from svx track
  SimEvaluator_hp::VertexStruct create_vertex( PHG4VtxPoint* vertex )
  {
    SimEvaluator_hp::VertexStruct vertexStruct;
    vertexStruct._id = vertex->get_id();
    vertexStruct._x = vertex->get_x();
    vertexStruct._y = vertex->get_y();
    vertexStruct._z = vertex->get_z();
    vertexStruct._t = vertex->get_t();
    return vertexStruct;
  }

  //! create track struct from struct from svx track
  SimEvaluator_hp::ParticleStruct create_particle( PHG4Particle* particle )
  {
    SimEvaluator_hp::ParticleStruct particleStruct;
    particleStruct._trkid = particle->get_track_id();
    particleStruct._parent_id = particle->get_parent_id();
    particleStruct._primary_id = particle->get_primary_id();
    particleStruct._pid = particle->get_pid();

    particleStruct._vtx_id = particle->get_vtx_id();

    // assign charge
    auto pdgParticle = TDatabasePDG::Instance()->GetParticle(particleStruct._pid);
    if( pdgParticle )
    { particleStruct._charge = pdgParticle->Charge()/3; }

    particleStruct._is_primary = is_primary( particle );
    particleStruct._px = particle->get_px();
    particleStruct._py = particle->get_py();
    particleStruct._pz = particle->get_pz();
    particleStruct._pt = get_pt( particle->get_px(), particle->get_py() );
    particleStruct._p = get_p( particle->get_px(), particle->get_py(), particle->get_pz() );
    particleStruct._e = particle->get_e();

    particleStruct._eta = get_eta( particleStruct._p, particleStruct._pz );
    return particleStruct;
  }

  //! create g4hit struct from G4Hit
  SimEvaluator_hp::G4HitStruct create_g4hit( PHG4Hit* g4hit )
  {
    SimEvaluator_hp::G4HitStruct g4hitstruct;
    g4hitstruct._layer = g4hit->get_layer();
    g4hitstruct._trkid = g4hit->get_trkid();
    g4hitstruct._t = g4hit->get_avg_t();
    g4hitstruct._x = g4hit->get_avg_x();
    g4hitstruct._y = g4hit->get_avg_y();
    g4hitstruct._z = g4hit->get_avg_z();
    g4hitstruct._r = get_r( g4hitstruct._x, g4hitstruct._y );
    g4hitstruct._phi = std::atan2( g4hitstruct._y, g4hitstruct._x );
    g4hitstruct._length = std::sqrt(
      square(g4hit->get_x(1) - g4hit->get_x(0))+
      square(g4hit->get_y(1) - g4hit->get_y(0))+
      square(g4hit->get_z(1) - g4hit->get_z(0)));
    return g4hitstruct;
  }

  std::ostream& operator << (std::ostream& out, const PHG4VtxPoint& vertex )
  {
    out << "( " << vertex.get_x() << ", " << vertex.get_y() << ", " << vertex.get_z() << ", " << vertex.get_t() << ")";
    return out;
  }

}

//_____________________________________________________________________
void SimEvaluator_hp::Container::Reset()
{
  _events.clear();
  _vertex_list.clear();
  _particle_list.clear();
  _g4hits.clear();
}

//_____________________________________________________________________
SimEvaluator_hp::SimEvaluator_hp( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int SimEvaluator_hp::Init(PHCompositeNode* topnode )
{

  // find DST node
  PHNodeIterator iter(topnode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "SimEvaluator_hp::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "SimEvaluator_hp::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  auto newNode = new PHIODataNode<PHObject>( new Container, "SimEvaluator_hp::Container", "PHObject" );
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SimEvaluator_hp::InitRun(PHCompositeNode*)
{
  // print_tpc( topnode );
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SimEvaluator_hp::process_event(PHCompositeNode* topnode)
{
  // load nodes
  auto res =  load_nodes(topnode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  // clear embedding id map
  m_g4embed_map.clear();

  // for debugging
  check_genevent();

  if( m_flags&EvalEvent) fill_event();
  if( m_flags&EvalVertices) fill_vertices();
  if( m_flags&EvalParticles)
  {
    fill_g4particle_map();
    fill_particles();
  }
  if( m_flags&EvalHits) fill_hits();
  if( m_flags&PrintVertices) print_vertices();

  m_g4particle_map.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int SimEvaluator_hp::End(PHCompositeNode*)
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int SimEvaluator_hp::load_nodes( PHCompositeNode* topnode )
{
  // local container
  m_container = findNode::getClass<Container>(topnode, "SimEvaluator_hp::Container");

  // hep mc
  m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topnode, "PHHepMCGenEventMap");

  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topnode, "G4TruthInfo");

  // g4hits
  m_g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topnode, "G4HIT_MVTX");
  m_g4hits_intt = findNode::getClass<PHG4HitContainer>(topnode, "G4HIT_INTT");
  m_g4hits_tpc = findNode::getClass<PHG4HitContainer>(topnode, "G4HIT_TPC");
  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topnode, "G4HIT_MICROMEGAS");
  m_g4hits_cemc = findNode::getClass<PHG4HitContainer>(topnode, "G4HIT_CEMC");
  m_g4hits_hcalin = findNode::getClass<PHG4HitContainer>(topnode, "G4HIT_HCALIN");
  m_g4hits_hcalout = findNode::getClass<PHG4HitContainer>(topnode, "G4HIT_HCALOUT");

  // event header
  m_eventheader = findNode::getClass<EventHeader>(topnode, "EventHeader");

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void SimEvaluator_hp::print_tpc( PHCompositeNode* topnode )
{
  // get relevant node
  const auto container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topnode, "CYLINDERCELLGEOM_SVTX");
  if( !container ) return;

  std::vector<int> phibins;
  std::vector<int> zbins;
  std::vector<float> thickness;
  std::vector<float> radius;
  const range_adaptor<PHG4TpcCylinderGeomContainer::ConstIterator> range(std::move( container->get_begin_end()));
  for( const auto&pair:range )
  {
    phibins.push_back( pair.second->get_phibins() );
    zbins.push_back( pair.second->get_zbins() );
    thickness.push_back( pair.second->get_thickness() );
    radius.push_back(pair.second->get_radius()-0.5*pair.second->get_thickness());
  }

  // add end radius for last layer
  radius.push_back( radius.back() + thickness.back() );

  print_vector( "phibins", phibins );
  print_vector( "zbins", zbins );
  print_vector( "thickness", thickness );
  print_vector( "radius", radius );

}

//_____________________________________________________________________
void SimEvaluator_hp::fill_g4particle_map()
{
  m_g4particle_map.clear();
  for( const auto& container: {m_g4hits_tpc, m_g4hits_intt, m_g4hits_mvtx, m_g4hits_micromegas} )
  {
    if( !container ) continue;

    // loop over hits
    const auto range = container->getHits();
    for( auto iter = range.first; iter != range.second; ++iter )
    {
      const auto map_iter = m_g4particle_map.lower_bound( iter->second->get_trkid() );
      if( map_iter != m_g4particle_map.end() && map_iter->first == iter->second->get_trkid() )
      {
        map_iter->second |= (1LL<<iter->second->get_layer());
      } else {
        m_g4particle_map.insert( map_iter, std::make_pair( iter->second->get_trkid(), 1LL<<iter->second->get_layer() ) );
      }
    }
  }

}

//_____________________________________________________________________
void SimEvaluator_hp::check_genevent()
{}

//_____________________________________________________________________
void SimEvaluator_hp::fill_event()
{
  if( !( m_container && m_geneventmap ) ) return;

  // clear vertices from previous event
  m_container->clearEventList();

  // create event and store pileup information
  auto event = create_event(m_geneventmap);

  // count number of primary particles with pt > 0.5
  if( m_g4truthinfo )
  {

    // get main vertex id
    event._main_vertex_id = m_g4truthinfo->GetPrimaryVertexIndex();

    // loop over particles
    for( const auto& [id, particle]:range_adaptor(m_g4truthinfo->GetPrimaryParticleRange()) )
    {
      if( !particle ) continue;

      // create convenient structure from particle
      auto pstruct = create_particle( particle );

      // check pt, charge and pseudo rapidity
      if( pstruct._pt > 0.5 && std::abs( pstruct._eta ) < 1 && pstruct._charge != 0 )
      { ++event._nparticles; }
    }
  }

  // add to container
  m_container->addEvent(event);

}

//_____________________________________________________________________
void SimEvaluator_hp::fill_vertices()
{
  if( !( m_container && m_g4truthinfo ) ) return;

  // clear vertices from previous event
  m_container->clearVertexList();

  // get main primary vertex id
  const auto main_vertex_id = m_g4truthinfo->GetPrimaryVertexIndex();

  auto range = m_g4truthinfo->GetVtxRange();
  // auto range = m_g4truthinfo->GetPrimaryVtxRange();
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    auto vertex = iter->second;
    if( vertex )
    {
      auto vertexStruct = create_vertex(vertex);
      vertexStruct._embed = m_g4truthinfo->isEmbededVtx(vertex->get_id());
      vertexStruct._is_main_vertex = (vertex->get_id() == main_vertex_id);
      m_container->addVertex(vertexStruct);
    }
  }
}

//_____________________________________________________________________
void SimEvaluator_hp::fill_particles()
{

  if( !( m_container && m_g4truthinfo ) ) return;

  // clear vertices from previous event
  m_container->clearParticleList();

//   // get the set of embeded ids from truth info
//   {
//     std::set<int> embed_ids;
//     for( const auto& [id,embed_id]:range_adaptor(m_g4truthinfo->GetEmbeddedTrkIds()) );
//     { embed_ids.insert( embed_id ); }
//   }

  for( const auto& [id,particle]:range_adaptor(m_g4truthinfo->GetParticleRange()))
  {
    if( !particle ) continue;
    auto particleStruct = create_particle( particle );

    // vertex
    auto vtx = m_g4truthinfo->GetVtx(particle->get_vtx_id());
    if( vtx )
    {
      particleStruct._x = vtx->get_x();
      particleStruct._y = vtx->get_y();
      particleStruct._z = vtx->get_z();
      particleStruct._t = vtx->get_t();
    }

    // embed index
    particleStruct._embed = get_embed( particle );

    // hit mask
    const auto iter( m_g4particle_map.find( particle->get_track_id() ) );
    if( iter !=  m_g4particle_map.cend() )
    { particleStruct._mask = iter->second; }

    m_container->addParticle( particleStruct );
  }
}

//_____________________________________________________________________
void SimEvaluator_hp::fill_hits()
{
  if( !m_container ) return;

  // clear container
  m_container->clearG4Hits();

  {
    // map tracker detector id to g4hit container
    std::map<TrkrDefs::TrkrId, PHG4HitContainer*> containers = {
      { TrkrDefs::mvtxId, m_g4hits_mvtx },
      { TrkrDefs::inttId, m_g4hits_intt },
      { TrkrDefs::tpcId, m_g4hits_tpc },
      { TrkrDefs::micromegasId, m_g4hits_micromegas }
    };

    // loop over tracker containers
    for( const auto& pair:containers )
    {
      const auto& container( pair.second );
      if( !container ) continue;

      // load g4hits
      const auto range = container->getHits();
      for( auto iter = range.first; iter != range.second; ++iter )
      {

        // create hit
        auto g4hit = create_g4hit( iter->second );

        // detector id
        g4hit._detid = static_cast<int16_t>( pair.first );

        // embed id
        if( m_g4truthinfo )
        {
          g4hit._embed = get_embed( iter->second );
          g4hit._pid = get_pid( iter->second );
        }

        // add
        m_container->addG4Hit( g4hit );
      }
    }
  }

  {
    // do the same with calorimeter
    // map tracker detector id to g4hit container
    std::map<SvtxTrack::CAL_LAYER, PHG4HitContainer*> containers = {
      { SvtxTrack::CEMC, m_g4hits_cemc },
      { SvtxTrack::HCALIN, m_g4hits_hcalin },
      { SvtxTrack::HCALOUT, m_g4hits_hcalout },
    };

    // loop over tracker containers
    for( const auto& pair:containers )
    {
      const auto& container( pair.second );
      if( !container ) continue;

      // load g4hits
      const auto range = container->getHits();
      for( auto iter = range.first; iter != range.second; ++iter )
      {

        // create hit
        auto g4hit = create_g4hit( iter->second );

        // detector id
        g4hit._caloid = static_cast<int16_t>( pair.first );

        // embed id
        if( m_g4truthinfo )
        {
          g4hit._embed = get_embed( iter->second );
          g4hit._pid = get_pid( iter->second );
        }

        // add
        m_container->addG4Hit( g4hit );
      }
    }
  }

}

//_____________________________________________________________________
void SimEvaluator_hp::print_vertices()
{
  if( !m_g4truthinfo ) return;

  // get main primary vertex id
  const auto main_vertex_id = m_g4truthinfo->GetPrimaryVertexIndex();
  const auto vertex = m_g4truthinfo->GetPrimaryVtx( main_vertex_id );
  if( vertex )
  {

    std::cout << "SimEvaluator_hp::print_vertices - main primary vertex: " << *vertex << std::endl;

  } else {

    std::cerr << "SimEvaluator_hp::print_vertices - no main primary vertex found." << std::endl;

  }

  auto range = m_g4truthinfo->GetPrimaryVtxRange();
  std::cout << "SimEvaluator_hp::print_vertices - primary vertex count: " << std::distance( range.first, range.second ) << std::endl;
  for( auto iter = range.first; iter != range.second; ++iter )
  {
    auto vertex = iter->second;
    if( vertex && vertex->get_id() != main_vertex_id )
    { std::cout << "SimEvaluator_hp::print_vertices - primary vertex: " << *vertex << std::endl; }
  }

}

//_____________________________________________________________________
int SimEvaluator_hp::get_embed( PHG4Hit* hit )
{
  if(!(m_g4truthinfo && hit)) return 0;

  // get trk id
  const auto trk_id = hit->get_trkid();
  if( trk_id > 0 ) return  m_g4truthinfo->isEmbeded( trk_id );
  else
  {
    // for secondary particles, check cache or get the corresponding PHParticle
    // check if already in map
    const auto iter = m_g4embed_map.lower_bound( trk_id );
    if( iter != m_g4embed_map.end() && iter->first == trk_id ) return iter->second;

    const auto particle = m_g4truthinfo->GetParticle( trk_id );
    const auto embed = get_embed( particle );
    return m_g4embed_map.insert( iter, std::make_pair( trk_id, embed ) )->second;
  }
}

//_____________________________________________________________________
int SimEvaluator_hp::get_embed( PHG4Particle* particle ) const
{ return (m_g4truthinfo && particle) ? m_g4truthinfo->isEmbeded( particle->get_primary_id() ):0; }

//_____________________________________________________________________
int SimEvaluator_hp::get_pid( PHG4Hit* hit )
{
  if(!(m_g4truthinfo && hit)) return 0;

  // get trk id
  const auto trk_id = hit->get_trkid();
  const auto iter = m_pid_map.lower_bound( trk_id );
  if( iter != m_pid_map.end() && iter->first == trk_id ) return iter->second;

  // get particle, pid, insert in map and return
  const auto particle = m_g4truthinfo->GetParticle( trk_id );
  const auto pid = particle->get_pid();
  return m_pid_map.insert( iter, std::make_pair( trk_id, pid ) )->second;
}
