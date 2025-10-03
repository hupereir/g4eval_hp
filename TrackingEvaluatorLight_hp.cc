#include "TrackingEvaluatorLight_hp.h"

#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawTowerDefs.h>

#include <ffarawobjects/Gl1RawHit.h>
#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/MicromegasRawHit.h>
#include <ffarawobjects/MicromegasRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/CMFlashCluster.h>
#include <trackbase/CMFlashClusterContainer.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/MvtxEventInfo.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrClusterCrossingAssoc.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <TVector3.h>

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <numeric>

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

  /// square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  //! radius
  template<class T> inline constexpr T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  //! pt
  template<class T> T get_pt( const T& px, const T& py ) { return std::sqrt( square(px) + square(py) ); }

  //! p
  template<class T> T get_p( const T& px, const T& py, const T& pz ) { return std::sqrt( square(px) + square(py) + square(pz) ); }

  //! eta
  template<class T> T get_eta( const T& p, const T& pz ) { return std::log( (p+pz)/(p-pz) )/2; }

  //! eta
  template<class T> T get_eta( const T& px, const T& py, const T& pz )
  {
    const auto p = get_p(px,py,pz);
    return get_eta(p,pz);
  }

  //! delta phi with proper rounding
  template<class T>
  T delta_phi( T phi )
  {
    while( phi >= M_PI ) {phi -= 2*M_PI; }
    while( phi <-M_PI ) {phi += 2*M_PI; }
    return phi;
  }

  //! map calorimeter names to layer type
  // using calo_names_map_t = std::map<SvtxTrack::CAL_LAYER, std::string>;
  enum
  {
    TOPO_HCAL = 7
  };

  using calo_names_map_t = std::map<int, std::string>;
  const calo_names_map_t m_calo_names = {
    { SvtxTrack::CEMC, "CEMC" },
    { SvtxTrack::HCALIN, "HCALIN" },
    { SvtxTrack::HCALOUT, "HCALOUT" },
    { SvtxTrack::OUTER_CEMC, "OUTER_CEMC" },
    { SvtxTrack::OUTER_HCALIN, "OUTER_HCALIN" },
    { SvtxTrack::OUTER_HCALOUT, "OUTER_HCALOUT" },
    { TOPO_HCAL, "TOPO_HCAL" }
  };

  //! needed for weighted linear interpolation
  struct interpolation_data_t
  {
    using list = std::vector<interpolation_data_t>;
    double x() const { return position.x(); }
    double y() const { return position.y(); }
    double z() const { return position.z(); }

    double px() const { return momentum.x(); }
    double py() const { return momentum.y(); }
    double pz() const { return momentum.z(); }

    TVector3 position;
    TVector3 momentum;
    double weight = 1;
  };


  //! calculate the interpolation of member function called on all members in collection to the provided y_extrap
  template<double (interpolation_data_t::*accessor)() const>
  double interpolate_y( const interpolation_data_t::list& hits, double y_extrap )
  {

    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double sw = 0;
    double swy = 0;
    double swy2 = 0;
    double swx = 0;
    double swyx = 0;

    bool valid( false );
    for( const auto& hit:hits )
    {

      const double x = (hit.*accessor)();
      const double w = hit.weight;
      if( w <= 0 ) continue;

      valid = true;
      const double y = hit.y();

      sw += w;
      swy += w*y;
      swy2 += w*square(y);
      swx += w*x;
      swyx += w*x*y;
    }

    if( !valid ) return NAN;

    const auto alpha = (sw*swyx - swy*swx);
    const auto beta = (swy2*swx - swy*swyx);
    const auto denom = (sw*swy2 - square(swy));

    return ( alpha*y_extrap + beta )/denom;
  }

  //! calculate the interpolation of member function called on all members in collection to the provided y_extrap
  template<double (interpolation_data_t::*accessor)() const>
  double interpolate_r( const interpolation_data_t::list& hits, double r_extrap )
  {

    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double sw = 0;
    double swr = 0;
    double swr2 = 0;
    double swx = 0;
    double swrx = 0;

    bool valid( false );
    for( const auto& hit:hits )
    {

      const double x = (hit.*accessor)();
      const double w = hit.weight;
      if( w <= 0 ) continue;

      valid = true;
      const double r = get_r(hit.x(), hit.y());

      sw += w;
      swr += w*r;
      swr2 += w*square(r);
      swx += w*x;
      swrx += w*x*r;
    }

    if( !valid ) return NAN;

    const auto alpha = (sw*swrx - swr*swx);
    const auto beta = (swr2*swx - swr*swrx);
    const auto denom = (sw*swr2 - square(swr));

    return ( alpha*r_extrap + beta )/denom;
  }

  //! calculate the average of member function called on all members in collection
  template<double (interpolation_data_t::*accessor)() const>
  double average( const interpolation_data_t::list& hits )
  {
    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double sw = 0;
    double swx = 0;

    bool valid( false );
    for( const auto& hit:hits )
    {

      const double x = (hit.*accessor)();
      if(std::isnan(x)) continue;

      const double w = hit.weight;
      if( w <= 0 ) continue;

      valid = true;
      sw += w;
      swx += w*x;
    }

    if( !valid ) return NAN;
    return swx/sw;
  }

  //! get cluster keys from a given seed
  std::vector<TrkrDefs::cluskey> get_cluster_keys( TrackSeed* seed )
  {
    std::vector<TrkrDefs::cluskey> out;
    if( seed )
    { std::copy( seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter( out ) ); }

    return out;
  }

  //! get cluster keys from a given track
  std::vector<TrkrDefs::cluskey> get_cluster_keys( SvtxTrack* track )
  {
    std::vector<TrkrDefs::cluskey> out;
    for( const auto& seed: { track->get_silicon_seed(), track->get_tpc_seed() } )
    {
      if( seed )
      { std::copy( seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter( out ) ); }
    }

    return out;
  }

  //! true if a track is a primary
  inline int is_primary( PHG4Particle* particle )
  { return particle->get_parent_id() == 0; }

  //! get mask from track clusters
  int64_t get_mask( SvtxTrack* track )
  {
    const auto cluster_keys = get_cluster_keys( track );
    return std::accumulate( cluster_keys.begin(), cluster_keys.end(), int64_t(0),
      []( int64_t value, const TrkrDefs::cluskey& key ) {
        return TrkrDefs::getLayer(key)<64 ? value|(1LL<<TrkrDefs::getLayer(key)) : value;
      } );
  }

  //! return number of clusters of a given type
  template<int type>
    int get_clusters( SvtxTrack* track )
  {
    const auto cluster_keys = get_cluster_keys( track );
    return std::count_if( cluster_keys.begin(), cluster_keys.end(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == type; } );
  }

  //! return number of clusters of a given type
  int get_clusters_micromegas_phi( SvtxTrack* track )
  {
    const auto cluster_keys = get_cluster_keys( track );
    return std::count_if( cluster_keys.begin(), cluster_keys.end(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getLayer(key) == 55; } );
  }

  //! return number of clusters of a given type
  int get_clusters_micromegas_z( SvtxTrack* track )
  {
    const auto cluster_keys = get_cluster_keys( track );
    return std::count_if( cluster_keys.begin(), cluster_keys.end(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getLayer(key) == 56; } );
  }

  //! return number of states of a given type
  template<int type>
    int get_states( SvtxTrack* track )
  {
    return std::count_if( track->begin_states(), track->end_states(),
      []( const std::pair<float, SvtxTrackState*>& state_pair ) { return TrkrDefs::getTrkrId(state_pair.second->get_cluskey()) == type; } );
  }

  //! fill basic information to track struct
  template<class T>
    void fill_track_struct( T& trackStruct, SvtxTrack* track )
  {
    trackStruct._charge = track->get_charge();
    trackStruct._nclusters = get_cluster_keys( track ).size();
    trackStruct._mask = get_mask( track );
    trackStruct._nclusters_mvtx = get_clusters<TrkrDefs::mvtxId>(track);
    trackStruct._nclusters_intt = get_clusters<TrkrDefs::inttId>(track);
    trackStruct._nclusters_tpc = get_clusters<TrkrDefs::tpcId>(track);
    trackStruct._nclusters_micromegas = get_clusters<TrkrDefs::micromegasId>( track );

    trackStruct._nstates = std::distance( track->begin_states(), track->end_states() );
    trackStruct._nstates_mvtx = get_states<TrkrDefs::mvtxId>(track);
    trackStruct._nstates_intt = get_states<TrkrDefs::inttId>(track);
    trackStruct._nstates_tpc = get_states<TrkrDefs::tpcId>(track);
    trackStruct._nstates_micromegas = get_states<TrkrDefs::micromegasId>(track);

    trackStruct._chisquare = track->get_chisq();
    trackStruct._ndf = track->get_ndf();

    trackStruct._px = track->get_px();
    trackStruct._py = track->get_py();
    trackStruct._pz = track->get_pz();
    trackStruct._p = get_p( trackStruct._px, trackStruct._py, trackStruct._pz );
    trackStruct._pt = get_pt( trackStruct._px, trackStruct._py );
    trackStruct._eta = get_eta( trackStruct._p, trackStruct._pz );
  }

  //! create central membrane cluster struct
  TrackingEvaluatorLight_hp::CaloClusterStruct create_calo_cluster( int layer, RawCluster* cluster )
  {
    TrackingEvaluatorLight_hp::CaloClusterStruct calo_cluster_struct;
    calo_cluster_struct._layer = layer;
    calo_cluster_struct._size = cluster->getNTowers();
    calo_cluster_struct._x = cluster->get_x();
    calo_cluster_struct._y = cluster->get_y();
    calo_cluster_struct._z = cluster->get_z();

    calo_cluster_struct._r = get_r( calo_cluster_struct._x, calo_cluster_struct._y );
    calo_cluster_struct._phi = std::atan2( calo_cluster_struct._y, calo_cluster_struct._x );
    calo_cluster_struct._eta = get_eta( calo_cluster_struct._x,calo_cluster_struct._y,calo_cluster_struct._z );

    calo_cluster_struct._e = cluster->get_energy();
    calo_cluster_struct._chisquare = cluster->get_chi2();

    // loop over towers
    for( const auto& [index, energy]:range_adaptor(cluster->get_towers()))
    {
      TrackingEvaluatorLight_hp::TowerStruct tower;

      // get indexes
      tower._ieta = RawTowerDefs::decode_index1(index);
      tower._iphi = RawTowerDefs::decode_index2(index);
      tower._e = energy;
      calo_cluster_struct._towers.emplace_back( tower );
    }

    return calo_cluster_struct;
  }

  //! create track struct from struct from svx track
  TrackingEvaluatorLight_hp::TrackStruct create_track( SvtxTrack* track )
  {
    TrackingEvaluatorLight_hp::TrackStruct trackStruct;

    // fill basic information, also used in TrackStruct_small
    fill_track_struct( trackStruct, track );

    // fill additional information
    trackStruct._nclusters_micromegas_phi = get_clusters_micromegas_phi( track );
    trackStruct._nclusters_micromegas_z = get_clusters_micromegas_z( track );

    trackStruct._x = track->get_x();
    trackStruct._y = track->get_y();
    trackStruct._z = track->get_z();
    trackStruct._r = get_r( trackStruct._x, trackStruct._y );
    trackStruct._phi = std::atan2( trackStruct._y, trackStruct._x );
    return trackStruct;
  }

  // add truth information
  void add_truth_information( TrackingEvaluatorLight_hp::TrackStruct& track, PHG4Particle* particle )
  {
    if( particle )
    {
      track._is_primary = is_primary( particle );
      track._pid = particle->get_pid();
      track._truth_px = particle->get_px();
      track._truth_py = particle->get_py();
      track._truth_pz = particle->get_pz();
      track._truth_pt = get_pt( track._truth_px, track._truth_py );
      track._truth_p = get_p( track._truth_px, track._truth_py, track._truth_pz );
      track._truth_eta = get_eta( track._truth_p, track._truth_pz );
    }
  }

  [[maybe_unused]] std::ostream& operator << (std::ostream& out, const TVector3& position)
  {
    out << "(" << position.x() << ", " << position.y() << ", " << position.z() << ")";
    return out;
  }

  /// streamer
  [[maybe_unused]] std::ostream& operator<<(std::ostream& out, const Acts::Vector3& v)
  {
    out << "(" << v.x() << "," << v.y() << "," << v.z() << ")";
    return out;
  }

  /// streamer
  template<class T>
  std::ostream& operator<<(std::ostream& out, std::vector<T> v)
  {
    if( v.empty() ) { out << "{}"; }
    else {
      out << "{";
      for( const auto& key:v )
      { out << " " << key; }
      out << "}";
    }
    return out;
  }

  template<class T>
  std::ostream& operator<<(std::ostream& out, std::set<T> v)
  {
    if( v.empty() ) { out << "{}"; }
    else {
      out << "{";
      for( const auto& key:v )
      { out << " " << key; }
      out << "}";
    }
    return out;
  }
}

//_____________________________________________________________________
void TrackingEvaluatorLight_hp::Container::Reset()
{
  _tracks.clear();
}

//_____________________________________________________________________
TrackingEvaluatorLight_hp::TrackingEvaluatorLight_hp( const std::string& name ):
  SubsysReco( name),
  m_calo_min_energy( {
    {SvtxTrack::CEMC, 0.15},
    {SvtxTrack::HCALIN, 0},
    {SvtxTrack::HCALOUT, 0},
    {TOPO_HCAL, 0}
  })
{}

//_____________________________________________________________________
int TrackingEvaluatorLight_hp::Init(PHCompositeNode* topNode )
{
  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "TrackingEvaluatorLight_hp::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "TrackingEvaluatorLight_hp::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  // add container to output tree
  auto newNode = new PHIODataNode<PHObject>( new Container, "TrackingEvaluatorLight_hp::Container","PHObject");

  // overwrite split level for easier offline browsing
  newNode->SplitLevel(99);
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluatorLight_hp::InitRun(PHCompositeNode* topNode)
{

  // TPC geometry (to check ADC clock frequency)
  auto geom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  assert(geom);

  const auto AdcClockPeriod = geom->GetFirstLayerCellGeom()->get_zstep();
  std::cout << "TrackingEvaluatorLight_hp::Init - AdcClockPeriod: " << AdcClockPeriod << std::endl;

  // print micromegas geometry
  load_nodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluatorLight_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  const auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  // cleanup output
  if( m_container ) m_container->Reset();

  fill_g4particle_map();
  evaluate_tracks();

  // clear internal maps
  m_g4hit_map.clear();
  m_g4particle_map.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
float TrackingEvaluatorLight_hp::get_dedx( TrackSeed* seed ) const
{
  if( !(seed && m_cluster_map && m_tpc_geom_container))
  { return 0; }

  // get clusters associated to track
  std::vector<float> dedxlist;
  for( const auto& cluster_key:get_cluster_keys(seed) )
  {
    // keep only TPC clusters
    if( TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::tpcId )
    { continue; }

    // cluster
    auto cluster = m_cluster_map->findCluster( cluster_key );
    if( !cluster )
    { continue; }

    // adc
    const auto adc = cluster->getAdc();

    // layer
    const auto layer = TrkrDefs::getLayer(cluster_key);

    // layer geometry
    const auto layergeom = m_tpc_geom_container->GetLayerCellGeom(layer);

    // layer thickness
    const auto thickness = layergeom->get_thickness();

    // layer radius
    const auto r = layergeom->get_radius();

    // rphi angle
    const auto alpha = std::abs(r*seed->get_qOverR()/2);
    const auto alphacorr = std::cos(alpha);

    // rz angle
    const auto beta = std::atan(seed->get_slope());
    const auto betacorr = std::cos(beta);

    // dedx
    dedxlist.emplace_back( adc*alphacorr*betacorr/thickness );
  }

  /*
   * question: would it not be more accurate to sum all the de,
   * all the dx and make the division in the end ?
   */

  // check dedx lists
  return dedxlist.empty() ? 0:std::accumulate( dedxlist.begin(), dedxlist.end(), 0. )/dedxlist.size();

}

//_____________________________________________________________________
float TrackingEvaluatorLight_hp::get_truth_dedx( TrackSeed* seed, int trkid ) const
{
  if( !(seed && m_cluster_map))
  { return 0; }

  // get clusters associated to track
  std::vector<float> dedxlist;
  for( const auto& cluster_key:get_cluster_keys(seed) )
  {
    // keep only TPC clusters
    if( TrkrDefs::getTrkrId(cluster_key) != TrkrDefs::tpcId )
    { continue; }

    // get associated g4hits
    for( const auto& hit:find_g4hits( cluster_key ) )
    {
      // check mc track id
      if( hit->get_trkid() != trkid ) continue;

      // calculate dedx and store
      const auto de = 1e6*hit->get_eion();
      const auto dx = std::sqrt(
        square(hit->get_x(1)-hit->get_x(0))+
        square(hit->get_y(1)-hit->get_y(0))+
        square(hit->get_z(1)-hit->get_z(0)));
      dedxlist.emplace_back(de/dx);
    }
  }

  const auto average = dedxlist.empty() ? 0:std::accumulate( dedxlist.begin(), dedxlist.end(), 0. )/dedxlist.size();

  // check dedx lists
  return average;
}

//_____________________________________________________________________
int TrackingEvaluatorLight_hp::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int TrackingEvaluatorLight_hp::load_nodes( PHCompositeNode* topNode )
{

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_tGeometry );

  // track map
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if( !m_cluster_map )
  { m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER"); }

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  // cluster hit association map
  m_hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");

  // local container
  m_container = findNode::getClass<Container>(topNode, "TrackingEvaluatorLight_hp::Container");

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  // g4hits
  m_g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  m_g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  m_g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  // tpc geometry
  m_tpc_geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  assert( m_tpc_geom_container );

  // micromegas geometry
  m_micromegas_geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL" );

  // load calocluster containers
  m_rawclustercontainermap.clear();
  for(const auto& [calo_layer, calo_name]:m_calo_names)
  {
    for( const std::string basename:{"CLUSTERINFO_", "CLUSTER_"} )
    {
      std::string clusterNodeName = basename + calo_name;
      auto clusterContainer = findNode::getClass<RawClusterContainer>(topNode, clusterNodeName.c_str());
      if( clusterContainer )
      {
        m_rawclustercontainermap.emplace( calo_layer, clusterContainer );
        break;
      }
    }
  }

  // also try loading HCAL topological clusters
  {
    std::string clusterNodeName = "TOPOCLUSTER_HCAL";
    auto clusterContainer = findNode::getClass<RawClusterContainer>(topNode, clusterNodeName.c_str());
    if( clusterContainer )
    {
      m_rawclustercontainermap.emplace( TOPO_HCAL, clusterContainer );
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
void TrackingEvaluatorLight_hp::evaluate_tracks()
{
  if( !( m_track_map && m_cluster_map && m_container ) ) return;

  // clear array
  m_container->clearTracks();

  for( const auto& [track_id,track]:*m_track_map )
  {
    // create track information
    auto track_struct = create_track( track );

    // truth information
    const auto [id,contributors] = get_max_contributor( track );
    track_struct._contributors = contributors;

    if( Verbosity() )
    {
      std::cout << "TrackingEvaluatorLight_hp::evaluate_tracks -"
        << " id: " << id
        << " contributors: " << contributors
        << std::endl;
    }

    // get associated particle and store relevant information
    if( m_g4truthinfo )
    {
      const auto particle = m_g4truthinfo->GetParticle(id);
      track_struct._embed = get_embed(particle);
      ::add_truth_information(track_struct, particle);
    }

    // get mask
    {
      const auto iter = m_g4particle_map.find( id );
      if( iter != m_g4particle_map.end() ) track_struct._truth_mask = iter->second;
      else if( Verbosity() )
      { std::cout << "TrackingEvaluatorLight_hp::evaluate_tracks - could not get mask for particle " << id << std::endl; }
    }

    // track crossing
    const auto crossing = track->get_crossing();
    if(crossing == SHRT_MAX)
    {
      std::cout << "TrackingEvaluatorLight_hp::evaluate_tracks - invalid crossing, track ignored." << std::endl;
      continue;
    }

    track_struct._crossing = crossing;

    // dedx
    track_struct._dedx = get_dedx( track->get_tpc_seed() );
    track_struct._truth_dedx = get_truth_dedx(track->get_tpc_seed(),id);

    // add matching calorimeter clusters
    {

      // emcal extrapolation
      if( const auto result = find_calo_cluster_emcal(track) )
      { track_struct._calo_clusters.push_back( result.value() ); }

      // calorimeter extrapolation
      for(const auto& calo_layer:std::initializer_list<int>{SvtxTrack::HCALIN, SvtxTrack::HCALOUT, TOPO_HCAL})
      {
        if( const auto result = find_calo_cluster_hcal(calo_layer,track) )
        { track_struct._calo_clusters.push_back( result.value() ); }
      }
    }

  }
}

//_____________________________________________________________________
TrackingEvaluatorLight_hp::G4HitSet TrackingEvaluatorLight_hp::find_g4hits( TrkrDefs::cluskey cluster_key ) const
{

  // check maps
  if( !( m_cluster_hit_map && m_hit_truth_map ) ) return G4HitSet();

  // check if in map
  auto map_iter = m_g4hit_map.lower_bound( cluster_key );
  if( map_iter != m_g4hit_map.end() && cluster_key == map_iter->first )
  { return map_iter->second; }

  // find hitset associated to cluster
  G4HitSet out;
  const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
  for(const auto& [first,hit_key]:range_adaptor( m_cluster_hit_map->getHits(cluster_key)))
  {

    // store hits to g4hit associations
    TrkrHitTruthAssoc::MMap g4hit_map;
    m_hit_truth_map->getG4Hits( hitset_key, hit_key, g4hit_map );

    // find corresponding g4 hist
    for( auto truth_iter = g4hit_map.begin(); truth_iter != g4hit_map.end(); ++truth_iter )
    {

      const auto g4hit_key = truth_iter->second.second;
      PHG4Hit* g4hit = nullptr;

      switch( TrkrDefs::getTrkrId( hitset_key ) )
      {
        case TrkrDefs::mvtxId:
        if( m_g4hits_mvtx ) g4hit = m_g4hits_mvtx->findHit( g4hit_key );
        break;

        case TrkrDefs::inttId:
        if( m_g4hits_intt ) g4hit = m_g4hits_intt->findHit( g4hit_key );
        break;

        case TrkrDefs::tpcId:
        if( m_g4hits_tpc ) g4hit = m_g4hits_tpc->findHit( g4hit_key );
        break;

        case TrkrDefs::micromegasId:
        if( m_g4hits_micromegas ) g4hit = m_g4hits_micromegas->findHit( g4hit_key );
        break;

        default: break;
      }

      if( g4hit ) out.insert( g4hit );
      else std::cout << "TrackingEvaluatorLight_hp::find_g4hits - g4hit not found " << g4hit_key << std::endl;

    }
  }

  // insert in map and return
  return m_g4hit_map.insert( map_iter, std::make_pair( cluster_key, std::move( out ) ) )->second;

}

//_____________________________________________________________________
std::pair<int,int> TrackingEvaluatorLight_hp::get_max_contributor( SvtxTrack* track ) const
{
  if(!(m_track_map && m_cluster_map && m_g4truthinfo)) return {0,0};

  // maps MC track id and number of matching g4hits
  using IdMap = std::map<int,int>;
  IdMap contributor_map;

  const auto cluster_keys = get_cluster_keys( track );

  if( Verbosity() )
  { std::cout << "TrackingEvaluatorLight_hp::get_max_contributor - clusters: " << cluster_keys.size() << std::endl; }

  // loop over clusters
  for( const auto& cluster_key:cluster_keys )
  {
    for( const auto& hit:find_g4hits( cluster_key ) )
    {
      const int trkid = hit->get_trkid();
      auto iter = contributor_map.lower_bound( trkid );
      if( iter == contributor_map.end() || iter->first != trkid )
      {
        contributor_map.insert(iter, std::make_pair(trkid,1));
      } else ++iter->second;
    }
  }

  if( contributor_map.empty() ) return {0,0};
  else return *std::max_element(
    contributor_map.cbegin(), contributor_map.cend(),
    []( const IdMap::value_type& first, const IdMap::value_type& second )
    { return first.second < second.second; } );

}

//_____________________________________________________________________
int TrackingEvaluatorLight_hp::get_embed( PHG4Particle* particle ) const
{ return (m_g4truthinfo && particle) ? m_g4truthinfo->isEmbeded( particle->get_primary_id() ):0; }

//_____________________________________________________________________
void TrackingEvaluatorLight_hp::fill_g4particle_map()
{
  m_g4particle_map.clear();

  // update all particle's masks for g4hits in TPC, intt and mvtx
  for( const auto& container: {m_g4hits_tpc, m_g4hits_intt, m_g4hits_mvtx} )
  {
    if( !container ) continue;

    // loop over hits
    const auto range = container->getHits();
    for( auto iter = range.first; iter != range.second; ++iter )
    {

      // get g4hit, track and layer
      const auto& g4hit = iter->second;
      const auto trkid = g4hit->get_trkid();
      const auto layer = g4hit->get_layer();

      // update relevant mask
      const auto map_iter = m_g4particle_map.lower_bound( trkid );
      if( map_iter != m_g4particle_map.end() && map_iter->first == trkid )
      {
        map_iter->second |= (1LL<<layer);
      } else {
        m_g4particle_map.insert( map_iter, std::make_pair( trkid, 1LL<<layer ) );
      }
    }
  }

  // special treatment for micromegas because one must check that the hits actually fires an existing tile
  if( m_g4hits_micromegas && m_micromegas_geom_container)
  {

    // loop over hits
    const auto range = m_g4hits_micromegas->getHits();
    for( auto iter = range.first; iter != range.second; ++iter )
    {

      // get g4hit, track and layer
      const auto& g4hit = iter->second;
      const auto trkid = g4hit->get_trkid();
      const auto layer = g4hit->get_layer();
      const auto tileid = g4hit->get_property_int( PHG4Hit::prop_index_i );
      if( tileid < 0 ) continue;

      // get relevant micromegas geometry
      const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(m_micromegas_geom_container->GetLayerGeom(layer));
      assert( layergeom );

      // get world coordinates
      TVector3 world( g4hit->get_avg_x(), g4hit->get_avg_y(), g4hit->get_avg_z() );

      // update relevant mask
      const auto map_iter = m_g4particle_map.lower_bound( trkid );
      if( map_iter != m_g4particle_map.end() && map_iter->first == trkid )
      {
        map_iter->second |= (1LL<<layer);
      } else {
        m_g4particle_map.insert( map_iter, std::make_pair( trkid, 1LL<<layer ) );
      }

    }

  }

}

//_____________________________________________________________
std::optional<TrackingEvaluatorLight_hp::CaloClusterStruct> TrackingEvaluatorLight_hp::find_calo_cluster_emcal( SvtxTrack* track ) const
{

  // select EMCAL layer
  const auto calo_layer = SvtxTrack::CEMC;
  if( m_rawclustercontainermap.find(calo_layer) == m_rawclustercontainermap.end() ) { return {}; }

  const auto container( m_rawclustercontainermap.at(calo_layer) );
  if( !container ) return {};

  // cut on cluster energy
  const auto min_energy = m_calo_min_energy.at(calo_layer);

  // calo cluster struct
  TrackingEvaluatorLight_hp::CaloClusterStruct calo_cluster_struct;
  float dmin = -1;

  for( const auto& [key,cluster]:range_adaptor( container->getClusters()))
  {

    if( cluster->get_energy() < min_energy ) continue;

    // cluster r
    const auto cluster_r = get_r(cluster->get_x(),cluster->get_y());

    // find matching state vector
    float dr_min = -1;
    auto state_iter = track->begin_states();
    for( auto iter = state_iter; iter != track->end_states(); ++iter )
    {
      const auto dr = std::abs( cluster_r - get_r( iter->second->get_x(), iter->second->get_y() ) );
      if( dr_min < 0 || dr < dr_min )
      {
        state_iter = iter;
        dr_min = dr;
      }
    }

    // no state found
    if( dr_min < 0 ) continue;

    // calculate distance between track state and cluster
    const auto& state = state_iter->second;

    // extrapolate track state to same r as cluster and calculate distance
    // need to extrapolate to the right r
    const auto trk_r = get_r( state->get_x(), state->get_y() );
    const auto dr = cluster_r - trk_r;

    const auto trk_drdt = (state->get_x()*state->get_px() + state->get_y()*state->get_py())/trk_r;
    const auto trk_dxdr = state->get_px()/trk_drdt;
    const auto trk_dydr = state->get_py()/trk_drdt;
    const auto trk_dzdr = state->get_pz()/trk_drdt;

    const auto trk_x = state->get_x() + dr*trk_dxdr;
    const auto trk_y = state->get_y() + dr*trk_dydr;
    const auto trk_z = state->get_z() + dr*trk_dzdr;

    const double d = square(trk_x - cluster->get_x()) + square(trk_y - cluster->get_y()) + square(trk_z- cluster->get_z());
    if( dmin < 0 || d < dmin )
    {
      dmin = d;
      calo_cluster_struct = create_calo_cluster(calo_layer, cluster);

      // add track information
      calo_cluster_struct._trk_x = trk_x;
      calo_cluster_struct._trk_y = trk_y;
      calo_cluster_struct._trk_z = trk_z;
      calo_cluster_struct._trk_r = get_r( trk_x, trk_y );
      calo_cluster_struct._trk_phi = std::atan2( trk_y, trk_x );

      calo_cluster_struct._trk_dr = dr_min;
      calo_cluster_struct._trk_eta = get_eta(calo_cluster_struct._trk_x,calo_cluster_struct._trk_y,calo_cluster_struct._trk_z);
    }
  }
  return dmin < 0 ? std::nullopt : std::optional(calo_cluster_struct);
}

//_____________________________________________________________
std::optional<TrackingEvaluatorLight_hp::CaloClusterStruct> TrackingEvaluatorLight_hp::find_calo_cluster_hcal( int calo_layer, SvtxTrack* track ) const
{

  std::cout << "TrackingEvaluatorLight_hp::find_calo_cluster_hcal- layer: " << calo_layer << std::endl;

  // check layer
  if( m_rawclustercontainermap.find(calo_layer) == m_rawclustercontainermap.end() ) { return {}; }

  // get container
  const auto container( m_rawclustercontainermap.at(calo_layer) );
  if( !container ) return {};

  // cut on cluster energy
  const auto min_energy = m_calo_min_energy.at(calo_layer);

  // calo cluster struct
  TrackingEvaluatorLight_hp::CaloClusterStruct calo_cluster_struct;
  float dmin = -1;

  for( const auto& [key,cluster]:range_adaptor( container->getClusters()))
  {

    if( cluster->get_energy() < min_energy ) continue;

    // cluster r
    const auto cluster_r = get_r(cluster->get_x(),cluster->get_y());

    // find iterators at smaller and larger radius than current cluster
    bool found = false;
    auto state_iter_before = track->begin_states();
    auto state_iter_after = track->begin_states();

    for( auto iter = state_iter_before; iter != track->end_states(); ++iter )
    {
      const auto& track_state( iter->second );

      // get track state radius and compare
      const auto state_r =  get_r( track_state->get_x(), track_state->get_y() );
      if( state_r < cluster_r )
      {
        state_iter_before = iter;
      } else if( state_r > cluster_r ) {
        state_iter_after = iter;
        found = true;
        break;
      }
    }

    // no state found
    if( !found ) continue;

    // calculate distance between track state and cluster
    const auto& state_before = state_iter_before->second;
    const auto& state_after = state_iter_after->second;

    const auto state_r_before = get_r( state_before->get_x(), state_before->get_y() );
    const auto state_r_after = get_r( state_after->get_x(), state_after->get_y() );

    const TVector3 position_before( state_before->get_x(),  state_before->get_y(),  state_before->get_z() );
    const TVector3 position_after( state_after->get_x(),  state_after->get_y(),  state_after->get_z() );
    const TVector3 position_cluster( cluster->get_x(),  cluster->get_y(),  cluster->get_z() );

    // find point closest to the cluster along the track
    const TVector3 track_vector = position_after-position_before;
    const TVector3 clus_vector = position_cluster-position_before;
    const auto lambda = clus_vector.Dot(track_vector)/track_vector.Mag2();

    std::cout << "TrackingEvaluatorLight_hp::evaluate_tracks -"
      << " calo_layer: " << calo_layer
      << " cluster_r: " << cluster_r
      << " state_r_before: " << state_r_before
      << " state_r_after: " << state_r_after
      << " lambda: " << lambda
      << std::endl;

    const TVector3 position_closest = position_before*(1.-lambda) + position_after*lambda;
    const auto d = (position_closest-position_cluster).Mag2();
    if( dmin < 0 || d < dmin )
    {
      dmin = d;
      calo_cluster_struct = create_calo_cluster(calo_layer, cluster);

      // add track information
      calo_cluster_struct._trk_x = position_closest.x();
      calo_cluster_struct._trk_y = position_closest.y();
      calo_cluster_struct._trk_z = position_closest.z();
      calo_cluster_struct._trk_r = get_r( calo_cluster_struct._trk_x, calo_cluster_struct._trk_y );
      calo_cluster_struct._trk_phi = std::atan2( calo_cluster_struct._trk_y, calo_cluster_struct._trk_x );
      calo_cluster_struct._trk_eta = get_eta(calo_cluster_struct._trk_x,calo_cluster_struct._trk_y,calo_cluster_struct._trk_z);
      calo_cluster_struct._trk_dr = 0;
    }
  }
  return dmin < 0 ? std::nullopt : std::optional(calo_cluster_struct);
}
