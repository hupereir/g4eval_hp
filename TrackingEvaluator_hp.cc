#include "TrackingEvaluator_hp.h"

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

  //! get cluster keys from a given track and detector
  template <int type>
  std::vector<TrkrDefs::cluskey> get_detector_cluster_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;
    for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (seed)
      {
        std::copy_if(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out),
          [](const TrkrDefs::cluskey& key)
        { return TrkrDefs::getTrkrId(key) == type;} );
      }
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
  TrackingEvaluator_hp::CaloClusterStruct create_calo_cluster( int layer, RawCluster* cluster )
  {
    TrackingEvaluator_hp::CaloClusterStruct calo_cluster_struct;
    calo_cluster_struct._layer = layer;
    calo_cluster_struct._size = cluster->getNTowers();
    calo_cluster_struct._x = cluster->get_x();
    calo_cluster_struct._y = cluster->get_y();
    calo_cluster_struct._z = cluster->get_z();

    calo_cluster_struct._r = get_r( calo_cluster_struct._x, calo_cluster_struct._y );
    calo_cluster_struct._phi = std::atan2( calo_cluster_struct._y, calo_cluster_struct._x );

    calo_cluster_struct._e = cluster->get_energy();
    calo_cluster_struct._chisquare = cluster->get_chi2();

    // loop over towers
    for( const auto& [index, energy]:range_adaptor(cluster->get_towers()))
    {
      TrackingEvaluator_hp::TowerStruct tower;

      // get indexes
      tower._ieta = RawTowerDefs::decode_index1(index);
      tower._iphi = RawTowerDefs::decode_index2(index);
      tower._e = energy;
      calo_cluster_struct._towers.emplace_back( tower );
    }

    return calo_cluster_struct;
  }

  //! create central membrane cluster struct
  TrackingEvaluator_hp::CMClusterStruct create_cm_cluster( unsigned int /*key*/, CMFlashCluster* cluster )
  {
    TrackingEvaluator_hp::CMClusterStruct cm_cluster_struct;
    cm_cluster_struct._nclusters = cluster->getNclusters();
    cm_cluster_struct._x = cluster->getX();
    cm_cluster_struct._y = cluster->getY();
    cm_cluster_struct._z = cluster->getZ();

    cm_cluster_struct._r = std::sqrt(square(cluster->getX()) + square(cluster->getY()));
    cm_cluster_struct._phi = std::atan2(cluster->getY(), cluster->getX());

    return cm_cluster_struct;
  }

  //! create track struct from struct from svx track
  TrackingEvaluator_hp::SeedStruct create_seed( TrackSeed* seed )
  {
    TrackingEvaluator_hp::SeedStruct seedStruct;

    const auto pos = TrackSeedHelper::get_xyz(seed);

    // fill additional information
    seedStruct._x = pos.x();
    seedStruct._y = pos.y();
    seedStruct._z = pos.z();
    seedStruct._r = get_r(seedStruct._x, seedStruct._y );
    seedStruct._phi = seed->get_phi();
    seedStruct._eta = seed->get_eta();
    return seedStruct;
  }

  //! create track struct from struct from svx track
  TrackingEvaluator_hp::TrackStruct create_track( SvtxTrack* track )
  {
    TrackingEvaluator_hp::TrackStruct trackStruct;

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

  //! create track struct from struct from svx track
  TrackingEvaluator_hp::TrackPairStruct create_track_pair( SvtxTrack* first, SvtxTrack* second )
  {
    TrackingEvaluator_hp::TrackPairStruct trackpair_struct;
    trackpair_struct._charge = first->get_charge() + second->get_charge();
    trackpair_struct._px = first->get_px() + second->get_px();
    trackpair_struct._py = first->get_py() + second->get_py();
    trackpair_struct._pz = first->get_pz() + second->get_pz();
    trackpair_struct._pt = get_pt( trackpair_struct._px, trackpair_struct._py );
    trackpair_struct._p = get_p( trackpair_struct._px, trackpair_struct._py, trackpair_struct._pz );
    trackpair_struct._eta = get_eta( trackpair_struct._p, trackpair_struct._pz );

    // electron mass
    static constexpr double electronMass = 0.511e-3;
    auto firstE = std::sqrt( square( electronMass ) + square( get_p( first->get_px(), first->get_py(), first->get_pz() ) ) );
    auto secondE = std::sqrt( square( electronMass ) + square( get_p( second->get_px(), second->get_py(), second->get_pz() ) ) );
    trackpair_struct._e = firstE + secondE;
    trackpair_struct._m = std::sqrt( square( trackpair_struct._e ) - square( trackpair_struct._p ) );

    return trackpair_struct;
  }

  //! number of hits associated to cluster
  void add_cluster_size( TrackingEvaluator_hp::ClusterStruct& cluster, TrkrDefs::cluskey clus_key, TrkrClusterHitAssoc* cluster_hit_map )
  {
    if( !cluster_hit_map ) return;
    const auto range = cluster_hit_map->getHits(clus_key);

    // store full size
    cluster._size =  std::distance( range.first, range.second );

    const auto detId = TrkrDefs::getTrkrId(clus_key);
    if(detId == TrkrDefs::micromegasId)
    {

      // for micromegas the directional cluster size depends on segmentation type
      auto segmentation_type = MicromegasDefs::getSegmentationType(clus_key);
      if( segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_Z ) cluster._z_size = cluster._size;
      else cluster._phi_size = cluster._size;

    } else {

      // for other detectors, one must loop over the constituting hits
      std::set<int> phibins;
      std::set<int> zbins;
      for(const auto& [first, hit_key]:range_adaptor(range))
      {
        switch( detId )
        {
          default: break;
          case TrkrDefs::mvtxId:
          {
            phibins.insert( MvtxDefs::getRow( hit_key ) );
            zbins.insert( MvtxDefs::getCol( hit_key ) );
            break;
          }
          case TrkrDefs::inttId:
          {
            phibins.insert( InttDefs::getRow( hit_key ) );
            zbins.insert( InttDefs::getCol( hit_key ) );
            break;
          }
          case TrkrDefs::tpcId:
          {
            phibins.insert( TpcDefs::getPad( hit_key ) );
            zbins.insert( TpcDefs::getTBin( hit_key ) );
            break;
          }
        }
      }
      cluster._phi_size = phibins.size();
      cluster._z_size = zbins.size();
    }
  }


  //! hit energy for a given cluster
  void add_cluster_energy( TrackingEvaluator_hp::ClusterStruct& cluster, TrkrDefs::cluskey clus_key,
    TrkrClusterHitAssoc* cluster_hit_map,
    TrkrHitSetContainer* hitsetcontainer )
  {

    // check container
    if(!(cluster_hit_map && hitsetcontainer)) return;

    // for now this is only filled for micromegas
    const auto detId = TrkrDefs::getTrkrId(clus_key);
    if(detId != TrkrDefs::micromegasId) return;

    const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(clus_key);
    const auto hitset = hitsetcontainer->findHitSet( hitset_key );
    if( !hitset ) return;

    const auto range = cluster_hit_map->getHits(clus_key);
    cluster._energy_max = 0;
    cluster._energy_sum = 0;

    for( const auto& pair:range_adaptor(range))
    {
      const auto hit = hitset->getHit( pair.second );
      if( hit )
      {
        const auto energy = hit->getAdc();
        cluster._energy_sum += energy;
        if( energy > cluster._energy_max ) cluster._energy_max = energy;
      }
    }

  }

  // add truth information
  void add_truth_information( TrackingEvaluator_hp::TrackStruct& track, PHG4Particle* particle )
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

  // calculate intersection between line and circle
  double line_circle_intersection( const TVector3& p0, const TVector3& p1, double radius )
  {
    const double A = square(p1.x() - p0.x()) + square(p1.y() - p0.y());
    const double B = 2*p0.x()*(p1.x()-p0.x()) + 2*p0.y()*(p1.y()-p0.y());
    const double C = square(p0.x()) + square(p0.y()) - square(radius);
    const double delta = square(B)-4*A*C;
    if( delta < 0 ) return -1;

    // check first intersection
    const double tup = (-B + std::sqrt(delta))/(2*A);
    if( tup >= 0 && tup < 1 ) return tup;

    // check second intersection
    const double tdn = (-B-sqrt(delta))/(2*A);
    if( tdn >= 0 && tdn < 1 ) return tdn;

    // no valid extrapolation
    return -1;
  }

  // print to stream
  [[maybe_unused]] std::ostream& operator << (std::ostream& out, const TrackingEvaluator_hp::ClusterStruct& cluster )
  {
    out << "ClusterStruct" << std::endl;
    out << "  cluster: (" << cluster._x << "," << cluster._y << "," << cluster._z << ")" << std::endl;
    out << "  track:   (" << cluster._trk_x << "," << cluster._trk_y << "," << cluster._trk_z << ")" << std::endl;
    out << "  truth:   (" << cluster._truth_x << "," << cluster._truth_y << "," << cluster._truth_z << ")" << std::endl;
    return out;
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
void TrackingEvaluator_hp::Container::Reset()
{
  _events.clear();
  _clusters.clear();
  _calo_clusters.clear();
  _cm_clusters.clear();
  _tracks.clear();
  _track_pairs.clear();
}

//_____________________________________________________________________
TrackingEvaluator_hp::TrackingEvaluator_hp( const std::string& name ):
  SubsysReco( name),
  m_calo_min_energy( {
    {SvtxTrack::CEMC, 0.15},
    {SvtxTrack::HCALIN, 0},
    {SvtxTrack::HCALOUT, 0},
    {TOPO_HCAL, 0}
  })
{}

//_____________________________________________________________________
int TrackingEvaluator_hp::Init(PHCompositeNode* topNode )
{
  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "TrackingEvaluator_hp::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "TrackingEvaluator_hp::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  // add container to output tree
  auto newNode = new PHIODataNode<PHObject>( new Container, "TrackingEvaluator_hp::Container","PHObject");

  // overwrite split level for easier offline browsing
  newNode->SplitLevel(99);
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::InitRun(PHCompositeNode* topNode)
{

  // TPC geometry (to check ADC clock frequency)
  auto geom = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  assert(geom);

  const auto AdcClockPeriod = geom->GetFirstLayerCellGeom()->get_zstep();
  std::cout << "TrackingEvaluator_hp::Init - AdcClockPeriod: " << AdcClockPeriod << std::endl;

  if( m_flags & UseTpcClusterMover )
  {
    // initialize cluster mover
    std::cout << "TrackingEvaluator_hp::Init - initializing TpcClusterMover" << std::endl;
    m_tpcClusterMover.initialize_geometry(geom);
  }

  // print micromegas geometry
  load_nodes(topNode);
  print_micromegas_geometry();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TrackingEvaluator_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  const auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  // cleanup output
  if( m_container ) m_container->Reset();

  if(m_flags&PrintClusters) print_clusters();
  if(m_flags&PrintSeeds) print_seeds(topNode);
  if(m_flags&PrintTracks) print_tracks();
  if(m_flags&EvalEvent) evaluate_event();
  if(m_flags&EvalClusters) evaluate_clusters();
  if(m_flags&EvalCaloClusters) evaluate_calo_clusters();
  if(m_flags&EvalCMClusters) evaluate_cm_clusters();
  if(m_flags&EvalTracks)
  {
    fill_g4particle_map();
    evaluate_tracks();
  }
  if(m_flags&EvalTrackPairs) evaluate_track_pairs();

  // clear maps
  m_g4hit_map.clear();
  m_g4particle_map.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
float TrackingEvaluator_hp::get_dedx( TrackSeed* seed ) const
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
float TrackingEvaluator_hp::get_truth_dedx( TrackSeed* seed, int trkid ) const
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
int TrackingEvaluator_hp::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int TrackingEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_tGeometry );

  // gl1 raw hit
  m_gl1rawhit = findNode::getClass<Gl1RawHit>(topNode,"GL1RAWHIT");

  // track map
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);

  // load raw hits container
  m_micromegas_rawhitcontainer = findNode::getClass<MicromegasRawHitContainer>(topNode, "MICROMEGASRAWHIT");

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if( !m_cluster_map )
  { m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER"); }

  // central membrane clusters
  m_cm_cluster_map = findNode::getClass<CMFlashClusterContainer>(topNode, "CORRECTED_CM_CLUSTER" );

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  // cluster hit association map
  m_hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");

  // local container
  m_container = findNode::getClass<Container>(topNode, "TrackingEvaluator_hp::Container");

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

  // tpc distortion corrections
  m_globalPositionWrapper.loadNodes(topNode);

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
void TrackingEvaluator_hp::evaluate_event()
{
  if(!m_container) return;

  // TPC threshold (adc)
  // from PHG4TpcDigitizer::InitRun - ADCThreshold * ADCNoiseConversionGain*1024/2200
  static constexpr float tpc_threshold = 78.4645;

  // create event struct
  EventStruct event;

  // GTM bco
  if( m_gl1rawhit )
  {
    event._gtm_bco = m_gl1rawhit->get_bco()&0xFFFFFFFFFF;

    if( Verbosity() )
    {
      std::cout << "TrackingEvaluator_hp::evaluate_event - bco: "
        << event._gtm_bco
        << std::endl;
    }
  }

  // number of micromegas raw hits
  if( m_micromegas_rawhitcontainer )
  { event._nrawhits_micromegas = m_micromegas_rawhitcontainer->get_nhits(); }

  if( m_hitsetcontainer )
  {
    // fill hit related information
    for(const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
    {
      const auto trkrId = TrkrDefs::getTrkrId(hitsetkey);
      const auto layer = TrkrDefs::getLayer(hitsetkey);
      assert(layer<EventStruct::max_layer);

      const auto hit_range = hitset->getHits();

      // nhits per layer
      if( trkrId == TrkrDefs::tpcId )
      {
        const auto accepted = std::count_if( hit_range.first, hit_range.second,
          []( const TrkrHitSet::Map::value_type pair )
          { return pair.second->getAdc() >= tpc_threshold; });
        event._nhits[layer] += accepted;
      } else {
        event._nhits[layer] += std::distance( hit_range.first, hit_range.second );
      }

      // nhits per layer, with no threshold cut for the TPC
      event._nhits_raw[layer] += std::distance( hit_range.first, hit_range.second );

      // also count corresponding clusters
      if(m_cluster_map)
      {

        // fill cluster related information
        for( const auto[key,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)) )
        {
          const auto trkrid = TrkrDefs::getTrkrId(key);
          switch( trkrid )
          {
            case TrkrDefs::mvtxId: ++event._nclusters_mvtx; break;
            case TrkrDefs::inttId: ++event._nclusters_intt; break;
            case TrkrDefs::tpcId: ++event._nclusters_tpc; break;
            case TrkrDefs::micromegasId: ++event._nclusters_micromegas; break;
          }

          const auto layer = static_cast<size_t>(TrkrDefs::getLayer(key));
          assert(layer<EventStruct::max_layer);
          ++event._nclusters[layer];
        }
      }
    }
  }

  auto count_g4hits = []( PHG4HitContainer* container )
  {  return container ? container->size():0; };

  event._ng4hits_mvtx = count_g4hits( m_g4hits_mvtx );
  event._ng4hits_intt = count_g4hits( m_g4hits_intt );
  event._ng4hits_tpc = count_g4hits( m_g4hits_tpc );
  event._ng4hits_micromegas = count_g4hits( m_g4hits_micromegas );

  // number of cm clusters
  if( m_cm_cluster_map )
  {
    const auto range = m_cm_cluster_map->getClusters();
    event._ncmclusters = std::distance( range.first, range.second );
  }

  // store
  m_container->addEvent(event);
}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_clusters()
{

  if(!(m_cluster_map&&m_hitsetcontainer&&m_container)) return;

  // clear array
  m_container->clearClusters();

  // print total hits in TPC
  if( false )
  {
    int n_tpc_hits = 0;
    int n_tpc_clusters = 0;
    for( const auto& [hitsetkey,hitset]:range_adaptor( m_hitsetcontainer->getHitSets( TrkrDefs::tpcId )))
    {
      n_tpc_hits += hitset->size();
      const auto range = m_cluster_map->getClusters(hitsetkey);
      n_tpc_clusters += std::distance( range.first, range.second );
    }

    std::cout << "TrackingEvaluator_hp::evaluate_clusters - n_tpc_hits: " << n_tpc_hits << std::endl;
    std::cout << "TrackingEvaluator_hp::evaluate_clusters - n_tpc_clusters: " << n_tpc_clusters << std::endl;
  }

  // first loop over hitsets
  for( const auto& [hitsetkey,hitset]:range_adaptor(m_hitsetcontainer->getHitSets()))
  {
    for( const auto& [key,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)))
    {
      // create cluster structure
      auto cluster_struct = create_cluster( key, cluster );
      add_cluster_size( cluster_struct, key, m_cluster_hit_map );
      add_cluster_energy( cluster_struct, key, m_cluster_hit_map, m_hitsetcontainer );

      // truth information
      const auto g4hits = find_g4hits( key );
      const bool is_micromegas( TrkrDefs::getTrkrId(key) == TrkrDefs::micromegasId );
      if( is_micromegas ) add_truth_information_micromegas( cluster_struct, g4hits );
      else add_truth_information( cluster_struct, g4hits );

      // add in array
      m_container->addCluster( cluster_struct );
    }
  }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_calo_clusters()
{
  if(!m_container) return;
  m_container->clearCaloClusters();

  // loop over calo types
  for(const auto& [calo_layer, clusterContainer]:m_rawclustercontainermap)
  {
    const auto clusters = clusterContainer->getClustersMap();
    for( const auto& [key,calo_cluster]:clusters )
    { m_container->addCaloCluster( create_calo_cluster( calo_layer, calo_cluster ) ); }

  }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_cm_clusters()
{
  if(!(m_cm_cluster_map&&m_container)) return;

  // clear
  m_container->clearCMClusters();

  // loop over cm clusters
  for( const auto& [key,cluster]:range_adaptor(m_cm_cluster_map->getClusters()))
  {
    auto cluster_struct = create_cm_cluster( key, cluster );
    m_container->addCMCluster( cluster_struct );
  }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_tracks()
{
  if( !( m_track_map && m_cluster_map && m_container ) ) return;

  // clear array
  m_container->clearTracks();

  for( const auto& [track_id,track]:*m_track_map )
  {
    // create track information
    auto track_struct = create_track( track );

    {
      // silicon seed information
      auto si_seed = track->get_silicon_seed();
      if(si_seed)
      {
        track_struct._has_si_seed = true;
        track_struct._si_seed = create_seed(si_seed);
      }
    }

    {
      // tpc seed information
      auto tpc_seed = track->get_tpc_seed();
      if(tpc_seed)
      {
        const auto tpc_seed_struct = create_seed(tpc_seed);

//         // create extrapolated tpc seed
//         auto tpc_seed_extrap = tpc_seed_struct;
//         if( track_struct._has_si_seed )
//         {
//           const auto dr = track_struct._si_seed._r - tpc_seed_struct._r;
// //           const auto seed_drdt = (tpc_seed_struct._x*tpc_seed->get_px() + tpc_seed_struct._y*tpc_seed->get_py())/tpc_seed_struct._r;
// //           const auto seed_dxdr = tpc_seed->get_px()/seed_drdt;
// //           const auto seed_dydr = tpc_seed->get_py()/seed_drdt;
// //           const auto seed_dzdr = tpc_seed->get_pz()/seed_drdt;
// //
// //           tpc_seed_extrap._x = tpc_seed_struct._x + dr*seed_dxdr;
// //           tpc_seed_extrap._y = tpc_seed_struct._y + dr*seed_dydr;
// //           tpc_seed_extrap._z = tpc_seed_struct._z + dr*seed_dzdr;
// //           tpc_seed_extrap._r = get_r( tpc_seed_extrap._x, tpc_seed_extrap._y);
// //           tpc_seed_extrap._z = tpc_seed_struct._z + dr/std::tan(tpc_seed->get_theta());
//           tpc_seed_extrap._z = tpc_seed_struct._z - dr/std::tan(tpc_seed->get_theta());
//         }

        // assign
        track_struct._has_tpc_seed = true;
        track_struct._tpc_seed = std::move(tpc_seed_struct);
//         track_struct._tpc_seed_extrap = std::move(tpc_seed_extrap);
      }
    }

    // truth information
    const auto [id,contributors] = get_max_contributor( track );
    track_struct._contributors = contributors;

    if( Verbosity() )
    {
      std::cout << "TrackingEvaluator_hp::evaluate_tracks -"
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
      { std::cout << "TrackingEvaluator_hp::evaluate_tracks - could not get mask for particle " << id << std::endl; }
    }

    // track crossing
    const auto crossing = track->get_crossing();
    if(crossing == SHRT_MAX)
    {
      std::cout << "TrackingEvaluator_hp::evaluate_tracks - invalid crossing, track ignored." << std::endl;
      continue;
    }

    track_struct._crossing = crossing;

    // dedx
    track_struct._dedx = get_dedx( track->get_tpc_seed() );
    track_struct._truth_dedx = get_truth_dedx(track->get_tpc_seed(),id);

    // store cluster and keys
    using cluster_map_t = std::map<TrkrDefs::cluskey,ClusterStruct>;
    cluster_map_t cluster_map;

    // loop over clusters
    for( const auto& cluster_key:get_cluster_keys( track ) )
    {
      auto cluster = m_cluster_map->findCluster( cluster_key );
      if( !cluster )
      {
        std::cout << "TrackingEvaluator_hp::evaluate_tracks -"
          << " unable to find cluster for key " << cluster_key
          << " detector: " << (int) TrkrDefs::getTrkrId( cluster_key )
          << std::endl;
        continue;
      }

      // create new cluster struct
      auto cluster_struct = create_cluster( cluster_key, cluster, crossing );
      add_cluster_size( cluster_struct, cluster_key, m_cluster_hit_map );
      add_cluster_energy( cluster_struct, cluster_key, m_cluster_hit_map, m_hitsetcontainer );

      cluster_map.emplace(cluster_key, cluster_struct);
    }

    if( m_flags & UseTpcClusterMover )
    {
      // move back clusters to surface using TPC Cluster mover
      using position_pair_t = std::pair<TrkrDefs::cluskey, Acts::Vector3>;
      using position_list_t = std::vector<position_pair_t>;
      position_list_t position_list_raw;
      std::transform( cluster_map.begin(), cluster_map.end(), std::back_inserter( position_list_raw ),
        []( const cluster_map_t::value_type& cluster_pair )
        { return std::make_pair( cluster_pair.first, Acts::Vector3{cluster_pair.second._x,cluster_pair.second._y,cluster_pair.second._z} ); } );

      auto position_list_moved = m_tpcClusterMover.processTrack(position_list_raw);

      // reassign to original cluster map
      for( const auto& [ckey,position]:position_list_moved )
      {
        auto&& cluster_struct = cluster_map.at(ckey);
        cluster_struct._x = position.x();
        cluster_struct._y = position.y();
        cluster_struct._z = position.z();
        cluster_struct._r = get_r( cluster_struct._x, cluster_struct._y );
        cluster_struct._phi = std::atan2( cluster_struct._y, cluster_struct._x );
      }
    }

    // truth information
    for( auto&& [cluster_key, cluster_struct]:cluster_map )
    {
      const auto g4hits = find_g4hits( cluster_key );
      const bool is_micromegas( TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::micromegasId );
      if( is_micromegas ) add_truth_information_micromegas( cluster_struct, g4hits );
      else add_truth_information( cluster_struct, g4hits );

      {
        // check g4hits associated particle id and update track _correct_mask
        const auto iter = std::find_if( g4hits.begin(), g4hits.end(), [id]( const auto& g4hit ){ return g4hit->get_trkid() == id; } );
        if( iter != g4hits.end() ) track_struct._correct_mask |= (1LL<<(*iter)->get_layer());

        // get track id energy as max contributor
        using eion_map_t = std::map<int, float>;
        using eion_pair_t = std::pair<int, float>;
        eion_map_t eion_map;
        for( const auto& g4hit:g4hits )
        {
          auto iter = eion_map.find( g4hit->get_trkid() );
          if( iter == eion_map.end() ) eion_map.insert( std::make_pair( g4hit->get_trkid(), g4hit->get_eion() ) );
          else iter->second += g4hit->get_eion();
        }

        if( !eion_map.empty() )
        {
          const auto max_element_iter = std::max_element( eion_map.begin(), eion_map.end(), []( const eion_pair_t& first, const eion_pair_t& second )
            { return first.second < second.second; } );
          if( max_element_iter->first == id ) track_struct._correct_mask_strict |= (1LL<<(*iter)->get_layer());
        }

      }
    }

    // track states
    for( auto&& [cluster_key, cluster_struct]:cluster_map )
    {

      // find track state that match cluster
      bool found = false;
      for( auto state_iter = track->begin_states(); state_iter != track->end_states(); ++state_iter )
      {
        const auto& [pathlengh, state] = *state_iter;
        if( state->get_cluskey() == cluster_key )
        {
          // store track state in cluster struct
          const bool is_micromegas( TrkrDefs::getTrkrId(cluster_key) == TrkrDefs::micromegasId );
          if( is_micromegas ) add_trk_information_micromegas( cluster_struct, state );
          else add_trk_information( cluster_struct, state );
          found = true;

          if( Verbosity())
          {
            std::cout << "TrackingEvaluator_hp::evaluate_tracks -"
              << " track id: " << track_id
              << " state vector found for layer: " << cluster_struct._layer
              << std::endl;
          }

          break;
        }
      }

      if((!found) && Verbosity())
      {
        std::cout << "TrackingEvaluator_hp::evaluate_tracks -"
          << " track id: " << track_id
          << " no state vector found for layer: " << cluster_struct._layer
          << std::endl;
      }

      // some printout
      if( Verbosity() )
      {
        std::cout << "TrackingEvaluator_hp::evaluate_tracks -"
          << " cluster key: " <<  cluster_key
          << " _trk_rphi_error: " << cluster_struct._trk_r*cluster_struct._trk_phi_error
          << " _trk_z_error: " << cluster_struct._trk_z_error << std::endl;
      }

    }

    // add clusters to track
    for( const auto&[key,cluster_struct]:cluster_map )
    { track_struct._clusters.push_back( std::move(cluster_struct) ); }

    // add matching calorimeter clusters
    if(m_flags&EvalCaloClusters)
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

    if(track_struct._nclusters_micromegas>0 || !(m_flags&MicromegasOnly))
    {m_container->addTrack( track_struct );}

  }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::evaluate_track_pairs()
{
  if( !( m_track_map && m_container ) ) return;

  // clear array
  m_container->clearTrackPairs();
  for( auto firstIter = m_track_map->begin(); firstIter != m_track_map->end(); ++firstIter )
  {
    const auto first = firstIter->second;

    TrackStruct_small first_struct;
    fill_track_struct( first_struct, first );
    for( auto  secondIter = m_track_map->begin(); secondIter != firstIter; ++secondIter )
    {
      const auto second = secondIter->second;
      auto trackpair_struct = create_track_pair( first, second );

      // copy first track information
      trackpair_struct._tracks[0] = first_struct;

      // fill second track information
      fill_track_struct( trackpair_struct._tracks[1], second );

      // add to container
      m_container->addTrackPair( trackpair_struct );
    }
  }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_clusters() const
{

  if(!m_cluster_map) return;

  for(const auto& hitsetkey:m_cluster_map->getHitSetKeys())
  {
    // get corresponding clusters
    for(const auto& [clusterkey,cluster]:range_adaptor(m_cluster_map->getClusters(hitsetkey)))
    {

      const auto trkrId = TrkrDefs::getTrkrId( clusterkey );
      if( trkrId==TrkrDefs::mvtxId )
      { print_cluster( clusterkey, cluster ); }
    }
  }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_tracks() const
{
  if( !m_track_map ) return;
  for(const auto& [id,track]:*m_track_map)
  { print_track( track ); }
}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_cluster( TrkrDefs::cluskey ckey, TrkrCluster* cluster ) const
{
  // get detector type
  const auto trkrId = TrkrDefs::getTrkrId( ckey );

  std::cout << "TrackingEvaluator_hp::print_cluster - MVTX cluster: " << ckey
    << " position: (" << cluster->getLocalX() << ", " << cluster->getLocalY() << ")"
    << " size: " << (int) cluster->getSize()
    << " stave: " << (int) MvtxDefs::getStaveId(ckey) << " chip: " << (int)MvtxDefs::getChipId(ckey) << " strobe: " << (int)MvtxDefs::getStrobeId(ckey) << std::endl;

  // get associated hits
  if( false )
  {

    // loop over hits associated to clusters
    const auto range = m_cluster_hit_map->getHits(ckey);
    std::cout
      << "TrackingEvaluator_hp::print_cluster -"
      << " hit_count: " << std::distance( range.first, range.second )
      << std::endl;

    for(const auto& [hitsetkey,hitkey]:range_adaptor(range))
    {

      switch( trkrId )
      {
        case TrkrDefs::tpcId:
        std::cout
          << "TrackingEvaluator_hp::print_cluster -"
          << " hit: " << hitkey
          << " sector: " << TpcDefs::getSectorId( ckey )
          << " pad: " << TpcDefs::getPad( hitkey )
          << " time bin: " << TpcDefs::getTBin( hitkey )
          << std::endl;
        break;

        default: break;

      }
    }
  }

  // get associated g4 hits
  if( false )
  {

    auto g4hits = find_g4hits( ckey );
    for( const auto& g4hit:g4hits )
    {

      std::cout
        << "TrackingEvaluator_hp::print_cluster -"
        << " layer: " << g4hit->get_layer()
        << " track: " << g4hit->get_trkid()
        << " in: (" << g4hit->get_x(0) << "," << g4hit->get_y(0) << "," << g4hit->get_z(0) << ")"
        << " out: (" << g4hit->get_x(1) << "," << g4hit->get_y(1) << "," << g4hit->get_z(1) << ")"
        << " polar in: (" << get_r( g4hit->get_x(0), g4hit->get_y(0) ) << "," << std::atan2( g4hit->get_y(0), g4hit->get_x(0) ) << "," << g4hit->get_z(0) << ")"
        << " polar out: (" << get_r( g4hit->get_x(1), g4hit->get_y(1) ) << "," << std::atan2( g4hit->get_y(1), g4hit->get_x(1) ) << "," << g4hit->get_z(1) << ")"
        << std::endl;
    }

    // convert hits to list of interpolation_data_t
    interpolation_data_t::list hits;
    for( const auto& g4hit:g4hits )
    {
      const auto weight = g4hit->get_edep();
      for( int i = 0; i < 2; ++i )
      {
        const TVector3 g4hit_world(g4hit->get_x(i), g4hit->get_y(i), g4hit->get_z(i));
        const TVector3 momentum_world(g4hit->get_px(i), g4hit->get_py(i), g4hit->get_pz(i));
        hits.push_back( {.position = g4hit_world, .momentum = momentum_world, .weight = weight } );
      }
    }

    // average g4hits positions at the same radius as the cluster to get resolution
    const auto xextrap = average<&interpolation_data_t::x>(hits);
    const auto yextrap = average<&interpolation_data_t::y>(hits);
    const auto zextrap = average<&interpolation_data_t::z>(hits);

    // print interpolation
    std::cout
      << "TrackingEvaluator_hp::print_cluster -"
      << " interpolation: (" << xextrap << "," << yextrap << "," << zextrap << ")"
      << " polar: (" << get_r( xextrap, yextrap ) << "," << std::atan2( yextrap, xextrap ) << "," << zextrap << ")"
      << std::endl;

  }

  // std::cout << std::endl;

}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_track(SvtxTrack* track) const
{

  if( !track ) return;
  if( !track->get_silicon_seed() ) return;
  if( !track->get_tpc_seed() ) return;

  // print track position and momentum
  std::cout << "TrackingEvaluator_hp::print_track - id: " << track->get_id() << std::endl;
  std::cout << "TrackingEvaluator_hp::print_track - position: (" << track->get_x() << ", " << track->get_y() << ", " << track->get_z() << ")" << std::endl;
  std::cout << "TrackingEvaluator_hp::print_track - momentum: (" << track->get_px() << ", " << track->get_py() << ", " << track->get_pz() << ")" << std::endl;
  std::cout << "TrackingEvaluator_hp::print_track - clusters: " << get_cluster_keys( track ).size() << ", states: " << track->size_states() << std::endl;

  std::cout << "TrackingEvaluator_hp::print_track - silicon seed id: " << track->get_silicon_seed() << std::endl;
  std::cout << "TrackingEvaluator_hp::print_track - tpc seed id: " << track->get_tpc_seed() << std::endl;
  std::cout << " MVTX cluster keys: " << get_detector_cluster_keys<TrkrDefs::mvtxId>(track) << std::endl;

  std::cout << " INTT cluster keys: " << get_detector_cluster_keys<TrkrDefs::inttId>(track) << std::endl;
  std::cout << " TPOT cluster keys: " << get_detector_cluster_keys<TrkrDefs::micromegasId>(track) << std::endl;
  std::cout << " TPC cluster keys: " << get_detector_cluster_keys<TrkrDefs::tpcId>(track) << std::endl;

  // print MVTX cluster keys
  for( const auto ckey:get_detector_cluster_keys<TrkrDefs::mvtxId>(track) )
  {
    // get cluster from map
    if( m_cluster_map )
    {
      auto cluster = m_cluster_map->findCluster( ckey );
      std::cout << " MVTX cluster: " << ckey
        << " position: (" << cluster->getLocalX() << ", " << cluster->getLocalY() << ")"
        << " size: " << (int)cluster->getSize()
        << " layer: " << (int)TrkrDefs::getLayer(ckey)
        << " stave: " << (int) MvtxDefs::getStaveId(ckey)
        << " chip: " << (int)MvtxDefs::getChipId(ckey)
        << " strobe: " << (int)MvtxDefs::getStrobeId(ckey)
        << " index: " << (int)TrkrDefs::getClusIndex(ckey)
        << std::endl;
    } else {
      std::cout << " MVTX cluster: " << ckey << " stave: " << (int) MvtxDefs::getStaveId(ckey) << " chip: " << (int)MvtxDefs::getChipId(ckey) << " strobe: " << (int)MvtxDefs::getStrobeId(ckey) << std::endl;
    }
  }

  // also print the TPC tracks states
  if( false )
  {
    for( auto iter = track->begin_states(); iter != track->end_states(); ++iter )
    {
      const auto& [pathlenght, state] = *iter;
      std::cout << " TrackState - "
        << "ckey: " << state->get_cluskey()
        << " layer: " << (int)TrkrDefs::getLayer(state->get_cluskey())
        << " position: (" << state->get_x() << "," << state->get_y() << "," << state->get_z() << ")"
        << " momentum: (" << state->get_px() << "," << state->get_py() << "," << state->get_px() << ")"
        << std::endl;
    }
  }

  std::cout << std::endl;

}

//_____________________________________________________________________
void TrackingEvaluator_hp::print_seeds( PHCompositeNode* topNode ) const
{

  const auto svtx_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  assert( svtx_seed_container );

  const auto silicon_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  assert( silicon_seed_container );

  const auto tpc_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  assert( tpc_seed_container );

  // print seed map sizes
  std::cout << "TrackingEvaluator_hp::print_seeds - svtx_seed_container size: " << svtx_seed_container->size() << std::endl;
  std::cout << "TrackingEvaluator_hp::print_seeds - silicon_seed_container size: " << silicon_seed_container->size() << std::endl;
  std::cout << "TrackingEvaluator_hp::print_seeds - tpc_seed_container size: " << tpc_seed_container->size() << std::endl;

  // loop over tpc seeds
  for( unsigned int index = 0; const auto& seed:*tpc_seed_container )
  {
    if( !seed ) continue;
    std::cout << "TrackingEvaluator_hp::print_seeds - seed index: " << index++ << std::endl;

    // loop over cluster keys
    std::set<TrkrDefs::cluskey> cluster_keys;
    std::copy( seed->begin_cluster_keys(), seed->end_cluster_keys(), std::inserter(cluster_keys, cluster_keys.end()) );
    std::cout << "TrackingEvaluator_hp::print_seeds - keys: " << cluster_keys << std::endl;
  }

//   const auto silicon_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
//   if( silicon_seed_container )
//   {
//     // loop over seeds
//     for( unsigned int index = 0; const auto& seed:*silicon_seed_container )
//     {
//       if( !seed ) continue;
//       if( seed->get_crossing() == SHRT_MAX ) continue;
//
//       std::cout << "TrackingEvaluator_hp::print_seeds - seed index: " << index++ << " crossing: " << seed->get_crossing() << std::endl;
//
//       // loop over cluster keys
//       std::set<TrkrDefs::cluskey> cluster_keys;
//       std::copy( seed->begin_cluster_keys(), seed->end_cluster_keys(), std::inserter(cluster_keys, cluster_keys.end()) );
//       std::cout << "TrackingEvaluator_hp::print_seeds - keys: " << cluster_keys << std::endl;
//
//       for( auto key_iter = seed->begin_cluster_keys(); key_iter != seed->end_cluster_keys(); ++key_iter )
//       {
//         const auto& ckey(*key_iter);
//         const auto trkrId = TrkrDefs::getTrkrId(ckey);
//         if( trkrId == TrkrDefs::mvtxId )
//         {
//           // get cluster from map
//           if( m_cluster_map )
//           {
//             auto cluster = m_cluster_map->findCluster( ckey );
//             std::cout << " MVTX cluster: " << ckey
//               << " position: (" << cluster->getLocalX() << ", " << cluster->getLocalY() << ")"
//               << " size: " << (int) cluster->getSize()
//               << " stave: " << (int) MvtxDefs::getStaveId(ckey) << " chip: " << (int)MvtxDefs::getChipId(ckey) << " strobe: " << (int)MvtxDefs::getStrobeId(ckey) << std::endl;
//           } else {
//             std::cout << " MVTX cluster: " << ckey << " stave: " << (int) MvtxDefs::getStaveId(ckey) << " chip: " << (int)MvtxDefs::getChipId(ckey) << " strobe: " << (int)MvtxDefs::getStrobeId(ckey) << std::endl;
//           }
//
//         }
//       }
//     }
//
// //   // GL1 BCO
// //   auto gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
// //   if (gl1)
// //   {
// //     std::cout << "TrackingEvaluator_hp::print_seeds - GL1 BCO" << gl1->lValue(0, "BCO") << std::endl;
// //   }
// //
// //   // print all cluster to crossing associations
// //   const auto clustercrossingassoc = findNode::getClass<TrkrClusterCrossingAssoc>(topNode, "TRKR_CLUSTERCROSSINGASSOC");
// //   if( false )
// //   {
// //     if( clustercrossingassoc )
// //     {
// //       for( const auto&[clus_key,crossing]:range_adaptor(clustercrossingassoc->getAll()))
// //       {
// //         std::cout << "TrackingEvaluator_hp::print_seeds - [" << clus_key << "," << crossing << "]" << std::endl;
// //       }
// //     } else {
// //       std::cout << "TrackingEvaluator_hp::print_seeds - TrkrClusterCrossingAssoc not found" << std::endl;
// //     }
// //   }
// //
// //   // print MVTX strobes and BCO
// //   const auto mvtx_event_info = findNode::getClass<MvtxEventInfo>(topNode, "MVTXEVENTHEADER");
// //   if( mvtx_event_info )
// //   {
// //     std::cout << "TrackingEvaluator_hp::print_seeds - l1 bco:" << mvtx_event_info->get_L1_BCOs() << std::endl;
// //     std::cout << "TrackingEvaluator_hp::print_seeds - strobe bco:" << mvtx_event_info->get_strobe_BCOs() << std::endl;
// //   } else {
// //     std::cout << "TrackingEvaluator_hp::print_seeds - MvtxEventInfo not found" << std::endl;
// //   }
// //
// //   // loop over seeds
// //   // get INTT BCO
// //   // get MVTX strobes, starting BCO and compare. I need to put that in a histogram
// //   const auto silicon_seed_container = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
// //   if( silicon_seed_container )
// //   {
// //     // loop over seeds
// //     for( unsigned int index = 0; const auto& seed:*silicon_seed_container )
// //     {
// //       if( !seed ) continue;
// //       if( seed->get_crossing() == SHRT_MAX ) continue;
// //
// //       std::cout << "TrackingEvaluator_hp::print_seeds - seed index: " << index++ << " crossing: " << seed->get_crossing() << std::endl;
// //
// //       // store crossing
// //       std::set<int> crossings;
// //       std::set<int> crossings_map;
// //       std::set<int> strobes;
// //
// //       // loop over cluster keys
// //       for( auto key_iter = seed->begin_cluster_keys(); key_iter != seed->end_cluster_keys(); ++key_iter )
// //       {
// //         const auto& clus_key(*key_iter);
// //         const auto trkrId = TrkrDefs::getTrkrId(clus_key);
// //         switch(trkrId)
// //         {
// //           case TrkrDefs::inttId:
// //           crossings.insert(InttDefs::getTimeBucketId(clus_key));
// // //           if( clustercrossingassoc )
// // //           {
// // //             for(const auto& [key,crossing]:range_adaptor(clustercrossingassoc->getCrossings(clus_key)))
// // //             { crossings_map.insert(crossing);}
// // //           }
// //           break;
// //
// //           case TrkrDefs::mvtxId:
// //           strobes.insert(MvtxDefs::getStrobeId(clus_key));
// //           break;
// //
// //           default: break;
// //         }
// //       }
// //
// //       // print
// //       std::cout << "TrackingEvaluator_hp::print_seeds - strobes: " << strobes << " crossings: " << crossings << " from map: " << crossings_map << std::endl;
// //
// //       //
// //
// //     }
//
//   } else {
//     std::cout << "TrackingEvaluator_hp::print_seeds - SiliconTrackSeedContainer not found" << std::endl;
//   }


}

//_____________________________________________________________________
TrackingEvaluator_hp::G4HitSet TrackingEvaluator_hp::find_g4hits( TrkrDefs::cluskey cluster_key ) const
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
      else std::cout << "TrackingEvaluator_hp::find_g4hits - g4hit not found " << g4hit_key << std::endl;

    }
  }

  // insert in map and return
  return m_g4hit_map.insert( map_iter, std::make_pair( cluster_key, std::move( out ) ) )->second;

}

//_____________________________________________________________________
std::pair<int,int> TrackingEvaluator_hp::get_max_contributor( SvtxTrack* track ) const
{
  if(!(m_track_map && m_cluster_map && m_g4truthinfo)) return {0,0};

  // maps MC track id and number of matching g4hits
  using IdMap = std::map<int,int>;
  IdMap contributor_map;

  const auto cluster_keys = get_cluster_keys( track );

  if( Verbosity() )
  { std::cout << "TrackingEvaluator_hp::get_max_contributor - clusters: " << cluster_keys.size() << std::endl; }

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
int TrackingEvaluator_hp::get_embed( PHG4Particle* particle ) const
{ return (m_g4truthinfo && particle) ? m_g4truthinfo->isEmbeded( particle->get_primary_id() ):0; }

//_____________________________________________________________________
TrackingEvaluator_hp::ClusterStruct TrackingEvaluator_hp::create_cluster( TrkrDefs::cluskey key, TrkrCluster* cluster, short int crossing ) const
{
  // get global coordinates
  const auto global = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(key, cluster, crossing);

  // apply distortion corrections
  TrackingEvaluator_hp::ClusterStruct cluster_struct;
  cluster_struct._layer = TrkrDefs::getLayer(key);
  cluster_struct._x = global.x();
  cluster_struct._y = global.y();
  cluster_struct._z = global.z();
  cluster_struct._r = get_r( cluster_struct._x, cluster_struct._y );
  cluster_struct._phi = std::atan2( cluster_struct._y, cluster_struct._x );
  cluster_struct._phi_error = cluster->getRPhiError()/cluster_struct._r;
  cluster_struct._z_error = cluster->getZError();

  if(TrkrDefs::getTrkrId(key) == TrkrDefs::micromegasId)
  { cluster_struct._tileid = MicromegasDefs::getTileId(key); }

  return cluster_struct;
}

//_____________________________________________________________________
void TrackingEvaluator_hp::add_trk_information( TrackingEvaluator_hp::ClusterStruct& cluster, SvtxTrackState* state ) const
{
  // need to extrapolate to the right r
  const auto trk_r = get_r( state->get_x(), state->get_y() );
  const auto dr = cluster._r - trk_r;
  const auto trk_drdt = (state->get_x()*state->get_px() + state->get_y()*state->get_py())/trk_r;
  const auto trk_dxdr = state->get_px()/trk_drdt;
  const auto trk_dydr = state->get_py()/trk_drdt;
  const auto trk_dzdr = state->get_pz()/trk_drdt;

  // store state position
  cluster._trk_x = state->get_x() + dr*trk_dxdr;
  cluster._trk_y = state->get_y() + dr*trk_dydr;
  cluster._trk_z = state->get_z() + dr*trk_dzdr;
  cluster._trk_r = get_r( cluster._trk_x, cluster._trk_y );
  cluster._trk_phi = std::atan2( cluster._trk_y, cluster._trk_x );

  // store extrapolation delta r
  cluster._trk_dr = dr;

  /* store local momentum information */
  cluster._trk_px = state->get_px();
  cluster._trk_py = state->get_py();
  cluster._trk_pz = state->get_pz();

  /*
  store state angles in (r,phi) and (r,z) plans
  they are needed to study space charge distortions
  */
  const auto cosphi( std::cos( cluster._trk_phi ) );
  const auto sinphi( std::sin( cluster._trk_phi ) );
  const auto trk_pphi = -state->get_px()*sinphi + state->get_py()*cosphi;
  const auto trk_pr = state->get_px()*cosphi + state->get_py()*sinphi;
  const auto trk_pz = state->get_pz();
  cluster._trk_alpha = std::atan2( trk_pphi, trk_pr );
  cluster._trk_beta = std::atan2( trk_pz, trk_pr );
  cluster._trk_phi_error = state->get_phi_error();
  cluster._trk_z_error = state->get_z_error();

}

//_____________________________________________________________________
void TrackingEvaluator_hp::add_trk_information_micromegas( TrackingEvaluator_hp::ClusterStruct& cluster, SvtxTrackState* state ) const
{

  // get geometry cylinder from layer
  const auto layer = cluster._layer;
  const auto tileid = cluster._tileid;
  const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(m_micromegas_geom_container->GetLayerGeom(layer));
  assert( layergeom );

  // convert cluster position to local tile coordinates
  const TVector3 cluster_world( cluster._x, cluster._y, cluster._z );
  const auto cluster_local = layergeom->get_local_from_world_coords( tileid, m_tGeometry, cluster_world );

  // convert track position to local tile coordinates
  TVector3 track_world( state->get_x(), state->get_y(), state->get_z() );
  auto track_local = layergeom->get_local_from_world_coords( tileid, m_tGeometry, track_world );

  // convert direction to local tile coordinates
  const TVector3 direction_world( state->get_px(), state->get_py(), state->get_pz() );
  const auto direction_local = layergeom->get_local_from_world_vect( tileid, m_tGeometry, direction_world );

  // extrapolate to same local z (should be zero) as cluster
  const auto delta_z = cluster_local.z() - track_local.z();
  track_local += TVector3(
    delta_z*direction_local.x()/direction_local.z(),
    delta_z*direction_local.y()/direction_local.z(),
    delta_z);

  // convert back to global coordinates
  track_world = layergeom->get_world_from_local_coords( tileid, m_tGeometry, track_local );

  // store state position
  cluster._trk_x = track_world.x();
  cluster._trk_y = track_world.y();
  cluster._trk_z = track_world.z();
  cluster._trk_r = get_r( cluster._trk_x, cluster._trk_y );
  cluster._trk_phi = std::atan2( cluster._trk_y, cluster._trk_x );

  cluster._trk_dr = delta_z;

  /* store local momentum information */
  cluster._trk_px = state->get_px();
  cluster._trk_py = state->get_py();
  cluster._trk_pz = state->get_pz();

  /*
  store state angles in (r,phi) and (r,z) plans
  they are needed to study space charge distortions
  */
  const auto cosphi( std::cos( cluster._trk_phi ) );
  const auto sinphi( std::sin( cluster._trk_phi ) );
  const auto trk_pphi = -state->get_px()*sinphi + state->get_py()*cosphi;
  const auto trk_pr = state->get_px()*cosphi + state->get_py()*sinphi;
  const auto trk_pz = state->get_pz();
  cluster._trk_alpha = std::atan2( trk_pphi, trk_pr );
  cluster._trk_beta = std::atan2( trk_pz, trk_pr );
  cluster._trk_phi_error = state->get_phi_error();
  cluster._trk_z_error = state->get_z_error();

}

//_____________________________________________________________________
void TrackingEvaluator_hp::add_truth_information( TrackingEvaluator_hp::ClusterStruct& cluster, std::set<PHG4Hit*> g4hits ) const
{
  // store number of contributing g4hits
  cluster._truth_size = g4hits.size();

  // get layer, tpc flag and corresponding layer geometry
  const auto layer = cluster._layer;
  const bool is_tpc( layer >= 7 && layer < 55 );
  const auto layergeom = is_tpc ? m_tpc_geom_container->GetLayerCellGeom(layer):nullptr;
  const auto rin = layergeom ? layergeom->get_radius()-layergeom->get_thickness()/2:0;
  const auto rout = layergeom ? layergeom->get_radius()+layergeom->get_thickness()/2:0;

  {
    // count number of truth track ids participating to this cluster
    std::set<int> ids;
    for( const auto& g4hit:g4hits ) { ids.insert( g4hit->get_trkid() ); }
    cluster._ncontributors = ids.size();
  }

  // convert hits to list of interpolation_data_t
  interpolation_data_t::list hits;
  for( const auto& g4hit:g4hits )
  {
    interpolation_data_t::list tmp_hits;
    const auto weight = g4hit->get_edep();
    for( int i = 0; i < 2; ++i )
    {
      const TVector3 g4hit_world(g4hit->get_x(i), g4hit->get_y(i), g4hit->get_z(i));
      const TVector3 momentum_world(g4hit->get_px(i), g4hit->get_py(i), g4hit->get_pz(i));
      tmp_hits.push_back( {.position = g4hit_world, .momentum = momentum_world, .weight = weight } );
    }

    if( false && is_tpc )
    {
      // add layer boundary checks
      // ensure first hit has lowest r
      auto r0 = get_r(tmp_hits[0].x(),tmp_hits[0].y());
      auto r1 = get_r(tmp_hits[1].x(),tmp_hits[1].y());
      if( r0 > r1 )
      {
        std::swap(tmp_hits[0],tmp_hits[1]);
        std::swap(r0, r1);
      }

      // do nothing if out of bound
      if( r1 <= rin || r0 >= rout ) continue;

      // keep track of original deltar
      const auto dr_old = r1-r0;

      // clamp r0 to rin
      if( r0 < rin )
      {
        const auto t = line_circle_intersection( tmp_hits[0].position, tmp_hits[1].position, rin );
        if( t<0 ) continue;

        tmp_hits[0].position = tmp_hits[0].position*(1.-t) + tmp_hits[1].position*t;
        tmp_hits[0].momentum = tmp_hits[0].momentum*(1.-t) + tmp_hits[1].momentum*t;
        r0 = rin;
      }

      if( r1 > rout )
      {
        const auto t = line_circle_intersection( tmp_hits[0].position, tmp_hits[1].position, rout );
        if( t<0 ) continue;

        tmp_hits[1].position = tmp_hits[0].position*(1.-t) + tmp_hits[1].position*t;
        tmp_hits[1].momentum = tmp_hits[0].momentum*(1.-t) + tmp_hits[1].momentum*t;
        r1 = rout;
      }

      // update weights, only if clamping occured
      const auto dr_new = r1-r0;
      tmp_hits[0].weight *= dr_new/dr_old;
      tmp_hits[1].weight *= dr_new/dr_old;
    }

    // store in global list
    hits.push_back(std::move(tmp_hits[0]));
    hits.push_back(std::move(tmp_hits[1]));
  }

  const auto rextrap = cluster._r;

  // add truth position
  cluster._truth_x = interpolate_r<&interpolation_data_t::x>( hits, rextrap );
  cluster._truth_y = interpolate_r<&interpolation_data_t::y>( hits, rextrap );
  cluster._truth_z = interpolate_r<&interpolation_data_t::z>( hits, rextrap );

  // add truth momentum information
  cluster._truth_px = interpolate_r<&interpolation_data_t::px>( hits, rextrap );
  cluster._truth_py = interpolate_r<&interpolation_data_t::py>( hits, rextrap );
  cluster._truth_pz = interpolate_r<&interpolation_data_t::pz>( hits, rextrap );

//   // add truth position
//   cluster._truth_x = average<&interpolation_data_t::x>( hits );
//   cluster._truth_y = average<&interpolation_data_t::y>( hits );
//   cluster._truth_z = average<&interpolation_data_t::z>( hits );
//
//   // add truth momentum information
//   cluster._truth_px = average<&interpolation_data_t::px>( hits );
//   cluster._truth_py = average<&interpolation_data_t::py>( hits );
//   cluster._truth_pz = average<&interpolation_data_t::pz>( hits );

  cluster._truth_r = get_r( cluster._truth_x, cluster._truth_y );
  cluster._truth_phi = std::atan2( cluster._truth_y, cluster._truth_x );

  /*
  store state angles in (r,phi) and (r,z) plans
  they are needed to study space charge distortions
  */
  const auto cosphi( std::cos( cluster._truth_phi ) );
  const auto sinphi( std::sin( cluster._truth_phi ) );
  const auto truth_pphi = -cluster._truth_px*sinphi + cluster._truth_py*cosphi;
  const auto truth_pr = cluster._truth_px*cosphi + cluster._truth_py*sinphi;

  cluster._truth_alpha = std::atan2( truth_pphi, truth_pr );
  cluster._truth_beta = std::atan2( cluster._truth_pz, truth_pr );

}

//_____________________________________________________________________
void TrackingEvaluator_hp::add_truth_information_micromegas( TrackingEvaluator_hp::ClusterStruct& cluster, std::set<PHG4Hit*> g4hits ) const
{
  // store number of contributing g4hits
  cluster._truth_size = g4hits.size();

  {
    // count number of truth track ids participating to this cluster
    std::set<int> ids;
    for( const auto& g4hit:g4hits ) { ids.insert( g4hit->get_trkid() ); }
    cluster._ncontributors = ids.size();
  }

  const auto layer = cluster._layer;
  const auto tileid = cluster._tileid;
  const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(m_micromegas_geom_container->GetLayerGeom(layer));
  assert( layergeom );

  // convert cluster position to local tile coordinates
  const TVector3 cluster_world( cluster._x, cluster._y, cluster._z );
  const auto cluster_local = layergeom->get_local_from_world_coords( tileid, m_tGeometry, cluster_world );

  // convert hits to list of interpolation_data_t
  interpolation_data_t::list hits;
  for( const auto& g4hit:g4hits )
  {
    const auto weight = g4hit->get_edep();
    for( int i = 0; i < 2; ++i )
    {

      // convert position to local
      TVector3 g4hit_world(g4hit->get_x(i), g4hit->get_y(i), g4hit->get_z(i));
      TVector3 g4hit_local = layergeom->get_local_from_world_coords( tileid, m_tGeometry, g4hit_world );

      // convert momentum to local
      TVector3 momentum_world(g4hit->get_px(i), g4hit->get_py(i), g4hit->get_pz(i));
      TVector3 momentum_local = layergeom->get_local_from_world_vect( tileid, m_tGeometry, momentum_world );

      hits.push_back( {.position = g4hit_local, .momentum = momentum_local, .weight = weight } );
    }
  }

  // do position interpolation
  const TVector3 interpolation_local(
    average<&interpolation_data_t::x>(hits),
    average<&interpolation_data_t::y>(hits),
    average<&interpolation_data_t::z>(hits) );

  // convert back to global
  const TVector3 interpolation_world = layergeom->get_world_from_local_coords( tileid, m_tGeometry, interpolation_local );

  // do momentum interpolation
  const TVector3 momentum_local(
    average<&interpolation_data_t::px>(hits),
    average<&interpolation_data_t::py>(hits),
    average<&interpolation_data_t::pz>(hits));

  // convert back to global
  const TVector3 momentum_world = layergeom->get_world_from_local_vect( tileid, m_tGeometry, momentum_local );

  cluster._truth_x = interpolation_world.x();
  cluster._truth_y = interpolation_world.y();
  cluster._truth_z = interpolation_world.z();
  cluster._truth_r = get_r( cluster._truth_x, cluster._truth_y );
  cluster._truth_phi = std::atan2( cluster._truth_y, cluster._truth_x );

  /* add truth momentum information */
  cluster._truth_px = momentum_world.x();
  cluster._truth_py = momentum_world.y();
  cluster._truth_pz = momentum_world.z();

  /*
  store state angles in (r,phi) and (r,z) plans
  they are needed to study space charge distortions
  */
  const auto cosphi( std::cos( cluster._truth_phi ) );
  const auto sinphi( std::sin( cluster._truth_phi ) );
  const auto truth_pphi = -cluster._truth_px*sinphi + cluster._truth_py*cosphi;
  const auto truth_pr = cluster._truth_px*cosphi + cluster._truth_py*sinphi;

  cluster._truth_alpha = std::atan2( truth_pphi, truth_pr );
  cluster._truth_beta = std::atan2( cluster._truth_pz, truth_pr );

}

//_____________________________________________________________________
void TrackingEvaluator_hp::fill_g4particle_map()
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
std::optional<TrackingEvaluator_hp::CaloClusterStruct> TrackingEvaluator_hp::find_calo_cluster_emcal( SvtxTrack* track ) const
{

  // select EMCAL layer
  const auto calo_layer = SvtxTrack::CEMC;
  if( m_rawclustercontainermap.find(calo_layer) == m_rawclustercontainermap.end() ) { return {}; }

  const auto container( m_rawclustercontainermap.at(calo_layer) );
  if( !container ) return {};

  // cut on cluster energy
  const auto min_energy = m_calo_min_energy.at(calo_layer);

  // calo cluster struct
  TrackingEvaluator_hp::CaloClusterStruct calo_cluster_struct;
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
      calo_cluster_struct._trk_dphi = delta_phi( calo_cluster_struct._phi-calo_cluster_struct._trk_phi);
      calo_cluster_struct._trk_deta =
        get_eta(calo_cluster_struct._x,calo_cluster_struct._y,calo_cluster_struct._z)-
        get_eta(calo_cluster_struct._trk_x,calo_cluster_struct._trk_y,calo_cluster_struct._trk_z);
    }
  }
  return dmin < 0 ? std::nullopt : std::optional(calo_cluster_struct);
}

//_____________________________________________________________
std::optional<TrackingEvaluator_hp::CaloClusterStruct> TrackingEvaluator_hp::find_calo_cluster_hcal( int calo_layer, SvtxTrack* track ) const
{

  std::cout << "TrackingEvaluator_hp::find_calo_cluster_hcal- layer: " << calo_layer << std::endl;

  // check layer
  if( m_rawclustercontainermap.find(calo_layer) == m_rawclustercontainermap.end() ) { return {}; }

  // get container
  const auto container( m_rawclustercontainermap.at(calo_layer) );
  if( !container ) return {};

  // cut on cluster energy
  const auto min_energy = m_calo_min_energy.at(calo_layer);

  // calo cluster struct
  TrackingEvaluator_hp::CaloClusterStruct calo_cluster_struct;
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

    std::cout << "TrackingEvaluator_hp::evaluate_tracks -"
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
      calo_cluster_struct._trk_dr = 0;
    }
  }
  return dmin < 0 ? std::nullopt : std::optional(calo_cluster_struct);
}

//_____________________________________________________________________-
void TrackingEvaluator_hp::print_micromegas_geometry()
{
   // loop over layers
  const auto [begin,end] = m_micromegas_geom_container->get_begin_end();
  for( auto iter = begin; iter != end; ++iter )
  {

    // get layer geometry object
    auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(iter->second);

    // get layer
    const unsigned int layer = layergeom->get_layer();

    // loop over tiles
    const unsigned int tile_count = layergeom->get_tiles_count();
    for( unsigned int tileid = 0; tileid < tile_count; ++tileid )
    {
      // get tile center and print
      const auto center = layergeom->get_world_from_local_coords( tileid, m_tGeometry, {0,0} );
      std::cout << "TrackingEvaluator_hp::print_micromegas_geometry - layer: " << layer << " tile: " << tileid << " center: {" << center.x() << "," << center.y() << "," << center.z() << "}" << std::endl;
    }
  }
}
