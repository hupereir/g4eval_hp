#include "MicromegasTrackEvaluator_hp.h"

#include "MicromegasGeometryContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>
#include <micromegas/MicromegasMapping.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TFile.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <gsl/gsl_randist.h>

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

  //! square
  template<class T> inline constexpr T square( T x ) { return x*x; }

  //! radius
  template<class T> inline constexpr T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  //! pt
  template<class T> T get_pt( T px, T py ) { return std::sqrt( square(px) + square(py) ); }

  //! p
  template<class T> T get_p( T px, T py, T pz ) { return std::sqrt( square(px) + square(py) + square(pz) ); }

  //! eta
  template<class T> T get_eta( T p, T pz ) { return std::log( (p+pz)/(p-pz) )/2; }

  /// calculate intersection from circle to line, in 2d. return true on success
  /**
   * circle is defined as (x-xc)**2 + (y-yc)**2 = r**2
   * line is defined as nx(x-x0) + ny(y-y0) = 0
   * to solve we substitute y by y0 - nx/ny*(x-x0) in the circle equation and solve the 2nd order polynom
   * there is the extra complication that ny can be 0 (vertical line) to prevent this, we multiply all terms of the polynom by ny**2
   * and account for this special case when calculating x from y
   */
  bool circle_line_intersection(
      double r, double xc, double yc,
      double x0, double y0, double nx, double ny,
      double& xplus, double& yplus, double& xminus, double& yminus)
  {
    if (ny == 0)
    {
      // vertical lines are defined by ny=0 and x = x0
      xplus = xminus = x0;

      // calculate y accordingly
      const double delta = square(r) - square(x0 - xc);
      if (delta < 0)
      {
        return false;
      }

      const double sqdelta = std::sqrt(delta);
      yplus = yc + sqdelta;
      yminus = yc - sqdelta;
    }
    else
    {
      const double a = square(nx) + square(ny);
      const double b = -2. * (square(ny) * xc + square(nx) * x0 + nx * ny * (y0 - yc));
      const double c = square(ny) * (square(xc) - square(r)) + square(ny * (y0 - yc) + nx * x0);
      const double delta = square(b) - 4. * a * c;
      if (delta < 0)
      {
        return false;
      }

      const double sqdelta = std::sqrt(delta);
      xplus = (-b + sqdelta) / (2. * a);
      xminus = (-b - sqdelta) / (2. * a);

      yplus = y0 - (nx / ny) * (xplus - x0);
      yminus = y0 - (nx / ny) * (xminus - x0);
    }

    return true;
  }

  // line line intersection
  bool line_line_intersection(
      double m, double b,
      double x0, double y0, double nx, double ny,
      double& xplus, double& yplus, double& xminus, double& yminus)
  {
    if (ny == 0)
    {
      // vertical lines are defined by ny=0 and x = x0
      xplus = xminus = x0;

      // calculate y accordingly
      yplus = yminus = m * x0 + b;
    }
    else
    {

      double denom = nx + ny*m;
      if(denom == 0) {
        return false; // lines are parallel and there is no intersection
      }

      double x = (nx*x0 + ny*y0 - ny*b)/denom;
      double y = m*x + b;
      // a straight line has a unique intersection point
      xplus = xminus = x;
      yplus = yminus = y;

    }

    return true;
  }


  //* converninece trait for underlying type
  template<class T>
    using underlying_type_t = typename std::underlying_type<T>::type;

  //* convert an strong type enum to integral type
  template<class T>
    constexpr underlying_type_t<T>
    to_underlying_type(T value) noexcept
  { return static_cast<underlying_type_t<T>>(value);}

  // TVector3 streamer
  inline std::ostream& operator << (std::ostream& out, const TVector3& v )
  {
    out << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
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

  //! create track struct from struct from svx track
  MicromegasTrackEvaluator_hp::TrackStruct create_track( SvtxTrack* track )
  {
    MicromegasTrackEvaluator_hp::TrackStruct trackStruct;

    trackStruct._charge = track->get_charge();
    trackStruct._nclusters = track->size_cluster_keys();
    trackStruct._mask = get_mask( track );
    trackStruct._nclusters_mvtx = get_clusters<TrkrDefs::mvtxId>( track );
    trackStruct._nclusters_intt = get_clusters<TrkrDefs::inttId>( track );
    trackStruct._nclusters_tpc = get_clusters<TrkrDefs::tpcId>( track );
    trackStruct._nclusters_micromegas = get_clusters<TrkrDefs::micromegasId>( track );

    trackStruct._chisquare = track->get_chisq();
    trackStruct._ndf = track->get_ndf();

    trackStruct._px = track->get_px();
    trackStruct._py = track->get_py();
    trackStruct._pz = track->get_pz();
    trackStruct._p = get_p( trackStruct._px, trackStruct._py, trackStruct._pz );
    trackStruct._pt = get_pt( trackStruct._px, trackStruct._py );
    trackStruct._eta = get_eta( trackStruct._p, trackStruct._pz );

    return trackStruct;
  }

}

//_____________________________________________________________________
void MicromegasTrackEvaluator_hp::Container::Reset()
{
  _tracks.clear();
}

//_____________________________________________________________________
MicromegasTrackEvaluator_hp::MicromegasTrackEvaluator_hp( const std::string& name ):
  SubsysReco( name)
  {}

//_____________________________________________________________________
int MicromegasTrackEvaluator_hp::Init(PHCompositeNode* topNode )
{
  // print configuration
  std::cout << "MicromegasTrackEvaluator_hp::Init - m_use_default_pedestal: " << m_use_default_pedestal << std::endl;
  std::cout << "MicromegasTrackEvaluator_hp::Init - m_default_pedestal: " << m_default_pedestal << std::endl;
  std::cout
    << "MicromegasTrackEvaluator_hp::Init -"
    << " m_calibration_filename: "
    << (m_calibration_filename.empty() ? "unspecified":m_calibration_filename )
    << std::endl;

  std::cout << "MicromegasTrackEvaluator_hp::Init - m_extrapolate_from_silicon: " << m_extrapolate_from_silicon << std::endl;
  std::cout << "MicromegasTrackEvaluator_hp::Init - m_zero_field: " << m_zero_field << std::endl;
  std::cout << "MicromegasTrackEvaluator_hp::Init - m_min_tpc_layer: " << m_min_tpc_layer << std::endl;
  std::cout << "MicromegasTrackEvaluator_hp::Init - m_max_tpc_layer: " << m_max_tpc_layer << std::endl;

  // read calibrations
  if( !m_calibration_filename.empty() )
  { m_calibration_data.read( m_calibration_filename ); }

  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "MicromegasTrackEvaluator_hp::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // get EVAL node
  iter = PHNodeIterator(dstNode);
  auto evalNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "EVAL"));
  if( !evalNode )
  {
    // create
    std::cout << "MicromegasTrackEvaluator_hp::Init - EVAL node missing - creating" << std::endl;
    evalNode = new PHCompositeNode( "EVAL" );
    dstNode->addNode(evalNode);
  }

  // add container to output tree
  auto newNode = new PHIODataNode<PHObject>( new Container, "MicromegasTrackEvaluator_hp::Container","PHObject");

  // overwrite split level for easier offline browsing
  newNode->SplitLevel(99);
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasTrackEvaluator_hp::InitRun(PHCompositeNode* topnode )
{
  std::cout << "MicromegasTrackEvaluator_hp::InitRun" << std::endl;
  load_nodes(topnode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasTrackEvaluator_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;
  if( m_container ) m_container->Reset();

  evaluate_tracks();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasTrackEvaluator_hp::End(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int MicromegasTrackEvaluator_hp::load_nodes( PHCompositeNode* topNode )
{

  // acts geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  assert( m_tGeometry );

  // micromegas geometry
  m_micromegas_geomcontainer = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");

  // track map
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_trackmapname);

  // hitset container
  m_hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  assert( m_cluster_map );

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  // local container
  m_container = findNode::getClass<Container>(topNode, "MicromegasTrackEvaluator_hp::Container");
  assert(m_container);

  // tpc distortion corrections
  m_globalPositionWrapper.loadNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void MicromegasTrackEvaluator_hp::evaluate_tracks()
{
  if( !( m_track_map && m_cluster_map && m_container ) ) return;

  // clear array
  m_container->clear_tracks();

  for( const auto& [track_id,track]:*m_track_map )
  {
    // crossing
    const auto crossing = track->get_crossing();
    if (crossing == SHRT_MAX) { continue; }

    // create track information
    auto track_struct = create_track( track );
    track_struct._crossing = crossing;

    std::vector<Acts::Vector3> global_positions;
    std::vector<Acts::Vector3> global_positions_silicon;
    std::vector<Acts::Vector3> global_positions_mvtx;

    // try extrapolate track to Micromegas and find corresponding tile
    /* this is all copied from PHMicromegasTpcTrackMatching */
    const auto cluster_keys = get_cluster_keys(track);
    for( const auto& cluster_key:cluster_keys )
    {
      // detector type
      const auto detid = TrkrDefs::getTrkrId(cluster_key);
      switch( detid )
      {
        case TrkrDefs::micromegasId:
        {
          const unsigned int layer = TrkrDefs::getLayer(cluster_key);
          const auto cluster = m_cluster_map->findCluster(cluster_key);
          if( layer == 55 ) track_struct._found_cluster_phi = create_cluster(cluster_key,cluster);
          if( layer == 56 ) track_struct._found_cluster_z = create_cluster(cluster_key,cluster);
          break;
        }

        case TrkrDefs::tpcId:
        {
          // layer
          const unsigned int layer = TrkrDefs::getLayer(cluster_key);

          // check layer range
          if(!( layer >= m_min_tpc_layer && layer < m_max_tpc_layer) ) continue;

          // get matching
          const auto cluster = m_cluster_map->findCluster(cluster_key);
          const auto global_position = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluster_key, cluster, crossing);
          global_positions.push_back( global_position );
          break;
        }

        case TrkrDefs::mvtxId:
        {
          // get matching
          const auto cluster = m_cluster_map->findCluster(cluster_key);
          const auto global_position = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluster_key, cluster, crossing);
          global_positions_silicon.push_back( global_position );
          global_positions_mvtx.push_back( global_position );
          break;
        }

        case TrkrDefs::inttId:
        {
          // get matching
          const auto cluster = m_cluster_map->findCluster(cluster_key);
          const auto global_position = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(cluster_key, cluster, crossing);
          global_positions_silicon.push_back( global_position );
          break;
        }
      }
    }

    // check global positions
    if( m_extrapolate_from_silicon )
    {

      if( global_positions_mvtx.size()<3 ) { continue; }

    } else {

      if( global_positions.size()<3 ) { continue; }

    }

    // r,z linear fit
    const auto [slope_rz, intersect_rz] = m_extrapolate_from_silicon ?
      TrackFitUtils::line_fit(global_positions_mvtx):
      TrackFitUtils::line_fit(global_positions);

//     std::cout << "MicromegasTrackEvaluator_hp::evaluate_tracks -"
//       << " slope_rz: " << slope_rz
//       << " intersect_rz: " << intersect_rz
//       << std::endl;

    // circle fit parameters
    double R = 0, X0 = 0, Y0 = 0;

    // x,y straight fit paramers
    double slope_xy = 0, intersect_xy = 0;

    if( m_zero_field )
    {

      // x,y straight fit
      std::tie( slope_xy, intersect_xy ) = m_extrapolate_from_silicon ?
        TrackFitUtils::line_fit_xy(global_positions_silicon):
        TrackFitUtils::line_fit_xy(global_positions);

//       std::cout << "MicromegasTrackEvaluator_hp::evaluate_tracks -"
//         << " slope_xy: " << slope_xy
//         << " intersect_xy: " << intersect_xy
//         << std::endl;

    } else {

      // x,y circle
      std::tie( R, X0, Y0 ) = m_extrapolate_from_silicon ?
        TrackFitUtils::circle_fit_by_taubin(global_positions_silicon):
        TrackFitUtils::circle_fit_by_taubin(global_positions);

      if(R < 40.0)
      {
        continue;
      }

    }

    // look over micromegas layers
    const auto range = m_micromegas_geomcontainer->get_begin_end();
    for( const auto& [layer,baselayergeom]:range_adaptor(range))
    {

      const auto layergeom =  static_cast<const CylinderGeomMicromegas*>(baselayergeom);
      assert(layergeom);

      // get layer radius
      const auto layer_radius = layergeom->get_radius();

      // get intersection to track
      auto [xplus, yplus, xminus, yminus] =
        m_zero_field ?
        TrackFitUtils::line_circle_intersection(layer_radius, slope_xy, intersect_xy):
        TrackFitUtils::circle_circle_intersection(layer_radius, R, X0, Y0);

//       std::cout << "MicromegasTrackEvaluator_hp::evaluate_tracks -"
//         << " xplus: " << xplus
//         << " yplus: " << yplus
//         << std::endl;

      // finds the intersection of the fitted circle with the micromegas layer
      if(!std::isfinite(xplus))
      {
        continue;
      }

      // we can figure out which solution is correct based on the last cluster position in the TPC
      const double last_clus_phi = m_extrapolate_from_silicon ?
        std::atan2(global_positions_silicon.back().y(), global_positions_silicon.back().x()):
        std::atan2(global_positions.back().y(), global_positions.back().x());
      double phi_plus = std::atan2(yplus, xplus);
      double phi_minus = std::atan2(yminus, xminus);

      // calculate
      double r = layer_radius;
      double z = intersect_rz + slope_rz * r;

      // select the angle that is the closest to last cluster
      // store phi, apply coarse space charge corrections in calibration mode
      double phi = std::abs(last_clus_phi - phi_plus) < std::abs(last_clus_phi - phi_minus) ? phi_plus : phi_minus;

      // create cylinder intersection point in world coordinates
      const TVector3 world_intersection_cylindrical(r*std::cos(phi), r*std::sin(phi), z);

      // find matching tile
      int tileid = layergeom->find_tile_cylindrical(world_intersection_cylindrical);
      if (tileid < 0)
      {
        continue;
      }

      // now perform planar intersection
      /* note: in principle we could use the ACTS fit to do that */

      // get tile center and norm vector
      const auto tile_center = layergeom->get_world_from_local_coords(tileid, m_tGeometry, {0, 0});
      const double x0 = tile_center.x();
      const double y0 = tile_center.y();

      const auto tile_norm = layergeom->get_world_from_local_vect(tileid, m_tGeometry, {0, 0, 1});
      const double nx = tile_norm.x();
      const double ny = tile_norm.y();

      // calculate intersection to tile, now defined as a straight line
      if( m_zero_field )
      {

        if (!line_line_intersection(slope_xy, intersect_xy, x0, y0, nx, ny, xplus, yplus, xminus, yminus))
        { continue; }

      } else {

        if (!circle_line_intersection(R, X0, Y0, x0, y0, nx, ny, xplus, yplus, xminus, yminus))
        { continue; }

      }

      // select again angle closest to last cluster
      phi_plus = std::atan2(yplus, xplus);
      phi_minus = std::atan2(yminus, xminus);
      const bool is_plus = (std::abs(last_clus_phi - phi_plus) < std::abs(last_clus_phi - phi_minus));

      // calculate global x, y and z
      const double x = (is_plus ? xplus : xminus);
      const double y = (is_plus ? yplus : yminus);
      r = get_r(x, y);
      z = intersect_rz + slope_rz * r;

      // create state
      TrackStateStruct trk_state;
      trk_state._layer = layer;
      trk_state._tile = tileid;
      trk_state._x = x;
      trk_state._y = y;
      trk_state._z = z;

      // get local coordinates
      const auto local_intersection_planar = layergeom->get_local_from_world_coords(tileid, m_tGeometry, {x,y,z});
      trk_state._x_local = local_intersection_planar.x();
      trk_state._y_local = local_intersection_planar.y();

//       std::cout << "MicromegasTrackEvaluator_hp::evaluate_tracks -"
//         << " local_intersection_planar: " << local_intersection_planar.x() << ", " << local_intersection_planar.y()
//         << std::endl;

      // store in track
      const auto segmentation_type = layergeom->get_segmentation_type();
      if( segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_PHI ) track_struct._trk_state_phi = trk_state;
      else  track_struct._trk_state_z = trk_state;

      // keep track of minimum cluster distance to track state
      double dmin = -1;

      // get relevant TPOT clusters
      // generate tilesetid and get corresponding clusters
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation_type, tileid);

      // get cluster from cluster range
      const auto mm_clusrange = m_cluster_map->getClusters(hitsetkey);

      // loop
      for(const auto& [cluster_key, cluster]:range_adaptor(mm_clusrange))
      {
        ClusterStruct cluster_struct = create_cluster(cluster_key, cluster);
        cluster_struct._accepted = accept_cluster( track_struct, cluster_struct );
        track_struct._clusters.push_back(cluster_struct);

        // get distance to track state
        if(segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_PHI)
        {
          const double d = std::abs(trk_state._x_local - cluster_struct._x_local);
          if( dmin<0 || d<dmin )
          {
            dmin = d;
            track_struct._best_cluster_phi = cluster_struct;
          }
        } else {
          const double d = std::abs(trk_state._y_local - cluster_struct._y_local);
          if( dmin<0 || d<dmin )
          {
            dmin = d;
            track_struct._best_cluster_z = cluster_struct;
          }
        }
      }
    }

    // check track state
    if(
      track_struct._trk_state_phi._layer == 0 ||
      track_struct._trk_state_z._layer == 0 ||
      track_struct._trk_state_phi._tile !=
      track_struct._trk_state_z._tile )
    { continue; }

    // print
    if( track_struct._best_cluster_phi._accepted )
    {
      std::cout << "MicromegasTrackEvaluator_hp::evaluate_tracks -"
        << " tpc seed: " << track->get_tpc_seed()
        << " si seed: " << track->get_silicon_seed()
        << " ckey_min: " << track_struct._best_cluster_phi._key
        << std::endl;
    }

    if( track_struct._best_cluster_z._accepted )
    {
      std::cout << "MicromegasTrackEvaluator_hp::evaluate_tracks -"
        << " tpc seed: " << track->get_tpc_seed()
        << " si seed: " << track->get_silicon_seed()
        << " ckey_min: " << track_struct._best_cluster_z._key
        << std::endl;
    }

    // update found cluster accepted flags
    track_struct._found_cluster_phi._accepted = accept_cluster( track_struct, track_struct._found_cluster_phi );
    track_struct._found_cluster_z._accepted = accept_cluster( track_struct, track_struct._found_cluster_z );

    m_container->add_track( track_struct );
  }
}

//______________________________________________________________________________________________________________________
MicromegasTrackEvaluator_hp::ClusterStruct MicromegasTrackEvaluator_hp::create_cluster( TrkrDefs::cluskey cluster_key, TrkrCluster* cluster) const
{
  MicromegasTrackEvaluator_hp::ClusterStruct cluster_struct;

  cluster_struct._layer = TrkrDefs::getLayer(cluster_key);
  cluster_struct._tile = MicromegasDefs::getTileId(cluster_key);
  cluster_struct._size = cluster->getSize();
  cluster_struct._key = cluster_key;

  if( m_cluster_hit_map && m_hitsetcontainer )
  {

    // get associated hits
    const auto range = m_cluster_hit_map->getHits(cluster_key);
    cluster_struct._size = std::distance( range.first, range.second );

    /// Get the upper 32 bits from cluster keys
    const auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);

    // get hits from hitset map
    const auto hitset = m_hitsetcontainer->findHitSet(hitsetkey);

    // loop over assiciated hits
    unsigned int adc_max = 0;
    for( auto iter = range.first; iter != range.second; ++iter )
    {
      const auto& [key,hitkey] = *iter;
      assert( key == cluster_key );

      // get strip
      const auto strip = MicromegasDefs::getStrip( hitkey );

      // get associated hit
      const auto hit = hitset->getHit( hitkey );
      assert( hit );

      const auto adc = hit->getAdc();
      if( adc > adc_max )
      {
        adc_max = adc;
        cluster_struct._strip = strip;
      }

      // get adc, remove pedestal, increment total charge
      const double pedestal = m_use_default_pedestal ?
        m_default_pedestal:
        m_calibration_data.get_pedestal_mapped(hitsetkey, strip);
      cluster_struct._charge += (adc-pedestal);
    }

    // cluster region
    cluster_struct._region = cluster_struct._strip/64;
  }

  // local position
  cluster_struct._x_local = cluster->getLocalX();
  cluster_struct._y_local = cluster->getLocalY();

  // global position
  const auto globalPosition = m_tGeometry->getGlobalPosition(cluster_key, cluster);
  cluster_struct._x = globalPosition.x();
  cluster_struct._y = globalPosition.y();
  cluster_struct._z = globalPosition.z();

  return cluster_struct;
}


//______________________________________________________________________________________________________________________
bool MicromegasTrackEvaluator_hp::accept_cluster(
  const MicromegasTrackEvaluator_hp::TrackStruct& track_struct,
  const MicromegasTrackEvaluator_hp::ClusterStruct& cluster_struct ) const
{

  if( cluster_struct._layer <= 0 ) return false;

  // select proper track state
  const auto& trk_state = (cluster_struct._layer == 55) ?
    track_struct._trk_state_phi:
    track_struct._trk_state_z;

  if( trk_state._layer <= 0 ) return false;

  const float drphi = trk_state._x_local - cluster_struct._x_local;
  const float dz = trk_state._y_local - cluster_struct._y_local;
  const bool accept_phi = std::abs(drphi) < _rphi_search_win[cluster_struct._layer-55];
  const bool accept_z = std::abs(dz) < _z_search_win[cluster_struct._layer-55];
  return accept_phi && accept_z;
}
