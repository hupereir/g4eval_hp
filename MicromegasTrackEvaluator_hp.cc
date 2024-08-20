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

#include <tpc/TpcDistortionCorrectionContainer.h>

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

  auto newNode = new PHIODataNode<PHObject>( new Container, "MicromegasTrackEvaluator_hp::Container","PHObject");
  evalNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasTrackEvaluator_hp::InitRun(PHCompositeNode* /*topnode*/ )
{
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
  assert(m_hitsetcontainer);

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  assert( m_cluster_map );

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  assert( m_cluster_hit_map );

  // local container
  m_container = findNode::getClass<Container>(topNode, "MicromegasTrackEvaluator_hp::Container");
  assert(m_container);

  // tpc distortion corrections
  m_dcc_module_edge = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerModuleEdge");
  m_dcc_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  m_dcc_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerAverage");
  m_dcc_fluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerFluctuation");

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

    // try extrapolate track to Micromegas and find corresponding tile
    /* this is all copied from PHMicromegasTpcTrackMatching */
    const auto cluster_keys = get_cluster_keys(track);

    for( const auto cluster_key:cluster_keys )
    {
      // detector type
      const auto detid = TrkrDefs::getTrkrId(cluster_key);
      if( detid != TrkrDefs::tpcId ) continue;

      // layer
      const unsigned int layer = TrkrDefs::getLayer(cluster_key);
      if (layer < _min_tpc_layer) continue;

      // get matching
      const auto cluster = m_cluster_map->findCluster(cluster_key);
      const auto global_position = get_global_position(cluster_key, cluster, crossing);
      global_positions.push_back( global_position );

    }

    // check global positions
    if( global_positions.size()<3 ) { continue; }

    // helical fit
    const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin(global_positions);
    const auto [A, B] = TrackFitUtils::line_fit(global_positions);
    if(R < 40.0)
    {
      continue;
    }

    // look over micromegas layers
    const auto range = m_micromegas_geomcontainer->get_begin_end();
    for( auto iter = range.first; iter != range.second; ++iter )
    {
      const auto layergeom =  static_cast<CylinderGeomMicromegas*>(iter->second);
      assert(layergeom);

      // get layer radius
      const auto layer_radius = layergeom->get_radius();

      // get intersection to track
      auto [xplus, yplus, xminus, yminus] = TrackFitUtils::circle_circle_intersection(layer_radius, R, X0, Y0);

      // finds the intersection of the fitted circle with the micromegas layer
      if(!std::isfinite(xplus))
      {
        continue;
      }

      // we can figure out which solution is correct based on the last cluster position in the TPC
      const double last_clus_phi = std::atan2(global_positions.back().y(), global_positions.back().x());
      double phi_plus = std::atan2(yplus, xplus);
      double phi_minus = std::atan2(yminus, xminus);

      // calculate
      double r = layer_radius;
      double z = B + A * r;

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

      // calculate intersection to tile
      if (!circle_line_intersection(R, X0, Y0, x0, y0, nx, ny, xplus, yplus, xminus, yminus))
      {
        continue;
      }

      // select again angle closest to last cluster
      phi_plus = std::atan2(yplus, xplus);
      phi_minus = std::atan2(yminus, xminus);
      const bool is_plus = (std::abs(last_clus_phi - phi_plus) < std::abs(last_clus_phi - phi_minus));

      // calculate global x, y and z
      const double x = (is_plus ? xplus : xminus);
      const double y = (is_plus ? yplus : yminus);
      r = get_r(x, y);
      z = B + A * r;

      // create state
      TrackStateStruct trk_state;
      trk_state._layer = layergeom->get_layer();
      trk_state._tile = tileid;
      trk_state._x = x;
      trk_state._y = y;
      trk_state._z = z;

      // std::cout << "MicromegasTrackEvaluator_hp::evaluate_tracks - layer: " << (int) layergeom->get_layer() << " tile: " << tileid << " position: " << x << ", " << y << ", " << z << std::endl;

      // get local coordinates
      const auto local_intersection_planar = layergeom->get_local_from_world_coords(tileid, m_tGeometry, {x,y,z});
      trk_state._x_local = local_intersection_planar.x();
      trk_state._y_local = local_intersection_planar.y();

      // store in track
      const auto segmentation_type = layergeom->get_segmentation_type();
      if( segmentation_type == MicromegasDefs::SegmentationType::SEGMENTATION_PHI ) track_struct._trk_state_phi = trk_state;
      else  track_struct._trk_state_z = trk_state;

    }

    // check track state
    if(
      track_struct._trk_state_phi._layer == 0 ||
      track_struct._trk_state_z._layer == 0 ||
      track_struct._trk_state_phi._tile !=
      track_struct._trk_state_z._tile )
    { continue; }

    m_container->add_track( track_struct );
  }


}

//_________________________________________________________________________________
Acts::Vector3 MicromegasTrackEvaluator_hp::get_global_position( TrkrDefs::cluskey key, TrkrCluster* cluster, short int crossing ) const
{
  // get global position from Acts transform
  auto globalPosition = m_tGeometry->getGlobalPosition(key, cluster);

  // for the TPC calculate the proper z based on crossing and side
  const auto trkrid = TrkrDefs::getTrkrId(key);
  if(trkrid ==  TrkrDefs::tpcId)
  {
    const auto side = TpcDefs::getSide(key);
    globalPosition.z() = m_clusterCrossingCorrection.correctZ(globalPosition.z(), side, crossing);

    // apply distortion corrections
    if(m_dcc_module_edge)
    {
      globalPosition = m_distortionCorrection.get_corrected_position( globalPosition, m_dcc_module_edge );
    }

    if(m_dcc_static)
    {
      globalPosition = m_distortionCorrection.get_corrected_position( globalPosition, m_dcc_static );
    }

    if(m_dcc_average)
    {
      globalPosition = m_distortionCorrection.get_corrected_position( globalPosition, m_dcc_average );
    }

    if(m_dcc_fluctuation)
    {
      globalPosition = m_distortionCorrection.get_corrected_position( globalPosition, m_dcc_fluctuation );
    }
  }

  return globalPosition;
}
