#ifndef G4EVAL_MicromegasTrackEvaluator_hp_H
#define G4EVAL_MicromegasTrackEvaluator_hp_H

#include <fun4all/SubsysReco.h>
#include <gsl/gsl_rng.h>
#include <micromegas/MicromegasCalibrationData.h>
#include <micromegas/MicromegasDefs.h>
#include <phool/PHObject.h>
#include <tpc/TpcGlobalPositionWrapper.h>
#include <trackbase/TrkrDefs.h>

#include <map>
#include <memory>

class ActsGeometry;
class PHG4CylinderGeomContainer;
class TrkrCluster;
class TrkrClusterHitAssoc;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class SvtxTrack;
class SvtxTrackMap;

class MicromegasTrackEvaluator_hp : public SubsysReco
{
  public:

  /// constructor
  MicromegasTrackEvaluator_hp( const std::string& = "MicromegasTrackEvaluator_hp" );

  /// global initialization
  virtual int Init(PHCompositeNode*);

  /// run initialization
  virtual int InitRun(PHCompositeNode*);

  /// event processing
  virtual int process_event(PHCompositeNode*);

  /// end of processing
  virtual int End(PHCompositeNode*);

  /// set default pedestal
  void set_default_pedestal( double value )
  { m_default_pedestal = value; }

  /// set whether default pedestal is used or not
  void set_use_default_pedestal( bool value )
  { m_use_default_pedestal = value; }

  /// calibration file
  void set_calibration_file( const std::string& value )
  { m_calibration_filename = value; }

  class TrackStateStruct
  {
    public:

    unsigned short _layer = 0;
    unsigned short _tile = 0;

    double _x_local = 0;
    double _y_local = 0;

    double _x = 0;
    double _y = 0;
    double _z = 0;
  };

  class ClusterStruct
  {
    public:

    unsigned short _layer = 0;
    unsigned short _tile = 0;
    unsigned short _size = 0;
    double _charge = 0;
    int _strip = 0;
    int _region = 0;

    double _x_local = 0;
    double _y_local = 0;

    double _x = 0;
    double _y = 0;
    double _z = 0;

    using List = std::vector<ClusterStruct>;
  };


  // track information to be stored in tree
  class TrackStruct
  {
    public:

    // constructor
    explicit TrackStruct()
    {}

    using List = std::vector<TrackStruct>;

    int _charge = 0;
    unsigned int _nclusters = 0;

    /// crossing
    unsigned int _crossing = 0;

    /// mask of layers for which there is a cluster in the track
    int64_t _mask = 0LL;

    unsigned int _nclusters_mvtx = 0;
    unsigned int _nclusters_intt = 0;
    unsigned int _nclusters_tpc = 0;
    unsigned int _nclusters_micromegas = 0;

    float _chisquare = 0;
    int _ndf = 0;

    //!@name momentum
    //@{
    float _px = 0;
    float _py = 0;
    float _pz = 0;
    float _pt = 0;
    float _p = 0;
    float _eta = 0;
    //@}

    // track state extrapolated to phi layer
    TrackStateStruct _trk_state_phi;

    // track state extrapolated to z layer
    TrackStateStruct _trk_state_z;

    // all micromegas clusters in track
    ClusterStruct::List _clusters;

    // best micromegas clustser
    ClusterStruct _best_cluster_phi;

    // best micromegas cluster
    ClusterStruct _best_cluster_z;

  };

  class Container: public PHObject
  {

    public:

    //! constructor
    explicit Container() = default;

    //! copy constructor
    explicit Container(const Container &) = delete;

    //! assignment operator
    Container& operator = ( const Container& ) = delete;

    //! reset
    virtual void Reset();

    //!@name accessors
    //@{

    const TrackStruct::List& tracks() const
    { return _tracks; }

    //@}

    //!@name modifiers
    //@{

    void add_track( const TrackStruct& track )
    { _tracks.push_back( track ); }

    void clear_tracks()
    { _tracks.clear(); }

    //@}
    private:

    //! tracks array
    TrackStruct::List _tracks;

    //! unused
    /* it is only added here so that the dictionary is available in root file */
    ClusterStruct _unused_cluster;

    /* it is only added here so that the dictionary is available in root file */
    TrackStateStruct _unused_trk_state;

    ClassDef(Container,1)

  };

  //! tracl map name
  void set_trackmapname( const std::string& value )
  { m_trackmapname = value; }

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// evaluate tracks
  void evaluate_tracks();

  //! create cluster structure from cluster
  ClusterStruct create_cluster( TrkrDefs::cluskey, TrkrCluster* ) const;

  /// cluster array
  Container* m_container = nullptr;

  /// Acts tracking geometry for surface lookup
  ActsGeometry *m_tGeometry = nullptr;

  //! global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  // micromegas geometry
  PHG4CylinderGeomContainer* m_micromegas_geomcontainer = nullptr;

  //! hits
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;

  //! clusters
  TrkrClusterContainer* m_cluster_map = nullptr;

  //! cluster to hit association
  TrkrClusterHitAssoc* m_cluster_hit_map = nullptr;

  //! track map name
  std::string m_trackmapname = "SvtxTrackMap";

  //! tracks
  SvtxTrackMap* m_track_map = nullptr;

  //!@name calibration filename
  //@{

  /// if true, use default pedestal to get hit charge. Relies on calibration data otherwise
  bool m_use_default_pedestal = true;

  /// default pedestal
  double m_default_pedestal = 74.6;

  /// calibration filename
  std::string m_calibration_filename;

  /// calibration data
  MicromegasCalibrationData m_calibration_data;

  //@}

  // range of TPC layers to use in projection to micromegas
  unsigned int _min_tpc_layer = 38;

};

#endif  // G4EVAL_MicromegasTrackEvaluator_hp_H
