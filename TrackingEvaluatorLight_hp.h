#ifndef G4EVAL_TrackingEvaluatorLight_hp_H
#define G4EVAL_TrackingEvaluatorLight_hp_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <tpc/TpcClusterMover.h>
#include <tpc/TpcGlobalPositionWrapper.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/ActsTransformations.h>

#include <map>
#include <set>
#include <string>
#include <vector>

class ActsGeometry;
class CMFlashClusterContainer;
class MicromegasRawHitContainer;
class PHG4CylinderGeomContainer;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TpcGeomContainer;
class PHG4TruthInfoContainer;
class RawClusterContainer;
class SvtxTrack;
class SvtxTrackMap;
class SvtxTrackState;
class TrackSeed;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class TrkrHitTruthAssoc;

class Gl1RawHit;

class TrackingEvaluatorLight_hp : public SubsysReco
{
  public:

  //! constructor
  TrackingEvaluatorLight_hp( const std::string& = "TrackingEvaluatorLight_hp" );

  //! global initialization
  virtual int Init(PHCompositeNode*);

  //! run initialization
  virtual int InitRun(PHCompositeNode*);

  //! event processing
  virtual int process_event(PHCompositeNode*);

  //! end of processing
  virtual int End(PHCompositeNode*);

  // tower information
  class TowerStruct
  {
    public:

    using List = std::vector<TowerStruct>;

    int _ieta = 0;
    int _iphi = 0;
    float _e = 0;
  };

  // cluster information to be stored in tree
  class CaloClusterStruct
  {
    public:

    using List = std::vector<CaloClusterStruct>;

    //! cluster layer
    int _layer = SvtxTrack::PRES;

    //! number of hits belonging to the cluster
    unsigned int _size = 0;

    //!@name cluster position and energy
    //@{
    float _x = 0;
    float _y = 0;
    float _z = 0;
    float _r = 0;

    // geometrical phi and eta
    float _phi = 0;
    float _eta = 0;

    // energy
    float _e = 0;
    float _chisquare = 0;
    //@}

    //!@name matching track position
    //@{
    float _trk_x = 0;
    float _trk_y = 0;
    float _trk_z = 0;
    float _trk_r = 0;

    // geometrical phi and eta
    float _trk_phi = 0;
    float _trk_eta = 0;

    // extrapolation length
    float _trk_dr = 0;
    //@}

    //! towers
    TowerStruct::List _towers;
  };

  // track information to be stored in tree
  class TrackStruct
  {
    public:

    // constructor
    explicit TrackStruct() = default;

    using List = std::vector<TrackStruct>;

    int _charge = 0;
    unsigned int _nclusters = 0;

    /// crossing
    short int _crossing = 0;

    /// mask of layers for which there is a cluster in the track
    int64_t _mask = 0LL;

    /// mask of layers for which there is a cluster in the track, with correct associated g4hit
    int64_t _correct_mask = 0LL;

    /// mask of layers for which there is a cluster in the track, with correct associated g4hit
    int64_t _correct_mask_strict = 0LL;

    /// maks of layers for which there is a g4hit in the track
    int64_t _truth_mask = 0LL;

    unsigned int _nclusters_mvtx = 0;
    unsigned int _nclusters_intt = 0;
    unsigned int _nclusters_tpc = 0;
    unsigned int _nclusters_micromegas = 0;
    unsigned int _nclusters_micromegas_phi = 0;
    unsigned int _nclusters_micromegas_z = 0;

    /// number of track states on track
    unsigned int _nstates = 0;
    unsigned int _nstates_mvtx = 0;
    unsigned int _nstates_intt = 0;
    unsigned int _nstates_tpc = 0;
    unsigned int _nstates_micromegas = 0;

    float _chisquare = 0;
    int _ndf = 0;

    //!@name position
    //@{
    float _x = 0;
    float _y = 0;
    float _z = 0;
    float _r = 0;
    float _phi = 0;
    //@}

    //!@name momentum
    //@{
    float _px = 0;
    float _py = 0;
    float _pz = 0;
    float _pt = 0;
    float _p = 0;
    float _eta = 0;
    //@}

    //! dedx in TPC using cluster information only
    float _dedx = 0;

    //!@name truth momentum
    //@{
    int _pid = 0;
    int _embed = 0;
    bool _is_primary = false;

    // number of g4hits from this MC track that match
    unsigned int _contributors = 0;

    float _truth_px = 0;
    float _truth_py = 0;
    float _truth_pz = 0;
    float _truth_pt = 0;
    float _truth_p = 0;
    float _truth_eta = 0;

    // dedx in TPC using truth information (not sure how to calculate)
    float _truth_dedx = 0;

    //@}

    CaloClusterStruct _calo_cluster_emcal;
    CaloClusterStruct _calo_cluster_ihcal;
    CaloClusterStruct _calo_cluster_ohcal;
    CaloClusterStruct _calo_cluster_topo;
  };

  //! track container
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
    void Reset() override;

    //!@name accessors
    //@{

    const TrackStruct::List& tracks() const
    { return _tracks; }

    //@}

    //!@name modifiers
    //@{

    void addTrack( const TrackStruct& track )
    { _tracks.push_back( track ); }

    void clearTracks()
    { _tracks.clear(); }
    //@}

    private:

    //! tracks array
    TrackStruct::List _tracks;

    ClassDefOverride(Container,1)

  };

  //! track map name
  void set_trackmapname( const std::string& value )
  { m_trackmapname = value; }

  //! calorimeter min enery
  void set_calo_min_energy( int layer, double value )
  { m_calo_min_energy[layer] = value; }

  //!@name utility functions
  //@{
  //! tells if a given layer is in a layer bitwise mask
  static bool has_layer( int64_t mask, int layer )
  { return mask & (1LL<<layer); }

  //! get number of clusters in given range
  static int get_nclusters( int64_t mask, int first, int last )
  {
    int out = 0;
    for( int layer = first; layer < last; ++layer )
    { out += (int) has_layer( mask, layer ); }

    return out;
  }

  //! get number of mvtx clusters from mask
  static int get_nclusters_mvtx( int64_t mask )
  { return get_nclusters( mask, 0, 3 ); }

  //! get number of intt clusters from mask
  static int get_nclusters_intt( int64_t mask )
  { return get_nclusters( mask, 3, 7 ); }

  //! get number of tpc clusters from mask
  static int get_nclusters_tpc( int64_t mask )
  { return get_nclusters( mask, 7, 55 ); }

  //! get number of micromegas clusters from mask
  static int get_nclusters_micromegas( int64_t mask )
  { return get_nclusters( mask, 55, 57 ); }
  //@}

  private:

  //! load nodes
  int load_nodes( PHCompositeNode* );

  //! evaluate tracks
  void evaluate_tracks();

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits( TrkrDefs::cluskey ) const;

  //! get G4Particle id of max contributor to a given track
  std::pair<int,int> get_max_contributor( SvtxTrack* ) const;

  //! get embedded id for given g4track
  int get_embed(PHG4Particle*) const;

  //! fill MC track map
  void fill_g4particle_map();

  //! calculate dedx from clusters only
  float get_dedx( TrackSeed* ) const;

  //! calculate dedx from g4hits only
  float get_truth_dedx( TrackSeed*, int) const;

  //! find calorimeter cluster matching track
  std::optional<CaloClusterStruct> find_calo_cluster_emcal( SvtxTrack* ) const;

  //! find calorimeter cluster matching track
  std::optional<CaloClusterStruct> find_calo_cluster_hcal( int /*layer*/, SvtxTrack* ) const;

  //! evaluation node
  Container* m_container = nullptr;

  /// Acts tracking geometry for surface lookup
  ActsGeometry *m_tGeometry = nullptr;

  //! hits
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;

  //! clusters
  TrkrClusterContainer* m_cluster_map = nullptr;

  //! cluster to hit association
  TrkrClusterHitAssoc* m_cluster_hit_map = nullptr;

  //! hit to truth association
  TrkrHitTruthAssoc* m_hit_truth_map = nullptr;

  //! track map name
  std::string m_trackmapname = "SvtxTrackMap";

  //! tracks
  SvtxTrackMap* m_track_map = nullptr;

  //! calorimeter cluster maps
  using calo_clusters_map_t = std::map<int, RawClusterContainer*>;
  calo_clusters_map_t m_rawclustercontainermap;

  //! calorimeter min energy cut
  using calo_min_energy_map_t = std::map<int, float>;
  calo_min_energy_map_t m_calo_min_energy;

  //!@name geant4 hits
  //@{
  PHG4HitContainer* m_g4hits_tpc = nullptr;
  PHG4HitContainer* m_g4hits_intt = nullptr;
  PHG4HitContainer* m_g4hits_mvtx = nullptr;
  PHG4HitContainer* m_g4hits_micromegas = nullptr;
  //@}

  //! truth information
  PHG4TruthInfoContainer* m_g4truthinfo = nullptr;

  //! tpc geometry
  PHG4TpcGeomContainer* m_tpc_geom_container = nullptr;

  //! micromegas geometry
  PHG4CylinderGeomContainer* m_micromegas_geom_container = nullptr;

  //! map cluster keys to g4hits
  using G4HitMap = std::map<TrkrDefs::cluskey,G4HitSet>;
  mutable G4HitMap m_g4hit_map;

  //! map trk_id to layer mask
  /** copied from SimEvaluator_hp */
  using G4ParticleMap = std::map<int,int64_t>;
  G4ParticleMap m_g4particle_map;

};

#endif  // G4EVAL_TrackingEvaluatorLight_hp_H
