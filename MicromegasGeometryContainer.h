#ifndef G4EVAL_HP_MICROMEGASGEOMETRYCONTAINER_H
#define G4EVAL_HP_MICROMEGASGEOMETRYCONTAINER_H

/**
 * @file tpccalib/MicromegasGeometryContainer.h
 * @author Hugo Pereira Da Costa
 * @date June 2018
 * @brief Contains micromegas strip positions
 */

#include <phool/PHObject.h>
#include <TVector3.h>
#include <map>

/**
 * @brief Cluster container object
 */
class MicromegasGeometryContainer : public PHObject
{
  public:

  /// constructor
  MicromegasGeometryContainer() = default;

  /// destructor
  ~MicromegasGeometryContainer() override = default;

  /// strip id
  class StripId
  {
    public:
    unsigned int layer = 0;
    unsigned int tile = 0;
    unsigned int strip = 0;

    bool operator == (const StripId& other ) const
    { return other.layer == layer && other.tile == tile && other.strip == strip; }

    bool operator < (const StripId& other ) const
    {
      if( layer != other.layer ) return layer < other.layer;
      else if( tile != other.tile ) return tile < other.tile;
      else return strip < other.strip;
    }
  };

  // strip geometry
  class StripGeometry
  {
    public:

    TVector3 local_begin;
    TVector3 local_end;

    TVector3 global_begin;
    TVector3 global_end;
  };

  ///@name accessors
  //@{

  /// identify object
  void identify(std::ostream &/*os*/ = std::cout) const override;

  /// get strip begin from layer, tile and strip number
  TVector3 get_strip_local_begin( unsigned int /*layer*/, unsigned int /*tile*/, unsigned int /*strip*/ ) const;

  /// get strip begin from layer, tile and strip number
  TVector3 get_strip_local_end( unsigned int /*layer*/, unsigned int /*tile*/, unsigned int /*strip*/ ) const;

  /// get strip begin from layer, tile and strip number
  TVector3 get_strip_begin( unsigned int /*layer*/, unsigned int /*tile*/, unsigned int /*strip*/ ) const;

  /// get strip begin from layer, tile and strip number
  TVector3 get_strip_end( unsigned int /*layer*/, unsigned int /*tile*/, unsigned int /*strip*/ ) const;

  //@}

  ///@name modifiers
  //@{

  /// reset method
  void Reset() override;

  /// add strip
  void add_strip( const StripId&, const StripGeometry& );

  //@}

  private:

  using strip_map_t = std::map<StripId, StripGeometry>;
  strip_map_t m_strips;

  ClassDefOverride(MicromegasGeometryContainer, 1)

};

#endif
