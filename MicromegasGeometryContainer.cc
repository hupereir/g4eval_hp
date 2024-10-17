
#include "MicromegasGeometryContainer.h"

//___________________________________________________________
void MicromegasGeometryContainer::identify( std::ostream& out ) const
{
  out << "MicromegasGeometryContainer" << std::endl;
  out << " m_strips entries: " << m_strips.size() << std::endl;
}

//___________________________________________________________
TVector3 MicromegasGeometryContainer::get_strip_local_begin( unsigned int layer, unsigned int tile, unsigned int strip ) const
{ return m_strips.at({.layer=layer, .tile=tile, .strip=strip}).local_begin; }

//___________________________________________________________
TVector3 MicromegasGeometryContainer::get_strip_local_end( unsigned int layer, unsigned int tile, unsigned int strip ) const
{ return m_strips.at({.layer=layer, .tile=tile, .strip=strip}).local_end; }

//___________________________________________________________
TVector3 MicromegasGeometryContainer::get_strip_begin( unsigned int layer, unsigned int tile, unsigned int strip ) const
{ return m_strips.at({.layer=layer, .tile=tile, .strip=strip}).global_begin; }

//___________________________________________________________
TVector3 MicromegasGeometryContainer::get_strip_end( unsigned int layer, unsigned int tile, unsigned int strip ) const
{ return m_strips.at({.layer=layer, .tile=tile, .strip=strip}).global_end; }

//___________________________________________________________
void MicromegasGeometryContainer::Reset()
{ m_strips.clear(); }

//___________________________________________________________
void MicromegasGeometryContainer::add_strip(
  const MicromegasGeometryContainer::StripId& strip,
  const MicromegasGeometryContainer::StripGeometry& geometry )
{ m_strips[strip]=geometry; }
