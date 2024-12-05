#include "EventCounter_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHTimer.h>

//_____________________________________________________________________
EventCounter_hp::EventCounter_hp( const std::string& name, unsigned int granularity ):
  SubsysReco( name),
  m_granularity( granularity )
{ std::cout << "EventCounter_hp::EventCounter_hp." << std::endl; }

//_____________________________________________________________________
int EventCounter_hp::Init(PHCompositeNode*)
{

  std::cout << "EventCounter_hp::Init." << std::endl;

  // initialize timer
  m_timer.restart();
  m_running_timer.restart();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int EventCounter_hp::InitRun(PHCompositeNode* )
{
  std::cout << "EventCounter_hp::InitRun." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int EventCounter_hp::process_event(PHCompositeNode*)
{
  // print event number
  if( m_granularity > 0 && (m_ievent % m_granularity) == 0 )
  {
    std::cout << "EventCounter_hp::process_event -"
      << " event = " << m_ievent
      << " time (ms): " << m_running_timer.elapsed()
      << std::endl;
    m_running_timer.restart();
  }
  ++m_ievent;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int EventCounter_hp::End(PHCompositeNode*)
{
  std::cout << "EventCounter_hp::End." << std::endl;

  // print timer information
  m_timer.stop();
  std::cout
    << "EventCounter_hp::End -"
    << " events: " << m_ievent
    << " time per event (ms):" << m_timer.get_accumulated_time()/m_ievent
    << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}
