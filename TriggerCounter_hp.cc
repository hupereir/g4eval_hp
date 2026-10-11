#include "TriggerCounter_hp.h"

#include <calotrigger/TriggerAnalyzer.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <qautils/QAHistManagerDef.h>

#include <TH1.h>

#include <cassert>
#include <iostream>

namespace
{
  constexpr int n_trigger_bits = 64;
}

//_____________________________________________________________________
TriggerCounter_hp::TriggerCounter_hp( const std::string& name ):
  SubsysReco( name),
  m_triggeranalyzer( new TriggerAnalyzer() )
{ std::cout << "TriggerCounter_hp::TriggerCounter_hp." << std::endl; }

//_____________________________________________________________________
int TriggerCounter_hp::Init(PHCompositeNode*)
{

  std::cout << "TriggerCounter_hp::Init." << std::endl;

  // histogram manager
  auto* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // trigger counts
  m_h_trigger_stat = new TH1L( "m_h_trigger_stat", "Trigger statistics", n_trigger_bits, 0, n_trigger_bits );
  m_h_trigger_stat->GetXaxis()->SetTitle( "trigger bit" );
  m_h_trigger_stat->SetFillStyle(1001);
  m_h_trigger_stat->SetFillColor(kYellow);

  // assign trigger names to bins
  m_triggeranalyzer->decodeTriggers(topNode);
  for (int i = 0; i < n_trigger_bits; ++i)
  { m_h_trigger_stat->SetBinLabel(i+1, m_triggeranalyzer->getTriggerName(i).c_str() ); }

  // register
  hm->registerHisto(m_h_trigger_stat);

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TriggerCounter_hp::InitRun(PHCompositeNode* )
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int TriggerCounter_hp::process_event(PHCompositeNode* topNode )
{
  // decode trigger information
  m_triggeranalyzer->decodeTriggers(topNode);

  // loop over bits fill histogram for fired triggers
  for (int i = 0; i < n_trigger_bits; ++i)
  {
    if( m_triggeranalyzer->didTriggerFire(i) )
    { m_h_trigger_stat->Fill(i); }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TriggerCounter_hp::End(PHCompositeNode*)
{
  std::cout << "TriggerCounter_hp::End." << std::endl;

  // print trigger information
  if( m_h_trigger_stat )
  {
    for (int i = 0; i < n_trigger_bits; ++i)
    {
      std::cout << "TriggerCounter_hp::End -"
        << " trigger bit: " << i
        << " name: " << m_h_trigger_stat->GetBinLabel(i+1)
        << " counts: " << m_h_trigger_stat->GetBinContent(i+1)
        << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
