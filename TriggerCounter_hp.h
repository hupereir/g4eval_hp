#ifndef G4EVAL_TriggerCounter_hp_H
#define G4EVAL_TriggerCounter_hp_H

#include <fun4all/SubsysReco.h>

#include <memory>

class TH1;
class TriggerAnalyzer;

class TriggerCounter_hp : public SubsysReco
{
  public:

  /// constructor
  TriggerCounter_hp( const std::string& = "TriggerCounter_hp" );

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// run initialization
  int InitRun(PHCompositeNode*)  override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  private:

  // trigger analizer
  std::unique_ptr<TriggerAnalyzer> m_triggeranalyzer{};

  // trigger counts histoframs
  TH1* m_h_trigger_stat{nullptr};


};

#endif
