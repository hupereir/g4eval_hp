#ifndef G4EVAL_CLEARTRACKSEEDS_HP_H
#define G4EVAL_CLEARTRACKSEEDS_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>


class ClearTrackSeeds_hp : public SubsysReco
{
  public:

  //! constructor
  ClearTrackSeeds_hp( const std::string& = "ClearTrackSeeds_hp" );

  //! global initialization
  int Init(PHCompositeNode*) override;

  //! run initialization
  int InitRun(PHCompositeNode*) override;

  //! event processing
  int process_event(PHCompositeNode*) override;

  //! end of processing
  int End(PHCompositeNode*) override;
};

#endif
