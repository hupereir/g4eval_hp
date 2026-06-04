#ifndef G4EVAL_MULTIPLICITYEVALUATOR_HP_H
#define G4EVAL_MULTIPLICITYEVALUATOR_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>

class MultiplicityEvaluator_hp: public SubsysReco
{
  public:

  //! constructor
  MultiplicityEvaluator_hp( const std::string& = "MultiplicityEvaluator_hp" );

  /// global initialization
  int Init(PHCompositeNode*) override;

  /// event processing
  int process_event(PHCompositeNode*) override;

  /// end of processing
  int End(PHCompositeNode*) override;

  class MultiplicityStruct
  {
    public:

    using List = std::vector<MultiplicityStruct>;

    uint _rawhits = 0;
    uint _hits = 0;
    uint _clusters = 0;

    uint64_t _gtm_bco = 0;

  };

  /// track container
  class Container: public PHObject
  {

    public:

    /// constructor
    explicit Container() = default;

    /// copy constructor
    explicit Container(const Container &) = delete;

    /// assignment operator
    Container& operator = ( const Container& ) = delete;

//     /// reset
//     virtual void Reset()
//     {}

    //!@name accessors
    //@{
    const MultiplicityStruct& current_multiplicity() const { return _mult; }
    const MultiplicityStruct& previous_multiplicity() const { return _prev_mult; }
    //}

    //!@name modifiers
    //@{
    void set_current_multiplicity( const MultiplicityStruct& mult ) { _mult = mult; }
    void set_previous_multiplicity( const MultiplicityStruct& mult ) { _prev_mult = mult; }
    //@}

    private:

    // current event multiplicities
    MultiplicityStruct _mult;

    // previous event multiplicities
    MultiplicityStruct _prev_mult;

    ClassDefOverride(Container,1)

  };

  private:

  //! evaluation node
  Container* m_container = nullptr;

};

#endif
