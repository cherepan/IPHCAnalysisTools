#ifndef MyExample_h
#define MyExample_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SVFitStorage.h"

class MyExample : public Selection {

 public:
  MyExample(TString Name_, TString id_);
  virtual ~MyExample();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,PrimeVtx,NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> NumVertices;

};
#endif
