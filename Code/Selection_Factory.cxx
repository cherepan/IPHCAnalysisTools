#include "Selection_Factory.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

#include "Example.h"
#include "TauSpinExample.h"
#ifdef USE_cherepanov
#include "cherepanov/MyTest.h"
#include "cherepanov/NtupleValidation.h"
#include "cherepanov/SkimmingNtuples.h"
#endif

#ifdef USE_goe

#endif

#ifdef USE_lebihan

#endif
// #ifdef USE_<username>

// #endif


Selection_Factory::Selection_Factory(){
}

Selection_Factory::~Selection_Factory(){
}

Selection_Base* Selection_Factory::Factory(TString Analysis, TString UncertType,int mode, int runtype, double lumi){
  Selection_Base* s;
  Analysis.ToLower();

  // ensuring code will compile independently of user code
  // WARNING: be aware of the consequences of "Contains". Make sure that Class "foo" is put after "foobar".
  if(Analysis.Contains("example"))s=new Example(Analysis,UncertType);
  else if(Analysis.Contains("tauspin"))s=new TauSpinExample(Analysis,UncertType);
#ifdef USE_cherepanov
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
  else if(Analysis.Contains("ntuplevalidation"))s=new NtupleValidation(Analysis,UncertType);
  else if(Analysis.Contains("skimmingntuples"))s=new SkimmingNtuples(Analysis,UncertType);
#endif
// #ifdef USE_goe
//   else if(Analysis.Contains("bla"))s=new Bla(Analysis,UncertType);
// #endif

// #ifdef USE_lebihan
//   else if(Analysis.Contains("bla"))s=new Bla(Analysis,UncertType);
// #endif


  else{
	Logger(Logger::Error)<< "Invalid Analysis type \"" << Analysis << "\". Using default <Example.h> " << std::endl;
    s=new Example(Analysis,UncertType);
  }
  s->SetMode(mode);
  s->SetRunType(runtype);
  s->SetLumi(lumi);
  s->Configure();
  return s;
}
