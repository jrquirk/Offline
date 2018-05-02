#ifndef __TrkPatRec_RobustHelixFinderDiag_hh__
#define __TrkPatRec_RobustHelixFinderDiag_hh__

#include "TrkPatRec/inc/RobustHelixFinder_types.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "TrkReco/inc/RobustHelixFit.hh"

#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrkData/inc/TrkStrawHit.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "TH1F.h"

namespace mu2e {

  using namespace RobustHelixFinderTypes;
  
  class RobustHelixFinderDiag : public mu2e::ModuleHistToolBase {
  public:

    enum {
      kNEventHistSets = 10,
      kNHelixHistSets = 10,
      kNHitHistSets   = 10
    };

    
    struct Hist_t {
      TH1F*  nTimePeaks;
      TH1F*  ntclhits[2];
      TH1F*  nhits   [2];           // number of hits on a helix  
      TH1F*  nseeds  [2];
      TH1F*  ntripl0 ;
      TH1F*  ntripl1 [2];
      TH1F*  xyniter [2];
      TH1F*  fzniter [2];
      TH1F*  niter   [2];
 
      TH1F*  nXYSh   [2];
      TH1F*  nZPhiSh [2];


      TH1F*  lambda0 [2];
      TH1F*  lambda1 [2];

      TH1F*  rinit   [2];
      TH1F*  radius  [2];   
      TH1F*  chi2XY  [2];
      TH1F*  chi2ZPhi[2];

      TH1F*  pT      [2];
      TH1F*  p       [2];
      TH1F*  dr[2];
      TH1F*  chi2d_helix[2];
    };

  protected:
    int                              _mcTruth;
    //    std::unique_ptr<McUtilsToolBase> _mcUtils;
    std::string                      _shDigiLabel;
    Hist_t                           _hist;            // owned
    Data_t*                          _data;            // cached

    //    const        mu2e::PtrStepPointMCVectorCollection* _listOfMCStrawHits;

  public:

    RobustHelixFinderDiag(const fhicl::ParameterSet& PSet);
    ~RobustHelixFinderDiag();

  private:
    
    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };


//-----------------------------------------------------------------------------
  RobustHelixFinderDiag::RobustHelixFinderDiag(const fhicl::ParameterSet& PSet) {
    printf(" RobustHelixFinderDiag::RobustHelixFinderDiag : HOORAY! \n");
  }

//-----------------------------------------------------------------------------
  RobustHelixFinderDiag::~RobustHelixFinderDiag() {
  }

    
//-----------------------------------------------------------------------------
  int RobustHelixFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {
  
    _hist.nTimePeaks    = Tfs->make<TH1F>("ntpeaks"  , "number of time peaks"                      , 11, -0.5, 10.5);
    _hist.nseeds[0]     = Tfs->make<TH1F>("nseeds0"  , "number of track candidates: all events"    , 21, -0.5, 20.5);
    _hist.nseeds[1]     = Tfs->make<TH1F>("nseeds1"  , "number of track candidates: nhits > 15;"    , 21, -0.5, 20.5);
    _hist.ntclhits[0]   = Tfs->make<TH1F>("ntclhits0" , "number of hits on a time peak - no delta"  , 201, -0.5, 200.5);
    _hist.ntclhits[1]   = Tfs->make<TH1F>("ntclhits1" , "number of hits on a time peak - no delta: nhits > 15"  , 201, -0.5, 200.5);
    _hist.nhits[0]       = Tfs->make<TH1F>("nhitsNeg"    , "number of hits on a track candidate"       , 401, -0.5, 800.5);
    _hist.nhits[1]       = Tfs->make<TH1F>("nhitsPos"    , "number of hits on a track candidate"       , 401, -0.5, 800.5);
    _hist.ntripl0       = Tfs->make<TH1F>("ntripl0"    , "number of triplets fit Circle"           , 401, -0.5, 800.5);

    _hist.lambda0[0]    = Tfs->make<TH1F>("lambda0Neg"    , "initFZ, Neg;#lambda [mm/rad] "   , 401, -0.5, 400.5);
    _hist.lambda0[1]    = Tfs->make<TH1F>("lambda0Pos"    , "initFZ, Pos;#lambda [mm/rad] "   , 401, -0.5, 400.5);
    _hist.lambda1[0]    = Tfs->make<TH1F>("lambda1Neg"    , "fitFZ, Neg;#lambda [mm/rad] "    , 401, -0.5, 400.5);
    _hist.lambda1[1]    = Tfs->make<TH1F>("lambda1Pos"    , "fitFZ, Pos;#lambda [mm/rad]"     , 401, -0.5, 400.5);

    _hist.ntripl1[0]    = Tfs->make<TH1F>("ntripl1Neg"    , "number of triplets fit Helix, Neg "   , 401, -0.5, 800.5);
    _hist.ntripl1[1]    = Tfs->make<TH1F>("ntripl1Pos"    , "number of triplets fit Helix, Pos "   , 401, -0.5, 800.5);

    _hist.xyniter[0]    = Tfs->make<TH1F>("xyniterNeg"   , "number of iterations XY fit, Neg "   , 401, -0.5, 800.5);
    _hist.xyniter[1]    = Tfs->make<TH1F>("xyniterPos"   , "number of iterations XY fit, Pos "   , 401, -0.5, 800.5);

    _hist.fzniter[0]    = Tfs->make<TH1F>("fzniterNeg"   , "number of iterations z#phi fit, Neg "   , 401, -0.5, 800.5);
    _hist.fzniter[1]    = Tfs->make<TH1F>("fzniterPos"   , "number of iterations z#phi fit, Pos "   , 401, -0.5, 800.5);

    _hist.niter  [0]    = Tfs->make<TH1F>("niterNeg"     , "number of iterations helix fit, Neg "   , 401, -0.5, 800.5);
    _hist.niter  [1]    = Tfs->make<TH1F>("niterPos"     , "number of iterations helix fit, Pos "   , 401, -0.5, 800.5);

    _hist.nXYSh  [0]    = Tfs->make<TH1F>("nXYShNeg"     , "number of strawHits from the circle fit, Neg "   , 201, -0.5, 200.5);
    _hist.nZPhiSh[0]    = Tfs->make<TH1F>("nZPhiShNeg"   , "number of strawHits from the z#phi fit, Neg "   , 201, -0.5, 200.5);
    _hist.nXYSh  [1]    = Tfs->make<TH1F>("nXYShPos"     , "number of strawHits from the circle fit, Pos "   , 201, -0.5, 200.5);
    _hist.nZPhiSh[1]    = Tfs->make<TH1F>("nZPhiShPos"   , "number of strawHits from the z#phi fit, Pos "   , 201, -0.5, 200.5);

    _hist.rinit  [0]    = Tfs->make<TH1F>("rinitNeg"       , "helix radius fitCircle, Neg; r [mm]"      , 401, -0.5, 400.5);

    _hist.rinit  [1]    = Tfs->make<TH1F>("rinitPos"       , "helix radius fitCircle, Pos; r [mm]"      , 401, -0.5, 400.5);

    _hist.radius[0]     = Tfs->make<TH1F>("radius0"  , "helix radius; r [mm]"                      , 401, -0.5, 400.5);
    _hist.radius[1]     = Tfs->make<TH1F>("radius1"  , "helix radius nhits > 15; r [mm]"           , 401, -0.5, 400.5);
    _hist.pT [0]        = Tfs->make<TH1F>("pT0"      , "transverse momentum; pT [MeV/c]"           , 400, -0.5, 200.5);
    _hist.p  [0]        = Tfs->make<TH1F>("p0"       , "momentum; p [MeV/c]"                       , 400, -0.5, 200.5);
    _hist.pT [1]        = Tfs->make<TH1F>("pT1"      , "transverse momentum nhits > 15; pT [MeV/c]", 400, -0.5, 200.5);
    _hist.p  [1]        = Tfs->make<TH1F>("p1"       , "momentum nhits > 15; p [MeV/c]"            , 400, -0.5, 200.5);
    _hist.chi2XY[0]     = Tfs->make<TH1F>("chi2XY0"  , "normalized chi2-XY"                        , 200, 0., 20.);
    _hist.chi2XY[1]     = Tfs->make<TH1F>("chi2XY1"  , "normalized chi2-XY: nhits>15"              , 200, 0., 20.);
    _hist.chi2ZPhi[0]   = Tfs->make<TH1F>("chi2ZPhi0", "normalized chi2-ZPhi"                      , 200, 0., 20.);
    _hist.chi2ZPhi[1]   = Tfs->make<TH1F>("chi2ZPhi1", "normalized chi2-ZPhi: nhits>15"            , 200, 0., 20.);
    _hist.dr  [0]       = Tfs->make<TH1F>("dr0"      , "dr; r - r_{no-target} [mm]"                , 800, -200, 200);
    _hist.dr  [1]       = Tfs->make<TH1F>("dr1"      , "dr: nhits>15; r - r_{no-target} [mm]"      , 800, -200, 200);
    // _hist.shmeanr  [0]  = Tfs->make<TH1F>("shmeanr0" , "straw hit mean radius; r_{sh} [mm]"          , 1800, 0, 900);
    // _hist.shmeanr  [1]  = Tfs->make<TH1F>("shmeanr1" , "straw hit mean radius: nhits>15; r_{sh} [mm]", 1800, 0, 900);
    _hist.chi2d_helix[0]= Tfs->make<TH1F>("chi2dhel0" , "global chi2d; #chi^{2}/ndof"                   , 100, 0, 10); 
    _hist.chi2d_helix[1]= Tfs->make<TH1F>("chi2dhel1" , "global chi2d: nhits>15; #chi^{2}/ndof"         , 100, 0, 10); 

    return 0;
  }


//-----------------------------------------------------------------------------
// Mode is not used
//-----------------------------------------------------------------------------
  int RobustHelixFinderDiag::fillHistograms(void* Data, int Mode) {

    _data = (Data_t*) Data;

//-----------------------------------------------------------------------------
// fill helix-level histograms
//-----------------------------------------------------------------------------
    
    _hist.nTimePeaks->Fill(_data->nTimePeaks);

    int   nhelicities(2);

    for (int k=0; k<nhelicities; ++k){
      _hist.nseeds[k]->Fill(_data->nseeds[k]);
      _hist.nseeds[k]->Fill(_data->nseeds[k]);

      for (int i=0; i<_data->nseeds[k]; i++) {
	_hist.ntclhits   [k]->Fill(_data->ntclhits[k][i]   );
	if (k==1){
	  _hist.ntripl0 ->Fill(_data->ntriplet0[k][i]      );
	}
	_hist.nhits      [k]->Fill(_data->nhits[k][i]      );
	_hist.rinit      [k]->Fill(_data->rinit[k][i]      );    

	_hist.lambda0    [k]->Fill(fabs(_data->lambda0[k][i]));
	_hist.lambda1    [k]->Fill(fabs(_data->lambda1[k][i]));

	_hist.ntripl1    [k]->Fill(_data->ntriplet1[k][i]  );
	_hist.xyniter    [k]->Fill(_data->xyniter[k][i]    );
	_hist.fzniter    [k]->Fill(_data->fzniter[k][i]    );

	_hist.nXYSh      [k]->Fill(_data->nXYSh[k][i]      );
	_hist.nZPhiSh    [k]->Fill(_data->nZPhiSh[k][i]   );

	_hist.p          [k]->Fill(_data->p[k][i]          );
	_hist.pT         [k]->Fill(_data->pT[k][i]         );
	_hist.radius     [k]->Fill(_data->radius[k][i]     );
	_hist.chi2XY     [k]->Fill(_data->chi2XY[k][i]     );
	_hist.chi2ZPhi   [k]->Fill(_data->chi2ZPhi[k][i]   );
	_hist.dr         [k]->Fill(_data->dr[k][i]         );
	_hist.chi2d_helix[k]->Fill(_data->chi2d_helix[k][i]);
      
      }
    }

    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::RobustHelixFinderDiag)

#endif
