//
// Class to perform Cosmic fit
// Original author: Sophie Middleton 
//
// $Id: CosmicTrackFit.hh 
// $Author: sophie $ 
// $Date: 2018/01/12 18:56:10 $
//
#ifndef TrkReco_CosmicTrackFit_HH
#define TrkReco_CosmicTrackFit_HH

// framework
#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

// data
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/StraightTrackSeed.hh"

// Mu2e objects
#include "TrkReco/inc/CosmicTrackFinderData.hh"
#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"

//C++
#include <vector>
#include <utility>
#include <string>
#include <math.h>
#include <cmath>
#include <algorithm>

//ROOT
#include "TMatrixD.h"
namespace mu2e 
{
    class Tracker;
    class CosmicTrackFit
    {
     public:
                explicit CosmicTrackFit(fhicl::ParameterSet const&);
    		virtual ~CosmicTrackFit();

                bool initCosmicTrack(CosmicTrackFinderData& TrackData);
                CosmicTrack* InitLine(const ComboHit *FirstP1, const ComboHit *LastP2, CosmicTrack* line);
                //Step 1: Begin Fit- initializes the fit routine:
                void BeginFit(CosmicTrackFinderData& TrackData);
                
                //Step 2: RunFitChi2-holds the functions to find initial line, update, refine and add in drift
                void RunFitChi2(CosmicTrackFinderData& trackData);
		//Step 3: Fit All - finds the chi-squared anf line information when all hits in track taken into account. This will be the initial chi-squared value.
		CosmicTrack* FitAll(CosmicTrackFinderData& trackData,CosmicTrack* track, int WeightMode);
		
		void MulitpleTrackResolver(CosmicTrackFinderData& trackData,CosmicTrack* track);

		
                void UpdateFitErrors(std::vector<double> x, std::vector<double> y,std::vector<double> z, std::vector<double> err, StraightTrack* track,TMatrixD cov_x, std::vector<XYZVec> maj,std::vector<XYZVec> min);

                bool goodTrack(CosmicTrack* track);
                float pointToLineDCA(CosmicTrack* track, StrawHit hit);
		void DriftCorrection(CosmicTrackFinderData& trackData);
		

                float evalWeightXY  (const ComboHit& Hit, const CosmicTrack& track);
                const Tracker*            _tracker;
    		void  setTracker    (const Tracker*    Tracker) { _tracker     = Tracker; }
                
		
	private:
		
                //list function:
               
  		bool use(ComboHit const&) const;
    		void setOutlier(ComboHit&) const;
                float hitWeight(ComboHit const& hhit) const;
                int _dim; //dimensions of fit
		int _diag;
    		int _debug;		    // debug level
                StrawHitFlag _useflag, _dontuseflag;
	        int _minresid;
    		int      _minnsh;  // minimum # of StrawHits
    		unsigned _minCHHits; // minimum # of CH hits - should be at least 2 for fit to work....
		unsigned _n_outliers; //number of significant outliers/number of hits in track....helps with multiplicity(?)
    		float _maxresid; // max allowed pull which a hit can have to be classed as "good"
    		unsigned _maxniter; // maxium # of iterations to global minimum   
                
    		float _maxd; // maximum distance in hits to begin fit
		float _maxDOCA; //max allowed DOCA to allow hit into fit
    		float _maxchi2; //maximum allowed chi2
		
    
  };//end Fit class

}
#endif
