#include "MakeCrvRecoPulses.hh"
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>

namespace mu2eCrv
{

MakeCrvRecoPulses::MakeCrvRecoPulses(double calibrationFactor, double pedestal, bool usePulseArea) : 
                                     _calibrationFactor(calibrationFactor), _pedestal(pedestal),
                                     _usePulseArea(usePulseArea)
{}

void MakeCrvRecoPulses::SetWaveform(const std::vector<unsigned int> &waveform, unsigned int startTDC, double binWidth)
{
  _pulseTimes.clear();
  _pulseHeights.clear();
  _pulseWidths.clear();
  _pulseFitChi2s.clear();
  _fitParams0.clear();
  _fitParams1.clear();
  _fitParams2.clear();
  _t1s.clear();
  _t2s.clear();
  _PEs.clear();
  _LEtimes.clear();

  //find the maxima
  int nBins = static_cast<int>(waveform.size());
  std::vector<int> peaks;
  for(int bin=0; bin<=nBins; bin++) 
  {
    if(bin>1 && bin<nBins-1)  //don't search for peaks too close to the sample start or end
    {
      if(waveform[bin-1]<waveform[bin] && waveform[bin+1]<waveform[bin]) peaks.push_back(bin);
    }
  }

  for(size_t i=0; i<peaks.size(); i++)
  {
  //select a range of up to 4 points before and after the maximum point
  //-find up to 5 points before and after the maximum point for which the waveform is stricly decreasing
  //-remove 1 point on each side. this removes potentially "bad points" belonging to a second pulse (i.e. in double pulses)
    int maxBin = peaks[i];
    int startBin=maxBin;
    int endBin=maxBin;
    for(int bin=maxBin-1; bin>=0 && bin>=maxBin-5; bin--)
    {
      if(waveform[bin]<=waveform[bin+1]) startBin=bin;
      else break;
    }
    for(int bin=maxBin+1; bin<nBins && bin<=maxBin+5; bin++)
    {
      if(waveform[bin]<=waveform[bin-1]) endBin=bin;
      else break;
    }
    if(maxBin-startBin>1) startBin++;
    if(endBin-maxBin>1) endBin--;

    double t1=(startTDC+startBin)*binWidth;
    double t2=(startTDC+endBin)*binWidth;

    //fill the graph
    TGraph g;
    for(int bin=startBin; bin<=endBin; bin++) 
    {
      double t=(startTDC+bin)*binWidth;
      double v=waveform[bin]-_pedestal;
      g.SetPoint(g.GetN(), t, v);
    }

    //set the fit function
    TF1 f("peakfitter","[0]*(TMath::Exp(-(x-[1])/[2]-TMath::Exp(-(x-[1])/[2])))");
    f.SetParameter(0, (waveform[maxBin]-_pedestal)*2.718);
    f.SetParameter(1, (startTDC+maxBin)*binWidth);
    f.SetParameter(2, 15.0);

    //do the fit
    TFitResultPtr fr = g.Fit(&f,"NQS");
    if(!fr->IsValid()) continue;

    double fitParam0 = fr->Parameter(0);
    double fitParam1 = fr->Parameter(1);
    double fitParam2 = fr->Parameter(2);
    if(fitParam0<=0 || fitParam2<=0) continue;

    int    PEs          = lrint(fitParam0*fitParam2 / _calibrationFactor);
    double pulseTime    = fitParam1;
    double pulseHeight  = fitParam0/2.718;    //=fitParam0/e
    double pulseWidth   = fitParam2*1.283;    //=fitParam2*pi/sqrt(6)  // =standard deviation of the Gumbel distribution
    double pulseFitChi2 = fr->Chi2();

    if(!_usePulseArea) PEs = lrint((waveform[maxBin]-_pedestal) / _calibrationFactor);

    double LEtime=f.GetX(0.2*pulseHeight,pulseTime-50,pulseTime);   //i.e. at 20% of pulse height

    _pulseTimes.push_back(pulseTime);
    _pulseHeights.push_back(pulseHeight);
    _pulseWidths.push_back(pulseWidth);
    _pulseFitChi2s.push_back(pulseFitChi2);
    _fitParams0.push_back(fitParam0);
    _fitParams1.push_back(fitParam1);
    _fitParams2.push_back(fitParam2);
    _t1s.push_back(t1);
    _t2s.push_back(t2);
    _PEs.push_back(PEs);
    _LEtimes.push_back(LEtime);
  }
}

unsigned int MakeCrvRecoPulses::GetNPulses()
{
  return _PEs.size();
}

int MakeCrvRecoPulses::GetPEs(int pulse)
{
  int n = _PEs.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _PEs[pulse];
}

double MakeCrvRecoPulses::GetPulseTime(int pulse)
{
  int n = _pulseTimes.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseTimes[pulse];
}

double MakeCrvRecoPulses::GetPulseHeight(int pulse)
{
  int n = _pulseHeights.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseHeights[pulse];
}

double MakeCrvRecoPulses::GetPulseWidth(int pulse)
{
  int n = _pulseWidths.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseWidths[pulse];
}

double MakeCrvRecoPulses::GetPulseFitChi2(int pulse)
{
  int n = _pulseFitChi2s.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _pulseFitChi2s[pulse];
}

double MakeCrvRecoPulses::GetFitParam0(int pulse)
{
  int n = _fitParams0.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _fitParams0[pulse];
}

double MakeCrvRecoPulses::GetFitParam1(int pulse)
{
  int n = _fitParams1.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _fitParams1[pulse];
}

double MakeCrvRecoPulses::GetFitParam2(int pulse)
{
  int n = _fitParams2.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _fitParams2[pulse];
}

double MakeCrvRecoPulses::GetT1(int pulse)
{
  int n = _t1s.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _t1s[pulse];
}

double MakeCrvRecoPulses::GetT2(int pulse)
{
  int n = _t2s.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _t2s[pulse];
}

double MakeCrvRecoPulses::GetLEtime(int pulse)
{
  int n = _LEtimes.size();
  if(pulse<0 || pulse>=n) throw std::logic_error("invalid pulse number");
  return _LEtimes[pulse];
}

}

