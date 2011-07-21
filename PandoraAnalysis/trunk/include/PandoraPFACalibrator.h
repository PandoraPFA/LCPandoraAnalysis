#ifndef PandoraPFACalibrator_h
#define PandoraPFACalibrator_h

#include <vector>

#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/LCIO.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace lcio;
using namespace marlin;

/** PandoraPFACalibrator is a simple processor to aid calibration of PandoraPFA.
<pre>
There are three sets of constants that need to be calibrated:
  i) MokkaCaloDigi raw hit to GeV
 ii) PandoraPFA energy in GeV -> MIP equivalent conversion
iii) PandoraPFA MIP equivalent -> PFO energy
One might expect that factor iii) is the reciprocal of factor ii). However,
due to the isolation cuts used in PandoraPFA this is not quite the case.


MC Samples:
===========
  o the calibration constants will depend on the Mokka ECAL and HCAL 
    drivers, the GEANT4 physics list, and potentially the version of 
    GEANT4 used to generate the events.
  o to calibrate PandoraPFA generate the following samples with the
    same version of Mokka/Geant4 used for the main event samples:
       o 10000 10 GeV photons at cosine(theta)=0.02 with random values of
         phi. (cos(theta) = 0.02 is used so as to avoid the possible 
         gap between ECAL modules at zero).
       o 5000 10 GeV KL at cosine(theta)=0.02 with random values of
         phi.
       o 1000 10 GeV mu+ at cosine(theta)=0.02 with random values of
         phi.

Calibration Procedure:
======================
i)  Calibration of Mokka hits
  o Using the photon sample, tune the calibration parameters in MokkaCaloDigi 
     <parameter name="CalibrECAL" type="FloatVec">27.3 54.6  </parameter>
    to so that the histogram fCalEnergy is centred on 10 GeV. It is 
    recommended that the calibration constants for the different layers
    should be in the ratio of the absorber thicknesses, i.e. 1:2 for
    LDC01 and 1:3 for LDC00.
  o Repeat with the 10 GeV KL sample to fix 
     <parameter name="CalibrHCAL" type="FloatVec">28.1  </parameter>
ii) PandoraPFA MIP calibration
  o Using the sample of muons tune the parameters:
   <parameter name="ECALMIPcalibration" type="float">209.0  </parameter>
   <parameter name="HCALMIPcalibration" type="float">34.5  </parameter>
    so that fEcalMIPcorr and fHcalMIPcorr peak at approximately one. This 
    set the scale for the energy deposited by a "minimum ionizing particle".
iii) PandoraPFA ECAL/HCAL calibration
  o Using the sample of photons, run PandoraPFA and PandoraPFACalibrator
    and determine
  <parameter name="ECALEMMIPToGeV" type="float">0.004785  </parameter>
  <parameter name="ECALHadMIPToGeV" type="float">0.004785  </parameter>
    so that fPFA peaks at 10 GeV. It is recommended to set the response to
    hadrons to be the same as that for EM showers, i.e. 
          ECALHadMIPToGeV = ECALEMMIPToGeV
  o Using the sample of KLs do that same to obtain 
   <parameter name="HCALEMMIPToGeV" type="float">0.023  </parameter>
   <parameter name="HCALHadMIPToGeV" type="float">0.023  </parameter>
</pre>    
*/
class PandoraPFACalibrator : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new PandoraPFACalibrator ; }
  
  
  PandoraPFACalibrator() ;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  
  
  

 private:
  
  std::vector<std::string> _inputMCParticleCollections;
  std::vector<std::string> _ecalBarrelCollections;
  std::vector<std::string> _ecalEndCapCollections;
  std::vector<std::string> _hcalCollections;
  std::vector<std::string> _lcalCollections;
  std::vector<std::string> _bcalCollections;
  std::vector<std::string> _lhcalCollections;
  std::string _particleCollectionName;
  std::string _rootFile;

  int _nRun ;
  int _nEvt ;
  int   _digitalHcal;
  float _ecalToMIP;
  float _hcalToMIP;
  float _ecalMIPThreshold;
  float _hcalMIPThreshold;
  float _ecalEMMIPToGeV;
  float _hcalEMMIPToGeV;
  float _ecalBarrelHadMIPToGeV;
  float _ecalEndCapHadMIPToGeV;
  float _hcalHadMIPToGeV;
  float _zOfEndCap;
  float _cosTheta;
  float _cosThetaR;
  float _cosThetaX;
  float _zCoG;
  float _x;
  float _y;
  float _zmean;
  std::string _detectorName;

  TH1F* fPFA;
  TH1F* fPFAB;
  TH1F* fCalEnergy;
  TH1F* fCalEnergyH;
  TH2F* fCalEnergyVsCosTheta;
  TH2F* fCalEnergyVsCosThetaR;
  TH2F* fPFAVsCosTheta;
  TH2F* fPFAVsCosThetaX;
  TH2F* fPFAVsCosThetaR;
  TH2F* fPFAVsZCoG;
  TH2F* fXvsY;
  TH1F* fPFAH;
  TH1F* fEcalEnergy;
  TH2F* fEcalBarrelHcalEnergy;
  TH2F* fEcalEndCapHcalEnergy;
  TH1F* fHcalEnergy;
  TH1F* fLcalEnergy;
  TH1F* fEcalBarrelMIP;
  TH1F* fEcalEndCapMIP;
  TH1F* fEcal1MIP;
  TH1F* fHcalMIP;
  TH1F* fCosT;
  TH1F* fPhotonCosT;
  TH1F* fEcalBarrelMIPcorr;
  TH1F* fEcalEndCapMIPcorr;
  TH1F* fHcalMIPcorr;
  TH1F* fEcalBarrelEnergyByLayer;
  TH1F* fEcalEndCapEnergyByLayer;

} ;

#endif
