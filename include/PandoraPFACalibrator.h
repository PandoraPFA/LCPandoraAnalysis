#ifndef PandoraPFACalibrator_h
#define PandoraPFACalibrator_h

#include "marlin/Processor.h"

#include <vector>

class TH1F;
class TH2F;

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
class PandoraPFACalibrator : public Processor
{
public:
    virtual Processor*  newProcessor()
    {
        return new PandoraPFACalibrator;
    }

    PandoraPFACalibrator() ;

    /**
     *   @brief Called at the begin of the job before anything is read. Use to initialize the processor, e.g. book histograms.
     */
    virtual void init() ;

    /**
     *  @brief  Called for every run.
     */
    virtual void processRunHeader( LCRunHeader* run ) ;

    /**
     *  @brief  Called for every event - the working horse.
     */
    virtual void processEvent( LCEvent * evt ) ; 

    /**
     *  @brief  Check
     */
    virtual void check( LCEvent * evt ) ; 

    /**
     *  @brief  Called after data processing for clean up.
     */
    virtual void end() ;

private:
    typedef std::vector<std::string> StringVector;                  ///< 
    StringVector    m_inputMCParticleCollections;                   ///< 
    StringVector    m_ecalBarrelCollections;                        ///< 
    StringVector    m_ecalEndCapCollections;                        ///< 
    StringVector    m_hcalCollections;                              ///< 
    StringVector    m_muonCollections;                              ///< 
    StringVector    m_lcalCollections;                              ///< 
    StringVector    m_bcalCollections;                              ///< 
    StringVector    m_lhcalCollections;                             ///< 
    std::string     m_particleCollectionName;                       ///< 
    std::string     m_rootFile;                                     ///< 

    int             m_nRun;                                         ///< 
    int             m_nEvt;                                         ///< 

    float           m_ecalToMIP;                                    ///< 
    float           m_hcalToMIP;                                    ///< 
    float           m_muonToMIP;                                    ///< 

    float           m_ecalToEMGeVCalibration;                       ///< 
    float           m_hcalToHadGeVCalibration;                      ///< 

    float           m_ecalToHadGeVCalibrationBarrel;                ///< 
    float           m_ecalToHadGeVCalibrationEndCap;                ///< 
    float           m_hcalToEMGeVCalibration;                       ///< 

    float           m_maxHCalHitHadronicEnergy;                     ///< 

    float           m_zOfEndCap;                                    ///< 
    float           m_cosTheta;                                     ///< 
    float           m_cosThetaR;                                    ///< 
    float           m_cosThetaX;                                    ///< 
    float           m_zCoG;                                         ///< 
    float           m_x;                                            ///< 
    float           m_y;                                            ///< 
    std::string     m_detectorName;                                 ///< 

    TH1F           *m_PFA;                                          ///< 
    TH1F           *m_PFAB;                                         ///< 
    TH2F           *m_PFAVsCosTheta;                                ///< 
    TH2F           *m_PFAVsCosThetaR;                               ///< 
    TH2F           *m_PFAVsZCoG;                                    ///< 
    TH2F           *m_PFAVsCosThetaX;                               ///< 

    TH1F           *m_PFAE;                                         ///< 
    TH1F           *m_PFAH;                                         ///< 
    TH1F           *m_PFAM;                                         ///< 

    TH2F           *m_XvsY;                                         ///< 

    TH1F           *m_EcalEnergy;                                   ///< 
    TH1F           *m_HcalEnergy;                                   ///< 
    TH1F           *m_MuonEnergy;                                   ///< 
    TH1F           *m_LcalEnergy;                                   ///< 
    TH1F           *m_CalEnergy;                                    ///< 

    TH2F           *m_EcalBarrelHcalEnergyEM;                       ///< 
    TH2F           *m_EcalEndCapHcalEnergyEM;                       ///< 
    TH2F           *m_EcalBarrelHcalEnergyHAD;                      ///< 
    TH2F           *m_EcalEndCapHcalEnergyHAD;                      ///< 

    TH1F           *m_CalEnergyE;                                   ///< 
    TH1F           *m_CalEnergyH;                                   ///< 
    TH1F           *m_CalEnergyM;                                   ///< 

    TH2F           *m_CalEnergyVsCosTheta;                          ///< 
    TH2F           *m_CalEnergyVsCosThetaR;                         ///< 

    TH1F           *m_EcalBarrelEnergyByLayer;                      ///< 
    TH1F           *m_EcalEndCapEnergyByLayer;                      ///< 

    TH1F           *m_EcalBarrelMIP;                                ///< 
    TH1F           *m_EcalEndCapMIP;                                ///< 
    TH1F           *m_HcalMIP;                                      ///< 
    TH1F           *m_MuonMIP;                                      ///< 

    TH1F           *m_EcalBarrelMIPcorr;                            ///< 
    TH1F           *m_EcalEndCapMIPcorr;                            ///< 
    TH1F           *m_HcalMIPcorr;                                  ///< 
    TH1F           *m_MuonMIPcorr;                                  ///< 

    TH1F           *m_CosT;                                         ///< 
    TH1F           *m_PhotonCosT;                                   ///< 
};

#endif
