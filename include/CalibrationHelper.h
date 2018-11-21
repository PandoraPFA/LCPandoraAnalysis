/**
 *  @file   PandoraAnalysis/include/CalibrationHelper.h
 * 
 *  @brief  Header file for the calibration helper class.
 */

#ifndef CALIBRATION_HELPER_H
#define CALIBRATION_HELPER_H 1

#include "EVENT/LCStrVec.h"

#include "marlin/Processor.h"

#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"


#include <set>
#include <vector>

class TFile;
class TH1F;
class TTree;

//------------------------------------------------------------------------------------------------------------------------------------------

namespace pandora_analysis
{

/**
 *  @brief  CalibrationHelper class
 */
class CalibrationHelper
{
public:
    /**
     *  @brief  Settings class
     */
    class Settings
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Settings();

        LCStrVec    m_eCalCollections{};                      ///< Input calorimeter hit collection names
        LCStrVec    m_hCalCollections{};                      ///< Input calorimeter hit collection names
        LCStrVec    m_muonCollections{};                      ///< Input calorimeter hit collection names
        LCStrVec    m_bCalCollections{};                      ///< Input calorimeter hit collection names
        LCStrVec    m_lHCalCollections{};                     ///< Input calorimeter hit collection names
        LCStrVec    m_lCalCollections{};                      ///< Input calorimeter hit collection names

        LCStrVec    m_eCalCollectionsSimCaloHit{};            ///< Input simcalorimeter hit collection names
        LCStrVec    m_hCalBarrelCollectionsSimCaloHit{};      ///< Input simcalorimeter hit collection names
        LCStrVec    m_hCalEndCapCollectionsSimCaloHit{};      ///< Input simcalorimeter hit collection names
        LCStrVec    m_hCalOtherCollectionsSimCaloHit{};       ///< Input simcalorimeter hit collection names
        LCStrVec    m_eCalBarrelCollectionsSimCaloHit{};      ///< Input simcalorimeter hit collection names
        LCStrVec    m_eCalEndCapCollectionsSimCaloHit{};      ///< Input simcalorimeter hit collection names
        LCStrVec    m_eCalOtherCollectionsSimCaloHit{};       ///< Input simcalorimeter hit collection names
        LCStrVec    m_muonCollectionsSimCaloHit{};            ///< Input simcalorimeter hit collection names
        LCStrVec    m_bCalCollectionsSimCaloHit{};            ///< Input simcalorimeter hit collection names
        LCStrVec    m_lHCalCollectionsSimCaloHit{};           ///< Input simcalorimeter hit collection names
        LCStrVec    m_lCalCollectionsSimCaloHit{};            ///< Input simcalorimeter hit collection names

        int         m_nBinsMuonCaloHitEnergyHist;             ///< Number of bins in MuonDirectionCorrectedCaloHitEnergy histogram 
        float       m_xUpperValueMuonCaloHitEnergyHist;       ///< Upper value of x-range of MuonDirectionCorrectedCaloHitEnergy histogram 
    };

    typedef std::vector<const EVENT::ReconstructedParticle *> ParticleVector;
    
    CalibrationHelper(const CalibrationHelper&) = delete;
    CalibrationHelper& operator=(const CalibrationHelper&) = delete;

    /**
     *  @brief  Constructor
     * 
     *  @param  settings the settings
     */
    CalibrationHelper(const Settings &settings);

    /**
     *  @brief  Destructor
     */
    ~CalibrationHelper();

    /**
     *  @brief  Produces the calibration data
     * 
     *  @param  pLCEvent the lc event
     *  @param  particleVector reconstructed particle vector for the lc event
     *  @param  nPfoTargetsTotal number of total targets
     *  @param  nPfoTargetsTracks number of target tracks
     *  @param  nPfoTargetsNeutralHadrons of target neutral hadrons
     *  @param  nPfoTargetsPhotons of target photons
     *  @param  pfoTargetsEnergyTotal total targets energy
     */
    void Calibrate(const EVENT::LCEvent *pLCEvent, const ParticleVector &particleVector, const int nPfoTargetsTotal, const int nPfoTargetsTracks, const int nPfoTargetsNeutralHadrons, const int nPfoTargetsPhotons, const float pfoTargetsEnergyTotal);

    /**
     *  @brief  Set branch addresses for calibration variables
     * 
     *  @param  pTTree tree to set branch addresses to
     */
    void SetBranchAddresses(TTree *pTTree);

    /**
     *  @brief  Create calibration histograms
     */
    void CreateHistograms();

    /**
     *  @brief  Set directory for calibration histograms
     * 
     *  @param  pTFile directory to set histograms to
     */
    void SetHistogramDirectories(TFile *pTFile);

    /**
     *  @brief  Write calibration histograms to tree
     */
    void WriteHistograms();

    /**
     *  @brief  Clear calibration member variables
     */
    void Clear();

private:
    /**
     *  @brief  Get the minimum number of HCal layers between any calo hit in the ParticleVector and the edge of the detector
     * 
     *  @param  pParticleVector to be examined
     */
    int GetMinNHCalLayersFromEdge(const ParticleVector &pParticleVector) const;
    
    /**
     *  @brief  Get the calo hit energy in the last hcal layers of the detector 
     * 
     *  @param  pLCEvent the lc event
     *  @param  collectionNames the collection to be read
     */
    float GetHCalCaloHitEnergyFromEdge(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames) const;

    /**
     *  @brief  Read and save the calorimeter hit information for a specific collection
     * 
     *  @param  pLCEvent the lc event
     *  @param  collectionNames the collection to be read
     *  @param  hitEnergySum the sum of the hit energy for the collection
     */
    void ReadCaloHitEnergies(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames, float &hitEnergySum) const;

    /**
     *  @brief  Read and save the simcalorimeter hit information for a specific collection
     * 
     *  @param  pLCEvent the lc event
     *  @param  collectionNames the collection to be read
     *  @param  hitEnergySum the sum of the hit energy for the collection
     */
    void ReadSimCaloHitEnergies(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames, float &hitEnergySum) const;

    /**
     *  @brief  Add the direction corrected SimCaloHits and the direction corrections for a particular collection to separate histograms
     * 
     *  @param  pLCEvent a collection in the lc event
     *  @param  collectionNames a collection in the lc event
     *  @param  Barrel, Endcap, Other for the HCal collection
     *  @param  pTH1F_Energy histogram to store direction corrected SimCaloHits
     *  @param  pTH1F_Direction_Correction histogram to store direction corrections
     */
    void AddSimCaloHitEntries(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames, const unsigned int SimCaloHit_DC, TH1F *pTH1F) const;

    /**
     *  @brief  Add the direction corrected calo hits for a particular collection to a histogram
     * 
     *  @param  pLCEvent a collection in the lc event
     *  @param  collectionNames a collection in the lc event
     *  @param  pTH1F histogram to store direction corrected calo hits
     */
    void AddDirectionCorrectedCaloHitEntries(const EVENT::LCEvent *pLCEvent, const EVENT::LCStrVec &collectionNames, TH1F *pTH1F) const;

    /**
     *  @brief  Get the minimum number of HCal layers between a calo hit and the edge of the detector
     * 
     *  @param  pCaloHit to be examined
     */
    int GetNHCalLayersFromEdge(const EVENT::CalorimeterHit *const pCaloHit) const;

    /**
     *  @brief  Get the radial distance between a calo hit in the HCal and the edge of the HCal
     * 
     *  @param  pCaloHit to be examined
     *  @param  symmetryOrder of the region of the HCal the hit is in i.e. Barrel, Endcap, Ring
     *  @param  phi0 of the region of the HCal the hit is in i.e. Barrel, Endcap, Ring
     */
    float GetMaximumRadius(const EVENT::CalorimeterHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const;

    /**
     *  @brief A helper function to access geometry information via DD4HEP 
     *
     *  @param includeFlag calorimeter propereties to include
     *  @param excludeFlag calorimeter propereties to exclude
     */
    dd4hep::rec::LayeredCalorimeterData *GetExtension(unsigned int includeFlag, unsigned int excludeFlag = 0) const;

    const Settings  m_settings;                                    ///< The calibration helper settings

    int             m_pfoMinHCalLayerToEdge;                       ///< GetNHCalLayersFromEdge
    int             m_hCalCaloHitEnergyFromEdge;                   ///< GetHCalCaloHitEnergyFromEdge

    float           m_totalCaloHitEnergy;                          ///< Sum outputs from ReadCaloHitEnergies
    float           m_eCalTotalCaloHitEnergy;                      ///< ReadCaloHitEnergies
    float           m_hCalTotalCaloHitEnergy;                      ///< ReadCaloHitEnergies
    float           m_muonTotalCaloHitEnergy;                      ///< ReadCaloHitEnergies
    float           m_bCalTotalCaloHitEnergy;                      ///< ReadCaloHitEnergies
    float           m_lHCalTotalCaloHitEnergy;                     ///< ReadCaloHitEnergies
    float           m_lCalTotalCaloHitEnergy;                      ///< ReadCaloHitEnergies

    float           m_totalSimCaloHitEnergy;                       ///< Sum outputs from ReadSimCaloHitEnergies
    float           m_eCalTotalSimCaloHitEnergy;                   ///< ReadSimCaloHitEnergies
    float           m_hCalTotalSimCaloHitEnergy;                   ///< ReadSimCaloHitEnergies
    float           m_muonTotalSimCaloHitEnergy;                   ///< ReadSimCaloHitEnergies
    float           m_bCalTotalSimCaloHitEnergy;                   ///< ReadSimCaloHitEnergies
    float           m_lHCalTotalSimCaloHitEnergy;                  ///< ReadSimCaloHitEnergies
    float           m_lCalTotalSimCaloHitEnergy;                   ///< ReadSimCaloHitEnergies

    TH1F           *m_hECalDirectionCorrectedCaloHitEnergy;        ///< AddDirectionCorrectedCaloHitEntries
    TH1F           *m_hHCalDirectionCorrectedCaloHitEnergy;        ///< AddDirectionCorrectedCaloHitEntries
    TH1F           *m_hMuonDirectionCorrectedCaloHitEnergy;        ///< AddDirectionCorrectedCaloHitEntries

    TH1F           *m_hHCalBarrelDirectionCorrectedSimCaloHit;     ///< AddSimCaloHitEntries setting 0, pass HCal Barrel sim calo hit collections
    TH1F           *m_hHCalEndCapDirectionCorrectedSimCaloHit;     ///< AddSimCaloHitEntries setting 0, pass HCal EndCap sim calo hit collections
    TH1F           *m_hHCalOtherDirectionCorrectedSimCaloHit;      ///< AddSimCaloHitEntries setting 0, pass HCal Other sim calo hit collections

    TH1F           *m_hECalDirectionCorrectedSimCaloHit;           ///< AddSimCaloHitEntries setting 0, pass all ECal sim calo hit collections

    TH1F           *m_hHCalBarrelDirectionCorrectionSimCaloHit;    ///< AddSimCaloHitEntries setting 1, pass HCal Barrel sim calo hit collections
    TH1F           *m_hHCalEndCapDirectionCorrectionSimCaloHit;    ///< AddSimCaloHitEntries setting 1, pass HCal EndCap sim calo hit collections
    TH1F           *m_hHCalOtherDirectionCorrectionSimCaloHit;     ///< AddSimCaloHitEntries setting 1, pass HCal Other sim calo hit collections

    TH1F           *m_hECalBarrelDirectionCorrectionSimCaloHit;    ///< AddSimCaloHitEntries setting 1, pass ECal Barrel sim calo hit collections
    TH1F           *m_hECalEndCapDirectionCorrectionSimCaloHit;    ///< AddSimCaloHitEntries setting 1, pass ECal EndCap sim calo hit collections
    TH1F           *m_hECalOtherDirectionCorrectionSimCaloHit;     ///< AddSimCaloHitEntries setting 1, pass ECal Other sim calo hit collections
};

} // namespace pandora_analysis

#endif // #ifndef CALIBRATION_HELPER_H
