/**
 *  @file   PandoraAnalysis/include/PfoAnalysis.h
 * 
 *  @brief  Header file for the pfo analysis class.
 * 
 *  $Log: $
 */

#ifndef PFO_ANALYSIS_H
#define PFO_ANALYSIS_H 1

#include "EVENT/LCStrVec.h"

#include "marlin/Processor.h"

#include "CalibrationHelper.h"

#include <set>
#include <vector>

class TFile;
class TH1F;
class TTree;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  PfoAnalysis class
 */
class PfoAnalysis : public marlin::Processor
{
public:
    /**
     *  @brief  Default constructor
     */
    PfoAnalysis();

    /**
     *  @brief  Create new processor
     */
    virtual Processor *newProcessor();

    /**
     *  @brief  Initialize, called at startup
     */
    virtual void init();

    /**
     *  @brief  Process run header
     *
     *  @param  pLCRunHeader the lc run header
     */
    virtual void processRunHeader(lcio::LCRunHeader *pLCRunHeader);

    /**
     *  @brief  Process event, main entry point
     *
     *  @param  pLCEvent the lc event
     */
    virtual void processEvent(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  Checks for event
     *
     *  @param  pLCEvent the lc event
     */
    virtual void check(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  End, called at shutdown
     */
    virtual void end();

private:
    typedef std::vector<const EVENT::ReconstructedParticle *> ParticleVector;
    typedef std::vector<const EVENT::MCParticle*> MCParticleVector;
    typedef std::set<const EVENT::MCParticle*> MCParticleList;
    typedef std::vector<std::string> StringVector;

    /**
     *  @brief  Clear current event details
     */
    void Clear();

    /**
     *  @brief  Extract lcio collections
     * 
     *  @param  pLCEvent the lc event
     */
    void ExtractCollections(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  Apply pfo selection rules, starting with root particles
     * 
     *  @param  pMCParticle the address of a mc particle (initially call this recursive function with a top-level, root particle)
     *  @param  mcPfoCandidates to collect the list of mc pfo candidates
     */
    void ApplyPfoSelectionRules(const EVENT::MCParticle *pMCParticle, MCParticleList &mcPfoCandidates) const;

    /**
     *  @brief  Make quark variables
     * 
     *  @param  pLCEvent the lc event
     */
    void MakeQuarkVariables(EVENT::LCEvent *pLCEvent);

    /**
     *  @brief  Perform pfo analysis
     */
    void PerformPfoAnalysis();

    /**
     *  @brief  Sort mc pfo targets by decreasing energy
     * 
     *  @param  pLhs the address of the first mc pfo target
     *  @param  pRhs the address of the second mc pfo target
     */
    static bool SortPfoTargetsByEnergy(const EVENT::MCParticle *const pLhs, const EVENT::MCParticle *const pRhs);

    int                 m_nRun;                                 ///< 
    int                 m_nEvt;                                 ///< 

    int                 m_nRunSum;                              ///< 
    int                 m_nEvtSum;                              ///< 

    std::string         m_inputPfoCollection;                   ///< 
    std::string         m_mcParticleCollection;                 ///< 

    int                 m_printing;                             ///< 
    std::string         m_rootFile;                             ///< 

    int                 m_lookForQuarksWithMotherZ;             ///< 

    float               m_mcPfoSelectionRadius;                 ///< 
    float               m_mcPfoSelectionMomentum;               ///< 
    float               m_mcPfoSelectionLowEnergyNPCutOff;      ///< 

    ParticleVector      m_pfoVector;                            ///< 
    MCParticleVector    m_pfoTargetVector;                      ///< 

    int                 m_nPfosTotal;                           ///< 
    int                 m_nPfosNeutralHadrons;                  ///< 
    int                 m_nPfosPhotons;                         ///< 
    int                 m_nPfosTracks;                          ///< 
    float               m_pfoEnergyTotal;                       ///< 
    float               m_pfoEnergyNeutralHadrons;              ///< 
    float               m_pfoEnergyPhotons;                     ///< 

    float               m_pfoEnergyTracks;                      ///< 
    float               m_pfoECalToEmEnergy;                    ///< 
    float               m_pfoECalToHadEnergy;                   ///< 
    float               m_pfoHCalToEmEnergy;                    ///< 
    float               m_pfoHCalToHadEnergy;                   ///< 
    float               m_pfoMuonToEnergy;                      ///< 
    float               m_pfoOtherEnergy;                       ///< 

    float               m_pfoMassTotal;                         ///< 

    typedef std::vector<float> FloatVector;
    FloatVector         m_pfoEnergies;                          ///< 
    FloatVector         m_pfoPx;                                ///< 
    FloatVector         m_pfoPy;                                ///< 
    FloatVector         m_pfoPz;                                ///< 
    FloatVector         m_pfoCosTheta;                          ///< 

    FloatVector         m_pfoTargetEnergies;                    ///< 
    FloatVector         m_pfoTargetPx;                          ///< 
    FloatVector         m_pfoTargetPy;                          ///< 
    FloatVector         m_pfoTargetPz;                          ///< 
    FloatVector         m_pfoTargetCosTheta;                    ///< 

    typedef std::vector<int> IntVector;
    IntVector           m_pfoPdgCodes;                          ///< 
    IntVector           m_pfoTargetPdgCodes;                    ///< 

    int                 m_nPfoTargetsTotal;                     ///< 
    int                 m_nPfoTargetsNeutralHadrons;            ///< 
    int                 m_nPfoTargetsPhotons;                   ///< 
    int                 m_nPfoTargetsTracks;                    ///< 

    float               m_pfoTargetsEnergyTotal;                ///< 
    float               m_pfoTargetsEnergyNeutralHadrons;       ///< 
    float               m_pfoTargetsEnergyPhotons;              ///< 
    float               m_pfoTargetsEnergyTracks;               ///< 

    float               m_mcEnergyENu;                          ///< 
    float               m_mcEnergyFwd;                          ///< 
    float               m_eQQ;                                  ///< 
    float               m_eQ1;                                  ///< 
    float               m_eQ2;                                  ///< 
    float               m_costQQ;                               ///< 
    float               m_costQ1;                               ///< 
    float               m_costQ2;                               ///< 
    float               m_mQQ;                                  ///< 
    float               m_thrust;                               ///< 
    int                 m_qPdg;                                 ///< 

    TFile              *m_pTFile;                               ///< 
    TTree              *m_pTTree;                               ///< 
    TH1F               *m_hPfoEnergySum;                        ///< 
    TH1F               *m_hPfoEnergySumL7A;                     ///< 

    int                                             m_collectCalibrationDetails;        ///< Whether to collect calibration details
    pandora_analysis::CalibrationHelper            *m_pCalibrationHelper;               ///< The address of the calibration helper
    pandora_analysis::CalibrationHelper::Settings   m_calibrationHelperSettings;        ///< The calibration helper settings
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline marlin::Processor *PfoAnalysis::newProcessor()
{
    return new PfoAnalysis;
}

#endif // #ifndef PFO_ANALYSIS_H
