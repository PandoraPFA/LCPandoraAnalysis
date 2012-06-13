/**
 *  @file   PandoraAnalysis/include/PfoAnalysis.h
 * 
 *  @brief  Header file for the pfo analysis class.
 * 
 *  $Log: $
 */

#ifndef PFO_ANALYSIS_H
#define PFO_ANALYSIS_H 1

#include "EVENT/ReconstructedParticle.h"

#include "marlin/Processor.h"

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
     *  @brief  Make quark variables
     */
    void MakeQuarkVariables();

    /**
     *  @brief  Perform pfo analysis
     */
    void PerformPfoAnalysis();

    typedef std::vector<ReconstructedParticle *> ParticleVector;
    typedef std::vector<std::string> StringVector;

    int                 m_nRun;                                 ///< 
    int                 m_nEvt;                                 ///< 

    StringVector        m_inputQuarkParticleCollections;        ///< 
    StringVector        m_inputMCParticleCollections;           ///< 
    StringVector        m_inputParticleCollections;             ///< 
    StringVector        m_inputReclusterMonitoringCollections;  ///< 

    int                 m_printing;                             ///< 
    std::string         m_rootFile;                             ///< 

    ParticleVector      m_pfovec;                               ///< 
    ParticleVector      m_mcpfovec;                             ///< 
    ParticleVector      m_quarkpfovec;                          ///< 

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

    float               m_mcEnergyTotal;                        ///< 
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
    float               m_netEnergyChange;                      ///< 
    float               m_sumModulusEnergyChanges;              ///< 
    float               m_sumSquaredEnergyChanges;              ///< 

    TTree              *m_tree;                                 ///< 
    TH1F               *m_hPfoEnergySum;                        ///< 
    TH1F               *m_hPfoEnergySumL7A;                     ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline marlin::Processor *PfoAnalysis::newProcessor()
{
    return new PfoAnalysis;
}

#endif // #ifndef PFO_ANALYSIS_H
