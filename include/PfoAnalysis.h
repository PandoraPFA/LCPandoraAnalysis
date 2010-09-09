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

class MCPfo;
class TH1F;

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
    virtual void processEvent(lcio::LCEvent *pLCEvent);

    /**
     *  @brief  Checks for event
     *
     *  @param  pLCEvent the lc event
     */
    virtual void check(lcio::LCEvent *pLCEvent);

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
     *  @brief  Make quark variables
     */
    void MakeQuarkVariables() ;

    /**
     *  @brief  Analyse pfa performance
     */
    void AnalysePFAPerformance();

    typedef std::vector<ReconstructedParticle *> ParticleVector;
    typedef std::vector<std::string> StringVector;

    StringVector        m_inputQuarkParticleCollections;        ///< 
    StringVector        m_inputMCParticleCollections;           ///< 
    StringVector        m_inputParticleCollections;             ///< 

    int                 m_printing;                             ///< 
    std::string         m_rootFile;                             ///< 

    int                 m_nRun;                                 ///< 
    int                 m_nEvt;                                 ///< 

    float               m_ez;                                   ///< 
    float               m_eq1;                                  ///< 
    float               m_eq2;                                  ///< 
    float               m_mz;                                   ///< 

    float               m_costz;                                ///< 
    float               m_costq1;                               ///< 
    float               m_costq2;                               ///< 
    float               m_thrust;                               ///< 
    int                 m_qpdg;                                 ///< 

    ParticleVector      m_pfovec;                               ///< 
    ParticleVector      m_mcpfovec;                             ///< 
    ParticleVector      m_quarkpfovec;                          ///< 
    ParticleVector      m_pMcPFOs;                              ///< 

    TH1F               *fNPFO;                                  ///< 
    TH1F               *fEq;                                    ///< 
    TH1F               *fPFA;                                   ///< 
    TH1F               *fPFAnu;                                 ///<
    TH1F               *fPFAnufwd;                              ///<
    TH1F               *fPFAuds;                                ///<
    TH1F               *fPFAudsHP20;                            ///<
    TH1F               *fPFAudsHP10;                            ///<
    TH1F               *fPFAudsHM10;                            ///<
    TH1F               *fPFAudsHM20;                            ///<
    TH1F               *fPFAFudsHP20;                           ///<
    TH1F               *fPFAFudsHP10;                           ///<
    TH1F               *fPFAFudsHM10;                           ///<
    TH1F               *fPFAFudsHM20;                           ///<
    TH1F               *fPFAudscb;                              ///<
    TH1F               *fPFAcb;                                 ///<
    TH1F               *fPFA1;                                  ///<
    TH1F               *fPFA2;                                  ///<
    TH1F               *fPFA3;                                  ///<
    TH1F               *fPFA4;                                  ///<
    TH1F               *fPFA5;                                  ///<
    TH1F               *fPFA6;                                  ///<
    TH1F               *fPFA7;                                  ///<
    TH1F               *fPFAL7A;                                ///<
    TH1F               *fPFAL7Aud;                              ///<
    TH1F               *fPFAL7As;                               ///<
    TH1F               *fPFAL7Ac;                               ///<
    TH1F               *fPFAL7Ab;                               ///<
    TH1F               *fPFAL7B;                                ///<
    TH1F               *fPFAFL7A;                               ///<
    TH1F               *fPFAFL7Aud;                             ///<
    TH1F               *fPFAFL7As;                              ///<
    TH1F               *fPFAFL7Ac;                              ///<
    TH1F               *fPFAFL7Ab;                              ///<
    TH1F               *fPFA8;                                  ///<
    TH1F               *fPFA9;                                  ///<
    TH1F               *fPFA10;                                 ///<
    TH1F               *fPFA11;                                 ///<
    TH1F               *fPFA12;                                 ///<
    TH1F               *fPFA13;                                 ///<
    TH1F               *fPFA14;                                 ///<
    TH1F               *fPFAMZ;                                 ///<
    TH1F               *fPFAMW;                                 ///<
    TH1F               *fPFAMZa;                                ///<
    TH1F               *fPFAMWa;                                ///<
    TH1F               *fPFAQQ;                                 ///<
    TH1F               *fPFAQQ8;                                ///<
    TH1F               *fPFADMZ;                                ///<
    TH1F               *fPFADMZOMZ;                             ///<
    TH1F               *fPFADMZOMZQQ8;                          ///<
    TH1F               *fPFADMZ8;                               ///<
    TH1F               *fPFADMZQQ8;                             ///<
    TH1F               *fPFADMZP8;                              ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline marlin::Processor *PfoAnalysis::newProcessor()
{
    return new PfoAnalysis;
}

#endif // #ifndef PFO_ANALYSIS_H
