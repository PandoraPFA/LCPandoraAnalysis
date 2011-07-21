/**
 *  @file   PandoraAnalysis/include/MCPfoMaker.h
 * 
 *  @brief  Header file for mc pfo maker class.
 * 
 *  $Log: $
 */

#ifndef MC_PFO_MAKER_H
#define MC_PFO_MAKER_H 1

#include "marlin/Processor.h"

namespace pandora_analysis {class AnalysisMCTree;}

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MCPfoMaker class
 */
class MCPfoMaker : public marlin::Processor
{
public:
    /**
     *  @brief  Default constructor
     */
    MCPfoMaker();

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
     *  @brief  Clear current mc tree details
     */
    void Clear();

    typedef std::vector<std::string> StringVector;

    StringVector                        m_inputMCParticleCollections;       ///< The names of the input mc particle collections
    std::string                         m_outputMCParticleCollection;       ///< The output mc particle collection name
    std::string                         m_outputQuarkParticleCollection;    ///< The output quark particle collection name
    bool                                m_lookForQuarksWithMotherZ;         ///< Flag to look for quarks with mother

    pandora_analysis::AnalysisMCTree   *m_pAnalysisMCTree;                  ///< Address of the analysis mc tree
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline marlin::Processor *MCPfoMaker::newProcessor()
{
    return new MCPfoMaker;
}

#endif // #ifndef MC_PFO_MAKER_H
