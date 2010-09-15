/**
 *  @file   PandoraAnalysis/include/MCTree.h
 * 
 *  @brief  Header file for the mc tree class.
 * 
 *  $Log: $
 */

#ifndef MC_TREE_H
#define MC_TREE_H 1

#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"

#include "lcio.h"

typedef std::vector<lcio::MCParticle *> MCParticleVector;

/**
 *  @brief  mc tree class
 */
class MCTree
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pLCCollection address of the lc collection
     */
    MCTree(const EVENT::LCCollection *const pLCCollection);

    /**
     *  @brief  Destructor
     */
    ~MCTree();

    /**
     *  @brief  Get mc quarks
     * 
     *  @return the mc quark vector
     */
    const MCParticleVector &GetMCQuarks() const;

    /**
     *  @brief  Get mc pfos
     * 
     *  @return the mc pfo vector
     */
    const MCParticleVector &GetMCPfos() const;

private:
    /**
     *  @brief  Store mc quarks
     * 
     *  @param  pLCCollection address of the lc collection
     */
    void StoreMCQuarks(const EVENT::LCCollection *const pLCCollection);

    /**
     *  @brief  Store mc pfos
     * 
     *  @param  pLCCollection address of the lc collection
     */
    void StoreMCPfos(const EVENT::LCCollection *const pLCCollection);

    MCParticleVector        m_mcQuarks;             ///< The mc quark vector
    MCParticleVector        m_mcPFOs;               ///< The mc pfo vector
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const MCParticleVector &MCTree::GetMCQuarks() const
{
    return m_mcQuarks;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const MCParticleVector &MCTree::GetMCPfos() const
{
    return m_mcPFOs;
}

#endif // #ifndef MC_TREE_H
