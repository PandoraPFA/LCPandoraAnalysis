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

namespace pandora_analysis
{

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
     *  @param  lookForQuarksWithMother flag to look for quarks with mother (by default, the primary quarks are the ones having
     *          no parents; but in the SLIC framework, there are no such quarks; enable this flag if you are working with files
     *          generated in the SLIC environment)
     */
    MCTree(const EVENT::LCCollection *const pLCCollection, const bool lookForQuarksWithMotherZ);

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
     *  @param  lookForQuarksWithMother flag to look for quarks with mother (by default, the primary quarks are the ones having
     *          no parents; but in the SLIC framework, there are no such quarks; enable this flag if you are working with files
     *          generated in the SLIC environment)
     */
    void StoreMCQuarks(const EVENT::LCCollection *const pLCCollection, const bool lookForQuarksWithMotherZ);

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

} // namespace pandora_analysis

#endif // #ifndef MC_TREE_H
