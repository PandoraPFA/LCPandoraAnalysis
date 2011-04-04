/**
 *  @file   PandoraAnalysis/src/AnalysisMCTree.cc
 * 
 *  @brief  Implementation of the analysis mc tree class.
 * 
 *  $Log: $
 */

#include "AnalysisMCTree.h"

#include <cmath>
#include <cstdlib>

namespace pandora_analysis
{

AnalysisMCTree::AnalysisMCTree(const EVENT::LCCollection *const pLCCollection, const bool lookForQuarksWithMotherZ)
{
    if (pLCCollection->getTypeName() != lcio::LCIO::MCPARTICLE)
        throw;

    this->StoreMCQuarks(pLCCollection, lookForQuarksWithMotherZ);
    this->StoreMCPfos(pLCCollection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

AnalysisMCTree::~AnalysisMCTree()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisMCTree::StoreMCQuarks(const EVENT::LCCollection *const pLCCollection, const bool lookForQuarksWithMotherZ)
{
    try
    {
        for (unsigned int i = 0, nParticles = pLCCollection->getNumberOfElements(); i < nParticles; ++i)
        {
            lcio::MCParticle *pMCParticle = dynamic_cast<lcio::MCParticle*>(pLCCollection->getElementAt(i));

            if (NULL == pMCParticle)
                throw EVENT::Exception("Collection type mismatch");

            const int absPdgCode(std::abs(pMCParticle->getPDG()));

            // By default, the primary quarks are the ones without any parents
            if (lookForQuarksWithMotherZ == false)
            {
                if ((absPdgCode >= 1) && (absPdgCode <= 6) && pMCParticle->getParents().empty())
                {
                    m_mcQuarks.push_back(pMCParticle);
                }
            }
            else if (lookForQuarksWithMotherZ == true)
            {
                // For MC files generated in the SLIC environment, the primary quarks have parents
                if ((absPdgCode >= 1) && (absPdgCode <= 6))
                {
                    // The mother should be the Z-boson
                    if ((pMCParticle->getParents().size() == 1) && ((pMCParticle->getParents())[0]->getPDG() == 23))
                    {
                        m_mcQuarks.push_back(pMCParticle);
                    }
                }
            }
        }
    }
    catch (EVENT::Exception &exception)
    {
        std::cout << "AnalysisMCTree::StoreMCQuarks - Failed to extract mc particle: " << exception.what() << std::endl;
        throw exception;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AnalysisMCTree::StoreMCPfos(const EVENT::LCCollection *const pLCCollection)
{
    try
    {
        MCParticleVector usedMCParticles;

        for (unsigned int i = 0, nParticles = pLCCollection->getNumberOfElements(); i < nParticles; ++i)
        {
            lcio::MCParticle *pMCParticle = dynamic_cast<lcio::MCParticle*>(pLCCollection->getElementAt(i));

            if (NULL == pMCParticle)
                throw EVENT::Exception("Collection type mismatch");

            const float rs(std::sqrt((pMCParticle->getVertex()[0] * pMCParticle->getVertex()[0]) +
                (pMCParticle->getVertex()[1] * pMCParticle->getVertex()[1]) +
                (pMCParticle->getVertex()[2] * pMCParticle->getVertex()[2]) ));

            const float re(std::sqrt((pMCParticle->getEndpoint()[0] * pMCParticle->getEndpoint()[0]) +
                (pMCParticle->getEndpoint()[1] * pMCParticle->getEndpoint()[1]) +
                (pMCParticle->getEndpoint()[2] * pMCParticle->getEndpoint()[2]) ));

            float p(0.f);
            const float energySquared(pMCParticle->getEnergy() * pMCParticle->getEnergy());
            const float massSquared(pMCParticle->getMass() * pMCParticle->getMass());

            if (energySquared > massSquared)
            {
                p = std::sqrt(energySquared - massSquared);
            }

            // Cut on momentum and inner/outer radii
            if ((p > 0.01) && (rs < 500) && (re > 500))
            {
                usedMCParticles.push_back(pMCParticle);

                // Look to see if parent mc particle has already been used
                bool found = false;
                lcio::MCParticle *pTemporaryMCParticle = pMCParticle;
                MCParticleVector parents = pMCParticle->getParents();

                while (!parents.empty())
                {
                    pTemporaryMCParticle = parents[0];

                    for (MCParticleVector::const_iterator iter = usedMCParticles.begin(), iterEnd = usedMCParticles.end(); iter != iterEnd; ++iter)
                    {
                        if (pTemporaryMCParticle == *iter)
                        {
                            found = true;
                            break;
                        }
                    }

                    parents = pTemporaryMCParticle->getParents();

                    if (found)
                        break;
                }

                // Veto low energy neutrons and protons
                bool lownp = false;

                if ((pMCParticle->getPDG() == 2112) || (pMCParticle->getPDG() == 2212))
                {
                    if ((rs > 10) && (pMCParticle->getEnergy() < 1.2))
                    {
                        lownp = true;
                    }
                }

                // Add mc particle to the list of mc pfos
                if (!found && !lownp)
                {
                    m_mcPFOs.push_back(pMCParticle);
                }
            }
        }
    }
    catch (EVENT::Exception &exception)
    {
        std::cout << "AnalysisMCTree::StoreMCPfos - Failed to extract mc particle: " << exception.what() << std::endl;
        throw exception;
    }
}

} // namespace pandora_analysis
