/**
 *  @file   PandoraAnalysis/src/MCPfoMaker.cc
 * 
 *  @brief  Implementation of mc pfo maker class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"

#include "IMPL/LCCollectionVec.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "AnalysisMCTree.h"
#include "MCPfoMaker.h"

#include <cmath>

MCPfoMaker mcPfoMaker;

//------------------------------------------------------------------------------------------------------------------------------------------

MCPfoMaker::MCPfoMaker() :
    Processor("MCPfoMaker"),
    m_pAnalysisMCTree(NULL)
{
    _description = "MCPfoMaker creates mc pfos for comparison with PandoraPFANew output";

    registerInputCollections(LCIO::MCPARTICLE,
                             "InputMCParticleCollections", 
                             "Names of input mc particle collections",
                             m_inputMCParticleCollections,
                             StringVector());

    registerProcessorParameter("OutputMCParticleCollection",
                               "Output mc particle collection name",
                               m_outputMCParticleCollection,
                               std::string());

    registerProcessorParameter("OutputQuarkParticleCollection",
                               "Output quark particle collection name",
                               m_outputQuarkParticleCollection,
                               std::string());

    registerProcessorParameter("LookForQuarksWithMotherZ",
                               "Flag to look for quarks with mother Z",
                               m_lookForQuarksWithMotherZ,
                               bool(false));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCPfoMaker::init()
{
    this->Clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCPfoMaker::processRunHeader(lcio::LCRunHeader */*pLCRunHeader*/)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCPfoMaker::processEvent(EVENT::LCEvent *pLCEvent)
{
    this->Clear();

    for (StringVector::const_iterator iter = m_inputMCParticleCollections.begin(), iterEnd = m_inputMCParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection* col = pLCEvent->getCollection(*iter);

            // Find the MC tree and form collection of PerfectPFOs
            if (col->getTypeName() != LCIO::MCPARTICLE)
                throw;

            m_pAnalysisMCTree = new pandora_analysis::AnalysisMCTree(col, m_lookForQuarksWithMotherZ);

            const pandora_analysis::MCParticleVector &mcPfoVector(m_pAnalysisMCTree->GetMCPfos());

            // Make a new vector of particles (MCPFOs)
            LCCollectionVec *pMCPfoCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

            for (unsigned int i = 0; i < mcPfoVector.size(); ++i)
            {
                MCParticle *pMCParticle = mcPfoVector[i];
                ReconstructedParticleImpl *pRecoParticle = new ReconstructedParticleImpl();

                const float energy(mcPfoVector[i]->getEnergy());
                const float mom[3] = {pMCParticle->getMomentum()[0], pMCParticle->getMomentum()[1], pMCParticle->getMomentum()[2]};
                const float mass(std::sqrt(std::max(0.f, energy * energy - mom[0] * mom[0] - mom[1] * mom[1] - mom[2] * mom[2])));

                pRecoParticle->setMomentum(mom);
                pRecoParticle->setEnergy(energy);
                pRecoParticle->setMass(mass);
                pRecoParticle->setType(mcPfoVector[i]->getPDG());
                pMCPfoCol->addElement(pRecoParticle);
            }

            pLCEvent->addCollection(pMCPfoCol, m_outputMCParticleCollection.c_str());

            // Make a new vector of quarks (MCPFOs)
            LCCollectionVec *pMCQuarkCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

            const pandora_analysis::MCParticleVector &mcQuarks(m_pAnalysisMCTree->GetMCQuarks());

            for (unsigned int i = 0; i < mcQuarks.size(); ++i)
            {
                MCParticle *pQuark = mcQuarks[i];
                ReconstructedParticleImpl *pRecoParticle = new ReconstructedParticleImpl();

                const float energy(pQuark->getEnergy());
                const float mom[3] = {pQuark->getMomentum()[0], pQuark->getMomentum()[1], pQuark->getMomentum()[2]};
                const float mass(std::sqrt(std::max(0.f, energy * energy - mom[0] * mom[0] - mom[1] * mom[1] - mom[2] * mom[2])));

                pRecoParticle->setMomentum(mom);
                pRecoParticle->setEnergy(energy);
                pRecoParticle->setMass(mass);
                pRecoParticle->setType(pQuark->getPDG());
                pMCQuarkCol->addElement(pRecoParticle);
            }

            pLCEvent->addCollection(pMCQuarkCol, m_outputQuarkParticleCollection.c_str());
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Could not extract input mc particle collection: " << *iter << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCPfoMaker::check(EVENT::LCEvent */*pLCEvent*/)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCPfoMaker::end()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MCPfoMaker::Clear()
{
    delete m_pAnalysisMCTree;
    m_pAnalysisMCTree = NULL;
}
