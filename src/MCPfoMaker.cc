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

#include "MCTree.h"
#include "MCPfoMaker.h"

#include <cmath>

MCPfoMaker mcPfoMaker;

//------------------------------------------------------------------------------------------------------------------------------------------

MCPfoMaker::MCPfoMaker() :
    Processor("MCPfoMaker"),
    m_pMCTree(NULL)
{
    _description = "MCPfoMaker analyses output of PandoraPFA" ;

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

            m_pMCTree = new pandora_analysis::MCTree(col, m_lookForQuarksWithMotherZ);

            const pandora_analysis::MCParticleVector &mcPfoVector(m_pMCTree->GetMCPfos());

            // Make a new vector of particles (MCPFOs)
            LCCollectionVec *mcPfoCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

            for (unsigned int i = 0; i < mcPfoVector.size(); ++i)
            {
                ReconstructedParticleImpl *recoPart = new ReconstructedParticleImpl();
                float Mom[3];
                Mom[0] = mcPfoVector[i]->getMomentum()[0];
                Mom[1] = mcPfoVector[i]->getMomentum()[1];
                Mom[2] = mcPfoVector[i]->getMomentum()[2];

                float energy = mcPfoVector[i]->getEnergy();
                float mass = (energy*energy-Mom[0]*Mom[0]-Mom[1]*Mom[1]-Mom[2]*Mom[2]);

                if (mass > 0)
                {
                    mass=std::sqrt(mass);
                }
                else
                {
                    mass = 0.;
                }

                recoPart->setMomentum(Mom);
                recoPart->setEnergy(energy);
                recoPart->setMass(mass);
                recoPart->setType(mcPfoVector[i]->getPDG());
                mcPfoCol->addElement(recoPart);
            }

            pLCEvent->addCollection(mcPfoCol, m_outputMCParticleCollection.c_str());

            // Make a new vector of quarks (MCPFOs)
            LCCollectionVec *mcQuarkCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

            const pandora_analysis::MCParticleVector &mcQuarks(m_pMCTree->GetMCQuarks());

            for (unsigned int i = 0; i < mcQuarks.size(); ++i)
            {
                MCParticle *quark = mcQuarks[i];
                ReconstructedParticleImpl *recoPart = new ReconstructedParticleImpl();

                float Mom[3];
                Mom[0] = quark->getMomentum()[0];
                Mom[1] = quark->getMomentum()[1];
                Mom[2] = quark->getMomentum()[2];

                float energy = quark->getEnergy();
                float mass = (energy*energy-Mom[0]*Mom[0]-Mom[1]*Mom[1]-Mom[2]*Mom[2]);

                if(mass>0)
                {
                    mass=std::sqrt(mass);
                }
                else
                {
                    mass = 0.;
                }

                recoPart->setMomentum(Mom);
                recoPart->setEnergy(energy);
                recoPart->setMass(mass);
                recoPart->setType(quark->getPDG());
                mcQuarkCol->addElement(recoPart);
            }

            pLCEvent->addCollection(mcQuarkCol, m_outputQuarkParticleCollection.c_str());
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
    delete m_pMCTree;
    m_pMCTree = NULL;
}
