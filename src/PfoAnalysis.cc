/**
 *  @file   PandoraAnalysis/src/PfoAnalysis.cc
 * 
 *  @brief  Implementation of the pfo analysis class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"
#include "EVENT/LCGenericObject.h"
#include "EVENT/MCParticle.h"

#include "CalorimeterHitType.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include "PfoAnalysis.h"

#include <cmath>
#include <set>

PfoAnalysis pfoAnalysis;

//------------------------------------------------------------------------------------------------------------------------------------------

PfoAnalysis::PfoAnalysis() :
    Processor("PfoAnalysis")
{
    _description = "PfoAnalysis analyses output of PandoraPFANew";

    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
                             "InputQuarkParticleCollections",
                             "Names of input quark particle collections",
                             m_inputQuarkParticleCollections,
                             StringVector());

    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
                             "InputMCParticleCollections",
                             "Names of input mc particle collections",
                             m_inputMCParticleCollections,
                             StringVector());

    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
                             "InputParticleCollections",
                             "Names of input reconstructed particle collections",
                             m_inputParticleCollections,
                             StringVector());

    registerInputCollections(LCIO::LCGENERICOBJECT,
                             "InputReclusterMonitoringCollections",
                             "Names of input recluster monitoring collections",
                             m_inputReclusterMonitoringCollections,
                             StringVector());

    registerInputCollections(LCIO::MCPARTICLE,
                             "MCParticleCollections", 
                             "Names of mc particle collections",
                             m_mcParticleCollections,
                             StringVector());

    std::string rootFile("PFOAnalysis.root");
    registerProcessorParameter( "RootFile",
                                "Name of the output root file",
                                m_rootFile,
                                rootFile);

    registerProcessorParameter( "Printing",
                                "Set the debug print level",
                                m_printing,
                                (int)0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::init()
{
    m_nRun = 0;
    m_nEvt = 0;
    m_nRunSum = 0;
    m_nEvtSum = 0;
    this->Clear();

    m_pTFile = new TFile(m_rootFile.c_str(), "recreate");

    m_tree = new TTree("PfoAnalysisTree", "PfoAnalysisTree");
    m_tree->SetDirectory(m_pTFile);
    m_tree->Branch("run", &m_nRun, "run/I");
    m_tree->Branch("event", &m_nEvt, "event/I");
    m_tree->Branch("nPfosTotal", &m_nPfosTotal, "nPfosTotal/I");
    m_tree->Branch("nPfosNeutralHadrons", &m_nPfosNeutralHadrons, "nPfosNeutralHadrons/I");
    m_tree->Branch("nPfosPhotons", &m_nPfosPhotons, "nPfosPhotons/I");
    m_tree->Branch("nPfosTracks", &m_nPfosTracks, "nPfosTracks/I");
    m_tree->Branch("pfoEnergyTotal", &m_pfoEnergyTotal, "pfoEnergyTotal/F");
    m_tree->Branch("pfoEnergyNeutralHadrons", &m_pfoEnergyNeutralHadrons, "pfoEnergyNeutralHadrons/F");
    m_tree->Branch("pfoEnergyPhotons", &m_pfoEnergyPhotons, "pfoEnergyPhotons/F");
    m_tree->Branch("pfoEnergyTracks", &m_pfoEnergyTracks, "pfoEnergyTracks/F");
    m_tree->Branch("pfoECalToEmEnergy", &m_pfoECalToEmEnergy, "pfoECalToEmEnergy/F");
    m_tree->Branch("pfoECalToHadEnergy", &m_pfoECalToHadEnergy, "pfoECalToHadEnergy/F");
    m_tree->Branch("pfoHCalToEmEnergy", &m_pfoHCalToEmEnergy, "pfoHCalToEmEnergy/F");
    m_tree->Branch("pfoHCalToHadEnergy", &m_pfoHCalToHadEnergy, "pfoHCalToHadEnergy/F");
    m_tree->Branch("pfoOtherEnergy", &m_pfoOtherEnergy, "pfoOtherEnergy/F");
    m_tree->Branch("pfoMuonToEnergy", &m_pfoMuonToEnergy, "pfoMuonToEnergy/F");
    m_tree->Branch("pfoMassTotal", &m_pfoMassTotal, "pfoMassTotal/F");
    m_tree->Branch("mcEnergyTotal", &m_mcEnergyTotal, "mcEnergyTotal/F");
    m_tree->Branch("mcEnergyENu", &m_mcEnergyENu, "mcEnergyENu/F");
    m_tree->Branch("mcEnergyFwd", &m_mcEnergyFwd, "mcEnergyFwd/F");
    m_tree->Branch("eQQ", &m_eQQ, "eQQ/F");
    m_tree->Branch("eQ1", &m_eQ1, "eQ1/F");
    m_tree->Branch("eQ2", &m_eQ2, "eQ2/F");
    m_tree->Branch("costQQ", &m_costQQ, "costQQ/F");
    m_tree->Branch("costQ1", &m_costQ1, "costQ1/F");
    m_tree->Branch("costQ2", &m_costQ2, "costQ2/F");
    m_tree->Branch("mQQ", &m_mQQ, "mQQ/F");
    m_tree->Branch("thrust", &m_thrust, "thrust/F");
    m_tree->Branch("qPdg", &m_qPdg, "qPdg/I");
    m_tree->Branch("netEnergyChange", &m_netEnergyChange, "netEnergyChange/F");
    m_tree->Branch("sumModulusEnergyChanges", &m_sumModulusEnergyChanges, "sumModulusEnergyChanges/F");
    m_tree->Branch("sumSquaredEnergyChanges", &m_sumSquaredEnergyChanges, "sumSquaredEnergyChanges/F");
    m_tree->Branch("pfoEnergies", &m_pfoEnergies);
    m_tree->Branch("pfoPx", &m_pfoPx);
    m_tree->Branch("pfoPy", &m_pfoPy);
    m_tree->Branch("pfoPz", &m_pfoPz);
    m_tree->Branch("pfoCosTheta", &m_pfoCosTheta);
    m_tree->Branch("pfoTargetEnergies", &m_pfoTargetEnergies);
    m_tree->Branch("pfoTargetPx", &m_pfoTargetPx);
    m_tree->Branch("pfoTargetPy", &m_pfoTargetPy);
    m_tree->Branch("pfoTargetPz", &m_pfoTargetPz);
    m_tree->Branch("pfoTargetCosTheta", &m_pfoTargetCosTheta);
    m_tree->Branch("pfoPdgCodes", &m_pfoPdgCodes);
    m_tree->Branch("pfoTargetPdgCodes", &m_pfoTargetPdgCodes);
    m_tree->Branch("nPfoTargetsTotal", &m_nPfoTargetsTotal, "nPfoTargetsTotal/I");
    m_tree->Branch("nPfoTargetsNeutralHadrons", &m_nPfoTargetsNeutralHadrons, "nPfoTargetsNeutralHadrons/I");
    m_tree->Branch("nPfoTargetsPhotons", &m_nPfoTargetsPhotons, "nPfoTargetsPhotons/I");
    m_tree->Branch("nPfoTargetsTracks", &m_nPfoTargetsTracks, "nPfoTargetsTracks/I");
    m_tree->Branch("pfoTargetsEnergyTotal", &m_pfoTargetsEnergyTotal, "pfoTargetsEnergyTotal/F");
    m_tree->Branch("pfoTargetsEnergyNeutralHadrons", &m_pfoTargetsEnergyNeutralHadrons, "pfoTargetsEnergyNeutralHadrons/F");
    m_tree->Branch("pfoTargetsEnergyPhotons", &m_pfoTargetsEnergyPhotons, "pfoTargetsEnergyPhotons/F");
    m_tree->Branch("pfoTargetsEnergyTracks", &m_pfoTargetsEnergyTracks, "pfoTargetsEnergyTracks/F");

    m_hPfoEnergySum = new TH1F("fPFA", "total pfo energy", 10000, 0., 5000.);
    m_hPfoEnergySum->SetDirectory(m_pTFile);
    m_hPfoEnergySumL7A = new TH1F("fPFA_L7A", "total pfo energy < 0.7 A", 10000, 0., 5000.);
    m_hPfoEnergySumL7A->SetDirectory(m_pTFile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::processRunHeader(lcio::LCRunHeader */*pLCRunHeader*/)
{
    m_nRun = 0;
    m_nEvt = 0;
    ++m_nRunSum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::processEvent(EVENT::LCEvent *pLCEvent)
{
    m_nRun = pLCEvent->getRunNumber();
    m_nEvt = pLCEvent->getEventNumber();
    ++m_nEvtSum;
    this->Clear();
    this->ExtractCollections(pLCEvent);
    this->MakeQuarkVariables();
    this->PerformPfoAnalysis();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::check(EVENT::LCEvent */*pLCEvent*/)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::end()
{
    if (m_printing > -1)
    {
        std::cout << "PfoAnalysis::end() " << this->name() << " processed " << m_nEvtSum << " events in " << m_nRunSum << " runs " << std::endl
                  << "Rootfile: " << m_rootFile.c_str() << std::endl;
    }

    m_tree->Write();
    m_hPfoEnergySum->Write();
    m_hPfoEnergySumL7A->Write();

    m_pTFile->Write();
    m_pTFile->Close();
    delete m_pTFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::Clear()
{
    m_pfoVector.clear();
    m_mcPfoVector.clear();
    m_quarkPfoVector.clear();
    m_pfoTargetVector.clear();

    m_nPfosTotal = 0;
    m_nPfosNeutralHadrons = 0;
    m_nPfosPhotons = 0;
    m_nPfosTracks = 0;
    m_pfoEnergyTotal = 0.f;
    m_pfoEnergyNeutralHadrons = 0.f;
    m_pfoEnergyPhotons = 0.f;

    m_pfoEnergyTracks = 0.f;
    m_pfoECalToEmEnergy = 0.f;
    m_pfoECalToHadEnergy = 0.f;
    m_pfoHCalToEmEnergy = 0.f;
    m_pfoHCalToHadEnergy = 0.f;
    m_pfoMuonToEnergy = 0.f;
    m_pfoOtherEnergy = 0.f;

    m_pfoMassTotal = 0.f;
    m_mcEnergyTotal = 0.f;
    m_mcEnergyENu = 0.f;
    m_mcEnergyFwd = 0.f;
    m_eQQ = -99.f;
    m_eQ1 = -99.f;
    m_eQ2 = -99.f;
    m_costQQ = -99.f;
    m_costQ1 = -99.f;
    m_costQ2 = -99.f;
    m_mQQ = -99.f;
    m_thrust = -99.f;
    m_qPdg = -99;
    m_netEnergyChange = 0.f;
    m_sumModulusEnergyChanges = 0.f;
    m_sumSquaredEnergyChanges = 0.f;

    m_pfoEnergies.clear();
    m_pfoPx.clear();
    m_pfoPy.clear();
    m_pfoPz.clear();
    m_pfoCosTheta.clear();

    m_pfoTargetEnergies.clear();
    m_pfoTargetPx.clear();
    m_pfoTargetPy.clear();
    m_pfoTargetPz.clear();
    m_pfoTargetCosTheta.clear();

    m_pfoPdgCodes.clear();
    m_pfoTargetPdgCodes.clear();

    m_nPfoTargetsTotal = 0;
    m_nPfoTargetsNeutralHadrons = 0;
    m_nPfoTargetsPhotons = 0;
    m_nPfoTargetsTracks = 0;

    m_pfoTargetsEnergyTotal = 0.f;
    m_pfoTargetsEnergyNeutralHadrons = 0.f;
    m_pfoTargetsEnergyPhotons = 0.f;
    m_pfoTargetsEnergyTracks = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::ExtractCollections(EVENT::LCEvent *pLCEvent)
{
    // Extract quark particle collection
    for (StringVector::const_iterator iter = m_inputQuarkParticleCollections.begin(), iterEnd = m_inputQuarkParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
            {
                ReconstructedParticle *pReconstructedParticle = dynamic_cast<ReconstructedParticle*>(pLCCollection->getElementAt(i));

                if (NULL == pReconstructedParticle)
                    throw EVENT::Exception("Collection type mismatch");

                m_quarkPfoVector.push_back(pReconstructedParticle);
            }
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Could not extract input quark particle collection: " << *iter << std::endl;
        }
    }

    // Extract mc pfo collection
    for (StringVector::const_iterator iter = m_inputMCParticleCollections.begin(), iterEnd = m_inputMCParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
            {
                ReconstructedParticle *pReconstructedParticle = dynamic_cast<ReconstructedParticle*>(pLCCollection->getElementAt(i));

                if (NULL == pReconstructedParticle)
                    throw EVENT::Exception("Collection type mismatch");

                m_mcPfoVector.push_back(pReconstructedParticle);
            }
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Could not extract input mc particle collection: " << *iter << std::endl;
        }
    }

    // Extract reconstructed particle collection
    for (StringVector::const_iterator iter = m_inputParticleCollections.begin(), iterEnd = m_inputParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
            {
                ReconstructedParticle *pReconstructedParticle = dynamic_cast<ReconstructedParticle*>(pLCCollection->getElementAt(i));

                if (NULL == pReconstructedParticle)
                    throw EVENT::Exception("Collection type mismatch");

                m_pfoVector.push_back(pReconstructedParticle);
            }
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Could not extract input particle collection: " << *iter << std::endl;
        }
    }

    // Extract recluster monitoring information, if present
    for (StringVector::const_iterator iter = m_inputReclusterMonitoringCollections.begin(), iterEnd = m_inputReclusterMonitoringCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
            {
                EVENT::LCGenericObject *pLCGenericObject = dynamic_cast<EVENT::LCGenericObject*>(pLCCollection->getElementAt(i));

                if (NULL == pLCGenericObject)
                    throw EVENT::Exception("Collection type mismatch");

                m_netEnergyChange += pLCGenericObject->getFloatVal(0);
                m_sumModulusEnergyChanges += pLCGenericObject->getFloatVal(1);
                m_sumSquaredEnergyChanges += pLCGenericObject->getFloatVal(2);
            }
        }
        catch (...)
        {
            streamlog_out(WARNING) << "Could not extract ReclusterMonitoring information" << std::endl;
        }
    }

    // Extract mc particle collection
    std::set<MCParticle*> pfoTargetList;

    for (StringVector::const_iterator iter = m_mcParticleCollections.begin(), iterEnd = m_mcParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(*iter);

            for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
            {
                MCParticle *pMCParticle = dynamic_cast<MCParticle*>(pLCCollection->getElementAt(i));

                if (NULL == pMCParticle)
                    throw EVENT::Exception("Collection type mismatch");

                const float innerRadius(std::sqrt(pMCParticle->getVertex()[0] * pMCParticle->getVertex()[0] + pMCParticle->getVertex()[1] * pMCParticle->getVertex()[1] + pMCParticle->getVertex()[2] * pMCParticle->getVertex()[2]));
                const float outerRadius(std::sqrt(pMCParticle->getEndpoint()[0] * pMCParticle->getEndpoint()[0] + pMCParticle->getEndpoint()[1] * pMCParticle->getEndpoint()[1] + pMCParticle->getEndpoint()[2] * pMCParticle->getEndpoint()[2]));
                const float momentum(std::sqrt(pMCParticle->getMomentum()[0] * pMCParticle->getMomentum()[0] + pMCParticle->getMomentum()[1] * pMCParticle->getMomentum()[1] + pMCParticle->getMomentum()[2] * pMCParticle->getMomentum()[2]));

                if ((pfoTargetList.find(pMCParticle) == pfoTargetList.end()) && (outerRadius > 500.f) && (innerRadius <= 500.f) &&
                    (momentum > 0.01f) && !((pMCParticle->getPDG() == 2212 || pMCParticle->getPDG() == 2112) && (pMCParticle->getEnergy() < 1.2f)))
                {
                    pfoTargetList.insert(pMCParticle);
                    m_pfoTargetVector.push_back(pMCParticle);
                }
            }
        }
        catch (...)
        {
            streamlog_out(WARNING) << "Could not extract mc particle information" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::MakeQuarkVariables()
{
    if (!m_quarkPfoVector.empty())
    {
        m_qPdg = std::abs(m_quarkPfoVector[0]->getType());
        float energyTot(0.f);
        float costTot(0.f);

        for (unsigned int i = 0; i < m_quarkPfoVector.size(); ++i)
        {
            const float px(m_quarkPfoVector[i]->getMomentum()[0]);
            const float py(m_quarkPfoVector[i]->getMomentum()[1]);
            const float pz(m_quarkPfoVector[i]->getMomentum()[2]);
            const float energy(m_quarkPfoVector[i]->getEnergy());
            const float p(std::sqrt(px * px + py * py + pz * pz));
            const float cost(std::fabs(pz) / p);
            energyTot += energy;
            costTot += cost * energy;
        }

        m_thrust = costTot / energyTot;
    }

    if (m_quarkPfoVector.size() == 2)
    {
        const float pQ1[3] = {m_quarkPfoVector[0]->getMomentum()[0], m_quarkPfoVector[0]->getMomentum()[1], m_quarkPfoVector[0]->getMomentum()[2]};
        const float pQ2[3] = {m_quarkPfoVector[1]->getMomentum()[0], m_quarkPfoVector[1]->getMomentum()[1], m_quarkPfoVector[1]->getMomentum()[2]};
        const float pQQ[3] = {pQ1[0] + pQ2[0], pQ1[1] + pQ2[1], pQ1[2] + pQ2[2]};

        TLorentzVector q1(pQ1[0], pQ1[1], pQ1[2], m_quarkPfoVector[0]->getEnergy());
        TLorentzVector q2(pQ2[0], pQ2[1], pQ2[2], m_quarkPfoVector[1]->getEnergy());
        TLorentzVector qq = q1 + q2;

        m_mQQ = qq.M();
        m_eQQ = m_quarkPfoVector[0]->getEnergy() + m_quarkPfoVector[1]->getEnergy();

        const float pQ1Tot(std::sqrt(pQ1[0] * pQ1[0] + pQ1[1] * pQ1[1] + pQ1[2] * pQ1[2]));
        const float pQ2Tot(std::sqrt(pQ2[0] * pQ2[0] + pQ2[1] * pQ2[1] + pQ2[2] * pQ2[2]));
        const float pQQTot(std::sqrt(pQQ[0] * pQQ[0] + pQQ[1] * pQQ[1] + pQQ[2] * pQQ[2]));

        if (0.f != pQQTot)
            m_costQQ = pQQ[2] / pQQTot;

        if (0.f != pQ1Tot)
            m_costQ1 = pQ1[2] / pQ1Tot;

        if (0.f != pQ2Tot)
            m_costQ2 = pQ2[2] / pQ2Tot;

        m_eQ1 = pQ1Tot;
        m_eQ2 = pQ2Tot;

        if (m_printing > 0)
        {
            std::cout << " eQQ    = " << m_eQQ << std::endl
                      << " eQ1    = " << m_eQ1 << std::endl
                      << " eQ2    = " << m_eQ2 << std::endl
                      << " costQQ = " << m_costQQ << std::endl
                      << " costQ1 = " << m_costQ1 << std::endl
                      << " costQ2 = " << m_costQ2 << std::endl
                      << " mQQ    = " << m_mQQ << std::endl
                      << " Thrust = " << m_thrust << std::endl
                      << " QPDG   = " << m_qPdg << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::PerformPfoAnalysis()
{
    float momTot[3] = {0.f, 0.f, 0.f};

    // Extract quantities relating to reconstructed pfos
    for (ParticleVector::const_iterator iter = m_pfoVector.begin(), iterEnd = m_pfoVector.end(); iter != iterEnd; ++iter)
    {
        ReconstructedParticle *pPfo = *iter;

        ++m_nPfosTotal;
        m_pfoEnergyTotal += pPfo->getEnergy();
        m_pfoPdgCodes.push_back(pPfo->getType());
        m_pfoEnergies.push_back(pPfo->getEnergy());

        m_pfoPx.push_back(pPfo->getMomentum()[0]);
        m_pfoPy.push_back(pPfo->getMomentum()[1]);
        m_pfoPz.push_back(pPfo->getMomentum()[2]);

        const float momentum(std::sqrt(pPfo->getMomentum()[0] * pPfo->getMomentum()[0] + pPfo->getMomentum()[1] * pPfo->getMomentum()[1] + pPfo->getMomentum()[2] * pPfo->getMomentum()[2]));
        const float cosTheta((momentum > std::numeric_limits<float>::epsilon()) ? pPfo->getMomentum()[2] / momentum : -999.f);
        m_pfoCosTheta.push_back(cosTheta);

        if (!pPfo->getTracks().empty())
        {
            // Charged pfos
            ++m_nPfosTracks;
            m_pfoEnergyTracks += pPfo->getEnergy();
        }
        else
        {
            // Neutral pfos
            float cellEnergySum(0.f);
            const ClusterVec &clusterVec = pPfo->getClusters();

            for (ClusterVec::const_iterator iter = clusterVec.begin(), iterEnd = clusterVec.end(); iter != iterEnd; ++iter)
            {
                const CalorimeterHitVec &calorimeterHitVec = (*iter)->getCalorimeterHits();

                for (CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
                {
                    EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;
                    cellEnergySum += pCalorimeterHit->getEnergy();
                }
            }

            //if (cellEnergySum < std::numeric_limits<float>::epsilon())
            //{
            //    std::cout << "Error, pfo found with neither tracks nor clusters... " << std::endl;
            //    throw;
            //}

            const float correctionFactor((cellEnergySum < std::numeric_limits<float>::epsilon()) ? 0.f : pPfo->getEnergy() / cellEnergySum);

            if (22 == pPfo->getType())
            {
                // Photons
                ++m_nPfosPhotons;
                m_pfoEnergyPhotons += pPfo->getEnergy();

                for (ClusterVec::const_iterator iter = clusterVec.begin(), iterEnd = clusterVec.end(); iter != iterEnd; ++iter)
                {
                    const CalorimeterHitVec &calorimeterHitVec = (*iter)->getCalorimeterHits();

                    for (CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
                    {
                        EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;

                        const float hitEnergy(correctionFactor * pCalorimeterHit->getEnergy());
                        CHT cht(pCalorimeterHit->getType());

                        //std::cout << "EM, type: " << pCalorimeterHit->getType() << ", energy: " << pCalorimeterHit->getEnergy() << " correction " << correctionFactor << std::endl;
                        //std::cout << " cht.is(CHT::ecal) " << cht.is(CHT::ecal) << " cht.is(CHT::hcal) " << cht.is(CHT::hcal) << " cht.is(CHT::muon) " << cht.is(CHT::muon) << std::endl;

                        if (cht.is(CHT::ecal))
                        {
                            m_pfoECalToEmEnergy += hitEnergy;
                        }
                        else if (cht.is(CHT::hcal))
                        {
                            m_pfoHCalToEmEnergy += hitEnergy;
                        }
                        else if (cht.is(CHT::muon))
                        {
                            m_pfoMuonToEnergy += hitEnergy;
                        }
                        else
                        {
                            m_pfoOtherEnergy += hitEnergy;
                        }
                    }
                }
            }
            else
            {
                // Neutral hadrons
                ++m_nPfosNeutralHadrons;
                m_pfoEnergyNeutralHadrons += pPfo->getEnergy();

                for (ClusterVec::const_iterator iter = clusterVec.begin(), iterEnd = clusterVec.end(); iter != iterEnd; ++iter)
                {
                    const CalorimeterHitVec &calorimeterHitVec = (*iter)->getCalorimeterHits();

                    for (CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
                    {
                        EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;

                        const float hitEnergy(correctionFactor * pCalorimeterHit->getEnergy());
                        CHT cht(pCalorimeterHit->getType());

                        //std::cout << "HAD, type: " << pCalorimeterHit->getType() << ", energy: " << pCalorimeterHit->getEnergy() << " correction " << correctionFactor << std::endl;
                        //std::cout << " cht.is(CHT::ecal) " << cht.is(CHT::ecal) << " cht.is(CHT::hcal) " << cht.is(CHT::hcal) << " cht.is(CHT::muon) " << cht.is(CHT::muon) << std::endl;

                        if (cht.is(CHT::ecal))
                        {
                            m_pfoECalToHadEnergy += hitEnergy;
                        }
                        else if (cht.is(CHT::hcal))
                        {
                            m_pfoHCalToHadEnergy += hitEnergy;
                        }
                        else if (cht.is(CHT::muon))
                        {
                            m_pfoMuonToEnergy += hitEnergy;
                        }
                        else
                        {
                            m_pfoOtherEnergy += hitEnergy;
                        }
                    }
                }
            }
        }

        momTot[0] += pPfo->getMomentum()[0];
        momTot[1] += pPfo->getMomentum()[1];
        momTot[2] += pPfo->getMomentum()[2];

        if (m_printing > 0)
        {
            std::cout << " RECO PFO, pdg: " << pPfo->getType() << ", E: " << pPfo->getEnergy() << ", nTracks: " << pPfo->getTracks().size()
                      << ", nClusters: " << pPfo->getClusters().size() << ", charge: " << pPfo->getCharge() << std::endl;
        }
    }

    m_pfoMassTotal = std::sqrt(m_pfoEnergyTotal * m_pfoEnergyTotal - momTot[0] * momTot[0] - momTot[1] * momTot[1] - momTot[2] * momTot[2]);

    // Extract quantities relating to mc pfos, including energy in "primary" neutrinos
    for (ParticleVector::const_iterator iter = m_mcPfoVector.begin(), iterEnd = m_mcPfoVector.end(); iter != iterEnd; ++iter)
    {
        ReconstructedParticle *pMCPfo = *iter;

        m_mcEnergyTotal += pMCPfo->getEnergy();
        const int pdgCode(pMCPfo->getType());

        if ((std::abs(pdgCode) == 12) || (std::abs(pdgCode) == 14) || (std::abs(pdgCode) == 16))
            m_mcEnergyENu += pMCPfo->getEnergy();

        // Find polar angle of MC PFO
        const float px(pMCPfo->getMomentum()[0]);
        const float py(pMCPfo->getMomentum()[1]);
        const float pz(pMCPfo->getMomentum()[2]);

        if (std::fabs(pz) / std::sqrt(px * px + py * py + pz * pz) > 0.98)
            m_mcEnergyFwd += pMCPfo->getEnergy();
    }

    // Extract quantities relating to pfo targets
    for (MCParticleVector::const_iterator iter = m_pfoTargetVector.begin(), iterEnd = m_pfoTargetVector.end(); iter != iterEnd; ++iter)
    {
        MCParticle *pMCParticle = *iter;

        ++m_nPfoTargetsTotal;
        m_pfoTargetsEnergyTotal += pMCParticle->getEnergy();
        m_pfoTargetPdgCodes.push_back(pMCParticle->getPDG());
        m_pfoTargetEnergies.push_back(pMCParticle->getEnergy());

        m_pfoTargetPx.push_back(pMCParticle->getMomentum()[0]);
        m_pfoTargetPy.push_back(pMCParticle->getMomentum()[1]);
        m_pfoTargetPz.push_back(pMCParticle->getMomentum()[2]);

        const float momentum(std::sqrt(pMCParticle->getMomentum()[0] * pMCParticle->getMomentum()[0] + pMCParticle->getMomentum()[1] * pMCParticle->getMomentum()[1] + pMCParticle->getMomentum()[2] * pMCParticle->getMomentum()[2]));
        const float cosTheta((momentum > std::numeric_limits<float>::epsilon()) ? pMCParticle->getMomentum()[2] / momentum : -999.f);
        m_pfoTargetCosTheta.push_back(cosTheta);

        if (22 == pMCParticle->getPDG())
        {
            ++m_nPfoTargetsPhotons;
            m_pfoTargetsEnergyPhotons += pMCParticle->getEnergy();
        }
        else if ((11 == std::abs(pMCParticle->getPDG())) || (13 == std::abs(pMCParticle->getPDG())) || (211 == std::abs(pMCParticle->getPDG()))) // TODO, more options here?
        {
            ++m_nPfoTargetsTracks;
            m_pfoTargetsEnergyTracks += pMCParticle->getEnergy();
        }
        else
        {
            ++m_nPfoTargetsNeutralHadrons;
            m_pfoTargetsEnergyNeutralHadrons += pMCParticle->getEnergy();
        }
    }

    if (m_printing > 0)
    {
        std::cout << " EVENT                : " << m_nEvt << std::endl
                  << " NPFOs                : " << m_nPfosTotal << " (" << m_nPfosTracks << " + " << m_nPfosPhotons + m_nPfosNeutralHadrons << ")" << std::endl
                  << " RECONSTRUCTED ENERGY : " << m_pfoEnergyTotal << std::endl
                  << " RECO ENERGY + eNu    : " << m_pfoEnergyTotal + m_mcEnergyENu << std::endl;
    }

    // Fill basic histograms
    m_hPfoEnergySum->Fill(m_pfoEnergyTotal, 1.);

    if ((m_qPdg >= 1) && (m_qPdg <= 3) && (m_thrust <= 0.7))
        m_hPfoEnergySumL7A->Fill(m_pfoEnergyTotal + m_mcEnergyENu, 1.);

    m_tree->Fill();
}
