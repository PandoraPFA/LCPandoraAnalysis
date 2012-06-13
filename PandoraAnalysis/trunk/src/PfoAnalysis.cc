/**
 *  @file   PandoraAnalysis/src/PfoAnalysis.cc
 * 
 *  @brief  Implementation of the pfo analysis class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"
#include "EVENT/LCGenericObject.h"

#include "CalorimeterHitType.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include "PfoAnalysis.h"

#include <cmath>

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

    this->Clear();

    m_tree = new TTree("PfoAnalysisTree", "PfoAnalysisTree");
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

    m_hPfoEnergySum = new TH1F("fPFA", "total pfo energy", 10000, 0., 5000.);
    m_hPfoEnergySumL7A = new TH1F("fPFA_L7A", "total pfo energy < 0.7 A", 10000, 0., 5000.);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::processRunHeader(lcio::LCRunHeader */*pLCRunHeader*/)
{
    ++m_nRun;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::processEvent(EVENT::LCEvent *pLCEvent)
{
    ++m_nEvt;
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
        std::cout << "PfoAnalysis::end() " << this->name() << " processed " << m_nEvt << " events in " << m_nRun << " runs " << std::endl
                  << "Rootfile: " << m_rootFile.c_str() << std::endl;
    }

    TFile *pTFile = new TFile(m_rootFile.c_str(), "recreate");
    m_tree->Write();
    m_hPfoEnergySum->Write();
    m_hPfoEnergySumL7A->Write();

    pTFile->Close();
    delete pTFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::Clear()
{
    m_pfovec.clear();
    m_mcpfovec.clear();
    m_quarkpfovec.clear();

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

                m_quarkpfovec.push_back(pReconstructedParticle);
            }
        }
        catch (...)
        {
            streamlog_out(ERROR) << "Could not extract input quark particle collection: " << *iter << std::endl;
        }
    }

    // Extract mc particle collection
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

                m_mcpfovec.push_back(pReconstructedParticle);
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

                m_pfovec.push_back(pReconstructedParticle);
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
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::MakeQuarkVariables()
{
    if (!m_quarkpfovec.empty())
    {
        m_qPdg = std::abs(m_quarkpfovec[0]->getType());
        float energyTot(0.f);
        float costTot(0.f);

        for (unsigned int i = 0; i < m_quarkpfovec.size(); ++i)
        {
            const float px(m_quarkpfovec[i]->getMomentum()[0]);
            const float py(m_quarkpfovec[i]->getMomentum()[1]);
            const float pz(m_quarkpfovec[i]->getMomentum()[2]);
            const float energy(m_quarkpfovec[i]->getEnergy());
            const float p(std::sqrt(px * px + py * py + pz * pz));
            const float cost(std::fabs(pz) / p);
            energyTot += energy;
            costTot += cost * energy;
        }

        m_thrust = costTot / energyTot;
    }

    if (m_quarkpfovec.size() == 2)
    {
        const float pQ1[3] = {m_quarkpfovec[0]->getMomentum()[0], m_quarkpfovec[0]->getMomentum()[1], m_quarkpfovec[0]->getMomentum()[2]};
        const float pQ2[3] = {m_quarkpfovec[1]->getMomentum()[0], m_quarkpfovec[1]->getMomentum()[1], m_quarkpfovec[1]->getMomentum()[2]};
        const float pQQ[3] = {pQ1[0] + pQ2[0], pQ1[1] + pQ2[1], pQ1[2] + pQ2[2]};

        TLorentzVector q1(pQ1[0], pQ1[1], pQ1[2], m_quarkpfovec[0]->getEnergy());
        TLorentzVector q2(pQ2[0], pQ2[1], pQ2[2], m_quarkpfovec[1]->getEnergy());
        TLorentzVector qq = q1 + q2;

        m_mQQ = qq.M();
        m_eQQ = m_quarkpfovec[0]->getEnergy() + m_quarkpfovec[1]->getEnergy();

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
    for (unsigned int i = 0; i < m_pfovec.size(); ++i)
    {
        ReconstructedParticle *pPfo = m_pfovec[i];
        ++m_nPfosTotal;
        m_pfoEnergyTotal += pPfo->getEnergy();

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
            std::cout << " RECO PFO : " << i << " " << pPfo->getType() << " " << pPfo->getEnergy() << "  " << pPfo->getTracks().size()
                      << ":" << pPfo->getClusters().size() << " charge : " << pPfo->getCharge() << std::endl;
        }
    }

    m_pfoMassTotal = std::sqrt(m_pfoEnergyTotal * m_pfoEnergyTotal - momTot[0] * momTot[0] - momTot[1] * momTot[1] - momTot[2] * momTot[2]);

    // Extract quantities relating to mc pfos, including energy in "primary" neutrinos
    for (unsigned int imc = 0; imc < m_mcpfovec.size(); ++imc)
    {
        ReconstructedParticle *pMCPfo = m_mcpfovec[imc];
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
