/**
 *  @file   PandoraAnalysis/src/PfoAnalysis.cc
 * 
 *  @brief  Implementation of the pfo analysis class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"

#include "CalorimeterHitType.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

#include "PfoAnalysis.h"

#include <cmath>
#include <cstdlib>

PfoAnalysis pfoAnalysis;

//------------------------------------------------------------------------------------------------------------------------------------------

PfoAnalysis::PfoAnalysis() :
    Processor("PfoAnalysis"),
    m_nRun(0),
    m_nEvt(0),
    m_nRunSum(0),
    m_nEvtSum(0),
    m_printing(0),
    m_lookForQuarksWithMotherZ(0),
    m_mcPfoSelectionRadius(500.f),
    m_mcPfoSelectionMomentum(0.01f),
    m_mcPfoSelectionLowEnergyNPCutOff(1.2f),
    m_nPfosTotal(0),
    m_nPfosNeutralHadrons(0),
    m_nPfosPhotons(0),
    m_nPfosTracks(0),
    m_pfoEnergyTotal(0.f),
    m_pfoEnergyNeutralHadrons(0.f),
    m_pfoEnergyPhotons(0.f),
    m_pfoEnergyTracks(0.f),
    m_pfoECalToEmEnergy(0.f),
    m_pfoECalToHadEnergy(0.f),
    m_pfoHCalToEmEnergy(0.f),
    m_pfoHCalToHadEnergy(0.f),
    m_pfoMuonToEnergy(0.f),
    m_pfoOtherEnergy(0.f),
    m_pfoMassTotal(0.f),
    m_nPfoTargetsTotal(0),
    m_nPfoTargetsNeutralHadrons(0),
    m_nPfoTargetsPhotons(0),
    m_nPfoTargetsTracks(0),
    m_pfoTargetsEnergyTotal(0.f),
    m_pfoTargetsEnergyNeutralHadrons(0.f),
    m_pfoTargetsEnergyPhotons(0.f),
    m_pfoTargetsEnergyTracks(0.f),
    m_mcEnergyENu(0.f),
    m_mcEnergyFwd(0.f),
    m_eQQ(-99.f),
    m_eQ1(-99.f),
    m_eQ2(-99.f),
    m_costQQ(-99.f),
    m_costQ1(-99.f),
    m_costQ2(-99.f),
    m_mQQ(-99.f),
    m_thrust(-99.f),
    m_qPdg(-99.f),
    m_pTFile(NULL),
    m_pTTree(NULL),
    m_hPfoEnergySum(NULL),
    m_hPfoEnergySumL7A(NULL),
    m_collectCalibrationDetails(0),
    m_pCalibrationHelper(NULL),
    m_calibrationHelperSettings()
{
    _description = "PfoAnalysis analyses output of PandoraPFANew";

    registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                            "PfoCollection",
                            "Name of input pfo collection",
                            m_inputPfoCollection,
                            std::string());

    registerInputCollection(LCIO::MCPARTICLE,
                            "MCParticleCollection",
                            "Name of mc particle collections",
                            m_mcParticleCollection,
                            std::string());

    registerProcessorParameter("LookForQuarksWithMotherZ",
                            "Flag to look for quarks with mother Z",
                            m_lookForQuarksWithMotherZ,
                            int(0));

    registerProcessorParameter("MCPfoSelectionRadius",
                            "MC pfo selection radius",
                            m_mcPfoSelectionRadius,
                            float(500.f));

    registerProcessorParameter("MCPfoSelectionMomentum",
                            "MC pfo selection momentum",
                            m_mcPfoSelectionMomentum,
                            float(0.01f));

    registerProcessorParameter("MCPfoSelectionLowEnergyNPCutOff",
                            "MC pfo selection neutron and proton low energy cut-off",
                            m_mcPfoSelectionLowEnergyNPCutOff,
                            float(1.2f));

    registerProcessorParameter("RootFile",
                            "Name of the output root file",
                            m_rootFile,
                            std::string("PFOAnalysis.root"));

    registerProcessorParameter("Printing",
                           "Set the debug print level",
                           m_printing,
                           int(0));

    registerProcessorParameter("CollectCalibrationDetails",
                           "Whether to collect calibration details",
                           m_collectCalibrationDetails,
                           int(0));

    registerProcessorParameter("HCalRingOuterSymmetryOrder",
                           "Set the HCalRingOuterSymmetryOrder",
                           m_calibrationHelperSettings.m_hCalRingOuterSymmetryOrder,
                           int(8));

    registerProcessorParameter("HCalRingOuterPhi0",
                           "Set the HCalRingOuterPhi0",
                           m_calibrationHelperSettings.m_hCalRingOuterPhi0,
                           float(0));

    registerInputCollections(LCIO::CALORIMETERHIT,
                           "ECalCollections", 
                           "Name of the ECal collection of calo hits used to form clusters",
                           m_calibrationHelperSettings.m_eCalCollections,
                           LCStrVec());

    registerInputCollections(LCIO::CALORIMETERHIT,
                           "HCalCollections", 
                           "Name of the HCal collection of calo hits used to form clusters",
                           m_calibrationHelperSettings.m_hCalCollections,
                           LCStrVec());

    registerInputCollections(LCIO::CALORIMETERHIT,
                           "MuonCollections", 
                           "Name of the Muon collection of calo hits used to form clusters",
                           m_calibrationHelperSettings.m_muonCollections,
                           LCStrVec());

    registerInputCollections(LCIO::CALORIMETERHIT,
                           "BCalCollections", 
                           "Name of the BCal collection of calo hits used to form clusters",
                           m_calibrationHelperSettings.m_bCalCollections,
                           LCStrVec());

    registerInputCollections(LCIO::CALORIMETERHIT,
                           "LHCalCollections", 
                           "Name of the LHCal collection of calo hits used to form clusters",
                           m_calibrationHelperSettings.m_lHCalCollections,
                           LCStrVec());

    registerInputCollections(LCIO::CALORIMETERHIT,
                           "LCalCollections", 
                           "Name of the LCal collection of calo hits used to form clusters",
                           m_calibrationHelperSettings.m_lCalCollections,
                           LCStrVec());

    registerInputCollections(LCIO::SIMCALORIMETERHIT, 
                           "ECalCollectionsSimCaloHit" , 
                           "Name of the ECal collection post Mokka, pre digitisation" , 
                           m_calibrationHelperSettings.m_eCalCollectionsSimCaloHit , 
                           LCStrVec());

    registerInputCollections(LCIO::SIMCALORIMETERHIT, 
                           "HCalBarrelCollectionsSimCaloHit" , 
                           "Name of the HCal Barrel collection post Mokka, pre digitisation" , 
                           m_calibrationHelperSettings.m_hCalBarrelCollectionsSimCaloHit , 
                           LCStrVec());

    registerInputCollections(LCIO::SIMCALORIMETERHIT, 
                           "HCalEndCapCollectionsSimCaloHit" , 
                           "Name of the HCal EndCap collection post Mokka, pre digitisation" , 
                           m_calibrationHelperSettings.m_hCalEndCapCollectionsSimCaloHit , 
                           LCStrVec());

    registerInputCollections(LCIO::SIMCALORIMETERHIT, 
                           "HCalOtherCollectionsSimCaloHit" , 
                           "Name of the HCal Other collection post Mokka, pre digitisation" , 
                           m_calibrationHelperSettings.m_hCalOtherCollectionsSimCaloHit , 
                           LCStrVec());
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

    m_pTTree = new TTree("PfoAnalysisTree", "PfoAnalysisTree");
    m_pTTree->SetDirectory(m_pTFile);
    m_pTTree->Branch("run", &m_nRun, "run/I");
    m_pTTree->Branch("event", &m_nEvt, "event/I");
    m_pTTree->Branch("nPfosTotal", &m_nPfosTotal, "nPfosTotal/I");
    m_pTTree->Branch("nPfosNeutralHadrons", &m_nPfosNeutralHadrons, "nPfosNeutralHadrons/I");
    m_pTTree->Branch("nPfosPhotons", &m_nPfosPhotons, "nPfosPhotons/I");
    m_pTTree->Branch("nPfosTracks", &m_nPfosTracks, "nPfosTracks/I");
    m_pTTree->Branch("pfoEnergyTotal", &m_pfoEnergyTotal, "pfoEnergyTotal/F");
    m_pTTree->Branch("pfoEnergyNeutralHadrons", &m_pfoEnergyNeutralHadrons, "pfoEnergyNeutralHadrons/F");
    m_pTTree->Branch("pfoEnergyPhotons", &m_pfoEnergyPhotons, "pfoEnergyPhotons/F");
    m_pTTree->Branch("pfoEnergyTracks", &m_pfoEnergyTracks, "pfoEnergyTracks/F");
    m_pTTree->Branch("pfoECalToEmEnergy", &m_pfoECalToEmEnergy, "pfoECalToEmEnergy/F");
    m_pTTree->Branch("pfoECalToHadEnergy", &m_pfoECalToHadEnergy, "pfoECalToHadEnergy/F");
    m_pTTree->Branch("pfoHCalToEmEnergy", &m_pfoHCalToEmEnergy, "pfoHCalToEmEnergy/F");
    m_pTTree->Branch("pfoHCalToHadEnergy", &m_pfoHCalToHadEnergy, "pfoHCalToHadEnergy/F");
    m_pTTree->Branch("pfoOtherEnergy", &m_pfoOtherEnergy, "pfoOtherEnergy/F");
    m_pTTree->Branch("pfoMuonToEnergy", &m_pfoMuonToEnergy, "pfoMuonToEnergy/F");
    m_pTTree->Branch("pfoMassTotal", &m_pfoMassTotal, "pfoMassTotal/F");
    m_pTTree->Branch("pfoEnergies", &m_pfoEnergies);
    m_pTTree->Branch("pfoPx", &m_pfoPx);
    m_pTTree->Branch("pfoPy", &m_pfoPy);
    m_pTTree->Branch("pfoPz", &m_pfoPz);
    m_pTTree->Branch("pfoCosTheta", &m_pfoCosTheta);
    m_pTTree->Branch("pfoTargetEnergies", &m_pfoTargetEnergies);
    m_pTTree->Branch("pfoTargetPx", &m_pfoTargetPx);
    m_pTTree->Branch("pfoTargetPy", &m_pfoTargetPy);
    m_pTTree->Branch("pfoTargetPz", &m_pfoTargetPz);
    m_pTTree->Branch("pfoTargetCosTheta", &m_pfoTargetCosTheta);
    m_pTTree->Branch("pfoPdgCodes", &m_pfoPdgCodes);
    m_pTTree->Branch("pfoTargetPdgCodes", &m_pfoTargetPdgCodes);
    m_pTTree->Branch("nPfoTargetsTotal", &m_nPfoTargetsTotal, "nPfoTargetsTotal/I");
    m_pTTree->Branch("nPfoTargetsNeutralHadrons", &m_nPfoTargetsNeutralHadrons, "nPfoTargetsNeutralHadrons/I");
    m_pTTree->Branch("nPfoTargetsPhotons", &m_nPfoTargetsPhotons, "nPfoTargetsPhotons/I");
    m_pTTree->Branch("nPfoTargetsTracks", &m_nPfoTargetsTracks, "nPfoTargetsTracks/I");
    m_pTTree->Branch("pfoTargetsEnergyTotal", &m_pfoTargetsEnergyTotal, "pfoTargetsEnergyTotal/F");
    m_pTTree->Branch("pfoTargetsEnergyNeutralHadrons", &m_pfoTargetsEnergyNeutralHadrons, "pfoTargetsEnergyNeutralHadrons/F");
    m_pTTree->Branch("pfoTargetsEnergyPhotons", &m_pfoTargetsEnergyPhotons, "pfoTargetsEnergyPhotons/F");
    m_pTTree->Branch("pfoTargetsEnergyTracks", &m_pfoTargetsEnergyTracks, "pfoTargetsEnergyTracks/F");
    m_pTTree->Branch("mcEnergyENu", &m_mcEnergyENu, "mcEnergyENu/F");
    m_pTTree->Branch("mcEnergyFwd", &m_mcEnergyFwd, "mcEnergyFwd/F");
    m_pTTree->Branch("eQQ", &m_eQQ, "eQQ/F");
    m_pTTree->Branch("eQ1", &m_eQ1, "eQ1/F");
    m_pTTree->Branch("eQ2", &m_eQ2, "eQ2/F");
    m_pTTree->Branch("costQQ", &m_costQQ, "costQQ/F");
    m_pTTree->Branch("costQ1", &m_costQ1, "costQ1/F");
    m_pTTree->Branch("costQ2", &m_costQ2, "costQ2/F");
    m_pTTree->Branch("mQQ", &m_mQQ, "mQQ/F");
    m_pTTree->Branch("thrust", &m_thrust, "thrust/F");
    m_pTTree->Branch("qPdg", &m_qPdg, "qPdg/I");

    if (m_collectCalibrationDetails)
        m_pCalibrationHelper = new pandora_analysis::CalibrationHelper(m_calibrationHelperSettings);

    if (m_pCalibrationHelper)
    {
        m_pCalibrationHelper->SetBranchAddresses(m_pTTree);
        m_pCalibrationHelper->CreateHistograms();
        m_pCalibrationHelper->SetHistogramDirectories(m_pTFile);
    }

    m_hPfoEnergySum = new TH1F("fPFA", "total pfo energy", 10000, 0., 5000.);
    m_hPfoEnergySumL7A = new TH1F("fPFA_L7A", "total pfo energy < 0.7 A", 10000, 0., 5000.);
    m_hPfoEnergySum->SetDirectory(m_pTFile);
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

    if ((m_nEvtSum % 100) == 0)
        std::cout << " processed events: " << m_nEvtSum << std::endl;

    this->Clear();
    this->ExtractCollections(pLCEvent);
    this->MakeQuarkVariables(pLCEvent);
    this->PerformPfoAnalysis();

    if (m_pCalibrationHelper)
        m_pCalibrationHelper->Calibrate(pLCEvent, m_pfoVector, m_nPfoTargetsTotal, m_nPfoTargetsTracks, m_nPfoTargetsNeutralHadrons, m_pfoTargetsEnergyTotal);

    m_pTTree->Fill();
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

    m_pTFile->cd();
    m_pTTree->Write();
    m_hPfoEnergySum->Write();
    m_hPfoEnergySumL7A->Write();

    if (m_pCalibrationHelper)
        m_pCalibrationHelper->WriteHistograms();

    m_pTFile->Close();
    delete m_pTFile;

    if (m_pCalibrationHelper)
        delete m_pCalibrationHelper;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::Clear()
{
    m_pfoVector.clear();
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

    if (m_pCalibrationHelper)
        m_pCalibrationHelper->Clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::ExtractCollections(EVENT::LCEvent *pLCEvent)
{
    // Extract reconstructed pfo collection
    try
    {
        const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(m_inputPfoCollection);

        for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
        {
            const EVENT::ReconstructedParticle *pReconstructedParticle = dynamic_cast<EVENT::ReconstructedParticle*>(pLCCollection->getElementAt(i));

            if (NULL == pReconstructedParticle)
                throw EVENT::Exception("Collection type mismatch");

            m_pfoVector.push_back(pReconstructedParticle);
        }
    }
    catch (...)
    {
        streamlog_out(ERROR) << "Could not extract input particle collection: " << m_inputPfoCollection << std::endl;
    }

    // Extract mc pfo collection
    MCParticleList mcPfoCandidates;

    try
    {
        const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(m_mcParticleCollection);

        for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
        {
            const EVENT::MCParticle *pMCParticle = dynamic_cast<EVENT::MCParticle*>(pLCCollection->getElementAt(i));

            if (NULL == pMCParticle)
                throw EVENT::Exception("Collection type mismatch");

            if (!pMCParticle->getParents().empty())
                continue;

            this->ApplyPfoSelectionRules(pMCParticle, mcPfoCandidates);
        }
    }
    catch (...)
    {
        streamlog_out(WARNING) << "Could not extract mc particle collection " << m_mcParticleCollection << std::endl;
    }

    m_pfoTargetVector.insert(m_pfoTargetVector.begin(), mcPfoCandidates.begin(), mcPfoCandidates.end());
    std::sort(m_pfoTargetVector.begin(), m_pfoTargetVector.end(), PfoAnalysis::SortPfoTargetsByEnergy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::ApplyPfoSelectionRules(const EVENT::MCParticle *pMCParticle, MCParticleList &mcPfoCandidates) const
{
    const float innerRadius(std::sqrt(pMCParticle->getVertex()[0] * pMCParticle->getVertex()[0] + pMCParticle->getVertex()[1] * pMCParticle->getVertex()[1] + pMCParticle->getVertex()[2] * pMCParticle->getVertex()[2]));
    const float outerRadius(std::sqrt(pMCParticle->getEndpoint()[0] * pMCParticle->getEndpoint()[0] + pMCParticle->getEndpoint()[1] * pMCParticle->getEndpoint()[1] + pMCParticle->getEndpoint()[2] * pMCParticle->getEndpoint()[2]));
    const float momentum(std::sqrt(pMCParticle->getMomentum()[0] * pMCParticle->getMomentum()[0] + pMCParticle->getMomentum()[1] * pMCParticle->getMomentum()[1] + pMCParticle->getMomentum()[2] * pMCParticle->getMomentum()[2]));

    if ((mcPfoCandidates.find(pMCParticle) == mcPfoCandidates.end()) &&
        (outerRadius > m_mcPfoSelectionRadius) && (innerRadius <= m_mcPfoSelectionRadius) && (momentum > m_mcPfoSelectionMomentum) &&
        !((pMCParticle->getPDG() == 2212 || pMCParticle->getPDG() == 2112) && (pMCParticle->getEnergy() < m_mcPfoSelectionLowEnergyNPCutOff)))
    {
        mcPfoCandidates.insert(pMCParticle);
    }
    else
    {
        const EVENT::MCParticleVec &daughterVector(pMCParticle->getDaughters());

        for (EVENT::MCParticleVec::const_iterator iter = daughterVector.begin(), iterEnd = daughterVector.end(); iter != iterEnd; ++iter)
        {
            this->ApplyPfoSelectionRules(*iter, mcPfoCandidates);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoAnalysis::MakeQuarkVariables(EVENT::LCEvent *pLCEvent)
{
    MCParticleVector mcQuarkVector;

    try
    {
        const EVENT::LCCollection *pLCCollection = pLCEvent->getCollection(m_mcParticleCollection);

        for (unsigned int i = 0, nElements = pLCCollection->getNumberOfElements(); i < nElements; ++i)
        {
            const EVENT::MCParticle *pMCParticle = dynamic_cast<EVENT::MCParticle*>(pLCCollection->getElementAt(i));

            if (NULL == pMCParticle)
                throw EVENT::Exception("Collection type mismatch");

            const int absPdgCode(std::abs(pMCParticle->getPDG()));

            // By default, the primary quarks are the ones without any parents
            if (!m_lookForQuarksWithMotherZ)
            {
                if ((absPdgCode >= 1) && (absPdgCode <= 6) && pMCParticle->getParents().empty())
                    mcQuarkVector.push_back(pMCParticle);
            }
            else
            {
                // For MC files generated in the SLIC environment, the primary quarks have parents; the mother should be the Z-boson
                if ((absPdgCode >= 1) && (absPdgCode <= 6))
                {
                    if ((pMCParticle->getParents().size() == 1) && ((pMCParticle->getParents())[0]->getPDG() == 23))
                        mcQuarkVector.push_back(pMCParticle);
                }
            }
        }
    }
    catch (...)
    {
        streamlog_out(WARNING) << "Could not extract mc quark information" << std::endl;
    }

    if (!mcQuarkVector.empty())
    {
        m_qPdg = std::abs(mcQuarkVector[0]->getPDG());
        float energyTot(0.f);
        float costTot(0.f);

        for (unsigned int i = 0; i < mcQuarkVector.size(); ++i)
        {
            const float px(mcQuarkVector[i]->getMomentum()[0]);
            const float py(mcQuarkVector[i]->getMomentum()[1]);
            const float pz(mcQuarkVector[i]->getMomentum()[2]);
            const float energy(mcQuarkVector[i]->getEnergy());
            const float p(std::sqrt(px * px + py * py + pz * pz));
            const float cost(std::fabs(pz) / p);
            energyTot += energy;
            costTot += cost * energy;
        }

        m_thrust = costTot / energyTot;
    }

    if (mcQuarkVector.size() == 2)
    {
        const float pQ1x = mcQuarkVector[0]->getMomentum()[0];
        const float pQ1y = mcQuarkVector[0]->getMomentum()[1];
        const float pQ1z = mcQuarkVector[0]->getMomentum()[2];

        const float pQ2x = mcQuarkVector[1]->getMomentum()[0];
        const float pQ2y = mcQuarkVector[1]->getMomentum()[1];
        const float pQ2z = mcQuarkVector[1]->getMomentum()[2];

        const float pQ1[3] = {pQ1x, pQ1y, pQ1z};
        const float pQ2[3] = {pQ2x, pQ2y, pQ2z};
        const float pQQ[3] = {pQ1[0] + pQ2[0], pQ1[1] + pQ2[1], pQ1[2] + pQ2[2]};

        const TLorentzVector q1(pQ1[0], pQ1[1], pQ1[2], mcQuarkVector[0]->getEnergy());
        const TLorentzVector q2(pQ2[0], pQ2[1], pQ2[2], mcQuarkVector[1]->getEnergy());
        const TLorentzVector qq = q1 + q2;

        m_mQQ = qq.M();
        m_eQQ = mcQuarkVector[0]->getEnergy() + mcQuarkVector[1]->getEnergy();

        const float pQ1Tot(std::sqrt(pQ1[0] * pQ1[0] + pQ1[1] * pQ1[1] + pQ1[2] * pQ1[2]));
        const float pQ2Tot(std::sqrt(pQ2[0] * pQ2[0] + pQ2[1] * pQ2[1] + pQ2[2] * pQ2[2]));
        const float pQQTot(std::sqrt(pQQ[0] * pQQ[0] + pQQ[1] * pQQ[1] + pQQ[2] * pQQ[2]));

        if (std::fabs(pQQTot) > std::numeric_limits<float>::epsilon())
            m_costQQ = pQQ[2] / pQQTot;

        if (std::fabs(pQ1Tot) > std::numeric_limits<float>::epsilon())
            m_costQ1 = pQ1[2] / pQ1Tot;

        if (std::fabs(pQ2Tot) > std::numeric_limits<float>::epsilon())
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
        const EVENT::ReconstructedParticle *pPfo = *iter;

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
            const EVENT::ClusterVec &clusterVec(pPfo->getClusters());

            for (EVENT::ClusterVec::const_iterator iter = clusterVec.begin(), iterEnd = clusterVec.end(); iter != iterEnd; ++iter)
            {
                const EVENT::CalorimeterHitVec &calorimeterHitVec((*iter)->getCalorimeterHits());

                for (EVENT::CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
                {
                    const EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;
                    cellEnergySum += pCalorimeterHit->getEnergy();
                }
            }

            const float correctionFactor((cellEnergySum < std::numeric_limits<float>::epsilon()) ? 0.f : pPfo->getEnergy() / cellEnergySum);

            if (22 == pPfo->getType())
            {
                // Photons
                ++m_nPfosPhotons;
                m_pfoEnergyPhotons += pPfo->getEnergy();

                for (EVENT::ClusterVec::const_iterator iter = clusterVec.begin(), iterEnd = clusterVec.end(); iter != iterEnd; ++iter)
                {
                    const EVENT::CalorimeterHitVec &calorimeterHitVec((*iter)->getCalorimeterHits());

                    for (EVENT::CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
                    {
                        const EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;

                        const float hitEnergy(correctionFactor * pCalorimeterHit->getEnergy());
                        const CHT cht(pCalorimeterHit->getType());

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

                for (EVENT::ClusterVec::const_iterator iter = clusterVec.begin(), iterEnd = clusterVec.end(); iter != iterEnd; ++iter)
                {
                    const EVENT::CalorimeterHitVec &calorimeterHitVec((*iter)->getCalorimeterHits());

                    for (EVENT::CalorimeterHitVec::const_iterator hitIter = calorimeterHitVec.begin(), hitIterEnd = calorimeterHitVec.end(); hitIter != hitIterEnd; ++hitIter)
                    {
                        const EVENT::CalorimeterHit *pCalorimeterHit = *hitIter;

                        const float hitEnergy(correctionFactor * pCalorimeterHit->getEnergy());
                        const CHT cht(pCalorimeterHit->getType());

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

    // Extract quantities relating to pfo targets
    for (MCParticleVector::const_iterator iter = m_pfoTargetVector.begin(), iterEnd = m_pfoTargetVector.end(); iter != iterEnd; ++iter)
    {
        const EVENT::MCParticle *pMCParticle = *iter;

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

        if (std::fabs(cosTheta) > 0.98f)
            m_mcEnergyFwd += pMCParticle->getEnergy();

        if ((std::abs(pMCParticle->getPDG()) == 12) || (std::abs(pMCParticle->getPDG()) == 14) || (std::abs(pMCParticle->getPDG()) == 16))
            m_mcEnergyENu += pMCParticle->getEnergy();

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

    if ((m_qPdg >= 1) && (m_qPdg <= 3) && (m_thrust <= 0.7f))
        m_hPfoEnergySumL7A->Fill(m_pfoEnergyTotal + m_mcEnergyENu, 1.);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PfoAnalysis::SortPfoTargetsByEnergy(const EVENT::MCParticle *const pLhs, const EVENT::MCParticle *const pRhs)
{
    return (pLhs->getEnergy() > pRhs->getEnergy());
}
