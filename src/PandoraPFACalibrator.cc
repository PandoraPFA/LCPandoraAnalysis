/**
 *  @file   PandoraAnalysis/src/PandoraPFACalibrator.cc
 * 
 *  @brief  Implementation of pandora pfa calibrator class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCStrVec.h"
#include "EVENT/ReconstructedParticle.h"
#include "UTIL/CellIDDecoder.h"

#include "gear/GEAR.h"
#include "gear/CalorimeterParameters.h"

#include "marlin/Global.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "PandoraPFACalibrator.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <set>

PandoraPFACalibrator aPandoraPFACalibrator;

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraPFACalibrator::PandoraPFACalibrator() :
    Processor("PandoraPFACalibrator"),
    m_nRun(0),
    m_nEvt(0),
    m_ecalToMIP(1.f),
    m_hcalToMIP(1.f),
    m_muonToMIP(1.f),
    m_ecalToEMGeVCalibration(1.f),
    m_hcalToHadGeVCalibration(1.f),
    m_ecalToHadGeVCalibrationBarrel(1.f),
    m_ecalToHadGeVCalibrationEndCap(1.f),
    m_hcalToEMGeVCalibration(1.f),
    m_maxHCalHitHadronicEnergy(10000),
    m_zOfEndCap(std::numeric_limits<float>::max()),
    m_pTFile(NULL),
    m_hPfoEnergy(NULL),
    m_hPfoEnergyBarrel(NULL),
    m_hPfoEnergy95ECal(NULL),
    m_hPfoEnergy95HCal(NULL),
    m_hPfoEnergy95Muon(NULL),
    m_hPfoEnergyVsCosTheta(NULL),
    m_hPfoEnergyVsCosThetaReco(NULL),
    m_hCaloEnergy(NULL),
    m_hCaloEnergyECal(NULL),
    m_hCaloEnergyHCal(NULL),
    m_hCaloEnergyMuon(NULL),
    m_hCaloEnergy95ECal(NULL),
    m_hCaloEnergy95HCal(NULL),
    m_hCaloEnergy95Muon(NULL),
    m_hEcalBarrelEnergyByLayer(NULL),
    m_hEcalEndCapEnergyByLayer(NULL),
    m_hECalHCalEnergyEM(NULL),
    m_hECalHcalEnergyHAD(NULL),
    m_hECalBarrelHCalEnergyEM(NULL),
    m_hECalEndCapHCalEnergyEM(NULL),
    m_hECalBarrelHCalEnergyHAD(NULL),
    m_hECalEndCapHCalEnergyHAD(NULL),
    m_hCaloEnergyVsCosTheta(NULL),
    m_hCaloEnergyVsCosThetaReco(NULL),
    m_hECalBarrelMIP(NULL),
    m_hECalEndCapMIP(NULL),
    m_hHCalMIP(NULL),
    m_hMuonMIP(NULL),
    m_hECalBarrelMIPCorr(NULL),
    m_hECalEndCapMIPCorr(NULL),
    m_hHCalMIPCorr(NULL),
    m_hMuonMIPCorr(NULL)
{
    _description = "PandoraPFACalibrator for calibration of PandoraPFA";

    registerProcessorParameter("RootFile",
        "Output root file name",
        m_rootFile,
        std::string("PandoraPFACalibrator.root"));

    LCStrVec mcPfoCollections;
    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
        "MCPfoCollections",
        "Names of mc pfo collections",
        m_mcPfoCollections,
        mcPfoCollections);

    LCStrVec reconstructedPfoCollections;
    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
        "ReconstructedPfoCollections",
        "Names of reconstructed particle collections",
        m_recoPfoCollections,
        reconstructedPfoCollections);

    LCStrVec inputMCParticleCollections;
    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
        "InputMCParticleCollections",
        "Names of mc pfo collections",
        m_inputMCParticleCollections,
        inputMCParticleCollections);

    registerProcessorParameter("InputParticleCollectionName",
        "Name of reconstructed particle collection",
        m_particleCollectionName,
        std::string("PandoraPFANewPFOs"));

    LCStrVec ecalBarrelCollections;
    registerInputCollections(LCIO::CALORIMETERHIT,
        "ECALBarrelcollections", 
        "Name of the ECAL barrel collection used to form clusters",
        m_ecalBarrelCollections,
        ecalBarrelCollections);

    LCStrVec ecalEndCapCollections;
    registerInputCollections(LCIO::CALORIMETERHIT,
        "ECALEndCapcollections", 
        "Name of the ECAL EndCap collection used to form clusters",
        m_ecalEndCapCollections,
        ecalEndCapCollections);

    LCStrVec hcalCollections;
    registerInputCollections(LCIO::CALORIMETERHIT,
        "HCALcollections", 
        "Name of the HCAL collection used to form clusters",
        m_hcalCollections,
        hcalCollections);

    LCStrVec muonCollections;
    registerInputCollections(LCIO::CALORIMETERHIT,
        "MUONcollections", 
        "Name of the MUON collection used to form clusters",
        m_muonCollections,
        muonCollections);

    LCStrVec bcalCollections;
    registerInputCollections(LCIO::CALORIMETERHIT,
        "BCALcollections", 
        "Name of the BCAL collection used to form clusters",
        m_bcalCollections,
        bcalCollections);

    LCStrVec lhcalCollections;
    registerInputCollections(LCIO::CALORIMETERHIT,
        "LHCALcollections", 
        "Name of the LHCAL collection used to form clusters",
        m_lhcalCollections,
        lhcalCollections);

    LCStrVec lcalCollections;
    registerInputCollections(LCIO::CALORIMETERHIT,
        "LCALcollections", 
        "Name of the LCAL collection used to form clusters",
        m_lcalCollections,
        lcalCollections);

    registerProcessorParameter("ECalToMipCalibration",
        "Calibration from deposited ECAL energy to MIP",
        m_ecalToMIP,
        1.f);

    registerProcessorParameter("HCalToMipCalibration",
        "Calibration from deposited HCAL energy to MIP",
        m_hcalToMIP,
        1.f);

    registerProcessorParameter("MuonToMipCalibration",
        "Calibration from deposited MUON energy to MIP",
        m_muonToMIP,
        1.f);

    registerProcessorParameter("ECalToEMGeVCalibration",
        "Calibration from deposited ECAL to EM energy",
        m_ecalToEMGeVCalibration,
        1.f);

    registerProcessorParameter("HCalToHadGeVCalibration",
        "Calibration from deposited HCAL to Hadronic energy",
        m_hcalToHadGeVCalibration,
        1.f);

    registerProcessorParameter("ECalToHadGeVCalibrationBarrel",
        "Calibration from deposited ECAL barrel to Hadronic energy",
        m_ecalToHadGeVCalibrationBarrel,
        1.f);

    registerProcessorParameter("ECalToHadGeVCalibrationEndCap",
        "Calibration from deposited ECAL endcap to Hadronic energy",
        m_ecalToHadGeVCalibrationEndCap,
        1.f);

    registerProcessorParameter("HCalToEMGeVCalibration",
        "Calibration from deposited HCAL to EM energy",
        m_hcalToEMGeVCalibration,
        1.f);

    registerProcessorParameter("MaxHCalHitHadronicEnergy",
        "The maximum hadronic energy allowed for a single hcal hit",
        m_maxHCalHitHadronicEnergy,
        10000.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::init() 
{
    // Finalise steering - allow use of old and new config for input pfo lists
    typedef std::set<std::string> StringSet;
    StringSet uniqueMCPfoCollections;
    uniqueMCPfoCollections.insert(m_mcPfoCollections.begin(), m_mcPfoCollections.end());
    uniqueMCPfoCollections.insert(m_inputMCParticleCollections.begin(), m_inputMCParticleCollections.end());
    m_mcPfoCollections.clear();
    m_mcPfoCollections.insert(m_mcPfoCollections.begin(), uniqueMCPfoCollections.begin(), uniqueMCPfoCollections.end());

    StringSet uniqueRecoPfoCollections;
    uniqueRecoPfoCollections.insert(m_recoPfoCollections.begin(), m_recoPfoCollections.end());
    uniqueRecoPfoCollections.insert(m_particleCollectionName);
    m_recoPfoCollections.clear();
    m_recoPfoCollections.insert(m_recoPfoCollections.begin(), uniqueRecoPfoCollections.begin(), uniqueRecoPfoCollections.end());

    CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");

    this->printParameters();

    // Member variable initialisation
    m_nRun = 0;
    m_nEvt = 0;

    m_pTFile = new TFile(m_rootFile.c_str(), "recreate");

    // PFA total energy
    m_hPfoEnergy = new TH1F("PfoEnergy", "Pfo energy total", 2000, 0., 250.);
    m_hPfoEnergyBarrel = new TH1F("PfoEnergyBarrel", "Pfo energy barrel", 2000, 0., 250.);
    m_hPfoEnergy95ECal = new TH1F("PfoEnergy95ECal", "Pfo energy, events with 95% energy in ECal", 2000, 0., 250.);
    m_hPfoEnergy95HCal = new TH1F("PfoEnergy95HCal", "Pfo energy, events with 95% energy in HCal", 2000, 0., 250.);
    m_hPfoEnergy95Muon = new TH1F("PfoEnergy95Muon", "Pfo energy, events with 95% energy in Muon", 2000, 0., 250.);
    m_hPfoEnergyVsCosTheta = new TH2F("PfoEnergyVsCosTheta", "Pfo energy vs cos theta", 500, 0., 1., 2000, 0., 250.);
    m_hPfoEnergyVsCosThetaReco = new TH2F("PfoEnergyVsCosThetaReco", "Pfo energy vs cos theta reco", 500, 0., 1., 2000, 0., 250.);

    // ECal total energy
    m_hCaloEnergy = new TH1F("CaloEnergy", "Calorimeter energy total", 2000, 0., 250.);
    m_hCaloEnergyECal = new TH1F("CaloEnergyECal", "Calorimeter energy total ECal", 2000, 0., 250.);
    m_hCaloEnergyHCal = new TH1F("CaloEnergyHCal", "Calorimeter energy total HCal", 2000, 0., 250.);
    m_hCaloEnergyMuon = new TH1F("CaloEnergyMuon", "Calorimeter energy total Muon", 2000, 0., 250.);
    m_hCaloEnergy95ECal = new TH1F("CaloEnergy95ECal", "Calorimeter energy, events with 95% energy in ECal", 2000, 0., 250.);
    m_hCaloEnergy95HCal = new TH1F("CaloEnergy95HCal", "Calorimeter energy, events with 95% energy in HCal", 2000, 0., 250.);
    m_hCaloEnergy95Muon = new TH1F("CaloEnergy95Muon", "Calorimeter energy, events with 95% energy in Muon", 2000, 0., 250.);
    m_hEcalBarrelEnergyByLayer = new TH1F("ECalBarrelEnergyByLayer", "ECal barrel energy profile", 100, 0., 100.);
    m_hEcalEndCapEnergyByLayer = new TH1F("ECalEndCapEnergyByLayer", "ECal endcap energy profile", 100, 0., 100.);
    m_hECalHCalEnergyEM = new TH2F("ECalHCalEnergyEM", "ECal vs HCal energy EM", 2000, 0., 250., 2000, 0., 250.);
    m_hECalHcalEnergyHAD = new TH2F("ECalHcalEnergyHAD", "ECal vs HCal energy HAD", 2000, 0., 250., 2000, 0., 250.);
    m_hECalBarrelHCalEnergyEM = new TH2F("ECalBarrelHCalEnergyEM", "ECal barrel vs HCal energy EM", 2000, 0., 250., 2000, 0., 250.);
    m_hECalEndCapHCalEnergyEM = new TH2F("ECalEndCapHCalEnergyEM", "ECal endcap vs HCal energy EM", 2000, 0., 250., 2000, 0., 250.);
    m_hECalBarrelHCalEnergyHAD = new TH2F("ECalBarrelHCalEnergyHAD", "ECal barrel vs HCal energy HAD", 2000, 0., 250., 2000, 0., 250.);
    m_hECalEndCapHCalEnergyHAD = new TH2F("ECalEndCapHCalEnergyHAD", "ECal endcap vs HCal energy HAD", 2000, 0., 250., 2000, 0., 250.);
    m_hCaloEnergyVsCosTheta = new TH2F("CaloEnergyVsCosTheta", "Calorimeter energy vs cos theta", 500, 0., 1., 2000, 0., 250.);
    m_hCaloEnergyVsCosThetaReco = new TH2F("CaloEnergyVsCosThetaReco", "Calorimeter energy vs cos theta reco", 500, 0., 1.,2000, 0., 250.);

    // MIP Calibration
    m_hECalBarrelMIP = new TH1F("ECalBarrelMIP", "ECal barrel MIP", 200, 0., 5.);
    m_hECalEndCapMIP = new TH1F("ECalEndCapMIP", "ECal endcap MIP", 200, 0., 5.);
    m_hHCalMIP = new TH1F("HCalMIP", "HCal MIP", 200, 0., 5.);
    m_hMuonMIP = new TH1F("MuonMIP", "Muon MIP", 200, 0., 5.);
    m_hECalBarrelMIPCorr = new TH1F("ECalBarrelMIPCorr", "ECal barrel MIP, direction corrected", 200, 0., 5.);
    m_hECalEndCapMIPCorr = new TH1F("ECalEndCapMIPCorr", "ECal endcap MIP, direction corrected", 200, 0., 5.);
    m_hHCalMIPCorr = new TH1F("HCalMIPCorr", "HCal MIP, direction corrected", 200, 0., 5.);
    m_hMuonMIPCorr = new TH1F("MuonMIPCorr", "Muon MIP, direction corrected", 200, 0., 5.);

    m_hPfoEnergy->SetDirectory(m_pTFile);
    m_hPfoEnergyBarrel->SetDirectory(m_pTFile);
    m_hPfoEnergy95ECal->SetDirectory(m_pTFile);
    m_hPfoEnergy95HCal->SetDirectory(m_pTFile);
    m_hPfoEnergy95Muon->SetDirectory(m_pTFile);
    m_hPfoEnergyVsCosTheta->SetDirectory(m_pTFile);
    m_hPfoEnergyVsCosThetaReco->SetDirectory(m_pTFile);
    m_hCaloEnergy->SetDirectory(m_pTFile);
    m_hCaloEnergyECal->SetDirectory(m_pTFile);
    m_hCaloEnergyHCal->SetDirectory(m_pTFile);
    m_hCaloEnergyMuon->SetDirectory(m_pTFile);
    m_hCaloEnergy95ECal->SetDirectory(m_pTFile);
    m_hCaloEnergy95HCal->SetDirectory(m_pTFile);
    m_hCaloEnergy95Muon->SetDirectory(m_pTFile);
    m_hECalHCalEnergyEM->SetDirectory(m_pTFile);
    m_hECalHcalEnergyHAD->SetDirectory(m_pTFile);
    m_hEcalBarrelEnergyByLayer->SetDirectory(m_pTFile);
    m_hEcalEndCapEnergyByLayer->SetDirectory(m_pTFile);
    m_hECalBarrelHCalEnergyEM->SetDirectory(m_pTFile);
    m_hECalEndCapHCalEnergyEM->SetDirectory(m_pTFile);
    m_hECalBarrelHCalEnergyHAD->SetDirectory(m_pTFile);
    m_hECalEndCapHCalEnergyHAD->SetDirectory(m_pTFile);
    m_hCaloEnergyVsCosTheta->SetDirectory(m_pTFile);
    m_hCaloEnergyVsCosThetaReco->SetDirectory(m_pTFile);
    m_hECalBarrelMIP->SetDirectory(m_pTFile);
    m_hECalEndCapMIP->SetDirectory(m_pTFile);
    m_hHCalMIP->SetDirectory(m_pTFile);
    m_hMuonMIP->SetDirectory(m_pTFile);
    m_hECalBarrelMIPCorr->SetDirectory(m_pTFile);
    m_hECalEndCapMIPCorr->SetDirectory(m_pTFile);
    m_hHCalMIPCorr->SetDirectory(m_pTFile);
    m_hMuonMIPCorr->SetDirectory(m_pTFile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::processRunHeader(LCRunHeader *pLCRunHeader)
{
    m_nRun++;
    m_out(DEBUG) << " DETECTOR : " << pLCRunHeader->getDetectorName() << std::endl;

    const gear::CalorimeterParameters &ecalEndCapParameters(marlin::Global::GEAR->getEcalEndcapParameters());
    m_zOfEndCap = static_cast<float>(ecalEndCapParameters.getExtent()[2]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::processEvent(LCEvent *pLCEvent)
{
    ++m_nEvt;

    // Read mc pfo collection
    float cosTheta(std::numeric_limits<float>::max());
    this->ReadMCParticles(pLCEvent, m_mcPfoCollections, cosTheta);

    // Read hit collections
    float ecalBarrelEnergy(0.f), ecalEndCapEnergy(0.f);
    this->ReadHitEnergies(pLCEvent, m_ecalBarrelCollections, ecalBarrelEnergy, m_ecalToMIP, m_hECalBarrelMIP, m_hECalBarrelMIPCorr, "K-1", m_hEcalBarrelEnergyByLayer);
    this->ReadHitEnergies(pLCEvent, m_ecalEndCapCollections, ecalEndCapEnergy, m_ecalToMIP, m_hECalEndCapMIP, m_hECalEndCapMIPCorr, "K-1", m_hEcalEndCapEnergyByLayer);

    float hcalEnergy(0.f), muonEnergy(0.f);
    this->ReadHitEnergies(pLCEvent, m_hcalCollections, hcalEnergy, m_hcalToMIP, m_hHCalMIP, m_hHCalMIPCorr);
    this->ReadHitEnergies(pLCEvent, m_muonCollections, muonEnergy, m_muonToMIP, m_hMuonMIP, m_hMuonMIPCorr);

    float lcalEnergy(0.f), lhcalEnergy(0.f), bcalEnergy(0.f);
    this->ReadHitEnergies(pLCEvent, m_lcalCollections, lcalEnergy);
    this->ReadHitEnergies(pLCEvent, m_lhcalCollections, lhcalEnergy);
    this->ReadHitEnergies(pLCEvent, m_bcalCollections, bcalEnergy);

    // Read reconstructed pfo collection
    float cosThetaReco(0.f), totalPfoEnergy(0.f);
    this->ReadPfoCollections(pLCEvent, m_recoPfoCollections, totalPfoEnergy, cosThetaReco);

    const float totalCalEnergy(ecalBarrelEnergy + ecalEndCapEnergy + hcalEnergy + muonEnergy + lcalEnergy + bcalEnergy + lhcalEnergy);

    // Calorimeter energy histograms
    if (cosTheta < 0.95f)
    {
        m_hCaloEnergy->Fill(totalCalEnergy);
        m_hCaloEnergyECal->Fill(ecalBarrelEnergy + ecalEndCapEnergy);
        m_hCaloEnergyHCal->Fill(hcalEnergy);
        m_hCaloEnergyMuon->Fill(muonEnergy);

        m_hECalHCalEnergyEM->Fill(ecalBarrelEnergy + ecalEndCapEnergy, hcalEnergy * m_hcalToEMGeVCalibration);
        m_hECalHcalEnergyHAD->Fill(ecalBarrelEnergy * m_ecalToHadGeVCalibrationBarrel + ecalEndCapEnergy * m_ecalToHadGeVCalibrationEndCap, hcalEnergy);
        m_hECalBarrelHCalEnergyEM->Fill(ecalBarrelEnergy, hcalEnergy * m_hcalToEMGeVCalibration);
        m_hECalEndCapHCalEnergyEM->Fill(ecalEndCapEnergy, hcalEnergy * m_hcalToEMGeVCalibration);
        m_hECalBarrelHCalEnergyHAD->Fill((ecalBarrelEnergy * m_ecalToHadGeVCalibrationBarrel), hcalEnergy);
        m_hECalEndCapHCalEnergyHAD->Fill((ecalEndCapEnergy * m_ecalToHadGeVCalibrationEndCap), hcalEnergy);

        if ((totalCalEnergy > 0.f) && (((ecalBarrelEnergy + ecalEndCapEnergy) / totalCalEnergy) > 0.95f))
            m_hCaloEnergy95ECal->Fill(totalCalEnergy);

        if ((totalCalEnergy > 0.f) && ((hcalEnergy / totalCalEnergy) > 0.95f))
            m_hCaloEnergy95HCal->Fill(totalCalEnergy);

        if ((totalCalEnergy > 0.f) && ((muonEnergy / totalCalEnergy) > 0.95f))
            m_hCaloEnergy95Muon->Fill(totalCalEnergy);
    }

    m_hCaloEnergyVsCosTheta->Fill(cosTheta, totalCalEnergy);
    m_hCaloEnergyVsCosThetaReco->Fill(cosThetaReco, totalCalEnergy);

    // Particle flow object energy histograms
    if (cosTheta < 0.95f)
    {
        m_hPfoEnergy->Fill(totalPfoEnergy);

        if (cosTheta < 0.7f)
            m_hPfoEnergyBarrel->Fill(totalPfoEnergy);

        if ((totalCalEnergy > 0.f) && (((ecalBarrelEnergy + ecalEndCapEnergy) / totalCalEnergy) > 0.95f))
            m_hPfoEnergy95ECal->Fill(totalPfoEnergy);

        if ((totalCalEnergy > 0.f) && ((hcalEnergy / totalCalEnergy) > 0.95f))
            m_hPfoEnergy95HCal->Fill(totalPfoEnergy);

        if ((totalCalEnergy > 0.f) && ((muonEnergy / totalCalEnergy) > 0.95f))
            m_hPfoEnergy95Muon->Fill(totalPfoEnergy);
    }

    m_hPfoEnergyVsCosTheta->Fill(cosTheta, totalPfoEnergy);
    m_hPfoEnergyVsCosThetaReco->Fill(cosThetaReco, totalPfoEnergy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::ReadMCParticles(LCEvent *pLCEvent, const LCStrVec &collectionNames, float &cosTheta) const
{
    for (LCStrVec::const_iterator iter = collectionNames.begin(), iterEnd = collectionNames.end(); iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection *pLCCollection = pLCEvent->getCollection(iter->c_str());

            if (NULL != pLCCollection)
            {
                const int nElements(pLCCollection->getNumberOfElements());

                // Consider first particle in collection only
                if (nElements > 0)
                {
                    ReconstructedParticle *pReconstructedParticle = dynamic_cast<ReconstructedParticle*>(pLCCollection->getElementAt(0));

                    if (NULL == pReconstructedParticle)
                    {
                        m_out(ERROR) << "Collection type mismatch " << (*iter) << " expected to contain objects of type ReconstructedParticle " << std::endl;
                        throw;
                    }

                    const float px(pReconstructedParticle->getMomentum()[0]);
                    const float py(pReconstructedParticle->getMomentum()[1]);
                    const float pz(pReconstructedParticle->getMomentum()[2]);
                    const float p(std::sqrt(px * px + py * py + pz * pz));
                    cosTheta = ((p > 0.f) ? std::fabs(pz / p) : std::numeric_limits<float>::max());
                }
            }
        }
        catch (DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " << (*iter) << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::ReadHitEnergies(LCEvent *pLCEvent, const LCStrVec &collectionNames, float &hitEnergySum, const float mipConstant,
    TH1F *const pMipPlot, TH1F *const pMipPlotCorrected, const char *const pLayerEncoding, TH1F *const pEnergyByLayerPlot) const
{
    for (LCStrVec::const_iterator iter = collectionNames.begin(), iterEnd = collectionNames.end(); iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection *pLCCollection = pLCEvent->getCollection(iter->c_str());

            if (NULL != pLCCollection)
            {
                CellIDDecoder<CalorimeterHit> *pCellIDDecoder = NULL;

                if (NULL != pEnergyByLayerPlot)
                    pCellIDDecoder = new CellIDDecoder<CalorimeterHit>(pLCCollection);

                const int nElements(pLCCollection->getNumberOfElements());

                for (int iHit = 0; iHit < nElements; ++iHit)
                {
                    CalorimeterHit *pCalorimeterHit = dynamic_cast<CalorimeterHit*>(pLCCollection->getElementAt(iHit));

                    if (NULL == pCalorimeterHit)
                    {
                        m_out(ERROR) << "Collection type mismatch " << (*iter) << " expected to contain objects of type CalorimeterHit " << std::endl;
                        throw;
                    }

                    const float hitEnergy(pCalorimeterHit->getEnergy());
                    hitEnergySum += hitEnergy;

                    const float x(pCalorimeterHit->getPosition()[0]);
                    const float y(pCalorimeterHit->getPosition()[1]);
                    const float z(pCalorimeterHit->getPosition()[2]);
                    const float r(std::sqrt(x * x + y * y + z * z));
                    const float correction((std::fabs(z) < m_zOfEndCap) ? r / std::sqrt(x * x + y * y) : r / std::fabs(z));

                    const float energyInMips(pCalorimeterHit->getEnergy() * mipConstant);

                    if (NULL != pMipPlot)
                        pMipPlot->Fill(energyInMips);

                    if (NULL != pMipPlotCorrected)
                        pMipPlotCorrected->Fill(energyInMips / correction);

                    if (NULL != pEnergyByLayerPlot)
                    {
                        const int layerNumber((*pCellIDDecoder)(pCalorimeterHit)[pLayerEncoding] + 1);
                        pEnergyByLayerPlot->Fill(layerNumber, hitEnergy);
                    }
                }

                delete pCellIDDecoder;
            }
        }
        catch (DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " << (*iter) << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::ReadPfoCollections(LCEvent *pLCEvent, const LCStrVec &collectionNames, float &pfoEnergySum, float &cosTheta) const
{
    for (LCStrVec::const_iterator iter = collectionNames.begin(), iterEnd = collectionNames.end(); iter != iterEnd; ++iter)
    {
        try
        {
            LCCollection *pLCCollection = pLCEvent->getCollection(m_particleCollectionName.c_str());

            if (NULL != pLCCollection)
            {
                float px(0.f), py(0.f), pz(0.f);
                const int nElements(pLCCollection->getNumberOfElements());

                for (int iPfo = 0; iPfo < nElements; ++iPfo)
                {
                    ReconstructedParticle *pReconstructedParticle = dynamic_cast<ReconstructedParticle*>(pLCCollection->getElementAt(iPfo));

                    if (NULL == pReconstructedParticle)
                    {
                        m_out(ERROR) << "Collection type mismatch " << m_particleCollectionName << " expected to contain objects of type ReconstructedParticle " << std::endl;
                        throw;
                    }

                    pfoEnergySum += pReconstructedParticle->getEnergy();
                    px += pReconstructedParticle->getMomentum()[0];
                    py += pReconstructedParticle->getMomentum()[1];
                    pz += pReconstructedParticle->getMomentum()[2];
                }

                const float p(std::sqrt(px * px + py * py + pz * pz));
                cosTheta = ((p > 0.f) ? std::fabs(pz / p) : std::numeric_limits<float>::max());
            }
        }
        catch (DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " << (*iter) << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::check(LCEvent *pLCEvent)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::end()
{
    m_out(DEBUG) << "PandoraPFACalibrator::end()  " << name() << " processed " << m_nEvt << " events in " << m_nRun << " runs " << std::endl;

    m_pTFile->cd();
    m_hPfoEnergy->Write();
    m_hPfoEnergyBarrel->Write();
    m_hPfoEnergy95ECal->Write();
    m_hPfoEnergy95HCal->Write();
    m_hPfoEnergy95Muon->Write();
    m_hPfoEnergyVsCosTheta->Write();
    m_hPfoEnergyVsCosThetaReco->Write();
    m_hCaloEnergy->Write();
    m_hCaloEnergyECal->Write();
    m_hCaloEnergyHCal->Write();
    m_hCaloEnergyMuon->Write();
    m_hCaloEnergy95ECal->Write();
    m_hCaloEnergy95HCal->Write();
    m_hCaloEnergy95Muon->Write();
    m_hEcalBarrelEnergyByLayer->Write();
    m_hEcalEndCapEnergyByLayer->Write();
    m_hECalHCalEnergyEM->Write();
    m_hECalHcalEnergyHAD->Write();
    m_hECalBarrelHCalEnergyEM->Write();
    m_hECalEndCapHCalEnergyEM->Write();
    m_hECalBarrelHCalEnergyHAD->Write();
    m_hECalEndCapHCalEnergyHAD->Write();
    m_hCaloEnergyVsCosTheta->Write();
    m_hCaloEnergyVsCosThetaReco->Write();
    m_hECalBarrelMIP->Write();
    m_hECalEndCapMIP->Write();
    m_hHCalMIP->Write();
    m_hMuonMIP->Write();
    m_hECalBarrelMIPCorr->Write();
    m_hECalEndCapMIPCorr->Write();
    m_hHCalMIPCorr->Write();
    m_hMuonMIPCorr->Write();
    m_pTFile->Write();
    m_pTFile->Close();

    delete m_pTFile;
}
