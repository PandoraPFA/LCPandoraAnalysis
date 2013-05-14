/**
 *  @file   PandoraAnalysis/src/PandoraPFACalibrator.cc
 * 
 *  @brief  Implementation of pandora pfa calibrator class.
 * 
 *  $Log: $
 */

#include "EVENT/LCCollection.h"
#include "EVENT/LCRelation.h"
#include "EVENT/MCParticle.h"
#include "EVENT/Track.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCStrVec.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/ClusterImpl.h"
#include "UTIL/CellIDDecoder.h"

#include "gear/GEAR.h"
#include "gear/CalorimeterParameters.h"

#include "marlin/Global.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

#include "PandoraPFACalibrator.h"
#include "MCPfoMaker.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <set>

PandoraPFACalibrator aPandoraPFACalibrator;

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraPFACalibrator::PandoraPFACalibrator() : Processor("PandoraPFACalibrator") 
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
        -1.f);

    registerProcessorParameter("HCalToMipCalibration",
        "Calibration from deposited HCAL energy to MIP",
        m_hcalToMIP,
        -1.f);

    registerProcessorParameter("MuonToMipCalibration",
        "Calibration from deposited MUON energy to MIP",
        m_muonToMIP,
        -1.f);

    registerProcessorParameter("ECalToEMGeVCalibration",
        "Calibration from deposited ECAL to EM energy",
        m_ecalToEMGeVCalibration,
        -1.f);

    registerProcessorParameter("HCalToHadGeVCalibration",
        "Calibration from deposited HCAL to Hadronic energy",
        m_hcalToHadGeVCalibration,
        -1.f);

    registerProcessorParameter("ECalToHadGeVCalibrationBarrel",
        "Calibration from deposited ECAL barrel to Hadronic energy",
        m_ecalToHadGeVCalibrationBarrel,
        -1.f);

    registerProcessorParameter("ECalToHadGeVCalibrationEndCap",
        "Calibration from deposited ECAL endcap to Hadronic energy",
        m_ecalToHadGeVCalibrationEndCap,
        -1.f);

    registerProcessorParameter("HCalToEMGeVCalibration",
        "Calibration from deposited HCAL to EM energy",
        m_hcalToEMGeVCalibration,
        -1.f);

    registerProcessorParameter("MaxHCalHitHadronicEnergy",
        "The maximum hadronic energy allowed for a single hcal hit",
        m_maxHCalHitHadronicEnergy,
        -1.f);
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
    m_PFA = new TH1F("fPFAtot", "pfo energy", 2000, 0., 250.);
    m_PFAB = new TH1F("fPFABarrel", "pfo energy barrel", 2000, 0., 250.);
    m_PFAVsCosTheta = new TH2F("fPFAVsCosTheta", "pfo energy vs cosTheta", 500, 0., 1., 2000, 0., 250.);
    m_PFAVsCosThetaR = new TH2F("fPFAVsCosThetaR", "pfo energy vs cosThetaReco", 500, 0., 1., 2000, 0., 250.);
    m_PFAVsCosThetaF = new TH2F("fPFAVsCosThetaFull", "pfo energy vs full cosTheta range", 200, -1., 1., 2000, 0., 250.);
    m_PFAE = new TH1F("fPFAECAL", "pfo energy ECAL only events", 2000, 0., 250.);
    m_PFAH = new TH1F("fPFAHCAL", "pfo energy HCAL only events", 2000, 0., 250.);
    m_PFAM = new TH1F("fPFAMUON", "pfo energy MUON only events", 2000, 0., 250.);

    // ECAL total energy
    m_EcalEnergy = new TH1F("fEcalEnergy", "total ECAL energy", 2000, 0., 250.);
    m_HcalEnergy = new TH1F("fHcalEnergy", "total HCAL energy", 2000, 0., 250.);
    m_MuonEnergy = new TH1F("fMuonEnergy", "total MUON energy", 2000, 0., 250.);
    m_LcalEnergy = new TH1F("fLcalEnergy", "total LCAL energy", 2000, 0., 250.);
    m_CalEnergy  = new TH1F("fCalEnergy", "total calorimeter energy", 2000, 0., 250.);
    m_EcalHcalEnergyEM = new TH2F("fEcalHcalEnergyEM", "ECAL vs HCAL energy EM", 2000, 0., 250., 2000, 0., 250.);
    m_EcalHcalEnergyHAD = new TH2F("fEcalHcalEnergyHAD", "ECAL vs HCAL energy HAD", 2000, 0., 250., 2000, 0., 250.);
    m_EcalBarrelHcalEnergyEM = new TH2F("fEcalBarrelHcalEnergyEM", "ECAL barrel vs HCAL energy EM", 2000, 0., 250., 2000, 0., 250.);
    m_EcalEndCapHcalEnergyEM = new TH2F("fEcalEndCapHcalEnergyEM", "ECAL endcap vs HCAL energy EM", 2000, 0., 250., 2000, 0., 250.);
    m_EcalBarrelHcalEnergyHAD = new TH2F("fEcalBarrelHcalEnergyHAD", "ECAL barrel vs HCAL energy HAD", 2000, 0., 250., 2000, 0., 250.);
    m_EcalEndCapHcalEnergyHAD = new TH2F("fEcalEndCapHcalEnergyHAD", "ECAL endcap vs HCAL energy HAD", 2000, 0., 250., 2000, 0., 250.);
    m_CalEnergyE = new TH1F("fCalEnergyECAL", "total calorimeter energy ECAL", 2000, 0., 250.);
    m_CalEnergyH = new TH1F("fCalEnergyHCAL", "total calorimeter energy HCAL", 2000, 0., 250.);
    m_CalEnergyM = new TH1F("fCalEnergyMUON", "total calorimeter energy MUON", 2000, 0., 250.);
    m_CalEnergyVsCosTheta = new TH2F("fCalEnergyVsCosTheta", "total calorimeter energy vs cos theta", 500, 0., 1., 2000, 0., 250.);
    m_CalEnergyVsCosThetaR = new TH2F("fCalEnergyVsCosThetaR", "total calorimeter energy vs cos theta reco", 500, 0., 1.,2000, 0., 250.);
    m_EcalBarrelEnergyByLayer = new TH1F("fEcalBarrelEnergyByLayer", "ECAL barrel energy profile", 100, 0., 100.);
    m_EcalEndCapEnergyByLayer = new TH1F("fEcalEndCapEnergyByLayer", "ECAL endcap energy profile", 100, 0., 100.);

    // MIP Calibration
    m_EcalBarrelMIP = new TH1F("fEcalBarrelMIP", "ECAL barrel MIP peak ", 200, 0., 5.);
    m_EcalEndCapMIP = new TH1F("fEcalEndCapMIP", "ECAL endcap MIP peak ", 200, 0., 5.);
    m_HcalMIP = new TH1F("fHcalMIP", "HCAL MIP peak ", 200, 0., 5.);
    m_MuonMIP = new TH1F("fMuonMIP", "MUON MIP peak ", 200, 0., 5.);
    m_EcalBarrelMIPcorr = new TH1F("fEcalBarrelMIPcorr", "ECAL barrel MIP peak ", 200, 0., 5.);
    m_EcalEndCapMIPcorr = new TH1F("fEcalEndCapMIPcorr", "ECAL endcap MIP peak ", 200, 0., 5.);
    m_HcalMIPcorr = new TH1F("fHcalMIPcorr", "HCAL MIP peak ", 200, 0., 5.);
    m_MuonMIPcorr = new TH1F("fMuonMIPcorr", "MUON MIP peak ", 200, 0., 5.);

    m_PFA->SetDirectory(m_pTFile);
    m_PFAB->SetDirectory(m_pTFile);
    m_PFAVsCosTheta->SetDirectory(m_pTFile);
    m_PFAVsCosThetaR->SetDirectory(m_pTFile);
    m_PFAVsCosThetaF->SetDirectory(m_pTFile);
    m_PFAE->SetDirectory(m_pTFile);
    m_PFAH->SetDirectory(m_pTFile);
    m_PFAM->SetDirectory(m_pTFile);
    m_EcalEnergy->SetDirectory(m_pTFile);
    m_HcalEnergy->SetDirectory(m_pTFile);
    m_MuonEnergy->SetDirectory(m_pTFile);
    m_LcalEnergy->SetDirectory(m_pTFile);
    m_CalEnergy->SetDirectory(m_pTFile);
    m_EcalHcalEnergyEM->SetDirectory(m_pTFile);
    m_EcalHcalEnergyHAD->SetDirectory(m_pTFile);
    m_EcalBarrelHcalEnergyEM->SetDirectory(m_pTFile);
    m_EcalEndCapHcalEnergyEM->SetDirectory(m_pTFile);
    m_EcalBarrelHcalEnergyHAD->SetDirectory(m_pTFile);
    m_EcalEndCapHcalEnergyHAD->SetDirectory(m_pTFile);
    m_CalEnergyE->SetDirectory(m_pTFile);
    m_CalEnergyH->SetDirectory(m_pTFile);
    m_CalEnergyM->SetDirectory(m_pTFile);
    m_CalEnergyVsCosTheta->SetDirectory(m_pTFile);
    m_CalEnergyVsCosThetaR->SetDirectory(m_pTFile);
    m_EcalBarrelEnergyByLayer->SetDirectory(m_pTFile);
    m_EcalEndCapEnergyByLayer->SetDirectory(m_pTFile);
    m_EcalBarrelMIP->SetDirectory(m_pTFile);
    m_EcalEndCapMIP->SetDirectory(m_pTFile);
    m_HcalMIP->SetDirectory(m_pTFile);
    m_MuonMIP->SetDirectory(m_pTFile);
    m_EcalBarrelMIPcorr->SetDirectory(m_pTFile);
    m_EcalEndCapMIPcorr->SetDirectory(m_pTFile);
    m_HcalMIPcorr->SetDirectory(m_pTFile);
    m_MuonMIPcorr->SetDirectory(m_pTFile);
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
    this->ReadHitEnergies(pLCEvent, m_ecalBarrelCollections, ecalBarrelEnergy, m_ecalToMIP, m_EcalBarrelMIP, m_EcalBarrelMIPcorr, "K-1", m_EcalBarrelEnergyByLayer);
    this->ReadHitEnergies(pLCEvent, m_ecalEndCapCollections, ecalEndCapEnergy, m_ecalToMIP, m_EcalEndCapMIP, m_EcalEndCapMIPcorr, "K-1", m_EcalEndCapEnergyByLayer);

    float hcalEnergy(0.f), muonEnergy(0.f);
    this->ReadHitEnergies(pLCEvent, m_hcalCollections, hcalEnergy, m_hcalToMIP, m_HcalMIP, m_HcalMIPcorr);
    this->ReadHitEnergies(pLCEvent, m_muonCollections, muonEnergy, m_muonToMIP, m_MuonMIP, m_MuonMIPcorr);

    float lcalEnergy(0.f), lhcalEnergy(0.f), bcalEnergy(0.f);
    this->ReadHitEnergies(pLCEvent, m_lcalCollections, lcalEnergy);
    this->ReadHitEnergies(pLCEvent, m_lhcalCollections, lhcalEnergy);
    this->ReadHitEnergies(pLCEvent, m_bcalCollections, bcalEnergy);

    // Read reconstructed pfo collection
    float cosThetaReco(0.f), totalPfoEnergy(0.f);
    this->ReadPfoCollections(pLCEvent, m_recoPfoCollections, totalPfoEnergy, cosThetaReco);

    // Fill remaining histograms
    const float totalCalEnergy(ecalBarrelEnergy + ecalEndCapEnergy + hcalEnergy + muonEnergy + lcalEnergy + bcalEnergy + lhcalEnergy);

    if (cosTheta < 0.95f)
    {
        m_EcalEnergy->Fill(ecalBarrelEnergy + ecalEndCapEnergy);

        m_EcalHcalEnergyEM->Fill(ecalBarrelEnergy + ecalEndCapEnergy, hcalEnergy * m_hcalToEMGeVCalibration);
        m_EcalHcalEnergyHAD->Fill(ecalBarrelEnergy * m_ecalToHadGeVCalibrationBarrel + ecalEndCapEnergy * m_ecalToHadGeVCalibrationEndCap, hcalEnergy);

        m_EcalBarrelHcalEnergyEM->Fill(ecalBarrelEnergy, hcalEnergy * m_hcalToEMGeVCalibration);
        m_EcalEndCapHcalEnergyEM->Fill(ecalEndCapEnergy, hcalEnergy * m_hcalToEMGeVCalibration);

        m_EcalBarrelHcalEnergyHAD->Fill((ecalBarrelEnergy * m_ecalToHadGeVCalibrationBarrel), hcalEnergy);
        m_EcalEndCapHcalEnergyHAD->Fill((ecalEndCapEnergy * m_ecalToHadGeVCalibrationEndCap), hcalEnergy);

        m_HcalEnergy->Fill(hcalEnergy);
        m_MuonEnergy->Fill(muonEnergy);
        m_LcalEnergy->Fill(lcalEnergy);
        m_CalEnergy->Fill(totalCalEnergy);

        if ((totalCalEnergy > 0.f) && (((ecalBarrelEnergy + ecalEndCapEnergy) / totalCalEnergy) > 0.95f))
            m_CalEnergyE->Fill(totalCalEnergy);

        if ((totalCalEnergy > 0.f) && ((hcalEnergy / totalCalEnergy) > 0.95f))
            m_CalEnergyH->Fill(totalCalEnergy);

        if ((totalCalEnergy > 0.f) && ((muonEnergy / totalCalEnergy) > 0.95f))
            m_CalEnergyM->Fill(totalCalEnergy);
    }

    m_CalEnergyVsCosTheta->Fill(cosTheta, totalCalEnergy);
    m_CalEnergyVsCosThetaR->Fill(cosThetaReco, totalCalEnergy);

    m_PFAVsCosTheta->Fill(cosTheta, totalPfoEnergy);
    m_PFAVsCosThetaR->Fill(cosThetaReco, totalPfoEnergy);

    if (cosTheta < 0.95f)
    {
        m_PFA->Fill(totalPfoEnergy);

        if (cosTheta < 0.7f)
            m_PFAB->Fill(totalPfoEnergy);

        if ((totalCalEnergy > 0.f) && (((ecalBarrelEnergy + ecalEndCapEnergy) / totalCalEnergy) > 0.95f))
            m_PFAE->Fill(totalPfoEnergy);

        if ((totalCalEnergy > 0.f) && ((hcalEnergy / totalCalEnergy) > 0.95f))
            m_PFAH->Fill(totalPfoEnergy);

        if ((totalCalEnergy > 0.f) && ((muonEnergy / totalCalEnergy) > 0.95f))
            m_PFAM->Fill(totalPfoEnergy);
    }
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
                    const float correction((std::fabs(z) < m_zOfEndCap) ? r / std::sqrt(x * x + y * y) : r / z);

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

    m_PFA->Write();
    m_PFAB->Write();
    m_PFAVsCosTheta->Write();
    m_PFAVsCosThetaR->Write();
    m_PFAVsCosThetaF->Write();
    m_PFAE->Write();
    m_PFAH->Write();
    m_PFAM->Write();
    m_EcalEnergy->Write();
    m_HcalEnergy->Write();
    m_MuonEnergy->Write();
    m_LcalEnergy->Write();
    m_CalEnergy->Write();
    m_EcalHcalEnergyEM->Write();
    m_EcalHcalEnergyHAD->Write();
    m_EcalBarrelHcalEnergyEM->Write();
    m_EcalEndCapHcalEnergyEM->Write();
    m_EcalBarrelHcalEnergyHAD->Write();
    m_EcalEndCapHcalEnergyHAD->Write();
    m_CalEnergyE->Write();
    m_CalEnergyH->Write();
    m_CalEnergyM->Write();
    m_CalEnergyVsCosTheta->Write();
    m_CalEnergyVsCosThetaR->Write();
    m_EcalBarrelEnergyByLayer->Write();
    m_EcalEndCapEnergyByLayer->Write();
    m_EcalBarrelMIP->Write();
    m_EcalEndCapMIP->Write();
    m_HcalMIP->Write();
    m_MuonMIP->Write();
    m_EcalBarrelMIPcorr->Write();
    m_EcalEndCapMIPcorr->Write();
    m_HcalMIPcorr->Write();
    m_MuonMIPcorr->Write();
    m_pTFile->Write();
    m_pTFile->Close();

    delete m_pTFile;
}
