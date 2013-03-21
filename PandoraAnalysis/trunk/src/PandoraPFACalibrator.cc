#include "PandoraPFACalibrator.h"
#include "MCPfoMaker.h"

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

#include <iostream>
#include <cmath>
#include <limits>

using namespace lcio;
using namespace marlin;

PandoraPFACalibrator aPandoraPFACalibrator;

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraPFACalibrator::PandoraPFACalibrator() : Processor("PandoraPFACalibrator") 
{
    _description = "PandoraPFACalibrator for calibration of PandoraPFA";

    registerProcessorParameter("RootFile",
        "Output root file name",
        m_rootFile,
        std::string("PandoraPFACalibrator.root"));

    LCStrVec inputMCParticleCollections;
    inputMCParticleCollections.push_back("MCPFOs");
    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
        "InputMCParticleCollections",
        "Names of input mc particle collections",
        m_inputMCParticleCollections,
        inputMCParticleCollections);

    registerProcessorParameter("InputParticleCollectionName",
        "Particle Collection Name ",
        m_particleCollectionName,
        std::string("PandoraPFANewPFOs"));

    LCStrVec ecalBarrelCollections;
    ecalBarrelCollections.push_back(std::string("ECALBarrel"));
    registerInputCollections(LCIO::CALORIMETERHIT,
        "ECALBarrelcollections", 
        "Name of the ECAL barrel collection used to form clusters",
        m_ecalBarrelCollections,
        ecalBarrelCollections);

    LCStrVec ecalEndCapCollections;
    ecalEndCapCollections.push_back(std::string("ECALEndCap"));
    registerInputCollections(LCIO::CALORIMETERHIT,
        "ECALEndCapcollections", 
        "Name of the ECAL EndCap collection used to form clusters",
        m_ecalEndCapCollections,
        ecalEndCapCollections);

    LCStrVec hcalCollections;
    hcalCollections.push_back(std::string("HCALBarrel"));
    hcalCollections.push_back(std::string("HCALEndcap"));
    hcalCollections.push_back(std::string("HCALOther"));
    registerInputCollections(LCIO::CALORIMETERHIT,
        "HCALcollections", 
        "Name of the HCAL collection used to form clusters",
        m_hcalCollections,
        hcalCollections);

    LCStrVec muonCollections;
    muonCollections.push_back(std::string("MUON"));
    registerInputCollections(LCIO::CALORIMETERHIT,
        "MUONcollections", 
        "Name of the MUON collection used to form clusters",
        m_muonCollections,
        muonCollections);

    LCStrVec bcalCollections;
    bcalCollections.push_back(std::string("BCAL"));
    registerInputCollections(LCIO::CALORIMETERHIT,
        "BCALcollections", 
        "Name of the BCAL collection used to form clusters",
        m_bcalCollections,
        bcalCollections);

    LCStrVec lhcalCollections;
    lhcalCollections.push_back(std::string("LHCAL"));
    registerInputCollections(LCIO::CALORIMETERHIT,
        "LHCALcollections", 
        "Name of the LHCAL collection used to form clusters",
        m_lhcalCollections,
        lhcalCollections);

    LCStrVec lcalCollections;
    lcalCollections.push_back(std::string("LCAL"));
    registerInputCollections(LCIO::CALORIMETERHIT,
        "LCALcollections", 
        "Name of the LCAL collection used to form clusters",
        m_lcalCollections,
        lcalCollections);

    //********************************************************************
    //**************************** Calibration ***************************
    //********************************************************************
    // Copy from PandoraPFANewProcessor
    // Calibration constants
    // NOTE THE CALIBRATION SCHEME USED HERE IS NON-STANDARD
    // First hits in ECAL/HCAL are converted to normal incident MIP equivalents
    // Thresholds are applied at the MIP level
    // MIP equivalent signals are converted to EM or Hadronic energy 
    //   i.e. different calibrations for EM and hadronic showers
    registerProcessorParameter("ECalToMipCalibration",
        "Calibration from deposited ECAL energy to MIP",
        m_ecalToMIP,
        (float)160.0);

    registerProcessorParameter("HCalToMipCalibration",
        "Calibration from deposited HCAL energy to MIP",
        m_hcalToMIP,
        (float)34.8);

    registerProcessorParameter("MuonToMipCalibration",
        "Calibration from deposited MUON energy to MIP",
        m_muonToMIP,
        (float)1.0);

    registerProcessorParameter("ECalToEMGeVCalibration",
        "Calibration from deposited ECAL to EM energy",
        m_ecalToEMGeVCalibration,
        (float)1.0);

    registerProcessorParameter("HCalToHadGeVCalibration",
        "Calibration from deposited HCAL to Hadronic energy",
        m_hcalToHadGeVCalibration,
        (float)1.0);

    registerProcessorParameter("ECalToHadGeVCalibrationBarrel",
        "Calibration from deposited ECAL barrel to Hadronic energy",
        m_ecalToHadGeVCalibrationBarrel,
        (float)1.03);

    registerProcessorParameter("ECalToHadGeVCalibrationEndCap",
        "Calibration from deposited ECAL endcap to Hadronic energy",
        m_ecalToHadGeVCalibrationEndCap,
        (float)1.16);

    registerProcessorParameter("HCalToEMGeVCalibration",
        "Calibration from deposited HCAL to EM energy",
        m_hcalToEMGeVCalibration,
        (float)1.0);

    registerProcessorParameter("MaxHCalHitHadronicEnergy",
        "The maximum hadronic energy allowed for a single hcal hit",
        m_maxHCalHitHadronicEnergy,
        (float)1.0);

    CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::init() 
{
    // usually a good idea to
    printParameters();

    m_nRun = 0;
    m_nEvt = 0;

    m_pTFile = new TFile(m_rootFile.c_str(), "recreate");

    // PFA total energy
    m_PFA = new TH1F("fPFAtot", "total energy", 1000, 0., 250.);
    m_PFAB = new TH1F("fPFAB", "total energy barrel", 1000, 0., 250.);
    m_PFAVsCosTheta = new TH2F("fPFAVsCosTheta", "total energy vs CosTheta", 500, 0., 1., 1000, 0., 250.);
    m_PFAVsCosThetaR = new TH2F("fPFAVsCosThetaR", "total energy vs CosThetaReco", 500, 0., 1., 1000, 0., 250.);
    m_PFAVsZCoG = new TH2F("fPFAVsZCoG", "total energy vs zCog", 600, 0., 3000., 1000, 0., 200.);
    m_PFAVsCosThetaX = new TH2F("fPFAVsCosThetaX", "total energy vs CosTheta", 200, -1., 1., 1000, 0., 250.);
    m_PFAE = new TH1F("fPFAE", "total energy ECAL only events", 1000, 0., 250.);
    m_PFAH = new TH1F("fPFAH", "total energy HCAL only events", 1000, 0., 250.);
    m_PFAM = new TH1F("fPFAM", "total energy MUON only events", 1000, 0., 250.);
    m_XvsY = new TH2F("fXvsY", "x vs y", 1000, -2500., 2500., 1000, -2500., 2500.);

    // ECAL total energy
    m_EcalEnergy = new TH1F("fEcalEnergy", "total ecal energy", 1000, 0., 250.);
    m_HcalEnergy = new TH1F("fHcalEnergy", "total hcal energy", 1000, 0., 250.);
    m_MuonEnergy = new TH1F("fMuonEnergy", "total muon energy", 1000, 0., 250.);
    m_LcalEnergy = new TH1F("fLcalEnergy", "total lcal energy", 1000, 0., 250.);
    m_CalEnergy  = new TH1F("fCalEnergy", "total cal energy", 1000, 0., 250.);
    m_EcalHcalEnergyEM = new TH2F("fEcalHcalEnergyEM", "ecal vs hcal energy EM", 1000, 0., 250., 1000, 0., 250.);
    m_EcalHcalEnergyHAD = new TH2F("fEcalHcalEnergyHAD", "ecal vs hcal energy HAD", 1000, 0., 250., 1000, 0., 250.);
    m_EcalBarrelHcalEnergyEM = new TH2F("fEcalBarrelHcalEnergyEM", "ecal barrel vs hcal energy EM", 1000, 0., 250., 1000, 0., 250.);
    m_EcalEndCapHcalEnergyEM = new TH2F("fEcalEndCapHcalEnergyEM", "ecal endcap vs hcal energy EM", 1000, 0., 250., 1000, 0., 250.);
    m_EcalBarrelHcalEnergyHAD = new TH2F("fEcalBarrelHcalEnergyHAD", "ecal barrel vs hcal energy HAD", 1000, 0., 250., 1000, 0., 250.);
    m_EcalEndCapHcalEnergyHAD = new TH2F("fEcalEndCapHcalEnergyHAD", "ecal endcap vs hcal energy HAD", 1000, 0., 250., 1000, 0., 250.);
    m_CalEnergyE = new TH1F("fCalEnergyE", "total cal energy E", 1000, 0., 250.);
    m_CalEnergyH = new TH1F("fCalEnergyH", "total cal energy H", 1000, 0., 250.);
    m_CalEnergyM = new TH1F("fCalEnergyM", "total cal energy M", 1000, 0., 250.);
    m_CalEnergyVsCosTheta = new TH2F("fCalEnergyVsCosTheta", "total cal energy vs cos theta", 500, 0., 1., 1000, 0., 250.);
    m_CalEnergyVsCosThetaR = new TH2F("fCalEnergyVsCosThetaR", "total cal energy vs cos theta reco", 500, 0., 1.,1000, 0., 250.);
    m_EcalBarrelEnergyByLayer = new TH1F("fEcalBarrelEnergyByLayer", "ecal barrel energy profile", 100, 0., 100.);
    m_EcalEndCapEnergyByLayer = new TH1F("fEcalEndCapEnergyByLayer", "ecal endcap energy profile", 100, 0., 100.);

    // MIP Calibration
    m_EcalBarrelMIP = new TH1F("fEcalBarrelMIP", "ecal barrel MIP peak ", 100, 0., 5.);
    m_EcalEndCapMIP = new TH1F("fEcalEndCapMIP", "ecal endcap MIP peak ", 100, 0., 5.);
    m_HcalMIP = new TH1F("fHcalMIP", "hcal MIP peak ", 100, 0., 5.);
    m_MuonMIP = new TH1F("fMuonMIP", "muon MIP peak ", 100, 0., 5.);
    m_EcalBarrelMIPcorr = new TH1F("fEcalBarrelMIPcorr", "ecal barrel MIP peak ", 100, 0., 5.);
    m_EcalEndCapMIPcorr = new TH1F("fEcalEndCapMIPcorr", "ecal endcap MIP peak ", 100, 0., 5.);
    m_HcalMIPcorr = new TH1F("fHcalMIPcorr", "hcal MIP peak ", 100, 0., 5.);
    m_MuonMIPcorr = new TH1F("fMuonMIPcorr", "muon MIP peak ", 100, 0., 5.);
    m_CosT = new TH1F("fCosT", "cosTheta ", 100, 0., 1.);
    m_PhotonCosT = new TH1F("fPhotonCosT", "cosTheta Photon", 100, 0., 1.0);

    m_PFA->SetDirectory(m_pTFile);
    m_PFAB->SetDirectory(m_pTFile);
    m_PFAVsCosTheta->SetDirectory(m_pTFile);
    m_PFAVsCosThetaR->SetDirectory(m_pTFile);
    m_PFAVsZCoG->SetDirectory(m_pTFile);
    m_PFAVsCosThetaX->SetDirectory(m_pTFile);
    m_PFAE->SetDirectory(m_pTFile);
    m_PFAH->SetDirectory(m_pTFile);
    m_PFAM->SetDirectory(m_pTFile);
    m_XvsY->SetDirectory(m_pTFile);
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
    m_CosT->SetDirectory(m_pTFile);
    m_PhotonCosT->SetDirectory(m_pTFile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::processRunHeader(LCRunHeader *run)
{
    m_nRun++;
    m_detectorName = run->getDetectorName();

    std::cout << " DETECTOR : " << m_detectorName << std::endl;

    const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
    m_zOfEndCap = static_cast<float>(pEcalEndcap.getExtent()[2]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::processEvent(LCEvent *evt)
{
    m_nEvt ++;
    float ecalBarrelEnergy = 0, ecalEndCapEnergy = 0;
    float hcalEnergy = 0.;
    float muonEnergy = 0.;
    float lcalEnergy = 0.;
    float bcalEnergy = 0.;
    float lhcalEnergy = 0.;
    float totEnergy = 0.;
    float pfoEnergy = 0.;

    // read MCParticle collection
    for (unsigned int i = 0; i < m_inputMCParticleCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection(m_inputMCParticleCollections.at(i).c_str());

            if (col != 0)
            {
                const int nelem(col->getNumberOfElements());
                m_out(DEBUG) << std::endl << "MCParticle collection " << m_inputMCParticleCollections.at(i) << " with "<< nelem << " elements" << std::endl;

                for (int iMC = 0; iMC < nelem; ++iMC)
                {
                    ReconstructedParticle *part =  dynamic_cast<ReconstructedParticle*>(col->getElementAt(0));
                    const float x(part->getMomentum()[0]);
                    const float y(part->getMomentum()[1]);
                    const float z(part->getMomentum()[2]);
                    const float p(std::sqrt(x * x + y * y + z * z));

                    m_cosTheta = 0;
                    m_x = 0.;
                    m_y = 0.;

                    if (p > 0)
                    {
                        m_cosThetaX = z / std::sqrt(x * x + y * y + z * z);
                        m_cosTheta = std::fabs(z) / std::sqrt(x * x + y * y + z * z);

                        if (std::fabs(z) > 0)
                        {
                            m_x = m_zOfEndCap * x / z;
                            m_y = m_zOfEndCap * y / z;
                        }
                    }

                    const float energy(part->getEnergy());
                    m_out(DEBUG) << " First MC Particle : " << energy << " GeV " << " at " << x << "," << y << "," << z << " (|cosTheta| = " << m_cosTheta << ")" << std::endl;
                }
            }
        }
        catch (DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " << m_inputMCParticleCollections.at(i) << std::endl;
        }
    }

    // read ECAL BARREL hits
    for (unsigned int i = 0; i < m_ecalBarrelCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection(m_ecalBarrelCollections.at(i).c_str());

            if (col != 0)
            {
                CellIDDecoder<CalorimeterHit> id(col);
                const int nelem(col->getNumberOfElements());

                for (int ihit = 0; ihit < nelem; ++ihit)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    float hitEnergy = hit->getEnergy();
                    ecalBarrelEnergy += hitEnergy;
                    int layerNumber = id(hit)["K-1"] + 1;

                    m_EcalBarrelEnergyByLayer->Fill(layerNumber,hitEnergy);

                    float energyInMips=0;
                    const float x(hit->getPosition()[0]);
                    const float y(hit->getPosition()[1]);
                    const float z(hit->getPosition()[2]);
                    const float r(std::sqrt(x*x+y*y+z*z));
                    const float correction((std::fabs(z) < m_zOfEndCap) ? r / std::sqrt(x*x+y*y) : r / z);

                    energyInMips = hit->getEnergy() * m_ecalToMIP;
                    m_EcalBarrelMIP->Fill(energyInMips);
                    m_EcalBarrelMIPcorr->Fill(energyInMips / correction);
                }
                m_out(DEBUG)<< " ECAL BARREL hits : " << nelem << " energy = " << ecalBarrelEnergy << std::endl;
            }
        }
        catch (DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " <<  m_ecalBarrelCollections[i] << std::endl;
        }
    }

    // read ECAL ENDCAP hits
    for (unsigned int i = 0; i < m_ecalEndCapCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection(m_ecalEndCapCollections.at(i).c_str());

            if (col != 0)
            {
                CellIDDecoder<CalorimeterHit> id(col);
                const int nelem(col->getNumberOfElements());

                for (int ihit = 0; ihit < nelem; ++ihit)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    float hitEnergy = hit->getEnergy();
                    ecalEndCapEnergy += hitEnergy;
                    int layerNumber = id(hit)["K-1"] + 1;

                    m_EcalEndCapEnergyByLayer->Fill(layerNumber,hitEnergy);

                    float energyInMips=0;
                    const float x(hit->getPosition()[0]);
                    const float y(hit->getPosition()[1]);
                    const float z(hit->getPosition()[2]);
                    const float r(std::sqrt(x * x + y * y + z * z));
                    const float correction((std::fabs(z) < m_zOfEndCap) ? r / std::sqrt(x*x+y*y) : r / z);

                    energyInMips = hit->getEnergy() * m_ecalToMIP;
                    m_EcalEndCapMIP->Fill(energyInMips);
                    m_EcalEndCapMIPcorr->Fill(energyInMips / correction);
                }
                m_out(DEBUG) << " ECAL ENDCAP hits : " << nelem << " energy = " << ecalEndCapEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " <<  m_ecalEndCapCollections[i] << std::endl;
        }
    }

    // read HCAL hits
    for (unsigned int i = 0; i < m_hcalCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection(m_hcalCollections.at(i).c_str());

            if (col != 0)
            {
                const int nelem(col->getNumberOfElements());

                for (int ihit = 0; ihit < nelem; ++ihit)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));

                    const float hitEnergy(std::min(hit->getEnergy(), m_maxHCalHitHadronicEnergy));
                    hcalEnergy += hitEnergy;
                    float energyInMips=0;
                    const float x(hit->getPosition()[0]);
                    const float y(hit->getPosition()[1]);
                    const float z(hit->getPosition()[2]);
                    const float r(std::sqrt(x * x + y * y + z * z));
                    const float correction((std::fabs(z) < m_zOfEndCap) ? r / std::sqrt(x*x+y*y) : r / z);

                    energyInMips = hitEnergy * m_hcalToMIP;
                    m_HcalMIP->Fill(energyInMips);
                    m_HcalMIPcorr->Fill(energyInMips / correction);
                }
                m_out(DEBUG) << " HCAL hits : " << nelem << " energy = " << hcalEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " <<  m_hcalCollections[i] << std::endl;
        }
    }

    // read MUON hits
    for (unsigned int i = 0; i < m_muonCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection(m_muonCollections.at(i).c_str());

            if (col != 0) 
            {
                const int nelem(col->getNumberOfElements());

                for (int ihit = 0; ihit < nelem; ++ihit)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    muonEnergy += hit->getEnergy();
                    float energyInMips=0;
                    const float x(hit->getPosition()[0]);
                    const float y(hit->getPosition()[1]);
                    const float z(hit->getPosition()[2]);
                    const float r(std::sqrt(x * x + y * y + z * z));
                    const float correction((std::fabs(z) < m_zOfEndCap) ? r / std::sqrt(x*x+y*y) : r / z);

                    energyInMips = hit->getEnergy() * m_muonToMIP;
                    m_MuonMIP->Fill(energyInMips);
                    m_MuonMIPcorr->Fill(energyInMips / correction);
                }
                m_out(DEBUG) << " MUON hits : " << nelem << " energy = " << muonEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " <<  m_muonCollections[i] << std::endl;
        }
    }

    // read LCAL hits
    for (unsigned int i = 0; i < m_lcalCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection(m_lcalCollections.at(i).c_str());

            if (col != 0)
            {
                const int nelem(col->getNumberOfElements());

                for (int ihit = 0; ihit < nelem; ++ihit)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    lcalEnergy += hit->getEnergy();
                }
                m_out(DEBUG) << " LCAL hits : " << nelem << " energy = " << lcalEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " <<  m_lcalCollections[i] << std::endl;
        }
    }

    // read BCAL hits
    for (unsigned int i = 0; i < m_bcalCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection(m_bcalCollections.at(i).c_str());

            if (col != 0)
            {
                const int nelem(col->getNumberOfElements());

                for (int ihit = 0; ihit < nelem; ++ihit)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    bcalEnergy += hit->getEnergy();
                }
                m_out(DEBUG) << " BCAL hits : " << nelem << " energy = " << bcalEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " <<  m_bcalCollections[i] << std::endl;
        }
    }

    // read LHCAL hits
    for (unsigned int i = 0 ; i < m_lhcalCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection(m_lhcalCollections.at(i).c_str());

            if (col != 0)
            {
                const int nelem(col->getNumberOfElements());

                for (int ihit = 0; ihit < nelem; ++ihit)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    lhcalEnergy += hit->getEnergy();
                }
                m_out(DEBUG) << " LHCAL hits : " << nelem << " energy = " << lhcalEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &)
        {
            m_out(DEBUG) << "No Collection : " <<  m_lhcalCollections[i] << std::endl;
        }
    }

    // Total energy
    totEnergy = ecalBarrelEnergy + ecalEndCapEnergy + hcalEnergy + muonEnergy + lcalEnergy + bcalEnergy + lhcalEnergy;
    m_out(DEBUG) << " Total Calorimetric Energy : " << totEnergy << std::endl;

    if (m_cosTheta < 0.95)
    {
        m_EcalEnergy->Fill(ecalBarrelEnergy + ecalEndCapEnergy,1.);

        m_EcalHcalEnergyEM->Fill(ecalBarrelEnergy + ecalEndCapEnergy, hcalEnergy * m_hcalToEMGeVCalibration, 1.);
        m_EcalHcalEnergyHAD->Fill(ecalBarrelEnergy * m_ecalToHadGeVCalibrationBarrel + ecalEndCapEnergy * m_ecalToHadGeVCalibrationEndCap, hcalEnergy, 1.);

        m_EcalBarrelHcalEnergyEM->Fill(ecalBarrelEnergy, hcalEnergy * m_hcalToEMGeVCalibration, 1.);
        m_EcalEndCapHcalEnergyEM->Fill(ecalEndCapEnergy, hcalEnergy * m_hcalToEMGeVCalibration, 1.);

        m_EcalBarrelHcalEnergyHAD->Fill((ecalBarrelEnergy * m_ecalToHadGeVCalibrationBarrel), hcalEnergy, 1.);
        m_EcalEndCapHcalEnergyHAD->Fill((ecalEndCapEnergy * m_ecalToHadGeVCalibrationEndCap), hcalEnergy, 1.);

        m_HcalEnergy->Fill(hcalEnergy,1.);
        m_MuonEnergy->Fill(muonEnergy,1.);
        m_LcalEnergy->Fill(lcalEnergy,1.);
        m_CalEnergy->Fill(totEnergy,1.);

        if (totEnergy > 0 && (ecalBarrelEnergy + ecalEndCapEnergy) / totEnergy > 0.95)
            m_CalEnergyE->Fill(totEnergy);

        if (totEnergy > 0 && hcalEnergy / totEnergy > 0.95)
            m_CalEnergyH->Fill(totEnergy);

        if (totEnergy > 0 && muonEnergy / totEnergy > 0.95)
            m_CalEnergyM->Fill(totEnergy);
    }

    m_CalEnergyVsCosTheta->Fill(m_cosTheta, totEnergy, 1.);

    // Read reconstructed PFOs
    try
    {
        LCCollection * col = evt->getCollection(m_particleCollectionName.c_str());

        if (col != 0) 
        {
            const int nelem(col->getNumberOfElements());
            bool check=true;
            float photonE = 0;
            float p[3] = {0.,0.,0.};
            m_zCoG=0.;
            float esum = 0;

            for (int j = 0; j < nelem; ++j)
            {
                ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(j));
                pfoEnergy+=recoPart->getEnergy();

                TrackVec tracks = recoPart->getTracks();

                if (tracks.size() != 0)
                    check = false;

                if (recoPart->getType() == 22)
                    photonE += recoPart->getEnergy();

                p[0] += recoPart->getMomentum()[0];
                p[1] += recoPart->getMomentum()[1];
                p[2] += recoPart->getMomentum()[2];

                const ClusterVec &clusters = recoPart->getClusters();

                for (ClusterVec::const_iterator iter = clusters.begin(), iterEnd = clusters.end(); iter != iterEnd; ++iter)
                {
                    m_zCoG += (*iter)->getPosition()[2] * (*iter)->getEnergy();
                    esum += (*iter)->getEnergy();
                }
            }

            if (esum > 0)
            {
                m_zCoG = m_zCoG / esum;
            }
            else
            {
                m_zCoG = std::numeric_limits<float>::max();
            }

            m_cosThetaR = std::fabs(p[2]) / std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
            m_out(DEBUG) << " PFO objects : " << nelem << " energy = " << pfoEnergy << std::endl;

            if (check)
            {
                m_CosT->Fill(m_cosTheta,1.);

                if (pfoEnergy > 0 && photonE > 0.5*pfoEnergy)
                    m_PhotonCosT->Fill(m_cosTheta,1.);
            }

            m_PFAVsCosTheta->Fill(m_cosTheta,pfoEnergy,1.);
            m_PFAVsCosThetaR->Fill(m_cosThetaR,pfoEnergy,1.);
            m_PFAVsZCoG->Fill(fabs(m_zCoG),pfoEnergy,1.);
            m_PFAVsCosThetaX->Fill(m_cosThetaX,pfoEnergy,1.);

            if (m_cosTheta < 0.95)
            {
                m_PFA->Fill(pfoEnergy,1.);

                if (m_cosTheta < 0.7)
                    m_PFAB->Fill(pfoEnergy,1.);

                if (totEnergy > 0 && (ecalBarrelEnergy + ecalEndCapEnergy) / totEnergy > 0.95)
                    m_PFAE->Fill(pfoEnergy,1.);

                if (totEnergy > 0 && hcalEnergy / totEnergy > 0.95)
                    m_PFAH->Fill(pfoEnergy,1.);

                if (totEnergy > 0 && muonEnergy / totEnergy > 0.95)
                    m_PFAM->Fill(pfoEnergy,1.);
            }

            m_XvsY->Fill(m_x, m_y, 1.);
        }
    }
    catch(DataNotAvailableException &)
    {
        m_out(DEBUG) << "No Collection : " <<  m_particleCollectionName << std::endl;
    }

    m_CalEnergyVsCosThetaR->Fill(m_cosThetaR, totEnergy,1.);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::check(LCEvent *pLCEvent)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::end()
{
    std::cout << "PandoraPFACalibrator::end()  " << name() << " processed " << m_nEvt << " events in " << m_nRun << " runs " << std::endl;

    m_PFA->Write();
    m_PFAB->Write();
    m_PFAVsCosTheta->Write();
    m_PFAVsCosThetaR->Write();
    m_PFAVsZCoG->Write();
    m_PFAVsCosThetaX->Write();
    m_PFAE->Write();
    m_PFAH->Write();
    m_PFAM->Write();
    m_XvsY->Write();
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
    m_CosT->Write();
    m_PhotonCosT->Write();
    m_pTFile->Write();
    m_pTFile->Close();

    delete m_pTFile;
}
