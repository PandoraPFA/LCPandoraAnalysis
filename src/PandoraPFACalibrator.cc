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

using namespace lcio ;
using namespace marlin ;

PandoraPFACalibrator aPandoraPFACalibrator;

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraPFACalibrator::PandoraPFACalibrator() : Processor("PandoraPFACalibrator") 
{
    // modify processor description
    _description = "PandoraPFACalibrator for calibration of PandoraPFA" ;

    // register steering parameters: name, description, class-variable, default value
    registerProcessorParameter( "RootFile" , 
                  "Name of the Track collection used for clustering"  ,
                  _rootFile,
                  std::string("PandoraPFACalibrator.root") ) ;

    LCStrVec inputMCParticleCollections;
    inputMCParticleCollections.push_back("MCPFOs");
    registerInputCollections(LCIO::RECONSTRUCTEDPARTICLE,
               "InputMCParticleCollections",
               "Names of input mc particle collections",
               _inputMCParticleCollections,
                            inputMCParticleCollections);

    registerProcessorParameter( "InputParticleCollectionName" , 
                  "Particle Collection Name "  ,
                  _particleCollectionName,
                  std::string("PandoraPFANewPFOs"));

    LCStrVec ecalBarrelCollections;
    ecalBarrelCollections.push_back(std::string("ECALBarrel"));
    registerInputCollections( LCIO::CALORIMETERHIT,
               "ECALBarrelcollections" , 
               "Name of the ECAL barrel collection used to form clusters"  ,
               _ecalBarrelCollections,
               ecalBarrelCollections ) ;

    LCStrVec ecalEndCapCollections;
    ecalEndCapCollections.push_back(std::string("ECALEndCap"));
    registerInputCollections( LCIO::CALORIMETERHIT,
               "ECALEndCapcollections" , 
               "Name of the ECAL EndCap collection used to form clusters"  ,
               _ecalEndCapCollections,
               ecalEndCapCollections ) ;

    LCStrVec hcalCollections;
    hcalCollections.push_back(std::string("HCALBarrel"));
    hcalCollections.push_back(std::string("HCALEndcap"));
    hcalCollections.push_back(std::string("HCALOther"));
    registerInputCollections( LCIO::CALORIMETERHIT,
               "HCALcollections" , 
               "Name of the HCAL collection used to form clusters"  ,
               _hcalCollections,
               hcalCollections ) ;

    LCStrVec muonCollections;
    muonCollections.push_back(std::string("MUON"));
    registerInputCollections( LCIO::CALORIMETERHIT,
               "MUONcollections" , 
               "Name of the MUON collection used to form clusters"  ,
               _muonCollections,
               muonCollections ) ;

    LCStrVec bcalCollections;
    bcalCollections.push_back(std::string("BCAL"));
    registerInputCollections( LCIO::CALORIMETERHIT,
               "BCALcollections" , 
               "Name of the BCAL collection used to form clusters"  ,
               _bcalCollections,
               bcalCollections ) ;


    LCStrVec lhcalCollections;
    lhcalCollections.push_back(std::string("LHCAL"));
    registerInputCollections( LCIO::CALORIMETERHIT,
               "LHCALcollections" , 
               "Name of the LHCAL collection used to form clusters"  ,
               _lhcalCollections,
               lhcalCollections ) ;


    LCStrVec lcalCollections;
    lcalCollections.push_back(std::string("LCAL"));
    registerInputCollections( LCIO::CALORIMETERHIT,
               "LCALcollections" , 
               "Name of the LCAL collection used to form clusters"  ,
               _lcalCollections,
               lcalCollections ) ;


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

    registerProcessorParameter( "ECalToMipCalibration" , 
                  "Calibration from deposited ECAL energy to MIP"  ,
                  _ecalToMIP,
                   (float)160.0 ) ;

    registerProcessorParameter( "HCalToMipCalibration" , 
                  "Calibration from deposited HCAL energy to MIP"  ,
                  _hcalToMIP,
                   (float)34.8 ) ;

    registerProcessorParameter( "MuonToMipCalibration" , 
                  "Calibration from deposited MUON energy to MIP"  ,
                  _muonToMIP,
                   (float)1.0 ) ;

    registerProcessorParameter( "ECalMipThreshold" , 
                  "ECAL Threshold in MIPS"  ,
                  _ecalMIPThreshold,
                   (float)0.5 ) ;

    registerProcessorParameter( "HCalMipThreshold" , 
                  "HCAL Threshold in MIPS"  ,
                  _hcalMIPThreshold,
                   (float)0.5 ) ;

    registerProcessorParameter( "HCalToEMGeVCalibration" , 
                  "Calibration from deposited HCAL MIP to EM energy"  ,
                  _hcalEMMIPToGeV,
                   (float)1.0 ) ;

    registerProcessorParameter( "HCalToHadGeVCalibration" , 
                  "Calibration from deposited HCAL MIP to Hadronic energy"  ,
                  _hcalHadMIPToGeV,
                   (float)1.0 ) ;

    registerProcessorParameter( "ECalToEMGeVCalibration" , 
                  "Calibration from deposited ECAL MIP to EM energy"  ,
                  _ecalEMMIPToGeV,
                   (float)1.0 ) ;

    registerProcessorParameter( "ECalToHadGeVCalibrationBarrel" , 
                  "Calibration from deposited ECAL barrel MIP to Hadronic energy"  ,
                  _ecalBarrelHadMIPToGeV,
                   (float)1.03 ) ;

    registerProcessorParameter( "ECalToHadGeVCalibrationEndCap" , 
                  "Calibration from deposited ECAL endcap MIP to Hadronic energy"  ,
                  _ecalEndCapHadMIPToGeV,
                   (float)1.16 ) ;

    CellIDDecoder<CalorimeterHit>::setDefaultEncoding("M:3,S-1:3,I:9,J:9,K-1:6");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::init() 
{
    // usually a good idea to
    printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

    // PFA total energy
    fPFA            = new TH1F("fPFAtot", "total energy",1000, 0., 250.0);
    fPFAB           = new TH1F("fPFAB", "total energy barrel",1000, 0., 250.0);
    fPFAVsCosTheta  = new TH2F("fPFAVsCosTheta", "total energy vs CosTheta",500,0.,1.,1000, 0., 250.0);
    fPFAVsCosThetaR = new TH2F("fPFAVsCosThetaR", "total energy vs CosThetaReco",500,0.,1.,1000, 0., 250.0);
    fPFAVsZCoG      = new TH2F("fPFAVsZCoG", "total energy vs zCog",600,0.,3000.,1000, 0., 200.0);
    fPFAVsCosThetaX = new TH2F("fPFAVsCosThetaX", "total energy vs CosTheta",200,-1.,1.,1000, 0., 250.0);

    fPFAE    = new TH1F("fPFAE", "total energy ECAL only events",1000, 0., 250.0);
    fPFAH    = new TH1F("fPFAH", "total energy HCAL only events",1000, 0., 250.0);
    fPFAM    = new TH1F("fPFAM", "total energy MUON only events",1000, 0., 250.0);

    fXvsY    = new TH2F("fXvsY", "x vs y",1000,-2500.,2500.,1000, -2500., 2500.0);

    // ECAL total energy
    fEcalEnergy           = new TH1F("fEcalEnergy", "total ecal energy",1000, 0., 250.0);
    fHcalEnergy           = new TH1F("fHcalEnergy", "total hcal energy",1000, 0., 250.0);
    fMuonEnergy           = new TH1F("fMuonEnergy", "total muon energy",1000, 0., 250.0);
    fLcalEnergy           = new TH1F("fLcalEnergy", "total lcal energy",1000, 0., 250.0);
    fCalEnergy            = new TH1F("fCalEnergy" , "total cal energy",1000, 0., 250.0);

    fEcalBarrelHcalEnergy = new TH2F("fEcalBarrelHcalEnergy", "ecal barrel vs hcal energy",1000, 0., 250.0,1000,0.,250.);
    fEcalEndCapHcalEnergy = new TH2F("fEcalEndCapHcalEnergy", "ecal endcap vs hcal energy",1000, 0., 250.0,1000,0.,250.);

    fCalEnergyE           = new TH1F("fCalEnergyE" , "total cal energy E",1000, 0., 250.0);
    fCalEnergyH           = new TH1F("fCalEnergyH" , "total cal energy H",1000, 0., 250.0);
    fCalEnergyM           = new TH1F("fCalEnergyM" , "total cal energy M",1000, 0., 250.0);

    fCalEnergyVsCosTheta  = new TH2F("fCalEnergyVsCosTheta" , "total cal energy vs cos theta",500,0.,1.,1000, 0., 250.0);
    fCalEnergyVsCosThetaR = new TH2F("fCalEnergyVsCosThetaR" , "total cal energy vs cos theta reco",500,0.,1.,1000, 0., 250.0);

    fEcalBarrelEnergyByLayer = new TH1F("fEcalBarrelEnergyByLayer", "ecal barrel energy profile",100, 0., 100.0);
    fEcalEndCapEnergyByLayer = new TH1F("fEcalEndCapEnergyByLayer", "ecal endcap energy profile",100, 0., 100.0);

    // MIP Calibration
    fEcalBarrelMIP = new TH1F("fEcalBarrelMIP", "ecal barrel MIP peak ",100, 0., 5.0);
    fEcalEndCapMIP = new TH1F("fEcalEndCapMIP", "ecal endcap MIP peak ",100, 0., 5.0);
    fHcalMIP       = new TH1F("fHcalMIP", "hcal MIP peak ",100, 0., 5.0);
    fMuonMIP       = new TH1F("fMuonMIP", "muon MIP peak ",1000, 0., 50.0);

    fEcalBarrelMIPcorr  = new TH1F("fEcalBarrelMIPcorr", "ecal barrel MIP peak ",100, 0., 5.0);
    fEcalEndCapMIPcorr  = new TH1F("fEcalEndCapMIPcorr", "ecal endcap MIP peak ",100, 0., 5.0);
    fHcalMIPcorr  = new TH1F("fHcalMIPcorr", "hcal MIP peak ",100, 0., 5.0);
    fMuonMIPcorr  = new TH1F("fMuonMIPcorr", "muon MIP peak ",1000, 0., 50.0);

    fCosT       = new TH1F("fCosT", "cosTheta ",100, 0., 1.0);
    fPhotonCosT = new TH1F("fPhotonCosT", "cosTheta Photon",100, 0., 1.0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::processRunHeader(LCRunHeader *run)
{
    _nRun++ ;
    _detectorName = run->getDetectorName();

    std::cout << " DETECTOR : " << _detectorName << std::endl;

    const gear::CalorimeterParameters& pEcalEndcap = Global::GEAR->getEcalEndcapParameters();
    _zOfEndCap = static_cast<float>(pEcalEndcap.getExtent()[2]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::processEvent( LCEvent * evt ) 
{
    _nEvt ++ ;
    float ecalBarrelEnergy = 0, ecalEndCapEnergy = 0;
    float hcalEnergy  = 0.;
    float muonEnergy  = 0.;
    float lcalEnergy  = 0.;
    float bcalEnergy  = 0.;
    float lhcalEnergy = 0.;
    float totEnergy   = 0.;
    float pfoEnergy   = 0.;

    // read MCParticle collection
    for (unsigned int i = 0; i < _inputMCParticleCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection( _inputMCParticleCollections.at(i).c_str() );
            int nelem = col->getNumberOfElements();

            m_out(DEBUG)<<"\n MCParticle collection "<<_inputMCParticleCollections.at(i)<<" with "<<nelem<<" elements"<<std::endl;

            for(int iMC = 0; iMC < nelem; iMC++)
            {
                ReconstructedParticle *part =  dynamic_cast<ReconstructedParticle*>( col->getElementAt(0) );
                float x =  part->getMomentum()[0];
                float y =  part->getMomentum()[1];
                float z =  part->getMomentum()[2];
                float p =  sqrt(x*x+y*y+z*z);

                _cosTheta = 0;
                _x = 0.;
                _y = 0.;

                if (p > 0)
                {
                    _cosThetaX = z / sqrt(x*x+y*y+z*z);
                    _cosTheta = std::fabs(z) / std::sqrt(x*x+y*y+z*z);

                    if (std::fabs(z) > 0)
                    {
                        _x = _zOfEndCap * x / z;
                        _y = _zOfEndCap * y / z;
                    }
                }

                float energy = part->getEnergy();
                m_out(DEBUG) << " First MC Particle : " << energy << " GeV " << " at " << x << "," << y << "," << z << " (|cosTheta| = " << _cosTheta << ")" << std::endl;
            }/*end loop over iQuark*/      
        }
        catch(DataNotAvailableException &e)
        {
            m_out(DEBUG) << "No Collection : " << _inputMCParticleCollections.at(i) << std::endl;
        }
    } /*end loop over i*/

    // read ECAL BARREL hits
    _zmean = 0;

    for (unsigned int i = 0; i < _ecalBarrelCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection( _ecalBarrelCollections.at(i).c_str() );

            if (col != 0) 
            {
                CellIDDecoder<CalorimeterHit> id( col ) ;
                int nelem = col->getNumberOfElements();

                for(int ihit = 0; ihit < nelem; ihit++)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    float hitEnergy = hit->getEnergy();
                    ecalBarrelEnergy += hitEnergy;
                    int layerNumber = id( hit )["K-1"] + 1 ;

                    fEcalBarrelEnergyByLayer->Fill(layerNumber,hitEnergy);

                    float scale = 1.0;

                    if(hit->getType() == 1) // Not sure this still works: that's what enums are for...
                        scale = 2.0;

                    float energyInMips=0;
                    float x = hit->getPosition()[0];
                    float y = hit->getPosition()[1];
                    float z = hit->getPosition()[2];
                    float r = sqrt(x*x+y*y+z*z);
                    _zmean += z;
                    //         nz    +=1.;
                    float correction = 1.;

                    if(fabs(z) < _zOfEndCap)
                    {
                        correction = r / std::sqrt(x*x+y*y);
                    }
                    else
                    {
                        correction = r / z;
                    }

                    energyInMips = hit->getEnergy() * _ecalToMIP / scale;
                    fEcalBarrelMIP->Fill(energyInMips);
                    fEcalBarrelMIPcorr->Fill(energyInMips / correction);
                }
                m_out(DEBUG)<< " ECAL BARREL hits : " << nelem << " energy = " << ecalBarrelEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &e)
        {
            m_out(DEBUG) << "No Collection : " <<  _ecalBarrelCollections[i] << std::endl;
        }
    }

    // read ECAL ENDCAP hits
    _zmean = 0;

    for (unsigned int i = 0; i < _ecalEndCapCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection( _ecalEndCapCollections.at(i).c_str() );

            if (col != 0) 
            {
                CellIDDecoder<CalorimeterHit> id( col ) ;
                int nelem = col->getNumberOfElements();

                for(int ihit = 0; ihit < nelem; ihit++)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    float hitEnergy = hit->getEnergy();
                    ecalEndCapEnergy += hitEnergy;
                    int layerNumber = id( hit )["K-1"] + 1 ;

                    fEcalEndCapEnergyByLayer->Fill(layerNumber,hitEnergy);

                    float scale = 1.0;
                    if (hit->getType() == 1) // Not sure this still works: that's what enums are for...
                        scale =2.0;

                    float energyInMips=0;
                    float x = hit->getPosition()[0];
                    float y = hit->getPosition()[1];
                    float z = hit->getPosition()[2];
                    float r = sqrt(x*x+y*y+z*z);
                    _zmean += z;
                    float correction = 1.;

                    if(fabs(z) < _zOfEndCap)
                    {
                        correction = r / std::sqrt(x*x+y*y);
                    }
                    else
                    {
                        correction = r / z;
                    }

                    energyInMips = hit->getEnergy()*_ecalToMIP / scale;
                    fEcalEndCapMIP->Fill(energyInMips);
                    fEcalEndCapMIPcorr->Fill(energyInMips / correction);
                }
                m_out(DEBUG) << " ECAL ENDCAP hits : " << nelem << " energy = " << ecalEndCapEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &e)
        {
            m_out(DEBUG) << "No Collection : " <<  _ecalEndCapCollections[i] << std::endl;
        }
    }

    // read HCAL hits
    for (unsigned int i = 0; i < _hcalCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection( _hcalCollections.at(i).c_str() );

            if (col != 0) 
            {
                int nelem = col->getNumberOfElements();

                for(int ihit = 0; ihit < nelem; ihit++)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    hcalEnergy += hit->getEnergy();
                    float energyInMips=0;
                    float x = hit->getPosition()[0];
                    float y = hit->getPosition()[1];
                    float z = hit->getPosition()[2];
                    float r = sqrt(x*x+y*y+z*z);
                    float correction = 1.;

                    if(fabs(z) < _zOfEndCap)
                    {
                        correction = r / std::sqrt(x*x+y*y);
                    }
                    else
                    {
                        correction = r / z;
                    }

                    energyInMips = hit->getEnergy() * _hcalToMIP;
                    fHcalMIP->Fill(energyInMips);
                    fHcalMIPcorr->Fill(energyInMips / correction);
                }
                m_out(DEBUG) << " HCAL hits : " << nelem << " energy = " << hcalEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &e)
        {
            m_out(DEBUG) << "No Collection : " <<  _hcalCollections[i] << std::endl;
        }
    }

    // read MUON hits
    for (unsigned int i = 0; i < _muonCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection( _muonCollections.at(i).c_str() );

            if (col != 0) 
            {
                int nelem = col->getNumberOfElements();

                for(int ihit = 0; ihit < nelem; ihit++)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    muonEnergy += hit->getEnergy();
                    float energyInMips=0;
                    float x = hit->getPosition()[0];
                    float y = hit->getPosition()[1];
                    float z = hit->getPosition()[2];
                    float r = sqrt(x*x+y*y+z*z);
                    float correction = 1.;

                    if(fabs(z) < _zOfEndCap)
                    {
                        correction = r / std::sqrt(x*x+y*y);
                    }
                    else
                    {
                        correction = r / z;
                    }
//std::cout << "muon hit x: " << hit->getPosition()[0] << ", y: " << hit->getPosition()[1] << ", z: " << hit->getPosition()[2] << ", E: " << hit->getEnergy() << ", t: " << hit->getTime() << std::endl;
                    energyInMips = hit->getEnergy() * _muonToMIP;
                    fMuonMIP->Fill(energyInMips);
                    fMuonMIPcorr->Fill(energyInMips / correction);
                }
                m_out(DEBUG) << " MUON hits : " << nelem << " energy = " << muonEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &e)
        {
            m_out(DEBUG) << "No Collection : " <<  _muonCollections[i] << std::endl;
        }
    }

    // read LCAL hits
    for (unsigned int i = 0; i < _lcalCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection( _lcalCollections.at(i).c_str() );

            if (col != 0) 
            {
                int nelem = col->getNumberOfElements();

                for(int ihit = 0; ihit < nelem; ihit++)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    lcalEnergy += hit->getEnergy();
                }
                m_out(DEBUG) << " LCAL hits : " << nelem << " energy = " << lcalEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &e)
        {
            m_out(DEBUG) << "No Collection : " <<  _lcalCollections[i] << std::endl;
        }
    }

    // read BCAL hits
    for (unsigned int i = 0; i < _bcalCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection( _bcalCollections.at(i).c_str() );

            if (col != 0) 
            {
                int nelem = col->getNumberOfElements();

                for(int ihit = 0; ihit < nelem; ihit++)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    bcalEnergy += hit->getEnergy();
                }
                m_out(DEBUG) << " BCAL hits : " << nelem << " energy = " << bcalEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &e)
        {
            m_out(DEBUG) << "No Collection : " <<  _bcalCollections[i] << std::endl;
        }
    }

    // read LHCAL hits
    for (unsigned int i=0; i < _lhcalCollections.size(); ++i) 
    {
        try
        {
            LCCollection * col = evt->getCollection( _lhcalCollections.at(i).c_str() );

            if (col != 0) 
            {
                int nelem = col->getNumberOfElements();

                for(int ihit = 0; ihit < nelem; ihit++)
                {
                    CalorimeterHit * hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(ihit));
                    lhcalEnergy += hit->getEnergy();
                }
                m_out(DEBUG) << " LHCAL hits : " << nelem << " energy = " << lhcalEnergy << std::endl;
            }
        }
        catch(DataNotAvailableException &e)
        {
            m_out(DEBUG) << "No Collection : " <<  _lhcalCollections[i] << std::endl;
        }
    }

    // Total energy
    totEnergy = ecalBarrelEnergy + ecalEndCapEnergy + hcalEnergy + muonEnergy + lcalEnergy + bcalEnergy + lhcalEnergy;
    m_out(DEBUG) << " Total Calorimetric Energy : " << totEnergy << std::endl;

    if (_cosTheta < 0.95)
    {
        fEcalEnergy->Fill(ecalBarrelEnergy + ecalEndCapEnergy,1.);
        fEcalBarrelHcalEnergy->Fill((ecalBarrelEnergy * _ecalBarrelHadMIPToGeV) / _ecalEMMIPToGeV, hcalEnergy, 1.);
        fEcalEndCapHcalEnergy->Fill((ecalEndCapEnergy * _ecalEndCapHadMIPToGeV) / _ecalEMMIPToGeV, hcalEnergy, 1.);
        fHcalEnergy->Fill(hcalEnergy,1.);
        fMuonEnergy->Fill(muonEnergy,1.);
        fLcalEnergy->Fill(lcalEnergy,1.);
        fCalEnergy->Fill(totEnergy,1.);
        if (totEnergy > 0 && (ecalBarrelEnergy + ecalEndCapEnergy) / totEnergy > 0.95)
            fCalEnergyE->Fill(totEnergy);
        if (totEnergy > 0 && hcalEnergy / totEnergy > 0.95)
            fCalEnergyH->Fill(totEnergy);
        if (totEnergy > 0 && muonEnergy / totEnergy > 0.95)
            fCalEnergyM->Fill(totEnergy);
    }
    fCalEnergyVsCosTheta->Fill(_cosTheta,totEnergy,1.);

    // Read reconstructed PFOs
    try
    {
        LCCollection * col = evt->getCollection( _particleCollectionName.c_str() );

        if (col != 0) 
        {
            int nelem = col->getNumberOfElements();
            bool check=true;
            float photonE = 0;
            float p[3] = {0.,0.,0.};
            _zCoG=0.;
            float esum = 0;

            for(int j = 0; j < nelem; j++)
            {
                ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(j));
                pfoEnergy+=recoPart->getEnergy();

                TrackVec tracks = recoPart->getTracks();

                if (tracks.size()!= 0)
                    check = false;

                if (recoPart->getType()== 22)
                    photonE += recoPart->getEnergy();

                p[0] += recoPart->getMomentum()[0];
                p[1] += recoPart->getMomentum()[1];
                p[2] += recoPart->getMomentum()[2];

                ClusterVec clusters = recoPart->getClusters();
                for(unsigned int icluster = 0; icluster < clusters.size(); icluster++)
                {
                    _zCoG+= clusters[icluster]->getPosition()[2]*clusters[icluster]->getEnergy();
                    esum += clusters[icluster]->getEnergy();
                }
            }/*end loop over j*/

            if (esum > 0)
            {
                _zCoG = _zCoG / esum;
            }
            else
            {
                _zCoG = std::numeric_limits<float>::max();
            }

            _cosThetaR = std::fabs(p[2]) / std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
            m_out(DEBUG) << " PFO objects : " << nelem << " energy = " << pfoEnergy << std::endl;

            if(check)
            {
                fCosT->Fill(_cosTheta,1.);

                if (pfoEnergy > 0 && photonE > 0.5*pfoEnergy)
                    fPhotonCosT->Fill(_cosTheta,1.);
            }

            fPFAVsCosTheta->Fill(_cosTheta,pfoEnergy,1.);
            fPFAVsCosThetaR->Fill(_cosThetaR,pfoEnergy,1.);
            fPFAVsZCoG->Fill(fabs(_zCoG),pfoEnergy,1.);
            fPFAVsCosThetaX->Fill(_cosThetaX,pfoEnergy,1.);

            if (_cosTheta < 0.95)
            {
                fPFA->Fill(pfoEnergy,1.);

                if (_cosTheta < 0.7)
                    fPFAB->Fill(pfoEnergy,1.);

                if (totEnergy > 0 && (ecalBarrelEnergy + ecalEndCapEnergy) / totEnergy > 0.95)
                    fPFAE->Fill(pfoEnergy,1.);

                if (totEnergy > 0 && hcalEnergy / totEnergy > 0.95)
                    fPFAH->Fill(pfoEnergy,1.);

                if (totEnergy > 0 && muonEnergy / totEnergy > 0.95)
                    fPFAM->Fill(pfoEnergy,1.);
            }

            fXvsY->Fill(_x, _y, 1.);
        }
    }
    catch(DataNotAvailableException &e)
    {
        m_out(DEBUG) << "No Collection : " <<  _particleCollectionName << std::endl;
    }

    fCalEnergyVsCosThetaR->Fill(_cosThetaR, totEnergy,1.);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::check(LCEvent *pLCEvent)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraPFACalibrator::end()
{
    std::cout << "PandoraPFACalibrator::end()  " << name() << " processed " << _nEvt << " events in " << _nRun << " runs " << std::endl ;

    TFile *hfile = new TFile(_rootFile.c_str(),"recreate");
    fPFA->Write();
    fPFAB->Write();
    fPFAVsCosTheta->Write();
    fPFAVsCosThetaR->Write();
    fPFAVsZCoG->Write();
    fPFAVsCosThetaX->Write();

    fPFAE->Write();
    fPFAH->Write();
    fPFAM->Write();

    fXvsY->Write();

    fEcalEnergy->Write();
    fHcalEnergy->Write();
    fMuonEnergy->Write();
    fLcalEnergy->Write();
    fCalEnergy->Write();

    fEcalBarrelHcalEnergy->Write();
    fEcalEndCapHcalEnergy->Write();

    fCalEnergyE->Write();
    fCalEnergyH->Write();
    fCalEnergyM->Write();

    fCalEnergyVsCosTheta->Write();
    fCalEnergyVsCosThetaR->Write();

    fEcalBarrelEnergyByLayer->Write();
    fEcalEndCapEnergyByLayer->Write();

    fEcalBarrelMIP->Write();
    fEcalEndCapMIP->Write();
    fHcalMIP->Write();
    fMuonMIP->Write();

    fEcalBarrelMIPcorr->Write();
    fEcalEndCapMIPcorr->Write();
    fHcalMIPcorr->Write();
    fMuonMIPcorr->Write();

    fCosT->Write();
    fPhotonCosT->Write();

    hfile->Close();
    delete hfile;
}
