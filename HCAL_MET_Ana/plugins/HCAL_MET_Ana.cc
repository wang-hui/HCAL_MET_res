// -*- C++ -*-
//
// Package:    HCALTest/HCAL_MET_Ana
// Class:      HCAL_MET_Ana
//
/**\class HCAL_MET_Ana HCAL_MET_Ana.cc HCALTest/HCAL_MET_Ana/plugins/HCAL_MET_Ana.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//DataFormats and Geometry
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

//ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include <Math/VectorUtil.h>

//STL headers
#include <vector>
#include <string>
#include <iostream>
#include <numeric>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class HCAL_MET_Ana : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit HCAL_MET_Ana(const edm::ParameterSet&);
        ~HCAL_MET_Ana();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        std::vector<reco::Muon> sel_muons(std::vector<reco::Muon> Muons);
        std::vector<reco::GsfElectron> sel_electrons(std::vector<reco::GsfElectron> Electrons);
        std::vector<math::XYZTLorentzVector> sel_z_mumu(std::vector<reco::Muon> SelMuons);
        std::vector<math::XYZTLorentzVector> sel_z_ee(std::vector<reco::GsfElectron> SelElectrons);
        float calc_ht(std::vector<reco::PFJet> PFJets);

        bool PrintChannel_;
        bool IsMC_;
        std::string RunMod_;
        bool PassZSel = false;

        edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitsToken_;
        edm::EDGetTokenT<CaloTowerCollection> CaloTowersToken_;
        edm::EDGetTokenT<std::vector<reco::CaloMET>> CaloMETToken_;
        edm::EDGetTokenT<std::vector<reco::CaloMET>> CaloMETBEToken_;
        edm::EDGetTokenT<std::vector<reco::Muon>> MuonToken_;
        edm::EDGetTokenT<std::vector<reco::GsfElectron>> ElectronToken_;
        edm::EDGetTokenT<std::vector<reco::CaloJet>> CaloJetToken_;
        edm::EDGetTokenT<std::vector<reco::PFJet>> PFJetToken_;
        edm::EDGetTokenT<std::vector<reco::Vertex>> VertexToken_;
        edm::EDGetTokenT<std::vector<reco::GenMET>> GenMETToken_;
        edm::EDGetTokenT<std::vector<reco::GenJet>> GenJetToken_;

        TH1F * BaselineTest_h;
        TH1F * PU_h, * HT_h;
        TH1F * nAllEle_h, *nSelEle_h;
        TH1F * GenMET_h, * GenMET_phi_h;
        TH1F * CaloMET_h, * CaloMET_phi_h;
        TH1F * CaloMETBE_h, * CaloMETBE_phi_h, * CaloMETBE20_phi_h;
        TH1F * myCaloMETBE_h, * myCaloMETBE_phi_h;
        TH1F * myCaloMETBE_HB_h, * myCaloMETBE_HE_h;
        TH1F * myCaloMETBE1_h, * myCaloMETBE1_phi_h;
        TH1F * myCaloMETBE_HPU_h, * myCaloMETBE_phi_HPU_h;
        TH1F * myCaloMETBE1_HPU_h, * myCaloMETBE1_phi_HPU_h;
        TH1F * myCaloMETBE_HHT_h, * myCaloMETBE_phi_HHT_h;
        TH1F * myCaloMETBE1_HHT_h, * myCaloMETBE1_phi_HHT_h;
        TH1F * myCaloETBE_h, * myCaloETBE1_h;

        const int METResArraySize = 9;
        const double METResArray [10] = {0.0, 5.0, 15.0, 30.0, 50.0, 75.0, 105.0, 140.0, 180.0, 225.0};

        TH2F * CaloMETBE_vs_PU_h, * CaloMETBE_phi_vs_PU_h;
        TH2F * UPara_ratio_vs_Zpt_h, * UPara_vs_Zpt_h, * UVert_vs_Zpt_h;
        TH2F * CaloJet_vs_GenJet_h, * CaloJet_vs_GenJet_etaL_h, * CaloJet_vs_GenJet_etaM_h, * CaloJet_vs_GenJet_etaH_h;
        TH2F * CaloJet_vs_GenJet_pull_h, * CaloJet_vs_GenJet_etaL_pull_h, * CaloJet_vs_GenJet_etaM_pull_h, * CaloJet_vs_GenJet_etaH_pull_h;
        TH2F * myCaloETBE_vs_eta_h;
        TH2F * CaloTowerET_vs_eta_h, * CaloTowerEMET_vs_eta_h, * CaloTowerHadET_vs_eta_h;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HCAL_MET_Ana::HCAL_MET_Ana(const edm::ParameterSet& iConfig):
    PrintChannel_(iConfig.getUntrackedParameter<bool>("PrintChannel")),
    IsMC_(iConfig.getUntrackedParameter<bool>("IsMC")),
    RunMod_(iConfig.getUntrackedParameter<std::string>("RunMod"))
{
    if(RunMod_ != "Zmumu" && RunMod_ != "Zee")
    {
        std::cout << "RunMod is not Zmumu or Zee, set PassZSel to true" << std::endl;
        PassZSel = true;
    }


    //now do what ever initialization is needed
    HBHERecHitsToken_ = consumes<HBHERecHitCollection>(edm::InputTag("hbhereco"));
    CaloTowersToken_ = consumes<CaloTowerCollection>(edm::InputTag("towerMaker"));
    CaloMETToken_ = consumes<std::vector<reco::CaloMET>>(edm::InputTag("caloMet"));
    CaloMETBEToken_ = consumes<std::vector<reco::CaloMET>>(edm::InputTag("caloMetBE"));
    MuonToken_ = consumes<std::vector<reco::Muon>>(edm::InputTag("muons"));
    ElectronToken_ = consumes<std::vector<reco::GsfElectron>>(edm::InputTag("gedGsfElectrons"));
    CaloJetToken_ = consumes<std::vector<reco::CaloJet>>(edm::InputTag("ak4CaloJets"));
    PFJetToken_ = consumes<std::vector<reco::PFJet>>(edm::InputTag("ak4PFJetsCHS"));
    VertexToken_ = consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"));
    if(IsMC_)
    {
        GenMETToken_ = consumes<std::vector<reco::GenMET>>(edm::InputTag("genMetTrue"));
        GenJetToken_ = consumes<std::vector<reco::GenJet>>(edm::InputTag("ak4GenJetsNoNu"));
    }

    edm::Service<TFileService> TFS;
    BaselineTest_h = TFS->make<TH1F>("BaselineTest_h", "0: all. 1: pass", 2, 0.0, 2.0);
    PU_h = TFS->make<TH1F>("PU_h", "PU_h", 100, 0.0, 100.0);
    HT_h = TFS->make<TH1F>("HT_h", "HT_h", 200, 0.0, 1000.0);
    nAllEle_h = TFS->make<TH1F>("nAllEle_h", "nAllEle_h", 10, 0.0, 10.0);
    nSelEle_h = TFS->make<TH1F>("nSelEle_h", "nSelEle_h", 10, 0.0, 10.0);
    CaloMET_h = TFS->make<TH1F>("CaloMET_h", "CaloMET_h", 100, 0.0, 200.0);
    CaloMET_phi_h = TFS->make<TH1F>("CaloMET_phi_h", "CaloMET_phi_h", 100, -3.2, 3.2);
    CaloMETBE_h = TFS->make<TH1F>("CaloMETBE_h", "CaloMETBE_h", 100, 0.0, 200.0);
    CaloMETBE_phi_h = TFS->make<TH1F>("CaloMETBE_phi_h", "CaloMETBE_phi_h", 100, -3.2, 3.2);
    CaloMETBE20_phi_h = TFS->make<TH1F>("CaloMETBE20_phi_h", "CaloMETBE20_phi_h", 100, -3.2, 3.2);

    myCaloMETBE_h = TFS->make<TH1F>("myCaloMETBE_h", "myCaloMETBE_h", 100, 0.0, 200.0);
    myCaloMETBE_phi_h = TFS->make<TH1F>("myCaloMETBE_phi_h", "myCaloMETBE_phi_h", 100, -3.2, 3.2);
    myCaloMETBE_HB_h = TFS->make<TH1F>("myCaloMETBE_HB_h", "myCaloMETBE_HB_h", 100, 0.0, 200.0);
    myCaloMETBE_HE_h = TFS->make<TH1F>("myCaloMETBE_HE_h", "myCaloMETBE_HE_h", 100, 0.0, 200.0);
    myCaloMETBE1_h = TFS->make<TH1F>("myCaloMETBE1_h", "myCaloMETBE1_h", 100, 0.0, 200.0);
    myCaloMETBE1_phi_h = TFS->make<TH1F>("myCaloMETBE1_phi_h", "myCaloMETBE1_phi_h", 100, -3.2, 3.2);

    myCaloMETBE_HPU_h = TFS->make<TH1F>("myCaloMETBE_HPU_h", "myCaloMETBE_HPU_h", 100, 0.0, 200.0);
    myCaloMETBE_phi_HPU_h = TFS->make<TH1F>("myCaloMETBE_phi_HPU_h", "myCaloMETBE_phi_HPU_h", 100, -3.2, 3.2);
    myCaloMETBE1_HPU_h = TFS->make<TH1F>("myCaloMETBE1_HPU_h", "myCaloMETBE1_HPU_h", 100, 0.0, 200.0);
    myCaloMETBE1_phi_HPU_h = TFS->make<TH1F>("myCaloMETBE1_phi_HPU_h", "myCaloMETBE1_phi_HPU_h", 100, -3.2, 3.2);

    myCaloMETBE_HHT_h = TFS->make<TH1F>("myCaloMETBE_HHT_h", "myCaloMETBE_HHT_h", 100, 0.0, 200.0);
    myCaloMETBE_phi_HHT_h = TFS->make<TH1F>("myCaloMETBE_phi_HHT_h", "myCaloMETBE_phi_HHT_h", 100, -3.2, 3.2);
    myCaloMETBE1_HHT_h = TFS->make<TH1F>("myCaloMETBE1_HHT_h", "myCaloMETBE1_HHT_h", 100, 0.0, 200.0);
    myCaloMETBE1_phi_HHT_h = TFS->make<TH1F>("myCaloMETBE1_phi_HHT_h", "myCaloMETBE1_phi_HHT_h", 100, -3.2, 3.2);

    myCaloETBE_h = TFS->make<TH1F>("myCaloETBE_h", "myCaloETBE_h", 200, 0.0, 1000.0);
    myCaloETBE1_h = TFS->make<TH1F>("myCaloETBE1_h", "myCaloETBE1_h", 200, 0.0, 1000.0);

    CaloTowerET_vs_eta_h = TFS->make<TH2F>("CaloTowerET_vs_eta_h", "CaloTowerET_vs_eta_h", 100, -5.0, 5.0, 200, 0.0, 5.0);
    CaloTowerEMET_vs_eta_h = TFS->make<TH2F>("CaloTowerEMET_vs_eta_h", "CaloTowerEMET_vs_eta_h", 100, -5.0, 5.0, 200, 0.0, 5.0);
    CaloTowerHadET_vs_eta_h = TFS->make<TH2F>("CaloTowerHadET_vs_eta_h", "CaloTowerHadET_vs_eta_h", 100, -5.0, 5.0, 200, 0.0, 5.0);
    myCaloETBE_vs_eta_h = TFS->make<TH2F>("myCaloETBE_vs_eta_h", "myCaloETBE_vs_eta_h", 60, -3.0, 3.0, 200, 0.0, 1.0);
    CaloMETBE_vs_PU_h = TFS->make<TH2F>("CaloMETBE_vs_PU_h", "CaloMETBE_vs_PU_h", 100, 0.0, 100.0, 100, 0.0, 200.0);
    CaloMETBE_phi_vs_PU_h = TFS->make<TH2F>("CaloMETBE_phi_vs_PU_h", "CaloMETBE_phi_vs_PU_h", 100, 0.0, 100.0, 100, -3.2, 3.2);

    UPara_ratio_vs_Zpt_h = TFS->make<TH2F>("UPara_ratio_vs_Zpt_h", "UPara_ratio_vs_Zpt_h", 20, 0.0, 200.0, 200, -10.0, 10.0);
    UPara_vs_Zpt_h = TFS->make<TH2F>("UPara_vs_Zpt_h", "UPara_vs_Zpt_h", METResArraySize, METResArray, 200, -300.0, 200.0);
    UVert_vs_Zpt_h = TFS->make<TH2F>("UVert_vs_Zpt_h", "UVert_vs_Zpt_h", METResArraySize, METResArray, 200, -150.0, 150.0);
    if(IsMC_)
    {
        GenMET_h = TFS->make<TH1F>("GenMET_h", "GenMET_h", 100, 0.0, 200.0);
        GenMET_phi_h = TFS->make<TH1F>("GenMET_phi_h", "GenMET_phi_h", 100, -3.2, 3.2);
        CaloJet_vs_GenJet_h = TFS->make<TH2F>("CaloJet_vs_GenJet_h", "CaloJet_vs_GenJet_h", 200, 0.0, 1000.0, 200, 0.0, 1000.0);
        CaloJet_vs_GenJet_etaL_h = TFS->make<TH2F>("CaloJet_vs_GenJet_etaL_h", "CaloJet_vs_GenJet_etaL_h", 200, 0.0, 1000.0, 200, 0.0, 1000.0);
        CaloJet_vs_GenJet_etaM_h = TFS->make<TH2F>("CaloJet_vs_GenJet_etaM_h", "CaloJet_vs_GenJet_etaM_h", 200, 0.0, 1000.0, 200, 0.0, 1000.0);
        CaloJet_vs_GenJet_etaH_h = TFS->make<TH2F>("CaloJet_vs_GenJet_etaH_h", "CaloJet_vs_GenJet_etaH_h", 200, 0.0, 1000.0, 200, 0.0, 1000.0);
        CaloJet_vs_GenJet_pull_h = TFS->make<TH2F>("CaloJet_vs_GenJet_pull_h", "CaloJet_vs_GenJet_pull_h", 200, 0.0, 1000.0, 200, 0.0, 10.0);
        CaloJet_vs_GenJet_etaL_pull_h = TFS->make<TH2F>("CaloJet_vs_GenJet_etaL_pull_h", "CaloJet_vs_GenJet_etaL_pull_h", 200, 0.0, 1000.0, 200, 0.0, 10.0);
        CaloJet_vs_GenJet_etaM_pull_h = TFS->make<TH2F>("CaloJet_vs_GenJet_etaM_pull_h", "CaloJet_vs_GenJet_etaM_pull_h", 200, 0.0, 1000.0, 200, 0.0, 10.0);
        CaloJet_vs_GenJet_etaH_pull_h = TFS->make<TH2F>("CaloJet_vs_GenJet_etaH_pull_h", "CaloJet_vs_GenJet_etaH_pull_h", 200, 0.0, 1000.0, 200, 0.0, 10.0);
    }
}


HCAL_MET_Ana::~HCAL_MET_Ana()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------

void HCAL_MET_Ana::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    math::XYZTLorentzVector SelZ;

    if(RunMod_ == "Zmumu")
    {
        edm::Handle<std::vector<reco::Muon>> MuonHandle;
        iEvent.getByToken(MuonToken_, MuonHandle);
        auto Muons = MuonHandle.product();
        auto SelMuons = sel_muons(*Muons);
        if(SelMuons.size() == 2)
        {
            auto SelZs = sel_z_mumu(SelMuons);
            if(SelZs.size() == 1)
            {
                SelZ = SelZs.at(0);
                PassZSel = true;
            }
        }
    }

    else if(RunMod_ == "Zee")
    {
        edm::Handle<std::vector<reco::GsfElectron>> ElectronHandle;
        iEvent.getByToken(ElectronToken_, ElectronHandle);
        auto Electrons = ElectronHandle.product();
        auto SelElectrons = sel_electrons(*Electrons);
        nAllEle_h->Fill(Electrons->size());
        nSelEle_h->Fill(SelElectrons.size());
        if(SelElectrons.size() == 2)
        {
            auto SelZs = sel_z_ee(SelElectrons);
            if(SelZs.size() == 1)
            {
                SelZ = SelZs.at(0);
                PassZSel = true;
            }
        }
    }

    BaselineTest_h->Fill(0);
    if(PassZSel)
    {
        BaselineTest_h->Fill(1);

        edm::Handle<std::vector<reco::Vertex>> VertexHandle;
        iEvent.getByToken(VertexToken_, VertexHandle);
        auto Vertex = VertexHandle.product();
        int nPU = Vertex->size();
        PU_h->Fill(nPU);

        edm::Handle<std::vector<reco::PFJet>> PFJetHandle;
        iEvent.getByToken(PFJetToken_, PFJetHandle);
        auto PFJets = PFJetHandle.product();
        auto HT = calc_ht(*PFJets);
        HT_h->Fill(HT);

        edm::ESHandle<CaloGeometry> CaloGeoHandle;
        iSetup.get<CaloGeometryRecord>().get(CaloGeoHandle);
        CaloGeometry CaloGeo = *CaloGeoHandle;

        edm::Handle<CaloTowerCollection> CaloTowersHandle;
        iEvent.getByToken(CaloTowersToken_, CaloTowersHandle);
        auto CaloTowers = CaloTowersHandle.product();
        for(auto CaloTower : *CaloTowers)
        {
            //std::cout << "ET " << CaloTower.et() << ", hadEt " << CaloTower.hadEt() << ", emEt " << CaloTower.emEt() << std::endl;
            //std::cout << CaloTower.ieta() << ", " << CaloTower.eta() << ", " << CaloTower.p4().Eta() << std::endl;
            CaloTowerET_vs_eta_h->Fill(CaloTower.eta(), CaloTower.et()/nPU);
            CaloTowerEMET_vs_eta_h->Fill(CaloTower.eta(), CaloTower.emEt()/nPU);
            CaloTowerHadET_vs_eta_h->Fill(CaloTower.eta(), CaloTower.hadEt()/nPU);
        }

        edm::Handle<HBHERecHitCollection> HBHERecHitsHandle;
        iEvent.getByToken(HBHERecHitsToken_, HBHERecHitsHandle);
        auto HBHERecHits = HBHERecHitsHandle.product();

        math::PtEtaPhiMLorentzVector HBTotLV, HETotLV, HBHETotLV, HBHE1TotLV;
        float HBHEET = 0;
        float HBHE1ET = 0;
        for(auto HBHERecHit : *HBHERecHits)
        {
            auto Hid = HBHERecHit.id();
            auto RawId = Hid.rawId();
            auto SubDet = Hid.subdet();
            auto Depth = Hid.depth();
            auto Ieta = Hid.ieta();
            auto Iphi = Hid.iphi();
            auto Energy = HBHERecHit.energy();

            auto Eta = CaloGeo.getPosition(Hid).eta();
            auto Phi = CaloGeo.getPosition(Hid).phi();

            auto Theta = 2 * TMath::ATan(exp(-1 * Eta));
            auto ET = Energy * TMath::Sin(Theta);

            myCaloETBE_vs_eta_h->Fill(Eta, ET/nPU);

            if(Energy > 0)
            {
                HBHETotLV += math::PtEtaPhiMLorentzVector(ET, Eta, Phi, 0);
                if(SubDet == 1) HBTotLV += math::PtEtaPhiMLorentzVector(ET, Eta, Phi, 0);
                else if(SubDet == 2) HETotLV += math::PtEtaPhiMLorentzVector(ET, Eta, Phi, 0);
                HBHEET += ET;
                if(ET > 1)
                {
                    HBHE1TotLV += math::PtEtaPhiMLorentzVector(ET, Eta, Phi, 0);
                    HBHE1ET += ET;
                }
            }

            if(PrintChannel_) std::cout << Hid << ", " << RawId << ", " << SubDet << ", " << Depth << ", " << Ieta << ", " << Eta << ", " << Iphi << ", " << Phi << ", " << Energy << std::endl;
        }

        auto myCaloMETBE = - HBHETotLV;
        auto myCaloMETBE_pt = myCaloMETBE.Pt();
        auto myCaloMETBE_phi = myCaloMETBE.Phi();
        myCaloMETBE_h->Fill(myCaloMETBE_pt);
        myCaloMETBE_phi_h->Fill(myCaloMETBE_phi);
        if(nPU > 30)
        {
            myCaloMETBE_HPU_h->Fill(myCaloMETBE_pt);
            myCaloMETBE_phi_HPU_h->Fill(myCaloMETBE_phi);
        }
        if(HT > 100)
        {
            myCaloMETBE_HHT_h->Fill(myCaloMETBE_pt);
            myCaloMETBE_phi_HHT_h->Fill(myCaloMETBE_phi);
        }
        myCaloETBE_h->Fill(HBHEET);

        myCaloMETBE_HB_h->Fill((- HBTotLV).Pt());
        myCaloMETBE_HE_h->Fill((- HETotLV).Pt());

        //=========== MET response and resolution start =========
        auto SelZ_pt = SelZ.Pt();
        auto dPhiZMET = myCaloMETBE_phi - SelZ.Phi();
        auto METPara = myCaloMETBE_pt * TMath::Cos(dPhiZMET);
        auto METVert = myCaloMETBE_pt * TMath::Sin(dPhiZMET);
        auto UPara = - METPara - SelZ_pt;
        auto UVert = - METVert;

        UPara_ratio_vs_Zpt_h->Fill(SelZ_pt, - UPara / SelZ_pt);
        UPara_vs_Zpt_h->Fill(SelZ_pt, UPara);
        UVert_vs_Zpt_h->Fill(SelZ_pt, UVert);
        //=========== MET response and resolution end ==========

        auto myCaloMETBE1 = - HBHE1TotLV;
        auto myCaloMETBE1_pt = myCaloMETBE1.Pt();
        auto myCaloMETBE1_phi = myCaloMETBE1.Phi();
        myCaloMETBE1_h->Fill(myCaloMETBE1_pt);
        myCaloMETBE1_phi_h->Fill(myCaloMETBE1_phi);
        if(nPU > 30)
        {
            myCaloMETBE1_HPU_h->Fill(myCaloMETBE1_pt);
            myCaloMETBE1_phi_HPU_h->Fill(myCaloMETBE1_phi);
        }
        if(HT > 100)
        {
            myCaloMETBE1_HHT_h->Fill(myCaloMETBE1_pt);
            myCaloMETBE1_phi_HHT_h->Fill(myCaloMETBE1_phi);
        }
        myCaloETBE1_h->Fill(HBHE1ET);

        if(IsMC_)
        {
            edm::Handle<std::vector<reco::GenMET>> GenMETHandle;
            iEvent.getByToken(GenMETToken_, GenMETHandle);
            auto GenMET = GenMETHandle.product();

            if(GenMET->size() == 1)
            {
                GenMET_h->Fill(GenMET->at(0).p4().Pt());
                GenMET_phi_h->Fill(GenMET->at(0).p4().Phi());
            }

            edm::Handle<std::vector<reco::GenJet>> GenJetHandle;
            iEvent.getByToken(GenJetToken_, GenJetHandle);
            auto GenJets = GenJetHandle.product();

            edm::Handle<std::vector<reco::CaloJet>> CaloJetHandle;
            iEvent.getByToken(CaloJetToken_, CaloJetHandle);
            auto CaloJets = CaloJetHandle.product();

            for(auto GenJet : *GenJets)
            {
                auto GenJetP4 = GenJet.p4();
                if(fabs(GenJetP4.Eta()) < 3.0)
                {
                    for(auto CaloJet : *CaloJets)
                    {
                        auto CaloJetP4 = CaloJet.p4();
                        if(ROOT::Math::VectorUtil::DeltaR(GenJetP4, CaloJetP4) < 0.2)
                        {
                            auto GenPt = GenJetP4.Pt();
                            auto RecoPt = CaloJetP4.Pt();
                            float PtPull = 0;
                            if(GenPt + RecoPt > 0) PtPull = fabs(GenPt - RecoPt) / sqrt(GenPt + RecoPt);

                            CaloJet_vs_GenJet_h->Fill(GenPt, RecoPt);
                            CaloJet_vs_GenJet_pull_h->Fill(GenPt, PtPull);
                            if(fabs(GenJetP4.Eta()) < 1.3)
                            {
                                CaloJet_vs_GenJet_etaL_h->Fill(GenPt, RecoPt);
                                CaloJet_vs_GenJet_etaL_pull_h->Fill(GenPt, PtPull);
                            }
                            else if(fabs(GenJetP4.Eta()) < 2.5)
                            {
                                CaloJet_vs_GenJet_etaM_h->Fill(GenPt, RecoPt);
                                CaloJet_vs_GenJet_etaM_pull_h->Fill(GenPt, PtPull);
                            }
                            else
                            {
                                CaloJet_vs_GenJet_etaH_h->Fill(GenPt, RecoPt);
                                CaloJet_vs_GenJet_etaH_pull_h->Fill(GenPt, PtPull);
                            }
                        }
                    }
                }
            }
        }

        edm::Handle<std::vector<reco::CaloMET>> CaloMETHandle;
        iEvent.getByToken(CaloMETToken_, CaloMETHandle);
        auto CaloMET = CaloMETHandle.product();

        if(CaloMET->size() == 1)
        {
            CaloMET_h->Fill(CaloMET->at(0).p4().Pt());
            CaloMET_phi_h->Fill(CaloMET->at(0).p4().Phi());
        }

        edm::Handle<std::vector<reco::CaloMET>> CaloMETBEHandle;
        iEvent.getByToken(CaloMETBEToken_, CaloMETBEHandle);
        auto CaloMETBE = CaloMETBEHandle.product();

        if(CaloMETBE->size() == 1)
        {
            auto CaloMETBE_pt = CaloMETBE->at(0).p4().Pt();
            auto CaloMETBE_phi = CaloMETBE->at(0).p4().Phi();

            CaloMETBE_h->Fill(CaloMETBE_pt);
            CaloMETBE_phi_h->Fill(CaloMETBE_phi);
            if(CaloMETBE_pt > 20){CaloMETBE20_phi_h->Fill(CaloMETBE_phi);}
            CaloMETBE_vs_PU_h->Fill(nPU, CaloMETBE_pt);
            CaloMETBE_phi_vs_PU_h->Fill(nPU, CaloMETBE_phi);
        }
    }
}

std::vector<reco::Muon> HCAL_MET_Ana::sel_muons(std::vector<reco::Muon> Muons)
{
    auto MuIdCut = reco::Muon::CutBasedIdMedium;
    //auto MuIsoCut = reco::Muon::MiniIsoMedium;
    auto LepIsoCut = reco::Muon::PFIsoMedium;
    float LepPtCut = 10;
    float EtaCut = 2.4;

    std::vector<reco::Muon> SelMuons = {};
    for(auto Muon : Muons)
    {
        if (Muon.p4().Pt() > LepPtCut && fabs(Muon.p4().Eta()) < EtaCut
                && Muon.passed(MuIdCut) && Muon.passed(LepIsoCut))
        {SelMuons.push_back(Muon);}
    }
    return SelMuons;
}

std::vector<reco::GsfElectron> HCAL_MET_Ana::sel_electrons(std::vector<reco::GsfElectron> Electrons)
{
    float LepPtCut = 10;
    float EtaCut = 2.4;

    std::vector<reco::GsfElectron> SelElectrons = {};
    for(auto Electron : Electrons)
    {
        if (Electron.p4().Pt() > LepPtCut && fabs(Electron.p4().Eta()) < EtaCut)
        {SelElectrons.push_back(Electron);}
    }
    return SelElectrons;
}

std::vector<math::XYZTLorentzVector> HCAL_MET_Ana::sel_z_mumu(std::vector<reco::Muon> SelMuons)
{
    std::vector<math::XYZTLorentzVector> SelZmumu = {};
    auto Mu1 = SelMuons.at(0);
    auto Mu2 = SelMuons.at(1);
    if(std::max(Mu1.p4().Pt(), Mu2.p4().Pt()) > 20 && (Mu1.charge() + Mu2.charge() == 0))
    {
        auto ZmumuCand = Mu1.p4() + Mu2.p4();
        if(ZmumuCand.M() > 81 && ZmumuCand.M() < 101) SelZmumu.push_back(ZmumuCand);
    }
    return SelZmumu;
}

std::vector<math::XYZTLorentzVector> HCAL_MET_Ana::sel_z_ee(std::vector<reco::GsfElectron> SelElectrons)
{
    std::vector<math::XYZTLorentzVector> SelZee = {};
    auto Ele1 = SelElectrons.at(0);
    auto Ele2 = SelElectrons.at(1);
    if(std::max(Ele1.p4().Pt(), Ele2.p4().Pt()) > 20 && (Ele1.charge() + Ele2.charge() == 0))
    {
        auto ZeeCand = Ele1.p4() + Ele2.p4();
        if(ZeeCand.M() > 81 && ZeeCand.M() < 101) SelZee.push_back(ZeeCand);
    }
    return SelZee;
}

float HCAL_MET_Ana::calc_ht(std::vector<reco::PFJet> PFJets)
{
    float JetPtCut = 20;
    float EtaCut = 2.4;

    float HT = 0;
    for(auto PFJet : PFJets)
    {
        if (PFJet.p4().Pt() > JetPtCut && fabs(PFJet.p4().Eta()) < EtaCut)
        {HT += PFJet.p4().Pt();}
    }
    return HT;
}

// ------------ method called once each job just before starting event loop  ------------
void HCAL_MET_Ana::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HCAL_MET_Ana::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HCAL_MET_Ana::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    //Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HCAL_MET_Ana);
