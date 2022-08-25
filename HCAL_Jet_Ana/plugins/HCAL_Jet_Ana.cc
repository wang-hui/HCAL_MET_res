// -*- C++ -*-
//
// Package:    HCALTest/HCAL_Jet_Ana
// Class:      HCAL_Jet_Ana
//
/**\class HCAL_Jet_Ana HCAL_Jet_Ana.cc HCALTest/HCAL_Jet_Ana/plugins/HCAL_Jet_Ana.cc

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
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Geometry/Records/interface/HcalRecNumberingRecord.h"

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

class HCAL_Jet_Ana : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit HCAL_Jet_Ana(const edm::ParameterSet&);
        ~HCAL_Jet_Ana();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        std::vector<reco::CaloJet> select_CaloJets(const std::vector<reco::GenJet> * GenJets, const std::vector<reco::CaloJet> * CaloJets);

        edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitsToken_;
        edm::EDGetTokenT<std::vector<reco::CaloJet>> CaloJetToken_;
        edm::EDGetTokenT<std::vector<reco::Vertex>> VertexToken_;
        edm::EDGetTokenT<std::vector<reco::GenJet>> GenJetToken_;

        TH1F * PU_h;

        TTree * GenJetTree;
        std::vector<math::XYZTLorentzVector> GenJetVec_p4;

        TTree * CaloJetTree;
        std::vector<float> CaloJetVec_Energy;
        std::vector<int> CaloJetVec_CaloConstituentsVec_Index;
        std::vector<int> CaloJetVec_CaloConstituentsVec_Ieta;
        std::vector<int> CaloJetVec_CaloConstituentsVec_Iphi;
        std::vector<math::XYZTLorentzVector> CaloJetVec_CaloConstituentsVec_p4;
        std::vector<float> CaloJetVec_CaloConstituentsVec_EmEnergy;
        std::vector<float> CaloJetVec_CaloConstituentsVec_HadEnergy;
        std::vector<float> CaloJetVec_CaloConstituentsVec_HBEnergy;
        std::vector<float> CaloJetVec_CaloConstituentsVec_HEEnergy;
        std::vector<float> CaloJetVec_CaloConstituentsVec_HFEnergy;
        std::vector<float> CaloJetVec_CaloConstituentsVec_HOEnergy;
        std::vector<int> CaloJetVec_CaloConstituentsVec_HCALChannelVec_Index;
        std::vector<int> CaloJetVec_CaloConstituentsVec_HCALChannelVec_Ieta;
        std::vector<int> CaloJetVec_CaloConstituentsVec_HCALChannelVec_Iphi;
        std::vector<int> CaloJetVec_CaloConstituentsVec_HCALChannelVec_Depth;
        std::vector<float> CaloJetVec_CaloConstituentsVec_HCALChannelVec_Energy;
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
HCAL_Jet_Ana::HCAL_Jet_Ana(const edm::ParameterSet& iConfig)
{
    //now do what ever initialization is needed
    HBHERecHitsToken_ = consumes<HBHERecHitCollection>(edm::InputTag("hbhereco"));
    CaloJetToken_ = consumes<std::vector<reco::CaloJet>>(edm::InputTag("ak4CaloJets"));
    VertexToken_ = consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"));
    GenJetToken_ = consumes<std::vector<reco::GenJet>>(edm::InputTag("ak4GenJetsNoNu"));

    edm::Service<TFileService> TFS;
    PU_h = TFS->make<TH1F>("PU_h", "PU_h", 100, 0.0, 100.0);

    GenJetTree = TFS->make<TTree>("GenJetTree", "GenJetTree");
    GenJetTree->Branch("GenJetVec_p4", &GenJetVec_p4);

    CaloJetTree = TFS->make<TTree>("CaloJetTree", "CaloJetTree");
    CaloJetTree->Branch("CaloJetVec_Energy", &CaloJetVec_Energy);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_Index", &CaloJetVec_CaloConstituentsVec_Index);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_Ieta", &CaloJetVec_CaloConstituentsVec_Ieta);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_Iphi", &CaloJetVec_CaloConstituentsVec_Iphi);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_p4", &CaloJetVec_CaloConstituentsVec_p4);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_EmEnergy", &CaloJetVec_CaloConstituentsVec_EmEnergy);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_HadEnergy", &CaloJetVec_CaloConstituentsVec_HadEnergy);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_HBEnergy", &CaloJetVec_CaloConstituentsVec_HBEnergy);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_HEEnergy", &CaloJetVec_CaloConstituentsVec_HEEnergy);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_HFEnergy", &CaloJetVec_CaloConstituentsVec_HFEnergy);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_HOEnergy", &CaloJetVec_CaloConstituentsVec_HOEnergy);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_HCALChannelVec_Index", &CaloJetVec_CaloConstituentsVec_HCALChannelVec_Index);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_HCALChannelVec_Ieta", &CaloJetVec_CaloConstituentsVec_HCALChannelVec_Ieta);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_HCALChannelVec_Iphi", &CaloJetVec_CaloConstituentsVec_HCALChannelVec_Iphi);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_HCALChannelVec_Depth", &CaloJetVec_CaloConstituentsVec_HCALChannelVec_Depth);
    CaloJetTree->Branch("CaloJetVec_CaloConstituentsVec_HCALChannelVec_Energy", &CaloJetVec_CaloConstituentsVec_HCALChannelVec_Energy);
}

HCAL_Jet_Ana::~HCAL_Jet_Ana()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------

void HCAL_Jet_Ana::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<std::vector<reco::Vertex>> VertexHandle;
    iEvent.getByToken(VertexToken_, VertexHandle);
    auto Vertex = VertexHandle.product();
    int nPU = Vertex->size();
    PU_h->Fill(nPU);

    edm::Handle<HBHERecHitCollection> HBHERecHitsHandle;
    iEvent.getByToken(HBHERecHitsToken_, HBHERecHitsHandle);
    auto HBHERecHits = HBHERecHitsHandle.product();

    edm::Handle<std::vector<reco::GenJet>> GenJetHandle;
    iEvent.getByToken(GenJetToken_, GenJetHandle);
    auto GenJets = GenJetHandle.product();

    edm::Handle<std::vector<reco::CaloJet>> CaloJetHandle;
    iEvent.getByToken(CaloJetToken_, CaloJetHandle);
    auto CaloJets = CaloJetHandle.product();

    GenJetVec_p4.clear();

    CaloJetVec_Energy.clear();
    CaloJetVec_CaloConstituentsVec_Index.clear();
    CaloJetVec_CaloConstituentsVec_Ieta.clear();
    CaloJetVec_CaloConstituentsVec_Iphi.clear();
    CaloJetVec_CaloConstituentsVec_p4.clear();
    CaloJetVec_CaloConstituentsVec_EmEnergy.clear();
    CaloJetVec_CaloConstituentsVec_HadEnergy.clear();
    CaloJetVec_CaloConstituentsVec_HBEnergy.clear();
    CaloJetVec_CaloConstituentsVec_HEEnergy.clear();
    CaloJetVec_CaloConstituentsVec_HFEnergy.clear();
    CaloJetVec_CaloConstituentsVec_HOEnergy.clear();
    CaloJetVec_CaloConstituentsVec_HCALChannelVec_Index.clear();
    CaloJetVec_CaloConstituentsVec_HCALChannelVec_Ieta.clear();
    CaloJetVec_CaloConstituentsVec_HCALChannelVec_Iphi.clear();
    CaloJetVec_CaloConstituentsVec_HCALChannelVec_Depth.clear();
    CaloJetVec_CaloConstituentsVec_HCALChannelVec_Energy.clear();

    std::vector<reco::CaloJet> SelectedCaloJets = select_CaloJets(GenJets, CaloJets);
    if (SelectedCaloJets.size() == 2)
    {
        GenJetVec_p4.push_back(GenJets->at(0).p4());
        GenJetVec_p4.push_back(GenJets->at(1).p4());

        int CaloJetCounter = 0;
        int CaloConsCounter = 0;

        for (auto SelectedCaloJet : SelectedCaloJets)
        {
            CaloJetVec_Energy.push_back(SelectedCaloJet.energy());
            std::vector<CaloTowerPtr> CaloConsPtrs = SelectedCaloJet.getCaloConstituents();
            for (auto CaloConsPtr : CaloConsPtrs)
            {
                auto CaloConsIeta = CaloConsPtr->ieta();
                auto CaloConsIphi = CaloConsPtr->iphi();
                CaloJetVec_CaloConstituentsVec_Index.push_back(CaloJetCounter);
                CaloJetVec_CaloConstituentsVec_Ieta.push_back(CaloConsIeta);
                CaloJetVec_CaloConstituentsVec_Iphi.push_back(CaloConsIphi);
                CaloJetVec_CaloConstituentsVec_p4.push_back(CaloConsPtr->p4());
                CaloJetVec_CaloConstituentsVec_EmEnergy.push_back(CaloConsPtr->emEnergy());
                CaloJetVec_CaloConstituentsVec_HadEnergy.push_back(CaloConsPtr->hadEnergy());
                CaloJetVec_CaloConstituentsVec_HBEnergy.push_back(CaloConsPtr->energyInHB());
                CaloJetVec_CaloConstituentsVec_HEEnergy.push_back(CaloConsPtr->energyInHE());
                CaloJetVec_CaloConstituentsVec_HFEnergy.push_back(CaloConsPtr->energyInHF());
                CaloJetVec_CaloConstituentsVec_HOEnergy.push_back(CaloConsPtr->energyInHO());

                for(auto HBHERecHit : *HBHERecHits)
                {
                    auto Hid = HBHERecHit.id();
                    auto Depth = Hid.depth();
                    auto Ieta = Hid.ieta();
                    auto Iphi = Hid.iphi();
                    auto Energy = HBHERecHit.energy();

                    if(Ieta == CaloConsIeta && Iphi == CaloConsIphi)
                    {
                        CaloJetVec_CaloConstituentsVec_HCALChannelVec_Index.push_back(CaloConsCounter);
                        CaloJetVec_CaloConstituentsVec_HCALChannelVec_Ieta.push_back(Ieta);
                        CaloJetVec_CaloConstituentsVec_HCALChannelVec_Iphi.push_back(Iphi);
                        CaloJetVec_CaloConstituentsVec_HCALChannelVec_Depth.push_back(Depth);
                        CaloJetVec_CaloConstituentsVec_HCALChannelVec_Energy.push_back(Energy);
                    }//match SimHits to CaloTower
                }//loop SimHits
                CaloConsCounter++;
            }//loop ColoTower
            CaloJetCounter++;
        }//loop CaloJet

        GenJetTree->Fill();
        CaloJetTree->Fill();
    }
}


// ------------ method called once each job just before starting event loop  ------------
void HCAL_Jet_Ana::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HCAL_Jet_Ana::endJob()
{
}

std::vector<reco::CaloJet> HCAL_Jet_Ana::select_CaloJets(const std::vector<reco::GenJet> * GenJets, const std::vector<reco::CaloJet> * CaloJets)
{
    std::vector<reco::CaloJet> SelectedCaloJets;
    if(GenJets->size() > 1 && CaloJets->size() > 1)
    {
        auto GenJet0P4 = GenJets->at(0).p4();
        auto GenJet1P4 = GenJets->at(1).p4();
        if(fabs(GenJet0P4.Eta()) < 3.0 && fabs(GenJet1P4.Eta()) < 3.0
                && ROOT::Math::VectorUtil::DeltaR(GenJet0P4, GenJet1P4) > 1)
        {
            auto CaloJet0P4 = CaloJets->at(0).p4();
            auto CaloJet1P4 = CaloJets->at(1).p4();

            if(ROOT::Math::VectorUtil::DeltaR(GenJet0P4, CaloJet0P4) < 0.2
                    && ROOT::Math::VectorUtil::DeltaR(GenJet1P4, CaloJet1P4) < 0.2)
            {
                SelectedCaloJets.push_back(CaloJets->at(0));
                SelectedCaloJets.push_back(CaloJets->at(1));
            }
            else if(ROOT::Math::VectorUtil::DeltaR(GenJet0P4, CaloJet1P4) < 0.2
                    && ROOT::Math::VectorUtil::DeltaR(GenJet1P4, CaloJet0P4) < 0.2)
            {
                SelectedCaloJets.push_back(CaloJets->at(1));
                SelectedCaloJets.push_back(CaloJets->at(0));
            }
        }
    }
    return SelectedCaloJets;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HCAL_Jet_Ana::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
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
DEFINE_FWK_MODULE(HCAL_Jet_Ana);
