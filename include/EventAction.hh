#ifndef EventAction_h
#define EventAction_h 1
//#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "G4UserEventAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "FMSlargeHit.hh"
#include "FMSsmallHit.hh"
#include "CellGeo.hh"
#include "StackingAction.hh"
#include "PostSim.hh"
#include "Qt.h"
#include "Trajectory.hh"
#include <stdio.h>
#include "G4QT.h"
#include <vector>
#include <math.h>
#include "TF2.h"
#include "TF1.h" //Technically included in tf2, but good practice

//#include "TMath.h"

class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
  EventAction(FILE* fp0,TFile* p0file,TTree* geotr, Int_t eventNum);
  ~EventAction();

  private:
  Bool_t DidSurviveRadiation(Double_t nDays_, Double_t cosTh_, Double_t depth_, Bool_t isLarge_);
  Int_t CreateTheDamFunc(Double_t dayChangeL_=0.005, Double_t dayChangeS_=0.015, Double_t normEtaL_=2.9, Double_t normEtaS_=3.7);
  Int_t GetEnergyHists(); 

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    FILE* fp;

    TFile* pFile;
    TTree* pgeotr;
    int evnum_;
    int rowg;
    int colg;
    int nstbg;
    float gaing;
    float corrg;

    TH1D *PhotonWavelengthS;
    TH1D *PhotonWavelengthL;
    TH1D *OrigPhotonWavelengthS;
    TH1D *OrigPhotonWavelengthL;
    TH2D *EnergyDepositShapeS;
    TH2D *EnergyDepositShapeL;

//    TH1D *GenPhotonZ, *DetectPhotonZ;
    TH1D *EdepZL, *EdepZS;
    TH1D *cosL, *cosS;

    PostSim pss;

    TTree *testtree;
    Int_t eventNumber;
    

    //BS
    Int_t nsaved_traj;
    UInt_t nhits;
    Int_t parentID[40];
    Int_t trackID[40];
    Float_t mass[40];
    Float_t energy[40];
    Float_t vx[40];
    Float_t vy[40];
    Float_t vz[40];
    Float_t px[40];
    Float_t py[40];
    Float_t pz[40];
    Float_t pt[40];
    Float_t eta[40];
    Int_t nstb[40];
    Int_t row[40];
    Int_t col[40];

    //BS Radiation Diagnostics
    UInt_t nPhotonsGeneratedL, nPhotonsGeneratedS;
    Int_t nPhotonsNormL, nPhotonsRadL,nPhotonsNormS, nPhotonsRadS ; //number of detected photons before and after radaition damage
    Double_t zPhotonsNormL, zPhotonsNormS, zPhotonsRadL, zPhotonsRadS; //depth from photocathode of photon, holds total z to later be divided by number of photons
    Double_t zAvgNormL, zAvgNormS, zAvgRadL, zAvgRadS;
    Double_t lPhotonsNormL, lPhotonsNormS, lPhotonsRadL, lPhotonsRadS; //Total Path Length traveled by Photons, holds total l to later be divided by number of photons
    Double_t lAvgNormL, lAvgNormS, lAvgRadL, lAvgRadS;
    Double_t EdeptotL, EdeptotS; //total deposited energy via energy deposition
    Double_t EdeptotModL, EdeptotModS; //total deposited energy via energy deposition, modified to account for random change in var for each
    Double_t EdeptotQuadL, EdeptotQuadS; //total deposited energy via energy deposition, modified to account for group changes in quadrants
    Double_t EdeptotLeadL, EdeptotLeadS; //total deposited energy via energy deposition, modified to account for group changes in quadrants



    Trajectory *trj;
    G4QT test;
    std::vector<G4QT> myvectest;
    std::vector<std::string> name;
    UInt_t nPrimaries;

//    UInt_t Qtenergy[3700];
//    UInt_t nqtEnergy;
//    UInt_t *pnqtEnergy;
    G4QT energyG4QT;
    std::vector<G4QT> energyQTvec;

    G4QT EnergyModG4QT;
    std::vector<G4QT> EnergyModVec;
//    UInt_t QtEnergyMod[3700];
//    UInt_t nqtEnergyMod;
//    UInt_t *pnqtEnergyMod;

    G4QT RadDamageG4QT;
    std::vector<G4QT> RadDamageVec;


    //BS
    G4QT CellGeo1G4QT;
    std::vector<G4QT> CellGeo1Vec;
    G4QT CellGeo2G4QT;
    std::vector<G4QT> CellGeo2Vec;
    G4QT CellGeo3G4QT;
    std::vector<G4QT> CellGeo3Vec;
    G4QT CellGeo4G4QT;
    std::vector<G4QT> CellGeo4Vec;
    G4QT CellGeo5G4QT;
    std::vector<G4QT> CellGeo5Vec;
    //End BS
    //

    //BS Damage Rearranging
    TF1 *pDamGaus;
    TF2 *pNormFuncL, *pNormFuncS;
    TFile *pEShapeFile;
    TH1D *pParamL[3], *pParamS[3];
//    Int_t CreateNormFunctions(Double_t normEtaL_=2.9, Double_t normEtaS_=3.7);
    Int_t CreateNormFunctions(Double_t normEtaL_=2.9, Double_t normEtaS_=3.5);
//    Bool_t DidSurviveRad(Double_t etaPhoton_, Double_t depthPhoton_, Double_t cosThPhoton_, Bool_t isLarge_, Double_t ratDayChange_, Double_t nDays_);
//    Bool_t DidSurviveRad(Double_t etaPhoton_, Double_t depthPhoton_, Double_t cosThPhoton_, Bool_t isLarge_, Double_t ratDayChange_, Double_t nDays_, Bool_t useWavelength_=false, Double_t wavelength_=465);
    Bool_t DidSurviveRad(Double_t etaPhoton_, Double_t depthPhoton_, Double_t cosThPhoton_, Bool_t isLarge_, Double_t ratDayChange_, Double_t nDays_, Bool_t useWavelength_=false, Double_t wavelength_=465, Bool_t useAvg_=true);

//    Int_t radWhichID[40][100]; //BS NOTE. I origionally tried to make this a longer array and then only choose the first nDaysSimulated, however this fails with the multi-dimensional array branches. Instead, it would grab the first 100 elements from [0][0-99] and then get the others
//
//    THEREFORE MAKE SURE THAT THE SECOND DIMENSION IS THE SAME as nDaysSimulated
//    BS. Might replace this then with TClones or tobject array
    Int_t nPhotonsRadNewL[40][52];
    Int_t nPhotonsRadNewS[40][52];
    Float_t dayNum[40][52];
    Float_t zAvgNewL[40][52];
    Float_t zAvgNewS[40][52];
    Int_t radWhichID[40][52];
    Float_t EdeptotNewL[40], EdeptotNewS[40]; //total deposited energy via energy deposition.
    Float_t EdeptotMrigankaNewL[40], EdeptotMrigankaNewS[40]; //total deposited energy via energy deposition, modified by mriganaka
    Int_t nDaysSimulated;
    Double_t daySpacing;
    //Wavelength Aditions
    TF1 *pFuncWaveDam;
    Int_t nPhotonsWaveNewL[40][52];
    Int_t nPhotonsWaveNewS[40][52];
    Float_t zAvgWaveNewL[40][52];
    Float_t zAvgWaveNewS[40][52];

    Int_t nPhotonsCerenNewL[40][52];
    Int_t nPhotonsCerenNewS[40][52];
    //End BS Damage Rearranging

  private:
	int max;
    G4int collectionIDL;
    G4int collectionIDS;
    G4int verboseLevel;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
