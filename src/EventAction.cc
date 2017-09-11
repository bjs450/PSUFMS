#include "EventAction.hh"
#include "FMSlargeHit.hh"
#include "FMSsmallHit.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "Trajectory.hh"
#include "G4ios.hh"
#include "stdio.h"
#include "TTree.h"
#include "TVector3.h"
#include "PostSim.hh"
#include <TROOT.h> //Sebastian
#include <vector>


#include "TRandom1.h" //BS
#include "TF1.h"
#include "TMath.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(FILE* fp0,TFile* tf0, TTree* geotr, Int_t eventNum)
{
  pFile=tf0;
  eventNumber=eventNum;
  pFile->cd();

  pss.eventnumber=eventNumber;
  trj=new Trajectory();

  evnum_=-1;
  nsaved_traj=0;
  memset(parentID,0,sizeof(parentID));
  memset(trackID,0,sizeof(trackID));
  memset(mass,0,sizeof(mass));
  memset(energy,0,sizeof(energy));
  memset(vx,0,sizeof(vx));
  memset(vy,0,sizeof(vy));
  memset(vz,0,sizeof(vz));
  memset(px,0,sizeof(px));
  memset(py,0,sizeof(py));
  memset(pz,0,sizeof(pz));
  memset(nstb,0,sizeof(nstb));
  memset(row,0,sizeof(row));
  memset(col,0,sizeof(col));

  //  EventAction::PhotonWavelengthL = new TH1D("PhotonWavelengthL", "Photon Wavelength that reach photocathode, Large", 500, 300, 800);
  // EventAction::PhotonWavelengthS = new TH1D("PhotonWavelengthS", "Photon Wavelength that reach photocathode, Small", 500, 300, 800);
  //EventAction::OrigPhotonWavelengthL = new TH1D("OrigPhotonWavelengthL", "Origional Photon Wavelength that reach photocathode, Large", 500, 300, 800);
  //EventAction::OrigPhotonWavelengthS = new TH1D("OrigPhotonWavelengthS", "Origional Photon Wavelength that reach photocathode, Small", 500, 300, 800);
  //  EventAction::EnergyDepositShapeS = new TH2D("EnergyDepositShapeS","Eta vs distance from Photocathode, Small. Energy Weighted",48,-1,46,101,3,4.4);
  // EventAction::EnergyDepositShapeL = new TH2D("EnergyDepositShapeL","Eta vs distance from Photocathode, Large. Energy Weighted",63,-1,61,101,2.5,3.5);

  //EdepZL= new TH1D("EdeptZL", "Energy Deposition", 63, -1, 61);
  //EdepZS= new TH1D("EdeptZS", "Energy Deposition", 48, -1, 46);


  EventAction::testtree= new TTree("testtree","testtree");

  testtree->Branch("event", &(EventAction::eventNumber), "event/I");
  testtree->Branch("nsaved_traj", &(EventAction::nsaved_traj), "nsaved_traj/I");
  testtree->Branch("parentID", EventAction::parentID, "parentID[nsaved_traj]/I");
  testtree->Branch("trackID", EventAction::trackID, "trackID[nsaved_traj]/I");
  testtree->Branch("mass", EventAction::mass, "mass[nsaved_traj]/F");
  testtree->Branch("energy",EventAction::energy,"energy[nsaved_traj]/F");
  testtree->Branch("vx", EventAction::vx, "vx[nsaved_traj]/F");
  testtree->Branch("vy", EventAction::vy, "vy[nsaved_traj]/F");
  testtree->Branch("vz", EventAction::vz, "vz[nsaved_traj]/F");
  testtree->Branch("px", EventAction::px, "px[nsaved_traj]/F");
  testtree->Branch("py", EventAction::py, "py[nsaved_traj]/F");
  testtree->Branch("pz", EventAction::pz, "pz[nsaved_traj]/F");
  testtree->Branch("pt",EventAction::pt, "pt[nsaved_traj]/F");
  testtree->Branch("eta",EventAction::eta, "eta[nsaved_traj]/F");
  testtree->Branch("name",&name);
  testtree->Branch("nPrimaries",&nPrimaries);
  testtree->Branch("nstb",EventAction::nstb, "nstb[nsaved_traj]/I");
  testtree->Branch("row",EventAction::row, "row[nsaved_traj]/I");
  testtree->Branch("col",EventAction::col, "col[nsaved_traj]/I");

  testtree->Branch("myvectest.",&myvectest);
  testtree->Branch("energyQT.",&energyQTvec);
  testtree->Branch("EnergyModVec.",&EnergyModVec);
  testtree->Branch("RadDamageVec.",&RadDamageVec);

  testtree->Branch("All1Gain.",&CellGeo1Vec);
  //  testtree->Branch("10PercentVar",&CellGeo2Vec);
  testtree->Branch("Day142Var.",&CellGeo2Vec);
  //  testtree->Branch("20PercentVar",&CellGeo3Vec);
  //  testtree->Branch("30PercentVar",&CellGeo4Vec);

  testtree->Branch("D142Plus10.",&CellGeo3Vec);
  testtree->Branch("D142PlusGaus.",&CellGeo4Vec);

  /*
  testtree->Branch("myvectest",&myvectest);
  testtree->Branch("energyQT",&energyQTvec);
  testtree->Branch("EnergyModVec",&EnergyModVec);
  testtree->Branch("RadDamageVec",&RadDamageVec);

  testtree->Branch("All1Gain",&CellGeo1Vec);
  //  testtree->Branch("10PercentVar",&CellGeo2Vec);
  testtree->Branch("Day142Var",&CellGeo2Vec);
  //  testtree->Branch("20PercentVar",&CellGeo3Vec);
  //  testtree->Branch("30PercentVar",&CellGeo4Vec);

  testtree->Branch("D142Plus10",&CellGeo3Vec);
  testtree->Branch("D142PlusGaus",&CellGeo4Vec);
  */

  testtree->Branch("TriggersAll1",pss.Trig1,"TrigersAll1[9]/I");
  //  testtree->Branch("TriggersVar10",pss.Trig2,"TriggersVar10[9]/I");
  testtree->Branch("TriggersVar142",pss.Trig2,"TriggersVar142[9]/I");

  testtree->Branch("TriggersVar142Plus10",pss.Trig3,"TriggersVar142Plus10[9]/I");
  testtree->Branch("TriggersVar142Gaus",pss.Trig4,"TriggersVar142Gaus[9]/I");

  testtree->Branch("TriggersVar142Lead",pss.Trig5,"TriggersVar142Lead[9]/I");
  //  testtree->Branch("TriggersVar080",pss.Trig4,"TriggersVar080[9]/I");
  //  testtree->Branch("TriggersVar20",pss.Trig3,"TriggersVar20[9]/I");
  //  testtree->Branch("TriggersVar30",pss.Trig4,"TriggersVar30[9]/I");
  // testtree->Branch("TriggersVar956",pss.Trig4,"TriggersVar956[9]/I");

  //BS testing Radiation Diagonstics
  testtree->Branch("nPhotonsNormL", &(EventAction::nPhotonsNormL) );
  testtree->Branch("nPhotonsRadL", &(EventAction::nPhotonsRadL) );
  testtree->Branch("nPhotonsNormS", &(EventAction::nPhotonsNormS) );
  testtree->Branch("nPhotonsRadS", &(EventAction::nPhotonsRadS) );
  //  testtree->Branch("nPhotonsGeneratedL", &(EventAction::nPhotonsGeneratedL) );
  //  testtree->Branch("nPhotonsGeneratedS", &(EventAction::nPhotonsGeneratedS) );
  testtree->Branch("nPhotonsRadL", &(EventAction::nPhotonsRadL) );
  /*
     testtree->Branch("zAvgNormL", &(EventAction::zAvgNormL) );
     testtree->Branch("zAvgRadL", &(EventAction::zAvgRadL) );
     testtree->Branch("zAvgNormS", &(EventAction::zAvgNormS) );
     testtree->Branch("zAvgRadS", &(EventAction::zAvgRadS) );
     testtree->Branch("lAvgNormL", &(EventAction::lAvgNormL) );
     testtree->Branch("lAvgRadL", &(EventAction::lAvgRadL) );
     testtree->Branch("lAvgNormS", &(EventAction::lAvgNormS) );
     testtree->Branch("lAvgRadS", &(EventAction::lAvgRadS) );
     testtree->Branch("zAvgNormL", &(EventAction::zAvgNormL) );
     testtree->Branch("zAvgRadL", &(EventAction::zAvgRadL) );
     testtree->Branch("zAvgNormS", &(EventAction::zAvgNormS) );
     testtree->Branch("zAvgRadS", &(EventAction::zAvgRadS) );
     testtree->Branch("lAvgNormL", &(EventAction::lAvgNormL) );
     testtree->Branch("lAvgRadL", &(EventAction::lAvgRadL) );
     testtree->Branch("lAvgNormS", &(EventAction::lAvgNormS) );
     testtree->Branch("lAvgRadS", &(EventAction::lAvgRadS) );
     */
  testtree->Branch("EdeptotL", &(EventAction::EdeptotL) );
  testtree->Branch("EdeptotS", &(EventAction::EdeptotS) );
  testtree->Branch("EdeptotModL", &(EventAction::EdeptotModL) );
  testtree->Branch("EdeptotModS", &(EventAction::EdeptotModS) );
  testtree->Branch("EdeptotGausL", &(EventAction::EdeptotQuadL) );
  testtree->Branch("EdeptotGausS", &(EventAction::EdeptotQuadS) );
  testtree->Branch("EdeptotLeadL", &(EventAction::EdeptotLeadL) );
  testtree->Branch("EdeptotLeadS", &(EventAction::EdeptotLeadS) );

  //BS Damage Additions
  pDamGaus=0;
  pNormFuncL=0; pNormFuncS=0;
  pEShapeFile = new TFile("HistFiles.root","READ");
  for(Int_t dumInit=0; dumInit<3; dumInit++) {pParamL[dumInit]=pParamS[dumInit]=0;}
  Int_t didSetupDamage=CreateNormFunctions(2.9,3.7);
  if(didSetupDamage) G4cout<<"ERROR, Could not Create Functions!!"<<G4endl;
  pFile->cd();

  Double_t daySpacing=7;
  nDaysSimulated=52; //Number of Days with entries, not necessarily x consecutive days;
  TString dumString; dumString.Form("[%i]",nDaysSimulated);

  testtree->Branch("nPhotonsRadNewL",nPhotonsRadNewL,"nPhotonsRadNewL[nsaved_traj]"+dumString+"/I");
  testtree->Branch("nPhotonsRadNewS",nPhotonsRadNewS,"nPhotonsRadNewS[nsaved_traj]"+dumString+"/I");
  testtree->Branch("radWhichID",radWhichID,"radWhichID[nsaved_traj]"+dumString+"/I");
  testtree->Branch("dayNum",dayNum,"dayNum[nsaved_traj]"+dumString+"/F");
  testtree->Branch("EdeptotNewL",EdeptotNewL,"EdeptotNewL[nsaved_traj]/F");
  testtree->Branch("EdeptotNewS",EdeptotNewS,"EdeptotNewS[nsaved_traj]/F");
  testtree->Branch("EdeptotMrigankaNewL",EdeptotMrigankaNewL,"EdeptotMrigankaNewL[nsaved_traj]/F");
  testtree->Branch("EdeptotMrigankaNewS",EdeptotMrigankaNewS,"EdeptotMrigankaNewS[nsaved_traj]/F");

  testtree->Branch("zAvgNewL",zAvgNewL,"zAvgNewL[nsaved_traj]"+dumString+"/F");
  testtree->Branch("zAvgNewS",zAvgNewS,"zAvgNewS[nsaved_traj]"+dumString+"/F");

  testtree->Branch("nPhotonsWaveNewL",nPhotonsWaveNewL,"nPhotonsWaveNewL[nsaved_traj]"+dumString+"/I");
  testtree->Branch("nPhotonsWaveNewS",nPhotonsWaveNewS,"nPhotonsWaveNewS[nsaved_traj]"+dumString+"/I");
  testtree->Branch("zAvgWaveNewL",zAvgWaveNewL,"zAvgWaveNewL[nsaved_traj]"+dumString+"/F");
  testtree->Branch("zAvgWaveNewS",zAvgWaveNewS,"zAvgWaveNewS[nsaved_traj]"+dumString+"/F");

  testtree->Branch("nPhotonsCerenNewL",nPhotonsCerenNewL,"nPhotonsCerenNewL[nsaved_traj]"+dumString+"/I");
  testtree->Branch("nPhotonsCerenNewS",nPhotonsCerenNewS,"nPhotonsCerenNewS[nsaved_traj]"+dumString+"/I");

  for(Int_t i=0; i<40;i++)
  {
    for(Int_t j=0; j<nDaysSimulated; j++)
    {
      radWhichID[i][j]=-1; //BS, Making clear to seperate out ones with no contribution to a 0
      nPhotonsRadNewL[i][j]=0;
      nPhotonsRadNewS[i][j]=0;
      dayNum[i][j]=(Float_t) j*daySpacing;
      zAvgNewL[i][j]= 0;
      zAvgNewS[i][j]= 0;

      nPhotonsWaveNewL[i][j]=0;
      nPhotonsWaveNewS[i][j]=0;
      zAvgWaveNewL[i][j]= 0;
      zAvgWaveNewS[i][j]= 0;

      nPhotonsCerenNewL[i][j]=0;
      nPhotonsCerenNewS[i][j]=0;
    }
  }

  pFuncWaveDam=new TF1("pFuncWaveDam","pol3",300,900);
  Double_t waveParams[4]={4.08392,-0.0122715,1.50105E-5,-6.1995E-9};
  pFuncWaveDam->SetParameters(waveParams);
  //End BS Damage Additions
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{
  G4cout << ">>> begin event "<< G4endl;
  pss.dummyFirstInEvent=0;

  //BS, Radiation Damage Diagonsistics
  nPhotonsGeneratedL=0; nPhotonsGeneratedS=0;
  nPhotonsNormL=0; nPhotonsNormS=0; nPhotonsRadL=0; nPhotonsRadS=0;
  zPhotonsNormL=0; zPhotonsNormS=0; zPhotonsRadL=0; zPhotonsRadS=0;
  lPhotonsNormL=0; lPhotonsNormS=0; lPhotonsRadL=0; lPhotonsRadS=0;
  zAvgNormL=0; zAvgNormS=0; zAvgRadL=0; zAvgRadS=0;
  lAvgNormL=0; lAvgNormS=0; lAvgRadL=0; lAvgRadS=0;
  EdeptotL=0; EdeptotS=0, EdeptotModL=0, EdeptotModS=0, EdeptotQuadL=0, EdeptotQuadS=0,EdeptotLeadL=0, EdeptotLeadS=0;

  //Filling the random array before the event
  srand(time(NULL));
  TRandom1 myrand(time(NULL));

  Int_t whichFactor=0, TDC=0;
  Float_t modFactor[4]={1.055, 1.0, 0.979, 0.85};

  Bool_t dummyFirst[4]={0};
  Int_t TDCFirst[4]={0};

  for(Int_t instb=1; instb<5; instb++)
  {
    for(Int_t irow=0; irow<34; irow++)
    {
      for(Int_t icol=0; icol<17; icol++)
      {
	Double_t varmax=0.10; //Maximum 10% variation on gain
	pss.gainRandVar[instb-1][irow][icol]= 1.0+ ( (rand()%2000-1000) * 0.001*varmax);

	//Correlated TDCS
	////BS, optimise this later
	if(instb==1)
	{
	  if(irow<17)
	  {
	    TDC=(Int_t) (myrand.Gaus(7.709,0.6302) + 0.5);
	    if (!(dummyFirst[0])) { dummyFirst[0]=true; TDCFirst[0]=TDC;}
	  }
	  else
	  {
	    TDC=(Int_t) (myrand.Gaus(7.339,0.6974) + 0.5);
	    if (!(dummyFirst[1])) { dummyFirst[1]=true; TDCFirst[1]=TDC;}
	  }
	}
	else if(instb==2)
	{
	  if(irow<17)
	  {
	    TDC=(Int_t) (myrand.Gaus(7.818,0.6571) + 0.5);
	    if (!(dummyFirst[2])) { dummyFirst[2]=true; TDCFirst[2]=TDC;}
	  }
	  else
	  {
	    TDC=(Int_t) (myrand.Gaus(6.848,0.7212) + 0.5);
	    if (!(dummyFirst[3])) { dummyFirst[3]=true; TDCFirst[3]=TDC;}
	  }
	}
	else if(instb==3)
	{
	  if(irow<12)
	  {
	    TDC=(Int_t) (myrand.Gaus(7.709,0.6302) + 0.5);
	  }
	  else
	  {
	    TDC=(Int_t) ( myrand.Gaus(7.339,0.6974) + 0.5);
	  }
	}
	else
	{
	  if(irow<12)
	  {
	    TDC=(Int_t) (myrand.Gaus(7.818,0.6571) +0.5);
	  }
	  else
	  {
	    TDC=(Int_t) (myrand.Gaus(6.848,0.7212) +0.5);
	  }
	}

	if(TDC<6) TDC=6;
	else if(TDC>9) TDC=9;
	pss.gainGroupVar[instb-1][irow][icol]=modFactor[TDC-6];
      }
    }
  }

  for(Int_t dumint=0; dumint<4; dumint++)
  {
    if(TDCFirst[dumint]<6) TDCFirst[dumint]=6;
    if(TDCFirst[dumint]>9) TDCFirst[dumint]=9;
    pss.leadClusterGain[dumint]=modFactor[TDC-6];
  }

  //BS Damage Addditions
  for(Int_t i=0; i<40; i++)
  {
    EdeptotNewL[i]=0; EdeptotNewS[i]=0; EdeptotMrigankaNewL[i]=0; EdeptotMrigankaNewS[i]=0;
    for(Int_t j=0; j<nDaysSimulated; j++)
    {
      radWhichID[i][j]=-1; //BS, Making clear to seperate out ones with no contribution to a 0
      nPhotonsRadNewL[i][j]=0;
      nPhotonsRadNewS[i][j]=0;
      zAvgNewL[i][j]= 0;
      zAvgNewS[i][j]= 0;

      nPhotonsWaveNewL[i][j]=0;
      nPhotonsWaveNewS[i][j]=0;
      zAvgWaveNewL[i][j]= 0;
      zAvgWaveNewS[i][j]= 0;

      nPhotonsCerenNewL[i][j]=0;
      nPhotonsCerenNewS[i][j]=0;
    }
  }
  //End BS Damage Addtions


  //BS I still worry there's a memory leak on the individual TMatrices, but this should stop for now the failure for multi-event
  pss.MatrixArray.SetOwner();
  pss.MatrixArray.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4cout << "<<< EndOfEventAction"<< G4endl;

  Int_t dummycount=0;
  pss.nprimaries = (int) evt->GetNumberOfPrimaryVertex();
  pFile->cd();

  Double_t RadParameters[3]={0.1146, 1.73, 43.4};
  TF1 *DamageFunc = new TF1("DamageFunc","gaus");
  DamageFunc->SetParameters(RadParameters);

  //BS fill in here
  //  Double_t WaveParameters[4];


  //More Detailed Radiation Damage Treatment
  TFile *EnergyShapeFile = new TFile("HistFiles.root","READ");

  TH1D *Par0L=(TH1D*) EnergyShapeFile->Get("EnergyDepositShapeL_0");
  TH1D *Par1L=(TH1D*) EnergyShapeFile->Get("EnergyDepositShapeL_1");
  TH1D *Par2L=(TH1D*) EnergyShapeFile->Get("EnergyDepositShapeL_2");

  TH1D *Par0S=(TH1D*) EnergyShapeFile->Get("EnergyDepositShapeS_0");
  TH1D *Par1S=(TH1D*) EnergyShapeFile->Get("EnergyDepositShapeS_1");
  TH1D *Par2S=(TH1D*) EnergyShapeFile->Get("EnergyDepositShapeS_2");

  Double_t RatioChangeL=0.005; //Change in percentage that occurs in Large Cell Signal per day. Default 0.5% = 0.005
  Double_t RatioChangeS=0.015; //Change in percentage that occurs in Small Cell Signal. Default 1.5% = 0.015

  //  Double_t CleanMFPL=115.;
  //  Double_t CleanMFPS=23.;
  Double_t Days=56;


  //Change to adjust for wavelenths later
  //
  Int_t NormalizingBinL=Par0L->FindBin(2.9);
  Int_t NormalizingBinS=Par0S->FindBin(3.7);

  Double_t RadParametersUpdatedL[3] = {(*Par0L)[NormalizingBinL],(*Par1L)[NormalizingBinL],(*Par2L)[NormalizingBinL]};
  Double_t RadParametersUpdatedS[3] = {(*Par0S)[NormalizingBinS],(*Par1S)[NormalizingBinS],(*Par2S)[NormalizingBinS]};

  TF1 *DamageFuncUpdated = new TF1("DamageFuncUpdated","gaus");

  DamageFuncUpdated->SetParameters(RadParametersUpdatedS);

  Double_t NormFactorS= -log( pow(1.-RatioChangeS,Days) ) / DamageFuncUpdated->Integral(0,45.);

  DamageFuncUpdated->SetParameters(RadParametersUpdatedL);

  Double_t NormFactorL= -log( pow(1.-RatioChangeL,Days) ) / DamageFuncUpdated->Integral(0,60.);


  evnum_++;
  nsaved_traj=0;
  nhits=0;
  G4QT test;

  // get number of stored trajectories
  //
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  for(G4int iv=0;iv<n_trajectories;iv++)
  {
    trj=(Trajectory*)((*(evt->GetTrajectoryContainer()))[iv]);
    Int_t parentID_ = (Int_t) trj->GetParentID();
    //BS. the primary particles are added in opposite order of what you would think. Also, added a storing cut on the decay mode in the PostUserTracking Action.
    max=40;
    if(nsaved_traj >= max )
    {
      G4cout<<"MAXIMUM NUMBER OF TRAJECTORIES PER EVENT REACHED, STOPPING EVENT"<<G4endl;
      break;
    }
    else
    {
      EventAction::parentID[nsaved_traj] = trj->GetParentID();
      EventAction::trackID[nsaved_traj] = trj->GetTrackID();
      EventAction::mass[nsaved_traj] = trj->GetMass()/GeV;
      EventAction::energy[nsaved_traj] = trj->GetEnergy()/GeV;
      EventAction::vx[nsaved_traj] = (Float_t) ( trj->GetVertexPosition() ).getX()/cm;
      EventAction::vy[nsaved_traj] = (Float_t) ( trj->GetVertexPosition() ).getY()/cm;
      EventAction::vz[nsaved_traj] = (Float_t) ( trj->GetVertexPosition() ).getZ()/cm;

      EventAction::px[nsaved_traj] = (Float_t) ( trj->GetInitialMomentum() ).getX()/GeV;
      EventAction::py[nsaved_traj] = (Float_t) ( trj->GetInitialMomentum() ).getY()/GeV;
      EventAction::pz[nsaved_traj] = (Float_t) ( trj->GetInitialMomentum() ).getZ()/GeV;
      EventAction::pt[nsaved_traj]=sqrt(pow(px[nsaved_traj],2)+pow(py[nsaved_traj],2));
      EventAction::eta[nsaved_traj]=-log(tan(.5*(acos( (pz[nsaved_traj]) /( sqrt(pow(px[nsaved_traj],2)+pow(py[nsaved_traj],2) + pow(pz[nsaved_traj],2) ) ) ) ) ));

      EventAction::nstb[nsaved_traj] =  trj->GetNstb();
      EventAction::row[nsaved_traj] =  trj->GetRow();
      EventAction::col[nsaved_traj] =  trj->GetCol();

      name.push_back( (std::string) trj->GetParticleName() );
      EventAction::nPrimaries=evt->GetNumberOfPrimaryVertex();

      //BS Damage Additions
      for(Int_t dumLoop=0; dumLoop<nDaysSimulated; dumLoop++)
      {
	EventAction::radWhichID[nsaved_traj][dumLoop]=trj->GetTrackID();
      }

      //
      EventAction::nsaved_traj++;
    };
  };


  pss.FillMatrixArray(nsaved_traj,EventAction::trackID);

  /////////////////////////////////////////////////

  G4String colName;
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  collectionIDL = SDman->GetCollectionID(colName="FMSlargeColl");
  collectionIDS = SDman->GetCollectionID(colName="FMSsmallColl");

  G4HCofThisEvent* HCofEvent = evt->GetHCofThisEvent();
  FMSlargeHitsCollection* myCollectionL = 0;
  FMSsmallHitsCollection* myCollectionS = 0;
  if(HCofEvent)
  {
    myCollectionL = (FMSlargeHitsCollection*)(HCofEvent->GetHC(collectionIDL));
    myCollectionS = (FMSsmallHitsCollection*)(HCofEvent->GetHC(collectionIDS));
  };

  G4int nCer[1000];
  G4double Esum,nCerenkov;
  G4int dummy = nsaved_traj; 



  /* large cell action */
  if(myCollectionL)
  {
    nhits=myCollectionL->entries();
    G4cout<<"Try to Print " << nhits <<" Large Cell Hits"<<G4endl;
    Esum=0; 
    nCerenkov=0;
    for(UInt_t jj=0;jj<nhits;jj++)
    {
      pss.zdepth=30.- ( (*myCollectionL)[jj]->GetPos().getZ() /cm ); //distance to end of cell

      pss.isHit=0;
      if((*myCollectionL)[jj]->IsHit())pss.isHit=1;
      pss.isOpt=0;
      //BS
      Bool_t DidSurvive=true;
      pss.isRealHit=false;
      if((*myCollectionL)[jj]->IsOptical())
      {
	pss.isOpt=1;     

	nPhotonsGeneratedL++;

	Double_t kinE= (Double_t) (*myCollectionL)[jj]->GetKinE()/GeV;
	Double_t Wave= TMath::HC()/(TMath::Qe() *kinE); //in nm

	//	OrigPhotonWavelengthL->Fill(Wave);

	if(pss.isHit)
	{
	  TVector3 posvec( (*myCollectionL)[jj]->GetPos().getX() /cm ,(*myCollectionL)[jj]->GetPos().getY() /cm , 700.);
	  Double_t HitEta= posvec.Eta();

	  Double_t CosTh= (*myCollectionL)[jj]->GetCosTh();
	  Double_t wavelength;

	  Int_t CurrentBin = Par0L->FindBin(HitEta);

	  RadParametersUpdatedL[0]= NormFactorL * (*Par0L)[CurrentBin];
	  RadParametersUpdatedL[1]=(*Par1L)[CurrentBin];
	  RadParametersUpdatedL[2]=(*Par2L)[CurrentBin];

	  //BS Damage Additions
	  for(Int_t dumAncestLoop=0;dumAncestLoop<((*myCollectionL)[jj]->GetAncestors()).size();dumAncestLoop++)
	  {
	    for(Int_t dumNSavedLoop=0; dumNSavedLoop<nsaved_traj; dumNSavedLoop++)
	    {
	      if(((*myCollectionL)[jj]->GetAncestors()).at(dumAncestLoop)==EventAction::trackID[dumNSavedLoop])
	      {
		for(Int_t dumNDays=0; dumNDays<nDaysSimulated; dumNDays++)
		{

		  Double_t nDaysActual= dayNum[dumNSavedLoop][dumNDays];
		  if( EventAction::DidSurviveRad( HitEta, pss.zdepth, CosTh, true, 0.005, nDaysActual) ) 
		  { 
		    nPhotonsRadNewL[dumNSavedLoop][dumNDays]+=1;
		    zAvgNewL[dumNSavedLoop][dumNDays]+= (Float_t) pss.zdepth;

		    if( (*myCollectionL)[jj]->GetIsCerenkov() ) nPhotonsCerenNewL[dumNSavedLoop][dumNDays]+=1;
		  }
		  if( EventAction::DidSurviveRad( HitEta, pss.zdepth, CosTh, true, 0.005, nDaysActual, true, Wave) ) 
		  { 
		    nPhotonsWaveNewL[dumNSavedLoop][dumNDays]+=1;
		    zAvgWaveNewL[dumNSavedLoop][dumNDays]+= (Float_t) pss.zdepth;
		  }
		}//dumNDays
	      }//if ancestor=traj
	    }//Saved Trajectories Loop
	  } //Ancestor Loop
	  //END BS Damage Addition


	  DamageFuncUpdated->SetParameters(RadParametersUpdatedL);

	  Double_t ProbabilitySurvive=100* TMath::Exp(-1.* DamageFunc->Integral(60-pss.zdepth,60) );

	  Double_t ProbabilitySurviveUpdated=0;
	  if( CosTh >0)
	  {
	    ProbabilitySurviveUpdated=100*TMath::Exp(-1./fabs(CosTh) * DamageFuncUpdated->Integral(0,pss.zdepth) );
	  }
	  else if(CosTh<0)
	  {
	    ProbabilitySurviveUpdated=100*TMath::Exp(-1./fabs(CosTh) *( DamageFuncUpdated->Integral(pss.zdepth,60.) + DamageFuncUpdated->Integral(0,60.) ) );
	  }
	  else
	  {
	    ProbabilitySurviveUpdated = 0;
	  }

	  if( rand()%100 < ProbabilitySurviveUpdated) pss.isRealHit=true;

	  nPhotonsNormL+=1;
	  zPhotonsNormL+=pss.zdepth;
	  lPhotonsNormL+= (*myCollectionL)[jj]->GetTotalTrackLength();
	  if(pss.isRealHit)
	  {
	    nPhotonsRadL+=1;
	    zPhotonsRadL+=pss.zdepth;
	    lPhotonsRadL+= (*myCollectionL)[jj]->GetTotalTrackLength();
	  }


	  //	  std::cout<<"Large Cell Depth ="<<pss.zdepth<<" and extra Probability to survive ="<<ProbabilitySurviveUpdated<<std::endl;


	  //	PhotonWavelengthL->Fill(Wave);
	}
      }
      //
      int jc=(*myCollectionL)[jj]->GetCellNb();
      pss.nstb = cellgeo->GetNstb(jc,1);
      pss.row = cellgeo->GetRow(jc,1);
      pss.col = cellgeo->GetCol(jc,1);
      pss.gain = cellgeo->GetGain(pss.nstb,pss.row,pss.col);
      pss.corr = cellgeo->GetCorr(pss.nstb,pss.row,pss.col);
      pss.en = (*myCollectionL)[jj]->GetEdep()/GeV;
      pss.tot_en = (*myCollectionL)[jj]->GetTdep()/GeV;
      pss.kin = (*myCollectionL)[jj]->GetKinE()/GeV;
      if(pss.isOpt==1) nCer[jc]+=1;
      Esum+=(*myCollectionL)[jj]->GetEdep()/GeV;
      if(pss.isOpt==1)nCerenkov+=1;
      pss.verx = (*myCollectionL)[jj]->GetPos().x()/cm;
      pss.very = (*myCollectionL)[jj]->GetPos().y()/cm;
      pss.verz = (*myCollectionL)[jj]->GetPos().z()/cm;
      pss.px = (*myCollectionL)[jj]->GetMom().x();
      pss.py = (*myCollectionL)[jj]->GetMom().y();
      pss.pz = (*myCollectionL)[jj]->GetMom().z();
      sprintf(pss.name,"%s",(*myCollectionL)[jj]->GetName().data());

      pss.trNum = (*myCollectionL)[jj]->GetTrackNo();
      pss.cosTh = (*myCollectionL)[jj]->GetCosTh();
      pss.NCer = (*myCollectionL)[jj]->GetnCerenkov();
      pss.gtime = (*myCollectionL)[jj]->GetGlobalTime()/ns;
      pss.ltime = (*myCollectionL)[jj]->GetLocalTime()/ns;
      pss.originalID = (*myCollectionL)[jj]->GetOriginalID();
      pss.evnum = evnum_;
      //pss.photid = (*myCollectionL)[jj]->GetPhotonID();
      //BS
      TVector3 myvec( pss.verx/100.,pss.very/100.,7.);

      //      if(! pss.isOpt && pss.isHit) EnergyDepositShapeL->Fill(pss.zdepth,myvec.Eta(),pss.en);

      //BS
      if(pss.isHit)
      {
	EdeptotL+=(*myCollectionL)[jj]->GetEdep()/GeV;

	EdeptotModL+=(*myCollectionL)[jj]->GetEdep()/GeV*pss.gainRandVar[pss.nstb-1][pss.row][pss.col];
	EdeptotQuadL+=(*myCollectionL)[jj]->GetEdep()/GeV*pss.gainGroupVar[pss.nstb-1][pss.row][pss.col];

	//North top
	if( (pss.nstb==1 && pss.row<17) || (pss.nstb==3 && pss.row<12) ) EdeptotLeadL+=(*myCollectionL)[jj]->GetEdep()/GeV*pss.leadClusterGain[0];
	//North bottom 
	if( (pss.nstb==1 && pss.row>16) || (pss.nstb==3 && pss.row>11) ) EdeptotLeadL+=(*myCollectionL)[jj]->GetEdep()/GeV*pss.leadClusterGain[1];
	//South top
	if( (pss.nstb==2 && pss.row<17) || (pss.nstb==4 && pss.row<12) ) EdeptotLeadL+=(*myCollectionL)[jj]->GetEdep()/GeV*pss.leadClusterGain[2];
	//South bottom
	if( (pss.nstb==2 && pss.row>16) || (pss.nstb==4 && pss.row>11) ) EdeptotLeadL+=(*myCollectionL)[jj]->GetEdep()/GeV*pss.leadClusterGain[3];
      }


      float EModPar0=1.09412 * pow(10,-8);
      float EModPar1=-0.026195;
      float EModTot= (*myCollectionL)[jj]->GetEdep()/GeV * exp(EModPar0 + EModPar1* (30.- ( (*myCollectionL)[jj]->GetPos().getZ() /cm ) ) );
      if(EModTot <0){std::cout<<"EMOD TOT LESS THAN 0!!!!!!!!!  = "<<EModTot<<std::endl; EModTot=0;}

      pss.Emod=EModTot;
      //BS NEW METHOD 08/11/17
      pss.FillMatrices(0);//Always fill total event matrices
      for(Int_t iterAncest=0; iterAncest < ((*myCollectionL)[jj]->GetAncestors()).size(); iterAncest++)
      {
	Bool_t foundOne=false;
	for(Int_t iterTrackID=0; iterTrackID<nsaved_traj; iterTrackID++)
	{
	  if( ((*myCollectionL)[jj]->GetAncestors()).at(iterAncest)==EventAction::trackID[iterTrackID])
	  {
	    foundOne=true;
	    pss.FillMatrices(iterTrackID+1); //+1 because 0 is full event
	    EdeptotNewL[iterTrackID] += (Float_t) pss.en;
	    EdeptotMrigankaNewL[iterTrackID] += (Float_t) pss.Emod;
	  }
	}
	if( !foundOne && iterAncest) std::cout<<"ERROR!!!! AncestorID of "<< ((*myCollectionL)[jj]->GetAncestors()).at(iterAncest)<<" NOT FOUND IN TRAJECTORIES"<<std::endl;
      }


      // 
      /*
	 for(int index=0;index<((*myCollectionL)[jj]->GetAncestors()).size();index++)
	 {
	 int inderx=-1;
	 bool needthis=true;
	 while(inderx<nsaved_traj){
	 inderx++;
	 if(((*myCollectionL)[jj]->GetAncestors()).at(index)==EventAction::trackID[inderx]){break;}
	 else{if(inderx==nsaved_traj-1){needthis=false;}}
	 }
	 if(needthis){
	 if(index==0){
	 pss.FillMatrices(inderx+1,true);
	 }
	 else{
	 if(((*myCollectionL)[jj]->GetAncestors()).at(index)!=1){
	 pss.FillMatrices(inderx+1,false);}

      //BS Damage Additions
      EdeptotNewL[inderx]+=(Float_t) pss.en;
      EdeptotMrigankaNewL[inderx]+=(Float_t) pss.Emod;
      //End BS Damage Additions
      }
      }
      else{ if(((*myCollectionL)[jj]->GetAncestors()).at(index)!=0) { std::cout<<"\n*****ERROR: Saved TrackID MisMatch*****\n"<<((*myCollectionL)[jj]->GetAncestors()).at(index)<<" is not in trackID[]\n";}}
      //	else{ std::cout<<"\n*****ERROR: Saved TrackID MisMatch*****\n"<<((*myCollectionL)[jj]->GetAncestors()).at(index)<<" is not in trackID[]\n";}
      }
      */

    };
    G4cout<< "Esum=->"<<Esum<<" Ncerenkov="<<nCerenkov<<"("<<1.*nCerenkov/Esum<<" photons/GeV)"<<G4endl;
    //BS DAMAGE ADDITIONS
    for(Int_t dumTraj=0; dumTraj<nsaved_traj; dumTraj++)
    {
      for(Int_t dumDays=0; dumDays<nDaysSimulated; dumDays++)
      {
	if(nPhotonsRadNewL[dumTraj][dumDays]==0) {zAvgNewL[dumTraj][dumDays]=-5;}
	else {zAvgNewL[dumTraj][dumDays] /= nPhotonsRadNewL[dumTraj][dumDays];}

	if(nPhotonsWaveNewL[dumTraj][dumDays]==0) {zAvgWaveNewL[dumTraj][dumDays]=-5;}
	else {zAvgWaveNewL[dumTraj][dumDays] /= nPhotonsWaveNewL[dumTraj][dumDays];}

      }
    } //dumTraj block

    //END BS DAMAGE ADDITIONS

  } //myCollectionL Block

  /* small cell action */
  if(myCollectionS)
  {
    nhits=myCollectionS->entries();
    G4cout<<"Try to Print " << nhits <<" Small Cell Hits"<<G4endl;
    Esum=0; 
    nCerenkov=0;
    for(UInt_t jj=0;jj<nhits;jj++)
    {
      //if( primary!=0 && primary != (*myCollectionS)[jj]->GetOriginalID() ) {continue;}
      pss.isOpt = 0;
      int jc=(*myCollectionS)[jj]->GetCellNb();
      pss.isHit = 0;
      if((*myCollectionS)[jj]->IsHit()) pss.isHit = 1;

      pss.nstb = cellgeo->GetNstb(jc,0);

      Bool_t DidSurvive=true;
      pss.isRealHit=false;
      pss.zdepth=15.5- ( (*myCollectionS)[jj]->GetPos().getZ() /cm ); //distance to end of cell
      if((*myCollectionS)[jj]->IsOptical())
      {
	pss.isOpt=1;     

	nPhotonsGeneratedS++;

	Double_t kinE = (Double_t) (*myCollectionS)[jj]->GetKinE()/GeV;
	Double_t Wave = TMath::HC()/(TMath::Qe() *kinE); //in nm

	//	OrigPhotonWavelengthS->Fill(Wave);

	if( pss.isHit)
	{
	  TVector3 posvec( (*myCollectionS)[jj]->GetPos().getX() /cm ,(*myCollectionS)[jj]->GetPos().getY() /cm , 700.);

	  Double_t HitEta = posvec.Eta();
	  Double_t CosTh= (*myCollectionS)[jj]->GetCosTh();
	  Int_t CurrentBin = Par0S->FindBin(HitEta);

	  RadParametersUpdatedS[0]= NormFactorS * (*Par0S)[CurrentBin];
	  RadParametersUpdatedS[1]=(*Par1S)[CurrentBin];
	  RadParametersUpdatedS[2]=(*Par2S)[CurrentBin];

	  //BS Damage Additions
	  for(Int_t dumAncestLoop=1;dumAncestLoop<((*myCollectionS)[jj]->GetAncestors()).size();dumAncestLoop++)
	  {
	    for(Int_t dumNSavedLoop=0; dumNSavedLoop<nsaved_traj; dumNSavedLoop++)
	    {
	      if(((*myCollectionS)[jj]->GetAncestors()).at(dumAncestLoop)==EventAction::trackID[dumNSavedLoop])
	      {
		for(Int_t dumNDays=0; dumNDays<nDaysSimulated; dumNDays++)
		{
		  Double_t nDaysActual=7*dumNDays;
		  if( EventAction::DidSurviveRad( HitEta, pss.zdepth, CosTh, false, 0.005, nDaysActual) )
		  {
		    nPhotonsRadNewS[dumNSavedLoop][dumNDays]+=1;
		    zAvgNewS[dumNSavedLoop][dumNDays]+= (Float_t) pss.zdepth;
		    if( (*myCollectionS)[jj]->GetIsCerenkov() ) nPhotonsCerenNewS[dumNSavedLoop][dumNDays]+=1;
		  }
		  if( EventAction::DidSurviveRad( HitEta, pss.zdepth, CosTh, true, 0.005, nDaysActual, true, Wave) ) 
		  { 
		    nPhotonsWaveNewS[dumNSavedLoop][dumNDays]+=1;
		    zAvgWaveNewS[dumNSavedLoop][dumNDays]+= (Float_t) pss.zdepth;
		  }
		}
	      }
	    }
	  }
	  //END BS DamageAddition

	  DamageFuncUpdated->SetParameters(RadParametersUpdatedS);

	  Double_t ProbabilitySurviveUpdated=0;
	  if(CosTh>0)
	  {
	    ProbabilitySurviveUpdated=100*TMath::Exp(-1./fabs(CosTh)* DamageFuncUpdated->Integral(0,pss.zdepth) );
	  }
	  else if (CosTh<0)
	  {
	    ProbabilitySurviveUpdated=100*TMath::Exp(-1./fabs(CosTh)* DamageFuncUpdated->Integral(pss.zdepth,45) + DamageFuncUpdated->Integral(0,45) );
	  }
	  else
	  {
	    ProbabilitySurviveUpdated=0;
	  }
	  Double_t ProbabilitySurvive=100* TMath::Exp(-1./fabs(CosTh)* DamageFunc->Integral(45-pss.zdepth,45) );
	  if( rand()%100 < ProbabilitySurviveUpdated) pss.isRealHit = true;

	  nPhotonsNormS+=1;
	  zPhotonsNormS+=pss.zdepth;
	  lPhotonsNormS+= (*myCollectionS)[jj]->GetTotalTrackLength();
	  if(pss.isRealHit)
	  {
	    nPhotonsRadS+=1;
	    zPhotonsRadS+=pss.zdepth;
	    lPhotonsRadS+= (*myCollectionS)[jj]->GetTotalTrackLength();
	  }
	  //	  std::cout<<"Small Cell Depth ="<<pss.zdepth<<" and extra Probability to survive ="<<ProbabilitySurviveUpdated<<std::endl;
	}
      }


      pss.row = cellgeo->GetRow(jc,0);
      pss.col = cellgeo->GetCol(jc,0);
      pss.gain = cellgeo->GetGain(pss.nstb,pss.row,pss.col);
      pss.corr = cellgeo->GetCorr(pss.nstb,pss.row,pss.col);
      pss.en = (*myCollectionS)[jj]->GetEdep()/GeV;
      pss.tot_en = (*myCollectionS)[jj]->GetTdep()/GeV;
      pss.kin = (*myCollectionS)[jj]->GetKinE()/GeV;
      if(pss.isOpt==1) nCer[jc]+=1;
      Esum+=(*myCollectionS)[jj]->GetEdep()/GeV;
      if(pss.isOpt==1)nCerenkov+=1;
      pss.verx = (*myCollectionS)[jj]->GetPos().x()/cm;
      pss.very = (*myCollectionS)[jj]->GetPos().y()/cm;
      pss.verz = (*myCollectionS)[jj]->GetPos().z()/cm;
      pss.px = (*myCollectionS)[jj]->GetMom().x();
      pss.py = (*myCollectionS)[jj]->GetMom().y();
      pss.pz = (*myCollectionS)[jj]->GetMom().z();
      sprintf(pss.name,"%s",(*myCollectionS)[jj]->GetName().data());

      pss.trNum = (*myCollectionS)[jj]->GetTrackNo();
      pss.cosTh = (*myCollectionS)[jj]->GetCosTh();
      pss.NCer = (*myCollectionS)[jj]->GetnCerenkov();
      pss.gtime = (*myCollectionS)[jj]->GetGlobalTime()/ns;
      pss.ltime = (*myCollectionS)[jj]->GetLocalTime()/ns;
      pss.originalID = (*myCollectionS)[jj]->GetOriginalID();

      pss.evnum = evnum_;

      //BS
      TVector3 myvec( pss.verx/100.,pss.very/100.,7.);

      //      if(! pss.isOpt && pss.isHit) EnergyDepositShapeS->Fill(pss.zdepth,myvec.Eta(),pss.en);
      if(pss.isHit)
      {
	EdeptotS+=(*myCollectionS)[jj]->GetEdep()/GeV;
	EdeptotModS+=(*myCollectionS)[jj]->GetEdep()/GeV*pss.gainRandVar[pss.nstb-1][pss.row][pss.col];
	EdeptotQuadS+=(*myCollectionS)[jj]->GetEdep()/GeV*pss.gainGroupVar[pss.nstb-1][pss.row][pss.col];

	//North top
	if( (pss.nstb==1 && pss.row<17) || (pss.nstb==3 && pss.row<12) ) EdeptotLeadS+=(*myCollectionS)[jj]->GetEdep()/GeV*pss.leadClusterGain[0];
	//North bottom 
	if( (pss.nstb==1 && pss.row>16) || (pss.nstb==3 && pss.row>11) ) EdeptotLeadS+=(*myCollectionS)[jj]->GetEdep()/GeV*pss.leadClusterGain[1];
	//South top
	if( (pss.nstb==2 && pss.row<17) || (pss.nstb==4 && pss.row<12) ) EdeptotLeadS+=(*myCollectionS)[jj]->GetEdep()/GeV*pss.leadClusterGain[2];
	//South bottom
	if( (pss.nstb==2 && pss.row>16) || (pss.nstb==4 && pss.row>11) ) EdeptotLeadS+=(*myCollectionS)[jj]->GetEdep()/GeV*pss.leadClusterGain[3];
      }

      float EModPar0S=-2.02547 * pow(10,-8);
      float EModPar1S=-0.039964;
      float EModTot= (*myCollectionS)[jj]->GetEdep()/GeV * exp(EModPar0S + EModPar1S * (22.5 - ( (*myCollectionS)[jj]->GetPos().getZ() /cm ) ) );
      pss.Emod=EModTot;

      //BS NEW ADDITION 08/11/2017
      pss.FillMatrices(0);//Always fill total event matrices
      for(Int_t iterAncest=0; iterAncest < ((*myCollectionS)[jj]->GetAncestors()).size(); iterAncest++)
      {
	Bool_t foundOne=false;
	for(Int_t iterTrackID=0; iterTrackID<nsaved_traj; iterTrackID++)
	{
	  if( ((*myCollectionS)[jj]->GetAncestors()).at(iterAncest)==EventAction::trackID[iterTrackID])
	  {
	    foundOne=true;
	    pss.FillMatrices(iterTrackID+1); //+1 because 0 is full event
	    EdeptotNewS[iterTrackID] += (Float_t) pss.en;
	    EdeptotMrigankaNewS[iterTrackID] += (Float_t) pss.Emod;
	  }
	}
	if( !foundOne && iterAncest) std::cout<<"ERROR!!!! AncestorID of "<< ((*myCollectionS)[jj]->GetAncestors()).at(iterAncest)<<" NOT FOUND IN TRAJECTORIES"<<std::endl;
      }

      /*
	 for(int index=0;index<((*myCollectionS)[jj]->GetAncestors()).size();index++){
	 int inderxx=-1;
	 bool needthis2=true;
	 while(inderxx<nsaved_traj){
	 inderxx++;
	 if(((*myCollectionS)[jj]->GetAncestors()).at(index)==EventAction::trackID[inderxx]){break;}
	 else{if(inderxx==nsaved_traj-1){needthis2=false;}}
	 }
	 if(needthis2){
	 if(index==0){
	 pss.FillMatrices(inderxx+1,true);
	 }
	 else{
	 if(((*myCollectionS)[jj]->GetAncestors()).at(index)!=1){
	 pss.FillMatrices(inderxx+1,false);}

      //BS Damage Additions
      EdeptotNewS[inderxx]+=(Float_t) pss.en;
      EdeptotMrigankaNewS[inderxx]+=(Float_t) pss.Emod;
      //End BS Damage Additions
      }
      }
      else{if(((*myCollectionS)[jj]->GetAncestors()).at(index)!=0) std::cout<<"\n*****ERROR: Saved TrackID MisMatch*****\n"<<((*myCollectionS)[jj]->GetAncestors()).at(index)<<" is not in trackID[]\n";}
      }
      */

    };
    G4cout<< "Esum=->"<<Esum<<" Ncerenkov="<<nCerenkov<<"("<<1.*nCerenkov/Esum<<" photons/GeV)"<<G4endl;
    //BS DAMAGE ADDITIONS
    for(Int_t dumTraj=0; dumTraj<nsaved_traj; dumTraj++)
    {
      for(Int_t dumDays=0; dumDays<nDaysSimulated; dumDays++)
      {
	if(nPhotonsRadNewS[dumTraj][dumDays]==0) {zAvgNewS[dumTraj][dumDays]=-5;}
	else { zAvgNewS[dumTraj][dumDays] /= nPhotonsRadNewS[dumTraj][dumDays];}
	if(nPhotonsWaveNewS[dumTraj][dumDays]==0) {zAvgWaveNewS[dumTraj][dumDays]=-5;}
	else { zAvgWaveNewS[dumTraj][dumDays] /= nPhotonsWaveNewS[dumTraj][dumDays];}
      }
    }
    //END BS DAMAGE ADDITIONS
  }


  //BS Radiation Check
  zAvgNormL = (Double_t)  zPhotonsNormL/nPhotonsNormL;
  zAvgRadL = (Double_t)  zPhotonsRadL/nPhotonsRadL;
  lAvgNormL = (Double_t)  lPhotonsNormL/nPhotonsNormL;
  lAvgRadL = (Double_t)  lPhotonsRadL/nPhotonsRadL;

  zAvgNormS = (Double_t)  zPhotonsNormS/nPhotonsNormS;
  zAvgRadS = (Double_t)  zPhotonsRadS/nPhotonsRadS;
  lAvgNormS = (Double_t)  lPhotonsNormS/nPhotonsNormS;
  lAvgRadS = (Double_t)  lPhotonsRadS/nPhotonsRadS;

  // get HT and digitise; fill ttr & hsimu; then zero matrices
  for(int k=0;k<=nsaved_traj;k++){
    pss.GetHT(k);
    printf("HT (n%d,r%d,c%d)\n",pss.nstb_ht,pss.row_ht,pss.col_ht);
    for(int n=0; n<4; n++) pss.Digitise(n,k);

    for(Int_t dumInt=0;dumInt<3700;dumInt++)
    {
      if(dumInt<pss.nQTdata){test.QTdata[dumInt]=(pss.QTdata)[dumInt];}
      else{test.QTdata[dumInt]=0;}
    }

    for(Int_t dumIntEnergy=0;dumIntEnergy<3700;dumIntEnergy++)
    {
      if(dumIntEnergy<pss.nqtEnergy){energyG4QT.QTdata[dumIntEnergy]=(pss.QtEnergy)[dumIntEnergy]; }
      else{energyG4QT.QTdata[dumIntEnergy]=0;}
    }

    for(Int_t dumIntMod=0;dumIntMod<3700;dumIntMod++)
    {
      if(dumIntMod<pss.nqtEnergyMod){EnergyModG4QT.QTdata[dumIntMod]=(pss.QtEnergyMod)[dumIntMod];}
      else{EnergyModG4QT.QTdata[dumIntMod]=0;}
    }

    for(Int_t dumIntRad=0;dumIntRad<3700;dumIntRad++)
    {
      if(dumIntRad<pss.nqtRadDamage){RadDamageG4QT.QTdata[dumIntRad]=(pss.QtRadDamage)[dumIntRad];}
      else{RadDamageG4QT.QTdata[dumIntRad]=0;}
    }

    for(Int_t dumIntGeo1=0;dumIntGeo1<3700;dumIntGeo1++)
    {
      if(dumIntGeo1<pss.nqtGeo1){CellGeo1G4QT.QTdata[dumIntGeo1]=(pss.QtGeo1)[dumIntGeo1];}
      else{CellGeo1G4QT.QTdata[dumIntGeo1]=0;}
    }

    for(Int_t dumIntGeo2=0;dumIntGeo2<3700;dumIntGeo2++)
    {
      if(dumIntGeo2<pss.nqtGeo2){CellGeo2G4QT.QTdata[dumIntGeo2]=(pss.QtGeo2)[dumIntGeo2];}
      else{CellGeo2G4QT.QTdata[dumIntGeo2]=0;}
    }

    for(Int_t dumIntGeo3=0;dumIntGeo3<3700;dumIntGeo3++)
    {
      if(dumIntGeo3<pss.nqtGeo3){CellGeo3G4QT.QTdata[dumIntGeo3]=(pss.QtGeo3)[dumIntGeo3];}
      else{CellGeo3G4QT.QTdata[dumIntGeo3]=0;}
    }

    for(Int_t dumIntGeo4=0;dumIntGeo4<3700;dumIntGeo4++)
    {
      if(dumIntGeo4<pss.nqtGeo4){CellGeo4G4QT.QTdata[dumIntGeo4]=(pss.QtGeo4)[dumIntGeo4];}
      else{CellGeo4G4QT.QTdata[dumIntGeo4]=0;}
    }

    for(Int_t dumIntGeo5=0;dumIntGeo5<3700;dumIntGeo5++)
    {
      if(dumIntGeo5<pss.nqtGeo5){CellGeo5G4QT.QTdata[dumIntGeo5]=(pss.QtGeo5)[dumIntGeo5];}
      else{CellGeo5G4QT.QTdata[dumIntGeo5]=0;}
    }

    test.SetEvent(eventNumber);
    test.SetNQTdata(pss.nQTdata);
    test.SetNPrimaries( (UInt_t) pss.nprimaries);
    test.SetNStored( (UInt_t) dummy);
    myvectest.push_back(test);

    energyG4QT.SetEvent(eventNumber);
    energyG4QT.SetNQTdata(pss.nqtEnergy);
    energyG4QT.SetNPrimaries( (UInt_t) pss.nprimaries);
    energyG4QT.SetNStored( (UInt_t) dummy);
    energyQTvec.push_back(energyG4QT);

    EnergyModG4QT.SetEvent(eventNumber);
    EnergyModG4QT.SetNQTdata(pss.nqtEnergyMod);
    EnergyModG4QT.SetNPrimaries( (UInt_t) pss.nprimaries);
    EnergyModG4QT.SetNStored( (UInt_t) dummy);
    EnergyModVec.push_back(EnergyModG4QT);

    RadDamageG4QT.SetEvent(eventNumber);
    RadDamageG4QT.SetNQTdata(pss.nqtRadDamage);
    RadDamageG4QT.SetNPrimaries( (UInt_t) pss.nprimaries);
    RadDamageG4QT.SetNStored( (UInt_t) dummy);
    RadDamageVec.push_back(RadDamageG4QT);

    CellGeo1G4QT.SetEvent(eventNumber);
    CellGeo1G4QT.SetNQTdata(pss.nqtGeo1);
    CellGeo1G4QT.SetNPrimaries( (UInt_t) pss.nprimaries);
    CellGeo1G4QT.SetNStored( (UInt_t) dummy);
    CellGeo1Vec.push_back(CellGeo1G4QT);

    CellGeo2G4QT.SetEvent(eventNumber);
    CellGeo2G4QT.SetNQTdata(pss.nqtGeo2);
    CellGeo2G4QT.SetNPrimaries( (UInt_t) pss.nprimaries);
    CellGeo2G4QT.SetNStored( (UInt_t) dummy);
    CellGeo2Vec.push_back(CellGeo2G4QT);

    CellGeo3G4QT.SetEvent(eventNumber);
    CellGeo3G4QT.SetNQTdata(pss.nqtGeo3);
    CellGeo3G4QT.SetNPrimaries( (UInt_t) pss.nprimaries);
    CellGeo3G4QT.SetNStored( (UInt_t) dummy);
    CellGeo3Vec.push_back(CellGeo3G4QT);

    CellGeo4G4QT.SetEvent(eventNumber);
    CellGeo4G4QT.SetNQTdata(pss.nqtGeo4);
    CellGeo4G4QT.SetNPrimaries( (UInt_t) pss.nprimaries);
    CellGeo4G4QT.SetNStored( (UInt_t) dummy);
    CellGeo4Vec.push_back(CellGeo4G4QT);

    CellGeo5G4QT.SetEvent(eventNumber);
    CellGeo5G4QT.SetNQTdata(pss.nqtGeo5);
    CellGeo5G4QT.SetNPrimaries( (UInt_t) pss.nprimaries);
    CellGeo5G4QT.SetNStored( (UInt_t) dummy);
    CellGeo5Vec.push_back(CellGeo5G4QT);

  }

  testtree->Fill();
  //Replace this change into main program
  pFile->cd();
}

//BS Damage Additions

Int_t EventAction::CreateNormFunctions(Double_t normEtaL_ /*=2.9*/, Double_t normEtaS_/*=3.7*/)
{
  if( pEShapeFile==0 ) { G4cout<<"ERROR: YOU HAVE NOT PROVIDED THE HISTOGRAM FILE"<<G4endl; return -1;}
  else
  {
    TString dumString;//Dummmy string to hold conversions from numbers
    //Set the Histograms and parameters
    for(Int_t whichParam=0; whichParam<3; whichParam++)
    {
      dumString.Form("%i",whichParam);
      TString nameTot= "EnergyDepositShapeL_" + dumString;
      pParamL[whichParam]=(TH1D*) pEShapeFile->Get(nameTot);
      nameTot= "EnergyDepositShapeS_" + dumString;
      pParamS[whichParam]=(TH1D*) pEShapeFile->Get(nameTot);
    }

    pDamGaus=new TF1("pDamGaus","gaus");

    //Grab the nomalization bins
    Int_t normBinL=pParamL[0]->FindBin(normEtaL_);
    Int_t normBinS=pParamS[0]->FindBin(normEtaS_);
    ////Integrate, Create large/small cell norm func
    //Large
    pDamGaus->SetParameters((*(pParamL[0]))[normBinL],(*(pParamL[1]))[normBinL],(*(pParamL[2]))[normBinL]);
    Double_t integralValue= pDamGaus->Integral(0,60);
    dumString.Form("%f",integralValue);
    TString funcString=" -log( pow( 1. - x, y) ) / " + dumString;
    pNormFuncL = new TF2("pNormFuncL",funcString,0,0.99999999,0,100000); //x is ratio Change per day (techincally per unit y), y is number of days. String of 0.9 is to prevent log(0), may be better way to avoid that
    //Small
    pDamGaus->SetParameters((*(pParamS[0]))[normBinS],(*(pParamS[1]))[normBinS],(*(pParamS[2]))[normBinS]);
    integralValue= pDamGaus->Integral(0,45);
    dumString.Form("%f",integralValue);
    funcString=" -log( pow( 1. - x, y) ) / " + dumString;
    pNormFuncS = new TF2("pNormFuncS",funcString,0,0.99999999,0,100000); //x is ratio Change per day (techincally per unit y), y is number of days. String of 0.9 is to prevent log(0), may be better way to avoid that
  }
  return 0;
}

//Bool_t EventAction::DidSurviveRad(Double_t etaPhoton_, Double_t depthPhoton_, Double_t cosThPhoton_, Bool_t isLarge_, Double_t ratDayChange_, Double_t nDays_)
Bool_t EventAction::DidSurviveRad(Double_t etaPhoton_, Double_t depthPhoton_, Double_t cosThPhoton_, Bool_t isLarge_, Double_t ratDayChange_, Double_t nDays_, Bool_t useWavelength_/*=false*/, Double_t wavelength_/*=465*/)
{
  if(cosThPhoton_ == 0) return false; //Prevents division by 0;
  if(nDays_== 0) return true; //Removes ~1/5000 rounding errors
  Double_t normFactor=0, probSurvive=0;
  Int_t whichBin;
  Double_t gausParam[3]={0};
  Double_t fullLength= isLarge_ ? 60. : 45;
  Double_t waveFactor = useWavelength_ ? pFuncWaveDam->Eval(wavelength_) : 1.;

  if(isLarge_)
  {
    normFactor=pNormFuncL->Eval(ratDayChange_, nDays_);
    if(nDays_==0) G4cout<<"NDays=0 && normFactor="<<normFactor<<G4endl;
    whichBin= pParamL[0]->FindBin(etaPhoton_);
    for(Int_t dumGausParam=0; dumGausParam<3; dumGausParam++) {  gausParam[dumGausParam]= (*(pParamL[dumGausParam]))[whichBin]; }
    gausParam[0]*=normFactor;
    pDamGaus->SetParameters(gausParam);
    if(cosThPhoton_ > 0)
    {
      probSurvive=100*TMath::Exp(-1./fabs(cosThPhoton_) * waveFactor * pDamGaus->Integral(0,depthPhoton_) );
    }
    else //Backwards facing photon that arrived at Photocathode
    {
      //Assume one reflection max from front of glass
      probSurvive=100.*TMath::Exp(-1./fabs(cosThPhoton_) * waveFactor* ( pDamGaus->Integral(depthPhoton_,fullLength)  + pDamGaus->Integral(0,fullLength)) );
    }
  }
  else //Small cell
  {
    normFactor=pNormFuncS->Eval(ratDayChange_, nDays_);
    whichBin= pParamS[0]->FindBin(etaPhoton_);
    for(Int_t dumGausParam=0; dumGausParam<3; dumGausParam++) {  gausParam[dumGausParam]= (*(pParamL[dumGausParam]))[whichBin]; }
    gausParam[0]*=normFactor;
    pDamGaus->SetParameters(gausParam);
    if(cosThPhoton_>0)
    {
      probSurvive=100.*TMath::Exp(-1./fabs(cosThPhoton_) * waveFactor * pDamGaus->Integral(0,depthPhoton_) );
    }
    else
    {
      //Assume one reflection max from front of glass
      probSurvive=100.*TMath::Exp(-1./fabs(cosThPhoton_) * waveFactor * ( pDamGaus->Integral(depthPhoton_,fullLength)  + pDamGaus->Integral(0,fullLength) ) );
    }
  }

  if( rand()%100 <= probSurvive) return true;
  else {return false;}
}
//END BS Damage Additions
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
