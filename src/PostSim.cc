#include "PostSim.hh"
#include "Rtypes.h"
#include "Trajectory.hh"

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>
#include "G4ThreeVector.hh"
#include "G4ios.hh"     
#include "globals.hh" 
#include "G4ParticleDefinition.hh" 
#include "G4TrajectoryPoint.hh"   
#include "G4Track.hh"
#include "G4Step.hh"
#include <vector>

PostSim::PostSim()
{
  MatrixArray.SetOwner((Bool_t) true);
  evpos.SetOwner((Bool_t) true);

  //  numberMatrices = 14; //Number of Matrices. Note, nMatrices also exists. Fix this
  numberMatrices = 15; //Number of Matrices. Note, nMatrices also exists. Fix this

  evnum=0;
  nstb=0;
  row=0;
  col=0;
  gain=0;
  corr=0;
  en=0;
  tot_en=0;
  kin=0;
  verx=0;
  very=0;
  verz=0;
  px=0;
  py=0;
  pz=0;
  strcpy(name,"");
  isOpt=0;
  isHit=0;
  trNum=0;
  cosTh=0;
  NCer=0;
  ltime=0;
  gtime=0;
  arr_size=0;
  Emod=0; //BS

  // instantiate p_files and QT
  //  fu = new G4FileUtilities();
  FMSGAIN = fu.GetEnv("FMSGAIN");
  FMSCORR = fu.GetEnv("FMSCORR");
  QTMAP = fu.GetEnv("QTMAP");
  QTMAP2PP = fu.GetEnv("QTMAP2PP");
  FMSTXT = fu.GetEnv("FMSTXT");
  p_files = new FilesSet((TString) FMSTXT.data(),
      (TString) "fpdped.txt",
      (TString) FMSGAIN.data(),
      (TString) FMSCORR.data(),
      (TString) "fill.txt",
      (TString) "Fake",
      (TString) "spinpat",
      (TString) "geom.txt",
      (TString) QTMAP.data(),
      (TString) QTMAP2PP.data());
  QT = new Qt(p_files);


  nQTdataPointer = &nQTdata;
  pnqtEnergy=&nqtEnergy;
  pnqtEnergyMod=&nqtEnergyMod;
  pnqtRadDamage=&nqtRadDamage;

  pnqtGeo1=&nqtGeo1;
  pnqtGeo2=&nqtGeo2;
  pnqtGeo3=&nqtGeo3;
  pnqtGeo4=&nqtGeo4;
  pnqtGeo5=&nqtGeo5;


  ///////////
  //SetGains(FMSTXT, FMSGAIN, FMSCORR,FMSBIT)
  geo1.SetGains("/home/branden/TriggerTesting","","FmsCorrAll1.txt");
  //  geo2.SetGains("/home/branden/TriggerTesting","","FmsCorrVar10Percent.txt");
  geo2.SetGains("/home/branden/TriggerTesting","","FmsCorrAll1.txt");
  //  geo3.SetGains("/home/branden/TriggerTesting","","FmsCorrVar20Percent.txt");
  //  geo4.SetGains("/home/branden/TriggerTesting","","FmsCorrVar30Percent.txt");
  //  geo4.SetGains("/home/branden/TriggerTesting","","FmsCorr956.txt");
//  geo3.SetGains("/home/branden/TriggerTesting","","FmsCorrQ142.txt");
  geo3.SetGains("/home/branden/TriggerTesting","","FmsCorrAll1.txt");
  //  geo4.SetGains("/home/branden/TriggerTesting","","FmsCorrQ142.txt");
  //  geo4.SetGains("/home/branden/TriggerTesting","","FmsCorrD82a.txt");

  geo4.SetGains("/home/branden/TriggerTesting","","FmsCorrQ142.txt");
  geo5.SetGains("/home/branden/TriggerTesting","","FmsCorrQ142.txt");

}

PostSim::~PostSim()
{
  if (QT) delete QT;
  if (p_files) delete p_files;
}


// adds current gtree entry to appropriate matrices/objarrays
//void PostSim::FillMatrices(int i, bool t)
void PostSim::FillMatrices(int i)
{
  /*
  TMatrix ** sfadc = (TMatrix**) MatrixArray.At(0);
  TMatrix ** sedep = (TMatrix**) MatrixArray.At(1);
  TMatrix ** shits = (TMatrix**) MatrixArray.At(2);
  TMatrix ** snopt = (TMatrix**) MatrixArray.At(3);
  TMatrix ** sncer = (TMatrix**) MatrixArray.At(4);
  TMatrix ** sadc = (TMatrix **) MatrixArray.At(5);
  TMatrix ** sfadcM = (TMatrix **) MatrixArray.At(6);
  TMatrix ** sedepM = (TMatrix **) MatrixArray.At(7);
  TMatrix ** radhits = (TMatrix **) MatrixArray.At(8);
  TMatrix ** radadc = (TMatrix **) MatrixArray.At(9);
  TMatrix ** sadcgeo1 = (TMatrix **) MatrixArray.At(10);
  TMatrix ** sadcgeo2 = (TMatrix **) MatrixArray.At(11);
  TMatrix ** sadcgeo3 = (TMatrix **) MatrixArray.At(12);
  TMatrix ** sadcgeo4 = (TMatrix **) MatrixArray.At(13);
  TMatrix ** sadcgeo5 = (TMatrix **) MatrixArray.At(14);
  */


  TMatrix ** fadc1 = (TMatrix**) MatrixArray.At(i*numberMatrices);
  TMatrix ** sedep1 = (TMatrix**) MatrixArray.At(i*numberMatrices+1);
  TMatrix ** hits1 = (TMatrix**) MatrixArray.At(i*numberMatrices+2);
  TMatrix ** nopt1 = (TMatrix**) MatrixArray.At(i*numberMatrices+3);
  TMatrix ** ncer1 = (TMatrix**) MatrixArray.At(i*numberMatrices+4);
  TMatrix ** adc1 = (TMatrix**) MatrixArray.At(i*numberMatrices+5);
  TMatrix ** fadc1M = (TMatrix**) MatrixArray.At(i*numberMatrices+6);
  TMatrix ** sedep1M = (TMatrix**) MatrixArray.At(i*numberMatrices+7);
  TMatrix ** radhits1 = (TMatrix**) MatrixArray.At(i*numberMatrices+8);
  TMatrix ** radadc1 = (TMatrix**) MatrixArray.At(i*numberMatrices+8);
  TMatrix ** adcgeo1 = (TMatrix**) MatrixArray.At(i*numberMatrices+9);
  TMatrix ** adcgeo2 = (TMatrix**) MatrixArray.At(i*numberMatrices+10);
  TMatrix ** adcgeo3 = (TMatrix**) MatrixArray.At(i*numberMatrices+11);
  TMatrix ** adcgeo4 = (TMatrix**) MatrixArray.At(i*numberMatrices+12);
  TMatrix ** adcgeo5 = (TMatrix**) MatrixArray.At(i*numberMatrices+13);

  /*
  if(t){
    (*(sedepM[nstb-1]))(row,col) += Emod;
    (*(sedep[nstb-1]))(row,col) += en;
    (*(sncer[nstb-1]))(row,col) += NCer;
    if(isOpt)
    {
      (*(snopt[nstb-1]))(row,col) += 1;
      if(isHit) (*(shits[nstb-1]))(row,col) += 1;
      if(isRealHit) (*(radhits[nstb-1])) (row,col) +=1;
    };
  }
  */
  (*(sedep1M[nstb-1]))(row,col) += Emod;
  (*(sedep1[nstb-1]))(row,col) += en;
  (*(ncer1[nstb-1]))(row,col) += NCer;
  if(isOpt)
  {
    (*(nopt1[nstb-1]))(row,col) += 1;
    if(isHit)  (*(hits1[nstb-1]))(row,col) += 1;
    if(isRealHit)  (*(radhits1[nstb-1]))(row,col) += 1;
  };
}

void PostSim::FillMatrixArray(int m, Int_t n[])
{
  for(int i=0;i<=m;i++)
  {
    for(Int_t j=0; j<numberMatrices; j++)
    {
      MatrixArray.AddLast((TObject*) new TMatrix[4]);
    }
  }
  for(int i=0;i<=m;i++){
    /*BS. Possible change to TMatrix*  fadc1[4] */
    TMatrix ** fadc1 = (TMatrix**) MatrixArray.At(i*numberMatrices);
    TMatrix ** sedep1 = (TMatrix**) MatrixArray.At(i*numberMatrices+1);
    TMatrix ** hits1 = (TMatrix**) MatrixArray.At(i*numberMatrices+2);
    TMatrix ** nopt1 = (TMatrix**) MatrixArray.At(i*numberMatrices+3);
    TMatrix ** ncer1 = (TMatrix**) MatrixArray.At(i*numberMatrices+4);
    TMatrix ** adc1 = (TMatrix**) MatrixArray.At(i*numberMatrices+5);
    TMatrix ** fadc1M = (TMatrix**) MatrixArray.At(i*numberMatrices+6);
    TMatrix ** sedep1M = (TMatrix**) MatrixArray.At(i*numberMatrices+7);
    TMatrix ** radhits = (TMatrix**) MatrixArray.At(i*numberMatrices+8);
    TMatrix ** radadc = (TMatrix**) MatrixArray.At(i*numberMatrices+9);
    TMatrix ** adcgeo1 = (TMatrix**) MatrixArray.At(i*numberMatrices+10);
    TMatrix ** adcgeo2 = (TMatrix**) MatrixArray.At(i*numberMatrices+11);
    TMatrix ** adcgeo3 = (TMatrix**) MatrixArray.At(i*numberMatrices+12);
    TMatrix ** adcgeo4 = (TMatrix**) MatrixArray.At(i*numberMatrices+13);
    TMatrix ** adcgeo5 = (TMatrix**) MatrixArray.At(i*numberMatrices+14);

    for(int n=0; n<4; n++)
    {
      if(n<2) { nnr=34; nnc=17; }
      else { nnr=24; nnc=12; };
      fadc1[n] = new TMatrix(nnr,nnc);
      sedep1[n] = new TMatrix(nnr,nnc);
      hits1[n] = new TMatrix(nnr,nnc);
      nopt1[n] = new TMatrix(nnr,nnc);
      ncer1[n] = new TMatrix(nnr,nnc);
      adc1[n] = new TMatrix(nnr,nnc);
      fadc1M[n] = new TMatrix(nnr,nnc);
      sedep1M[n] = new TMatrix(nnr,nnc);
      radhits[n] = new TMatrix(nnr,nnc);
      radadc[n] = new TMatrix(nnr,nnc);
      adcgeo1[n] = new TMatrix(nnr,nnc);
      adcgeo2[n] = new TMatrix(nnr,nnc);
      adcgeo3[n] = new TMatrix(nnr,nnc);
      adcgeo4[n] = new TMatrix(nnr,nnc);
      adcgeo5[n] = new TMatrix(nnr,nnc);
    };
  }
  ZeroMatrices(m); //BS. Method in PostSim. Sets each element of the above list of matrices to 0 with a for loop.
}

// sets all elements of all private matrices to zero
void PostSim::ZeroMatrices(int m)
{
  int nrows,ncols;
  for(int i=0;i<=m;i++){
    TMatrix ** fadc0 = (TMatrix**) MatrixArray.At(i*numberMatrices);
    TMatrix ** sedep0 = (TMatrix**) MatrixArray.At(i*numberMatrices+1);
    TMatrix ** hits0 = (TMatrix**) MatrixArray.At(i*numberMatrices+2);
    TMatrix ** nopt0 = (TMatrix**) MatrixArray.At(i*numberMatrices+3);
    TMatrix ** ncer0 = (TMatrix**) MatrixArray.At(i*numberMatrices+4);
    TMatrix ** adc0 = (TMatrix**) MatrixArray.At(i*numberMatrices+5);
    TMatrix ** fadc0M = (TMatrix**) MatrixArray.At(i*numberMatrices+6);
    TMatrix ** sedep0M = (TMatrix**) MatrixArray.At(i*numberMatrices+7);
    TMatrix ** radhits0 = (TMatrix**) MatrixArray.At(i*numberMatrices+8);
    TMatrix ** radadc0 = (TMatrix**) MatrixArray.At(i*numberMatrices+9);
    TMatrix ** adcgeo1_0 = (TMatrix**) MatrixArray.At(i*numberMatrices+10);
    TMatrix ** adcgeo2_0 = (TMatrix**) MatrixArray.At(i*numberMatrices+11);
    TMatrix ** adcgeo3_0 = (TMatrix**) MatrixArray.At(i*numberMatrices+12);
    TMatrix ** adcgeo4_0 = (TMatrix**) MatrixArray.At(i*numberMatrices+13);
    TMatrix ** adcgeo5_0 = (TMatrix**) MatrixArray.At(i*numberMatrices+14);
    for(int n=0; n<4; n++)
    {
      nrows=34; ncols=17;
      if(n>1) { nrows=24; ncols=12; };
      for(int nr=0; nr<nrows; nr++)
      {
	for(int nc=0; nc<ncols; nc++)
	{
	  (*(fadc0[n]))(nr,nc)=0;
	  (*(sedep0[n]))(nr,nc)=0;
	  (*(hits0[n]))(nr,nc)=0;
	  (*(nopt0[n]))(nr,nc)=0;
	  (*(ncer0[n]))(nr,nc)=0;
	  (*(adc0[n]))(nr,nc)=0;
	  (*(fadc0M[n]))(nr,nc)=0;
	  (*(sedep0M[n]))(nr,nc)=0;
	  (*(radhits0[n]))(nr,nc)=0;
	  (*(radadc0[n]))(nr,nc)=0;
	  (*(adcgeo1_0[n]))(nr,nc)=0;
	  (*(adcgeo2_0[n]))(nr,nc)=0;
	  (*(adcgeo3_0[n]))(nr,nc)=0;
	  (*(adcgeo4_0[n]))(nr,nc)=0;
	  (*(adcgeo5_0[n]))(nr,nc)=0;
	}
      }
    }
  }
} 

// note: need to call "GetHT" before calling this
// digitise (hits matrix --> adc matrix)
// fills ttr
// encodes qt block and fills hsimu
// (see drupal entry on May 13 2013 for details)
void PostSim::Digitise(Int_t nn, int k)
{
  //  StackingAction *teststack = new StackingAction();
  //  TTree *mytree = teststack->GetCoordTr();
  //  Int_t nevents = (Int_t) mytree->GetEntries();
  //  TH1F *myhist= new TH1F("myhist","myhist",10000,0,100);
  //  mytree->Draw("energy>>myhist","generation == 0");
  // digitisation fit parameters
  /*
     ppL[0] = 0;
     ppL[1] = 0.00224282; //BS, replace this with just the lin slope
     ppL[2] = -4.46158e-08;
     ppL[3] = 1.30951e-12;
     ppL[4] = -1.34152e-17;
     ppS[0] = 0;
     ppS[1] = 0.00219184;
     ppS[2] = -4.57127e-08;
     ppS[3] = 1.25689e-12;
     ppS[4] = -1.21213e-17;

     //BS Initial Linear approximation
  ppL[0] = 0;
  ppL[1] = 0.00221733; //BS, linear slope. Average of values
  ppL[2] = 0;
  ppL[3] = 0;
  ppL[4] = 0; 
  ppS[0] = 0;
  ppS[1] = 0.00221733;
  ppS[2] = 0; 
  ppS[3] = 0;
  ppS[4] = 0; 
     */

  Double_t edeptot=0;
  Double_t efittot=0;
  Double_t eradtot=0;
  Double_t radhitstot=0;

  //BS. For clean cells, based on profiles of the random distribuiton I see an average of 547 photons/Gev around 40 GeV pions in Small Cells, and ~461 photons/GeV in large cells around 25 geV pions. Comparision in Large cells to gaussian fit of number around 25 gives peak at 460
  ppL[0] = 0;
  ppL[1] = 0.0021692; 
  ppL[2] = 0;
  ppL[3] = 0;
  ppL[4] = 0; 
  ppS[0] = 0;
  ppS[1] = 0.0018281;
  ppS[2] = 0; 
  ppS[3] = 0;
  ppS[4] = 0; 

  //  Int_t nMatrices=14;
  Int_t nMatrices=15;
  TMatrix ** sfadc = (TMatrix**) MatrixArray.At(k*nMatrices);
  TMatrix ** sedep = (TMatrix**) MatrixArray.At(k*nMatrices+1);
  TMatrix ** shits = (TMatrix**) MatrixArray.At(k*nMatrices+2);
  TMatrix ** snopt = (TMatrix**) MatrixArray.At(k*nMatrices+3);
  TMatrix ** sncer = (TMatrix**) MatrixArray.At(k*nMatrices+4);
  TMatrix ** sadc = (TMatrix**) MatrixArray.At(k*nMatrices+5);
  TMatrix ** sfadcM = (TMatrix**) MatrixArray.At(k*nMatrices+6);
  TMatrix ** sedepM = (TMatrix**) MatrixArray.At(k*nMatrices+7);
  TMatrix ** radhits = (TMatrix**) MatrixArray.At(k*nMatrices+8);
  TMatrix ** radadc = (TMatrix**) MatrixArray.At(k*nMatrices+9);
  TMatrix ** adcgeo1 = (TMatrix**) MatrixArray.At(k*nMatrices+10);
  TMatrix ** adcgeo2 = (TMatrix**) MatrixArray.At(k*nMatrices+11);
  TMatrix ** adcgeo3 = (TMatrix**) MatrixArray.At(k*nMatrices+12);
  TMatrix ** adcgeo4 = (TMatrix**) MatrixArray.At(k*nMatrices+13);
  TMatrix ** adcgeo5 = (TMatrix**) MatrixArray.At(k*nMatrices+14);


  for(Int_t r=0; r<shits[nn]->GetNrows(); r++)
  {
    for(Int_t c=0; c<shits[nn]->GetNcols(); c++)
    {
      hits = (*(shits[nn]))(r,c);
      edep = (*(sedep[nn]))(r,c);
      edepMod= (*(sedepM[nn]))(r,c);
      Float_t tempRadHits = (*(radhits[nn]))(r,c);
      radhitstot+= tempRadHits;

      if(hits>0 || edep>0)
      {
	gain_ = cellgeo->GetGain(nn+1,r,c);
	corr_ = cellgeo->GetCorr(nn+1,r,c);
	gaincorr = gain_*corr_;

	gain1=geo1.GetGain(nn+1,r,c);
	corr1=geo1.GetCorr(nn+1,r,c);
	gaincorr1=gain1*corr1;

	gain2=geo2.GetGain(nn+1,r,c);
	corr2=geo2.GetCorr(nn+1,r,c);
	gaincorr2=gain2*corr2;

	gain3=geo3.GetGain(nn+1,r,c);
	corr3=geo3.GetCorr(nn+1,r,c);
	gaincorr3=gain3*corr3;
	//BS ADDING RANDOM NUMBER TO GAIN
//	gaincorr3*=gainRandVar[nn][r][c];

	gain4=geo4.GetGain(nn+1,r,c);
	corr4=geo4.GetCorr(nn+1,r,c);
	gaincorr4=gain4*corr4;
//	gaincorr4*=gainGroupVar[nn][r][c];

	gain5=geo5.GetGain(nn+1,r,c);
	corr5=geo5.GetCorr(nn+1,r,c);
	gaincorr5=gain5*corr5;
	/*
	if(r< (shits[nn]->GetNrows()/2)) //top
	{
	  if(nn%2==0) { gaincorr5*=leadClusterGain[0]; } //north
	  else { gaincorr5*=leadClusterGain[2]; } //south
	}
	else //bottom
	{
	  if(nn%2==0) { gaincorr5*=leadClusterGain[1]; } //north
	  else { gaincorr5*=leadClusterGain[3]; } //south
	}
	*/


	if(nn<2)
	{
	  fac=17;
	  fit_edep = ppL[0]+
	    ppL[1]*hits+
	    ppL[2]*pow(hits,2)+
	    ppL[3]*pow(hits,3)+
	    ppL[4]*pow(hits,4);

	  fitRadEdep = ppL[1] * tempRadHits;
	}
	else
	{
	  fac=12;
	  fit_edep = ppS[0]+
	    ppS[1]*hits+
	    ppS[2]*pow(hits,2)+
	    ppS[3]*pow(hits,3)+
	    ppS[4]*pow(hits,4);

	  fitRadEdep = ppS[1] * tempRadHits;
	};
	adc = (Int_t)(fit_edep/gaincorr + 0.5); // round to nearest integer
	adcrad = (Int_t) (fitRadEdep/gaincorr + 0.5);
	fadc= (Int_t) ((edep)/(gaincorr) + 0.5);
	fadcMod=(Int_t) ( edepMod/(gaincorr) + 0.5);
	chan = c+r*fac+1;


	//BS
	adcGeoValue1= (Int_t) ((edep)/(gaincorr1) + 0.5);
	adcGeoValue2= (Int_t) ((edep)/(gaincorr2) + 0.5);
	adcGeoValue3= (Int_t) ((edep)/(gaincorr3) + 0.5);
	adcGeoValue4= (Int_t) ((edep)/(gaincorr4) + 0.5);
	adcGeoValue5= (Int_t) ((edep)/(gaincorr5) + 0.5);
	//End BS
	

	//Setting maximum adc value to 255. NOTE. need to do this for other adc as well

	adcGeoValue1 = adcGeoValue1>4095 ? 4095 : adcGeoValue1;
	adcGeoValue2 = adcGeoValue2>4095 ? 4095 : adcGeoValue2;
	adcGeoValue3 = adcGeoValue3>4095 ? 4095 : adcGeoValue3;
	adcGeoValue4 = adcGeoValue4>4095 ? 4095 : adcGeoValue4;
	adcGeoValue5 = adcGeoValue5>4095 ? 4095 : adcGeoValue5;

	Int_t origADC3=adcGeoValue3;
	Int_t origADC1=adcGeoValue1;
	// BS Oct 25
	// Adding in bitshift and subtraction of 1 (Run 15)
	// Right Now, assumes bitshift is positive
	Bool_t subtractOne=true;
	Int_t bitShift2=geo2.GetBitShift(nn+1,r,c);
	Int_t masker2 = 4095-pow(2,bitShift2)+1;

	adcGeoValue2= adcGeoValue2&masker2; //just bitshift
	adcGeoValue3= adcGeoValue3&masker2; //bitshift + subtraction
	if(adcGeoValue3>0 && subtractOne) adcGeoValue3 -= pow(2,bitShift2); //bitshift + subtraction

	/*
	if(origADC1>0 && k==0)
	{
	  std::cout<<"BITSHIFT="<<bitShift2<<std::endl;
	  std::cout<<"origADC1="<<origADC1<<std::endl;
	  std::cout<<"Final ADC2="<<adcGeoValue2<<std::endl;
	  std::cout<<"Final ADC3="<<adcGeoValue3<<std::endl;
	}
	*/


	Int_t bitShift5=geo4.GetBitShift(nn+1,r,c);
	Int_t masker5 = 4095-pow(2,bitShift5)+1;
	adcGeoValue5= adcGeoValue5&masker5;
	if(adcGeoValue5>0 && subtractOne) adcGeoValue5 -= pow(2,bitShift5); //bitshift + subtraction
	//

	(*(sadc[nn]))(r,c) = adc;
	(*(sfadc[nn]))(r,c)= fadc;
	(*(sfadcM[nn]))(r,c)=fadcMod;
	(*(radadc[nn]))(r,c) = adcrad;
	(*(adcgeo1[nn]))(r,c) = adcGeoValue1;
	(*(adcgeo2[nn]))(r,c) = adcGeoValue2;
	(*(adcgeo3[nn]))(r,c) = adcGeoValue3;
	(*(adcgeo4[nn]))(r,c) = adcGeoValue4;
	(*(adcgeo5[nn]))(r,c) = adcGeoValue5;

	nstb_tr = nn+1;
	row_tr = r;
	col_tr = c;

	edeptot+= edep;
	efittot+= fit_edep;
	eradtot+= fitRadEdep;
      }
    }
  }

  if(nn == 3)
  {
    for(int g=0; g<3700; g++) {QTdata[g]=0; QtEnergy[g]=0; QtEnergyMod[g]=0; QtGeo1[g]=0; QtGeo2[g]=0; QtGeo3[g]=0; QtGeo4[g]=0; QtGeo5[g]=0;}
    success = QT->encodeQT(sadc,nQTdataPointer,QTdata);
    Bool_t DidEncode= QT->encodeQT(sfadc,pnqtEnergy,QtEnergy);
    Bool_t DidEncodeMod= QT->encodeQT(sfadcM,pnqtEnergyMod,QtEnergyMod);
    Bool_t DidEncodeRadDamage= QT->encodeQT(radadc,pnqtRadDamage,QtRadDamage);
    Bool_t DidEncodeGeo1= QT->encodeQT(adcgeo1,pnqtGeo1,QtGeo1);
    Bool_t DidEncodeGeo2= QT->encodeQT(adcgeo2,pnqtGeo2,QtGeo2);
    Bool_t DidEncodeGeo3= QT->encodeQT(adcgeo3,pnqtGeo3,QtGeo3);
    Bool_t DidEncodeGeo4= QT->encodeQT(adcgeo4,pnqtGeo4,QtGeo4);
    //
    Bool_t DidEncodeGeo5= QT->encodeQT(adcgeo5,pnqtGeo5,QtGeo5);




    nQTdata = *nQTdataPointer;
    nqtEnergy=*pnqtEnergy;
    nqtEnergyMod=*pnqtEnergyMod;
    nqtRadDamage=*pnqtRadDamage;
    nqtGeo1=*pnqtGeo1;
    nqtGeo2=*pnqtGeo2;
    nqtGeo3=*pnqtGeo3;
    nqtGeo4=*pnqtGeo4;
    nqtGeo5=*pnqtGeo5;
    if(dummyFirstInEvent==0)
    {

      //BS CHANGE
      /*
	 G4cout<<"MATRIX 0="<<G4endl;
	 sedep[0]->Print();
	 G4cout<<"MATRIX 1="<<G4endl;
	 sedep[1]->Print();
	 G4cout<<"MATRIX 2="<<G4endl;
	 sedep[2]->Print();
	 G4cout<<"MATRIX 3="<<G4endl;
	 sedep[3]->Print();
	 */


      for(Int_t i=0; i<9; i++)
      {
	Trig1[i]=0;
	Trig2[i]=0;
	Trig3[i]=0;
	Trig4[i]=0;
	Trig5[i]=0;
      }

      QT->decodeQT(nqtGeo1,QtGeo1,0);
      QT->TriggersTest();
      for(Int_t qtdum1=0; qtdum1<9; qtdum1++) 
      {
	Trig1[qtdum1] = (Int_t) (QT->TriggerPass)[qtdum1];
//	std::cout<<"QT 1 Trig["<<qtdum1<<"]="<<Trig1[qtdum1]<<" AND TriggerPass="<<(QT->TriggerPass)[qtdum1]<<std::endl;
      }

      QT->decodeQT(nqtGeo2,QtGeo2,0);
      QT->TriggersTest();
      for(Int_t qtdum2=0; qtdum2<9; qtdum2++)
      {
	Trig2[qtdum2] = (Int_t) (QT->TriggerPass)[qtdum2];
//	std::cout<<"QT 2 Trig["<<qtdum2<<"]="<<Trig2[qtdum2]<<" AND TriggerPass="<<(QT->TriggerPass)[qtdum2]<<std::endl;
      }

      QT->decodeQT(nqtGeo3,QtGeo3,0);
      QT->TriggersTest();
      for(Int_t qtdum3=0; qtdum3<9; qtdum3++)
      {
	Trig3[qtdum3] = (Int_t) (QT->TriggerPass)[qtdum3];
//	std::cout<<"QT 3 Trig["<<qtdum3<<"]="<<Trig3[qtdum3]<<" AND TriggerPass="<<(QT->TriggerPass)[qtdum3]<<std::endl;
      }

      QT->decodeQT(nqtGeo4,QtGeo4,0);
      QT->TriggersTest();
      for(Int_t qtdum4=0; qtdum4<9; qtdum4++) Trig4[qtdum4] = (Int_t) (QT->TriggerPass)[qtdum4];

      QT->decodeQT(nqtGeo5,QtGeo5,0);
      QT->TriggersTest();
      for(Int_t qtdum5=0; qtdum5<9; qtdum5++) Trig5[qtdum5] = (Int_t) (QT->TriggerPass)[qtdum5];

      dummyFirstInEvent=1;
    }
  }
}

// loops through edep matrices and sets nstb_ht,row_ht,col_ht to
// the highest tower corrdinates
void PostSim::GetHT(int i)
{
  TMatrix ** sedep = (TMatrix**) MatrixArray.At(i*numberMatrices+1);
  max_edep=0;
  nstb_ht=row_ht=col_ht=0;

  for(int n=0; n<4; n++)
  {
    for(int r=0; r<sedep[n]->GetNrows(); r++)
    {
      for(int c=0; c<sedep[n]->GetNcols(); c++)
      {
	max_edep_tmp = (*(sedep[n]))(r,c);
	if(max_edep_tmp > max_edep)
	{
	  max_edep = max_edep_tmp;
	  nstb_ht = n+1;
	  row_ht = r;
	  col_ht = c;
	};
      };
    };
  };

}
