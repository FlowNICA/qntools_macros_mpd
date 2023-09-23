TChain* makeChain(string& fileName, const char* treename) {
  cout << "Adding files to chain:" << endl;
  TChain *chain = new TChain(treename);
  if (fileName.rfind(".root") < fileName.size())
    chain->Add(fileName.data());
  else {
    TFileCollection fc("fc", "", fileName.c_str());
    chain->AddFileInfoList((TCollection*)fc.GetList());
  }
  chain->ls();
  return chain;
}

double getMass(int pdg)
{
  if(pdg<1000000) 
    return TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  else 
    return 0.931*(pdg/10%1000);
}

//--------------------------------------------------------------------------------------------------------------------------//
//this function simply connects the gain values read in to the BBC azimuthal distribution
//since tiles 7 and 9 (+ 13 and 15) share a gain value it is ambiguous how to assign the geometry here
//I prefer assigning the angle between the tiles thus "greying out" the adcs. 
//Others have assigned all of the adc to one (exclusive) or the the other. 
Float_t BBC_GetPhi(const Int_t eastWest, const Int_t tileId)
{
  //float GetPhiInBBC(int eastWest, int bbcN) { //tileId=0 to 23
  const float Pi = TMath::Pi();
  const float Phi_div = Pi / 6;
  float bbc_phi = Phi_div;
  switch(tileId) {
    case 0: bbc_phi = 3.*Phi_div;
  break;
    case 1: bbc_phi = Phi_div;
  break;
    case 2: bbc_phi = -1.*Phi_div;
  break;
    case 3: bbc_phi = -3.*Phi_div;
  break;
    case 4: bbc_phi = -5.*Phi_div;
  break;
    case 5: bbc_phi = 5.*Phi_div;
  break;
    //case 6: bbc_phi= (mRndm.Rndm() > 0.5) ? 2.*Phi_div:4.*Phi_div;	//tiles 7 and 9 are gained together we randomly assign the gain to one XOR the other
    case 6: bbc_phi = 3.*Phi_div;
  break;
    case 7: bbc_phi = 3.*Phi_div;
  break;
    case 8: bbc_phi = Phi_div;
  break;
    case 9: bbc_phi = 0.;
  break;
    case 10: bbc_phi = -1.*Phi_div;
  break;
    //case 11: bbc_phi = (mRndm.Rndm() > 0.5) ? -2.*Phi_div:-4.*Phi_div;	//tiles 13 and 15 are gained together
    case 11: bbc_phi = -3.*Phi_div;
  break;
    case 12: bbc_phi = -3.*Phi_div;
  break;
    case 13: bbc_phi = -5.*Phi_div;
  break;
    case 14: bbc_phi = Pi;
  break;
    case 15: bbc_phi = 5.*Phi_div;
  break;
  }

  //if we're looking at the east BBC we need to flip around x in the STAR coordinates, 
  //a line parallel to the beam would go through tile 1 on the W BBC and tile 3 on the 
  if(0 == eastWest){
    if (bbc_phi > -0.001){ //this is not a >= since we are talking about finite adcs -- not to important
      bbc_phi = Pi - bbc_phi;
    }
    else {
      bbc_phi= -Pi - bbc_phi;
    }
  }

  if(bbc_phi < 0.0) bbc_phi += 2.*Pi;
  if(bbc_phi > 2.*Pi) bbc_phi -= 2.*Pi;

  return bbc_phi;
}

Double_t GetZDCPosition(Int_t eastwest, Int_t verthori, Int_t strip)
// Get position of each slat;strip starts from 0
{

  std::vector<Double_t> zdcsmd_x = {0.5,2,3.5,5,6.5,8,9.5};
  std::vector<Double_t> zdcsmd_y = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

  Double_t mZDCSMDCenterex = 4.72466;
  Double_t mZDCSMDCenterey = 5.53629;
  Double_t mZDCSMDCenterwx = 4.39604;
  Double_t mZDCSMDCenterwy = 5.19968;

  if(eastwest==0 && verthori==0) return zdcsmd_x.at(strip)-mZDCSMDCenterex;
  if(eastwest==1 && verthori==0) return mZDCSMDCenterwx-zdcsmd_x.at(strip);
  if(eastwest==0 && verthori==1) return zdcsmd_y.at(strip)/sqrt(2.)-mZDCSMDCenterey;
  if(eastwest==1 && verthori==1) return zdcsmd_y.at(strip)/sqrt(2.)-mZDCSMDCenterwy;

  return -999.;
}

Double_t GetZDCPhi(Int_t eastwest, Int_t verthori, Int_t strip)
{
  double position = GetZDCPosition(eastwest, verthori, strip);
  if (position == -999.) return -999.;
  
  double cos, sin, phi;
  cos = -999.;
  sin = -999.;
  if(eastwest==0 && verthori==0)
  {
    cos = position;
    sin = 0.;
  }
  if(eastwest==1 && verthori==0)
  {
    cos = position;
    sin = 0.;
  }
  if(eastwest==0 && verthori==1)
  {
    cos = 0.;
    sin = position;
  }
  if(eastwest==1 && verthori==1)
  {
    cos = 0.;
    sin = position;
  }

  if (cos == -999.) return -999.;
  if (sin == -999.) return -999.;

  phi = atan2(sin,cos);

  return phi;
}

pair<double, double> GetXY(double _val_m2, double _val_sig, 
  double _fit_mean_m2_pi, double _fit_mean_m2_ka,
  double _fit_mean_sig_pi, double _fit_mean_sig_ka,
  double _fit_sig_m2_pi,
  double _fit_sig_sig_pi)
{
  double f = _fit_sig_sig_pi / _fit_sig_m2_pi;
  double alpha = -1.*atan2( (_fit_mean_m2_ka - _fit_mean_m2_pi), ((_fit_mean_sig_ka - _fit_mean_sig_pi)/f) );

  double x0 = (_val_sig - _fit_mean_sig_pi)/f;
  double y0 = _val_m2 - _fit_mean_m2_pi;

  double x = cos(alpha)*x0 - sin(alpha)*y0;
  double y = sin(alpha)*x0 + cos(alpha)*y0;
  return {x,y};
}

inline std::function<bool(double)> range(double lo, double hi) 
  {return [=](double x) { return lo <= x && x <= hi; };}

inline std::function<bool(double)> rangeStrict(double lo, double hi) 
  {return [=](double x) { return lo < x && x < hi; };}

inline std::function<bool(double)> equal(int val) 
  {return [=](double _val) { return TMath::Nint(_val) == val; };}

inline bool is (double x) 
  {return fabs(x - 1) < 1e-5;}

inline bool isNot (double x) 
  {return fabs(x) < 1e-5;}

inline bool positive (double x) 
  {return x > 0;}

inline bool negative (double x) 
  {return x < 0;}
