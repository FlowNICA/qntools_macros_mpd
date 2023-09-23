#include "utils.h"
#include "makeQvectors.h"
#include <algorithm>
#include <random>

filteredDF defineVariables(definedDF &d, string collision_energy, string str_nucleus_mass);
void setupQvectors();

void makeQvectors(string inputFiles = "", string calibFilePath = "qa.root", string outFilePath = "qn.root",
                  string collision_energy = "0.", string str_nucleus_mass = "209")
{
  TStopwatch timer;
  timer.Start();

  ROOT::RDataFrame d(*makeChain(inputFiles, "t"));
  auto dd = defineVariables(d, collision_energy, str_nucleus_mass);
  init(dd, outFilePath, calibFilePath);
  setupQvectors();
  timer.Stop();
  timer.Print();
  timer.Start();
  run(dd);
  cout << "Done!\n";
  timer.Stop();
  timer.Print();
}

filteredDF defineVariables(definedDF &d, string collision_energy, string str_nucleus_mass)
{
  const double T = std::stod( collision_energy );
  const double M = 0.938;
  const double E = T + M;
  const double P = sqrt( E*E - M*M );
  const double Y_BEAM = 0.25 * log( (E + P) / (E - P) );
  const double nucleus_mass = std::stod(str_nucleus_mass);
  const double NUCLEUS_RADIUS = 1.25 * pow( nucleus_mass, 1.0 / 3.0 );

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Define function to get p from px, py, pz
  auto getP = [](const RVec<fourVector> &tracks)
  {
    RVecD mom;
    for (auto &track:tracks)
      mom.push_back(sqrt(track.Pt()*track.Pt()+track.Pz()*track.Pz()));
    return mom;
  };

  auto getNhits = [](const RVecI &hits)
  {
    RVecF nhits;
    for (auto &hit:hits)
      nhits.push_back((float)hit);
    return nhits;
  };

  auto getPdg = [](const RVecI &pids)
  {
    RVecF pdgs;
    for (auto &pid:pids)
      pdgs.push_back((float)pid);
    return pdgs;
  };

  auto getCharge = [](const RVecS &charge)
  {
    RVecF ch;
    for (auto &q:charge)
      ch.push_back((float)q);
    return ch;
  };

  auto recoRapidityCM = [Y_BEAM](const RVec<fourVector> &tracks, const RVecI &pdg_codes)
  {
    RVecD yCM;
    for (int i=0; i<(int)tracks.size(); i++) {
      double mass = getMass(pdg_codes.at(i));
      double ene = sqrt(mass*mass + tracks.at(i).P()*tracks.at(i).P());
      yCM.push_back( 0.5*log((ene+tracks.at(i).Pz())/(ene-tracks.at(i).Pz())) - Y_BEAM );
    }
    return yCM;
  };

  auto simRapidityCM = [Y_BEAM](const RVec<fourVector> &tracks)
  {
    RVecD yCM;
    for (auto &track:tracks) {
      yCM.push_back( 0.5*log((track.E()+track.Pz())/(track.E()-track.Pz())) - Y_BEAM );
    }
    return yCM;
  };

  auto getFhcalModId = [](const RVecI &modIds)
  {
    RVecD modules;
    for (auto &id:modIds)
      modules.push_back((double)id);
    return modules;
  };
  
  auto fhcalModPhi = [](const RVec<XYZVector> &modules)
  {
    RVecF phi;
    for (auto &mod:modules)
      phi.push_back(mod.Phi());
    return phi;
  };
  /////////////////////////////////////////////////////////////////////////////////////////////////

  auto dd = d.Define("evB", "mcB")
    .Define("evBnorm", [NUCLEUS_RADIUS](double b){ return b/NUCLEUS_RADIUS; }, {"evB"})
    .Define("evPhiRP", "mcRP")
    .Define("evVtxX", "mcVtxX")
    .Define("evVtxY", "mcVtxY")
    .Define("evVtxZ", "mcVtxZ")
    .Define("trPt", "recoGlobalMom.fCoordinates.fPt")
    .Define("trP", getP, {"recoGlobalMom"})
    .Define("trEta", "recoGlobalMom.fCoordinates.fEta")
    .Define("trYcm", recoRapidityCM, {"recoGlobalMom", "recoGlobalSimPdg"})
    .Define("trPhi", "recoGlobalMom.fCoordinates.fPhi")
    .Define("trCh", getCharge, {"recoGlobalCharge"})
    .Define("trNhits", getNhits, {"recoGlobalNhits"})
    .Define("trMotherId", "recoGlobalSimMotherId")
    .Define("trIsTof", "recoGlobalTofFlag")
    .Define("trM2", "recoGlobalTofMass2")
    .Define("trPID", getPdg, {"recoGlobalSimPdg"})
    .Define("simPt", "simMom.fCoordinates.fPt")
    .Define("simP", getP, {"simMom"})
    .Define("simEta", "simMom.fCoordinates.fEta")
    .Define("simYcm", simRapidityCM, {"simMom"})
    .Define("simPhi", "simMom.fCoordinates.fPhi")
    .Define("simCh", getCharge, {"simCharge"})
    .Define("simPID", getPdg, {"simPdg"})
    .Define("simHasFHCalHit", "simHasHitFHCal")
    .Define("fhcalModID", getFhcalModId, {"fhcalModId"})
    .Define("fhcalModPhi", fhcalModPhi, {"fhcalModPos"})
    .Define("fhcalModPosX", "fhcalModPos.fCoordinates.fX")
    .Define("fhcalModPosY", "fhcalModPos.fCoordinates.fY")
    .Filter("evB>=0.")
    ;

  varPatterns =
    {
      "ev(B|Bnorm|PhiRP)", // kEvent
      "fhcalMod(Id|ID|E|Phi|PosX|PosY)", // kChannel
      "tr(P|Pt|Eta|Ycm|Phi|Ch|Nhits|MotherId|M2|PID)", // kRecParticle
      "sim(P|Pt|Eta|Ycm|Phi|Ch|MotherId|PID|HasFHCalHit)" // kSimParticle
    };

  return dd;
}

void setupQvectors()
{
  std::vector<int> f1_modules = { 14, 15, 16, 21, 23, 28, 29, 30 };
  std::vector<int> f2_modules = { 6, 7, 8, 9, 10, 13, 17, 20, 24, 27, 31, 34, 35, 36, 37, 38 };
  std::vector<int> f3_modules = { 0, 1, 2, 3, 4, 5, 11, 12, 18, 19, 25, 26, 32, 33, 39, 40, 41, 42, 43, 44 };
  
  vector<Qn::AxisD> corrAxesEvent =
    {
      {"evBnorm", 10, 0., 2.}
    };

  vector<Qn::AxisD> corrAxesParticle =
    {
      {"trPt", 30, 0., 3.},
      {"trYcm", 15, -1.5, 1.5},

    };

  vector<Qn::AxisD> corrAxesSimParticle =
    {
      {"simPt", 30, 0., 3.},
      {"simYcm", 15, -1.5, 1.5},

    };

  for (auto &axis : corrAxesEvent)
    man.AddCorrectionAxis(axis);

  Qn::Recentering recentering;
  recentering.SetApplyWidthEqualization(false);
  Qn::TwistAndRescale twistRescale;
  twistRescale.SetApplyRescale(true);
  twistRescale.SetTwistAndRescaleMethod(Qn::TwistAndRescale::Method::DOUBLE_HARMONIC);

  auto sumW = Qn::QVector::Normalization::M;
  auto track = Qn::DetectorType::TRACK;
  auto channel = Qn::DetectorType::CHANNEL;
  auto plain = Qn::QVector::CorrectionStep::PLAIN;
  auto recentered = Qn::QVector::CorrectionStep::RECENTERED;
  auto twisted = Qn::QVector::CorrectionStep::TWIST;
  auto rescaled = Qn::QVector::CorrectionStep::RESCALED;

  std::string name;

  name = "u_TPC_proton";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1, 2}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trNhits"}, [](float nhits){ return (nhits > 16); }, "nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trPID"}, [](float pid){ auto pid_code = std::round(pid); return (pid_code == 2212); }, "pid_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"}, [](float charge){ return (charge>0); }, "charge_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto2D(name.c_str(), {{"trYcm", 400, -2., 2.}, {"trPt", 300, 0.0, 3.0}}, "Ones");

  name = "u_TPC_pip";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1, 2}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trNhits"}, [](float nhits){ return (nhits > 16); }, "nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trPID"}, [](float pid){ auto pid_code = std::round(pid); return (pid_code == 211); }, "pid_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"}, [](float charge){ return (charge>0); }, "charge_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto2D(name.c_str(), {{"trYcm", 400, -2., 2.}, {"trPt", 300, 0.0, 3.0}}, "Ones");

  name = "u_TPC_pim";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", corrAxesParticle, {1, 2}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trNhits"}, [](float nhits){ return (nhits > 16); }, "nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trPID"}, [](float pid){ auto pid_code = std::round(pid); return (pid_code == -211); }, "pid_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"}, [](float charge){ return (charge<0); }, "charge_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto2D(name.c_str(), {{"trYcm", 400, -2., 2.}, {"trPt", 300, 0.0, 3.0}}, "Ones");

  name = "Q_TPC_Tp";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1, 2}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trNhits"}, [](float nhits){ return (nhits > 16); }, "nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trPID"}, [](float pid){ auto pid_code = std::round(pid); return (pid_code == 2212); }, "pid_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"}, [](float charge){ return (charge>0); }, "charge_cut");
  man.AddCutOnDetector(name.c_str(), {"trYcm"}, [](float ycm){ return (ycm > -1.2 && ycm < -0.6); }, "yCM_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt > 0.2 && pt < 2.); }, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trYcm", 400, -2., 2.}, {"trPt", 300, 0.0, 3.0}}, "Ones");

  name = "Q_TPC_Tpi";
  man.AddDetector(name.c_str(), track, "trPhi", "Ones", {}, {1, 2}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kRecParticle), "recParticle");
  man.AddCutOnDetector(name.c_str(), {"trNhits"}, [](float nhits){ return (nhits > 16); }, "nhits_cut");
  man.AddCutOnDetector(name.c_str(), {"trPID"}, [](float pid){ auto pid_code = std::round(pid); return (pid_code == -211); }, "pid_cut");
  man.AddCutOnDetector(name.c_str(), {"trCh"}, [](float charge){ return (charge<0); }, "charge_cut");
  man.AddCutOnDetector(name.c_str(), {"trYcm"}, [](float ycm){ return (ycm > -1.5 && ycm < -0.2); }, "yCM_cut");
  man.AddCutOnDetector(name.c_str(), {"trPt"}, [](float pt){ return (pt > 0. && pt < 0.35); }, "pt_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"trPhi", 100, -3.15, 3.15}, "Ones");
  man.AddHisto2D(name.c_str(), {{"trYcm", 400, -2., 2.}, {"trPt", 300, 0.0, 3.0}}, "Ones");

  name = "Q_FHCal_F1";
  man.AddDetector(name.c_str(), channel, "fhcalModPhi", "fhcalModE", {}, {1}, sumW);
  man.AddCutOnDetector(name.c_str(), {"fhcalModID"}, [f1_modules](double mod_id){ auto id = std::round(mod_id); return (std::find(f1_modules.begin(), f1_modules.end(), id) != f1_modules.end()); }, "Subevent_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"fhcalModId", 100, 0., 100}, "fhcalModE");
  man.AddHisto1D(name.c_str(), {"fhcalModPhi", 100, -3.15, 3.15}, "fhcalModE");
  man.AddHisto2D(name.c_str(), {{"fhcalModPosX", 100, -100., 100.}, {"fhcalModPosY", 100, -100., 100.}}, "fhcalModE");

  name = "Q_FHCal_F2";
  man.AddDetector(name.c_str(), channel, "fhcalModPhi", "fhcalModE", {}, {1}, sumW);
    man.AddCutOnDetector(name.c_str(), {"fhcalModID"}, [f2_modules](double mod_id){ auto id = std::round(mod_id); return (std::find(f2_modules.begin(), f2_modules.end(), id) != f2_modules.end()); }, "Subevent_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"fhcalModId", 100, 0., 100}, "fhcalModE");
  man.AddHisto1D(name.c_str(), {"fhcalModPhi", 100, -3.15, 3.15}, "fhcalModE");
  man.AddHisto2D(name.c_str(), {{"fhcalModPosX", 100, -100., 100.}, {"fhcalModPosY", 100, -100., 100.}}, "fhcalModE");

  name = "Q_FHCal_F3";
  man.AddDetector(name.c_str(), channel, "fhcalModPhi", "fhcalModE", {}, {1}, sumW);
  man.AddCutOnDetector(name.c_str(), {"fhcalModID"}, [f3_modules](double mod_id){ auto id = std::round(mod_id); return (std::find(f3_modules.begin(), f3_modules.end(), id) != f3_modules.end()); }, "Subevent_cut");
  man.AddCorrectionOnQnVector(name.c_str(), recentering);
  man.AddCorrectionOnQnVector(name.c_str(), twistRescale);
  man.SetOutputQVectors(name.c_str(), {plain, recentered, twisted, rescaled});
  man.AddHisto1D(name.c_str(), {"fhcalModId", 100, 0., 100}, "fhcalModE");
  man.AddHisto1D(name.c_str(), {"fhcalModPhi", 100, -3.15, 3.15}, "fhcalModE");
  man.AddHisto2D(name.c_str(), {{"fhcalModPosX", 100, -100., 100.}, {"fhcalModPosY", 100, -100., 100.}}, "fhcalModE");

  name = "u_sim_proton";
  man.AddDetector(name.c_str(), track, "simPhi", "Ones", corrAxesSimParticle, {1, 2}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"simPID"}, [](float pid){ auto pid_code = std::round(pid); return (pid_code == 2212); }, "pid_cut");
  man.AddCutOnDetector(name.c_str(), {"simCh"}, [](float charge){ return (charge>0); }, "charge_cut");
  man.SetOutputQVectors(name.c_str(), {plain});
  man.AddHisto2D(name.c_str(), {{"simYcm", 400, -2., 2.}, {"simPt", 300, 0.0, 3.0}}, "Ones");

  name = "u_sim_pip";
  man.AddDetector(name.c_str(), track, "simPhi", "Ones", corrAxesSimParticle, {1, 2}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"simPID"}, [](float pid){ auto pid_code = std::round(pid); return (pid_code == 211); }, "pid_cut");
  man.AddCutOnDetector(name.c_str(), {"simCh"}, [](float charge){ return (charge>0); }, "charge_cut");
  man.SetOutputQVectors(name.c_str(), {plain});
  man.AddHisto2D(name.c_str(), {{"simYcm", 400, -2., 2.}, {"simPt", 300, 0.0, 3.0}}, "Ones");

  name = "u_sim_pim";
  man.AddDetector(name.c_str(), track, "simPhi", "Ones", corrAxesSimParticle, {1, 2}, sumW);
  man.AddCutOnDetector(name.c_str(), {"particleType"}, equal(kSimParticle), "simParticle");
  man.AddCutOnDetector(name.c_str(), {"simPID"}, [](float pid){ auto pid_code = std::round(pid); return (pid_code == -211); }, "pid_cut");
  man.AddCutOnDetector(name.c_str(), {"simCh"}, [](float charge){ return (charge<0); }, "charge_cut");
  man.SetOutputQVectors(name.c_str(), {plain});
  man.AddHisto2D(name.c_str(), {{"simYcm", 400, -2., 2.}, {"simPt", 300, 0.0, 3.0}}, "Ones");

  name = "Q_RP";
  man.AddDetector(name.c_str(), track, "evPhiRP", "Ones", {}, {1, 2}, sumW);
  man.SetOutputQVectors(name.c_str(), {plain});
  man.AddHisto1D(name.c_str(), {"evPhiRP", 100, 0., 6.3}, "Ones");
}
