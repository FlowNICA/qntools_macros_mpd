//#include "/mnt/pool/nica/7/parfenovpeter/Soft/QnTools/install/include/QnTools/QnDataFrame.hpp"
#include "QnDataFrame.hpp"
#include "utils.h"

vector <vector<string>> QvQv=
{
  {"Q_TPC_Tp_RESCALED", "Q_FHCal_F1_RESCALED"},
  {"Q_TPC_Tp_RESCALED", "Q_FHCal_F2_RESCALED"},
  {"Q_TPC_Tp_RESCALED", "Q_FHCal_F3_RESCALED"},
  {"Q_FHCal_F1_RESCALED", "Q_FHCal_F3_RESCALED"},
  {"Q_TPC_Tpi_RESCALED", "Q_FHCal_F1_RESCALED"},
  {"Q_TPC_Tpi_RESCALED", "Q_FHCal_F2_RESCALED"},
  {"Q_TPC_Tpi_RESCALED", "Q_FHCal_F3_RESCALED"}
};

vector <vector<string>> uvQv=
{
  {"u_TPC_proton_RESCALED", "Q_FHCal_F1_RESCALED"},
  {"u_TPC_proton_RESCALED", "Q_FHCal_F2_RESCALED"},
  {"u_TPC_proton_RESCALED", "Q_FHCal_F3_RESCALED"},
  {"u_TPC_pip_RESCALED", "Q_FHCal_F1_RESCALED"},
  {"u_TPC_pip_RESCALED", "Q_FHCal_F2_RESCALED"},
  {"u_TPC_pip_RESCALED", "Q_FHCal_F3_RESCALED"},
  {"u_TPC_pim_RESCALED", "Q_FHCal_F1_RESCALED"},
  {"u_TPC_pim_RESCALED", "Q_FHCal_F2_RESCALED"},
  {"u_TPC_pim_RESCALED", "Q_FHCal_F3_RESCALED"},
  {"u_sim_proton_PLAIN", "Q_RP_PLAIN"},
  {"u_sim_pip_PLAIN", "Q_RP_PLAIN"},
  {"u_sim_pim_PLAIN", "Q_RP_PLAIN"}
};

vector <vector<string>> uvQvQv=
{
  {"u_TPC_proton_RESCALED", "Q_FHCal_F1_RESCALED", "Q_FHCal_F3_RESCALED"},
  {"u_TPC_pip_RESCALED", "Q_FHCal_F1_RESCALED", "Q_FHCal_F3_RESCALED"},
  {"u_TPC_pim_RESCALED", "Q_FHCal_F1_RESCALED", "Q_FHCal_F3_RESCALED"},
  {"u_sim_proton_PLAIN", "Q_RP_PLAIN", "Q_RP_PLAIN"},
  {"u_sim_pip_PLAIN", "Q_RP_PLAIN", "Q_RP_PLAIN"},
  {"u_sim_pim_PLAIN", "Q_RP_PLAIN", "Q_RP_PLAIN"}
};

void correlate(string inputFiles="qn.root", string outputFile="corr.root")
{
	TStopwatch timer;
	timer.Start();
  int nSamples = 50;
  Qn::AxisD centAxis({"evBnorm", 10, 0., 2.});
  auto axes_correlation = Qn::MakeAxes(centAxis);//, vtxAxis);
  TChain *c=makeChain(inputFiles, "tree");
  ROOT::RDataFrame d(*c);
  auto d_samples = Qn::Correlation::Resample(d, nSamples);

  namespace P2 = Qn::Correlation::TwoParticle;
  namespace P3 = Qn::Correlation::MixedHarmonics;
  namespace P4 = Qn::Correlation::FourParticle;
  auto wn = Qn::Correlation::UseWeights::No;
  auto wy = Qn::Correlation::UseWeights::Yes;
  auto wUnity = [](const Qn::QVector &a, const Qn::QVector &b) { return 1; };
  auto wUnity1P = [](const Qn::QVector &a) { return 1; };
  auto wSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights(); };
  auto wSumWu1P = [](const Qn::QVector &a) { return a.sumweights(); };
  auto wSumWu3P = [](const Qn::QVector &a, const Qn::QVector &b, const Qn::QVector &c) { return a.sumweights(); };
  auto wSumWuEP = [](const Qn::QVector &a, const Qn::QVector &b) { return a.sumweights()/b.sumweights(); };
  auto wDenomASumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return 1./a.sumweights(); };
  auto wDenomBSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return 1./b.sumweights(); };
  auto wDenomABSumWu = [](const Qn::QVector &a, const Qn::QVector &b) { return 1./(a.sumweights()*b.sumweights()); };

  auto corrBuilder =
    Qn::Correlation::CorrelationBuilder{&d_samples, nSamples, axes_correlation};
  for (auto &corr:QvQv)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1);

    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X1X1_", P2::xx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y1Y1_", P2::yy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X1Y1_", P2::xy(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y1X1_", P2::yx(1, 1), wUnity, wn, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP1_", P2::ScalarProduct(1, 1), wUnity, wn, qn, qn);  
  }
  for (auto &corr:uvQv)
  {
    std::array<std::string, 2> qn{corr.at(0), corr.at(1)};
    string corrName=corr.at(0)+"_"+corr.at(1); 

    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X1X1_", P2::xx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y1Y1_", P2::yy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X1Y1_", P2::xy(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y1X1_", P2::yx(1, 1), wSumWu, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP1_", P2::ScalarProduct(1, 1), wSumWu, wy, qn, qn);
  }

  for (auto &corr:uvQvQv)
  {
    std::array<std::string, 3> qn{corr.at(0), corr.at(1), corr.at(2)};
    string corrName=corr.at(0)+"_"+corr.at(1)+"_"+corr.at(2);

    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2X1X1_", P3::xxx(2, 1, 1), wSumWu3P, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2Y1Y1_", P3::xyy(2, 1, 1), wSumWu3P, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2X1Y1_", P3::yxy(2, 1, 1), wSumWu3P, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2Y1X1_", P3::yyx(2, 1, 1), wSumWu3P, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2Y1Y1_", P3::yyy(2, 1, 1), wSumWu3P, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2Y1X1_", P3::xyx(2, 1, 1), wSumWu3P, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_Y2X1X1_", P3::yxx(2, 1, 1), wSumWu3P, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_X2X1Y1_", P3::xxy(2, 1, 1), wSumWu3P, wy, qn, qn);
    corrBuilder.AddCorrelationWithInternalReader(corrName+"_SP2_", P3::ScalarProduct(2, 1, 1), wSumWu3P, wy, qn, qn);
  }

  // ---------------- //
  // saving to output //
  // ---------------- //
  auto corrFile = TFile::Open(outputFile.c_str(), "RECREATE");
  corrFile->cd();
  auto results = corrBuilder.GetResults();
  for (auto &res : results) {
    res->Write();
  }
  corrFile->Close();
	timer.Stop();
	timer.Print();
}
