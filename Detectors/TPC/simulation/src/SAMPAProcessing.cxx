// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file SAMPAProcessing.cxx
/// \brief Implementation of the SAMPA response
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include "TPCSimulation/SAMPAProcessing.h"
#include "TPCBase/CDBInterface.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "FairLogger.h"

using namespace o2::TPC;

SAMPAProcessing::SAMPAProcessing() : mSaturationSpline(), mRandomNoiseRing(RandomRing::RandomType::Gaus)
{
  updateParameters();
  importSaturationCurve("SAMPA_saturation.dat");
}

SAMPAProcessing::~SAMPAProcessing() = default;

void SAMPAProcessing::updateParameters()
{
  auto& cdb = CDBInterface::instance();
  mGasParam = &(cdb.getParameterGas());
  mDetParam = &(cdb.getParameterDetector());
  mEleParam = &(cdb.getParameterElectronics());
  mPedestalMap = &(cdb.getPedestals());
  mNoiseMap = &(cdb.getNoise());
}

bool SAMPAProcessing::importSaturationCurve(std::string file)
{
  std::string inputDir;
  const char* aliceO2env = std::getenv("O2_ROOT");
  if (aliceO2env) {
    inputDir = aliceO2env;
  }
  inputDir += "/share/Detectors/TPC/files/";

  std::ifstream saturationFile(inputDir + file, std::ifstream::in);
  if (!saturationFile) {
    LOG(FATAL) << "TPC::SAMPAProcessing - Input file '" << inputDir + file
               << "' does not exist! No SAMPA saturation curve loaded!" << FairLogger::endl;
    return false;
  }
  std::vector<std::pair<float, float>> saturation;
  for (std::string line; std::getline(saturationFile, line);) {
    float x, y;
    std::istringstream is(line);
    while (is >> x >> y) {
      saturation.emplace_back(x, y);
    }
  }
  const int nPoints = saturation.size();
  double xSat[nPoints], ySat[nPoints];
  for (int i = 0; i < nPoints; ++i) {
    xSat[i] = saturation[i].first;
    ySat[i] = saturation[i].second;
  }

  mSaturationSpline = std::unique_ptr<TSpline3>(new TSpline3("SAMPA saturation", xSat, ySat, nPoints, "b2e2", 0, 0));
  return true;
}

void SAMPAProcessing::getShapedSignal(float ADCsignal, float driftTime, std::vector<float>& signalArray) const
{
  const float timeBinTime = getTimeBinTime(driftTime);
  const float offset = driftTime - timeBinTime;
  for (float bin = 0; bin < mEleParam->getNShapedPoints(); bin += Vc::float_v::Size) {
    Vc::float_v binvector;
    for (int i = 0; i < Vc::float_v::Size; ++i) {
      binvector[i] = bin + i;
    }
    Vc::float_v time = timeBinTime + binvector * mEleParam->getZBinWidth();
    Vc::float_v signal = getGamma4(time, Vc::float_v(timeBinTime + offset), Vc::float_v(ADCsignal));
    for (int i = 0; i < Vc::float_v::Size; ++i) {
      signalArray[bin + i] = signal[i];
    }
  }
}
