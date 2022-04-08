// Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include "algorithms/iuwt/iuwt_mask.h"

#include <boost/numeric/conversion/bounds.hpp>

#include "algorithms/iuwt/iuwt_decomposition.h"

namespace radler::algorithms::iuwt {

std::string IUWTMask::Summary(const IUWTDecomposition& iuwt) const {
  std::ostringstream str;
  str << "IUWTMask with " << _masks.size()
      << " scale masks (iuwt: " << iuwt.Summary() << ")\n";
  for (size_t i = 0; i != _masks.size(); ++i) {
    double maxVal = boost::numeric::bounds<double>::lowest();
    double minVal = std::numeric_limits<double>::max();
    size_t count = 0;
    for (size_t j = 0; j != _masks[i].size(); ++j) {
      if (_masks[i][j]) {
        ++count;
        if (iuwt[i][j] > maxVal) maxVal = iuwt[i][j];
        if (iuwt[i][j] < minVal) minVal = iuwt[i][j];
      }
    }
    if (maxVal == boost::numeric::bounds<double>::lowest()) {
      maxVal = std::numeric_limits<double>::quiet_NaN();
      minVal = std::numeric_limits<double>::quiet_NaN();
    }
    str << "Scale " << i << ": " << count << " (" << minVal << " - " << maxVal
        << ")\n";
  }
  return str.str();
}
}  // namespace radler::algorithms::iuwt
