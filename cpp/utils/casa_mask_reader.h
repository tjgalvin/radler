// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_UTILS_CASA_MASK_READER_H_
#define RADLER_UTILS_CASA_MASK_READER_H_

#include <string>

namespace radler::utils {
class CasaMaskReader {
 public:
  explicit CasaMaskReader(const std::string& path);

  void Read(bool* mask);

  size_t Width() const { return _width; }
  size_t Height() const { return _height; }
  size_t NPolarizations() const { return _nPolarizations; }
  size_t NChannels() const { return _nChannels; }

 private:
  std::string _path;
  size_t _width, _height, _nPolarizations, _nChannels;
};
}  // namespace radler::utils
#endif
