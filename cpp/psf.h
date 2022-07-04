// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_PSF_H_
#define RADLER_PSF_H_

#include <iomanip>
#include <ostream>

#include <aocommon/imageaccessor.h>

namespace radler {

class Psf {
 public:
  explicit Psf(std::unique_ptr<aocommon::ImageAccessor> accessor)
      : accessor(std::move(accessor)) {}

  explicit Psf(int x, int y, std::unique_ptr<aocommon::ImageAccessor> accessor)
      : x(x), y(y), accessor(std::move(accessor)) {}

  Psf() = default;
  ~Psf() = default;

  Psf(const Psf&) = delete;
  Psf& operator=(const Psf&) = delete;

  Psf(Psf&&) = default;
  Psf& operator=(Psf&&) = default;

  /** The x-offset in pixels from the corner position. */
  int x{0};

  /** The y-offset in pixels from the corner position. */
  int y{0};

  /** The PSF image. */
  std::unique_ptr<aocommon::ImageAccessor> accessor{nullptr};

  friend std::ostream& operator<<(std::ostream& output, const Psf& psf) {
    return output << "[x: " << psf.x << ", y: " << psf.y
                  << ", accessor: " << (psf.accessor ? "set" : "nullptr")
                  << ']';
  }
};

}  // namespace radler

#endif  // RADLER_PSF_H_
