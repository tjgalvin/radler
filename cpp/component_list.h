// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_COMPONENT_LIST_H_
#define RADLER_COMPONENT_LIST_H_

#include <vector>

#include <aocommon/image.h>
#include <aocommon/uvector.h>

#include <schaapcommon/fitters/spectralfitter.h>

#include "image_set.h"

namespace radler {
class Radler;  // Forward declared to avoid a circular dependency.

namespace algorithms {
// Forward declared since the class isn't part of Radler's public interface.
class DeconvolutionAlgorithm;
class MultiScaleAlgorithm;
}  // namespace algorithms

class ComponentList {
 public:
  ComponentList()
      : width_(0),
        height_(0),
        n_frequencies_(0),
        components_added_since_last_merge_(0),
        max_components_before_merge_(0),
        list_per_scale_() {}

  /**
   * Constructor for single-scale clean
   */
  ComponentList(size_t width, size_t height, ImageSet& image_set)
      : width_(width),
        height_(height),
        n_frequencies_(image_set.Size()),
        components_added_since_last_merge_(0),
        max_components_before_merge_(100000),
        list_per_scale_(1) {
    LoadFromImageSet(image_set, 0);
  }

  /**
   * Constructor for multi-scale clean
   */
  ComponentList(size_t width, size_t height, size_t n_scales,
                size_t n_frequencies)
      : width_(width),
        height_(height),
        n_frequencies_(n_frequencies),
        components_added_since_last_merge_(0),
        max_components_before_merge_(100000),
        list_per_scale_(n_scales) {}

  ComponentList(const ComponentList&) = default;
  ComponentList(ComponentList&&) = default;
  ComponentList& operator=(const ComponentList&) = default;
  ComponentList& operator=(ComponentList&&) = default;

  struct Position {
    Position(size_t x, size_t y) : x(x), y(y) {}
    size_t x, y;
  };

  void Add(size_t x, size_t y, size_t scale_index, const float* values) {
    list_per_scale_[scale_index].values.push_back(values,
                                                  values + n_frequencies_);
    list_per_scale_[scale_index].positions.emplace_back(x, y);
    ++components_added_since_last_merge_;
    if (components_added_since_last_merge_ >= max_components_before_merge_)
      MergeDuplicates();
  }

  void Add(const ComponentList& other, int offset_x, int offset_y) {
    assert(other.n_frequencies_ == n_frequencies_);
    if (other.NScales() > NScales()) SetNScales(other.NScales());
    for (size_t scale = 0; scale != other.NScales(); ++scale) {
      const ScaleList& list = other.list_per_scale_[scale];
      for (size_t i = 0; i != list.positions.size(); ++i) {
        Add(list.positions[i].x + offset_x, list.positions[i].y + offset_y,
            scale, &list.values[i * n_frequencies_]);
      }
    }
  }

  void WriteSources(const Radler& radler, const std::string& filename,
                    long double pixel_scale_x, long double pixel_scale_y,
                    long double phase_centre_ra, long double phase_centre_dec,
                    long double l_shift, long double m_shift) const;

  /**
   * @brief Write component lists over all scales, typically
   * used for writing components of a multiscale clean.
   */
  void Write(const std::string& filename,
             const algorithms::MultiScaleAlgorithm& multiscale,
             long double pixel_scale_x, long double pixel_scale_y,
             long double phase_centre_ra, long double phase_centre_dec,
             long double l_shift = 0.0, long double m_shift = 0.0) const;

  void WriteSingleScale(const std::string& filename,
                        const algorithms::DeconvolutionAlgorithm& algorithm,
                        long double pixel_scale_x, long double pixel_scale_y,
                        long double phase_centre_ra,
                        long double phase_centre_dec, long double l_shift,
                        long double m_shift) const;

  void MergeDuplicates() {
    if (components_added_since_last_merge_ != 0) {
      for (size_t scale_index = 0; scale_index != list_per_scale_.size();
           ++scale_index) {
        MergeDuplicates(scale_index);
      }
      components_added_since_last_merge_ = 0;
    }
  }

  void Clear() {
    for (ScaleList& list : list_per_scale_) {
      list.positions.clear();
      list.values.clear();
    }
  }

  size_t Width() const { return width_; }
  size_t Height() const { return height_; }

  size_t ComponentCount(size_t scale_index) const {
    return list_per_scale_[scale_index].positions.size();
  }

  void GetComponent(size_t scale_index, size_t index, size_t& x, size_t& y,
                    float* values) const {
    assert(scale_index < list_per_scale_.size());
    assert(index < list_per_scale_[scale_index].positions.size());
    x = list_per_scale_[scale_index].positions[index].x;
    y = list_per_scale_[scale_index].positions[index].y;
    for (size_t f = 0; f != n_frequencies_; ++f)
      values[f] =
          list_per_scale_[scale_index].values[index * n_frequencies_ + f];
  }

  /**
   * @brief Multiply the components for a given scale index, position index and
   * channel index with corresponding (primary beam) correction factors.
   */
  void MultiplyScaleComponent(size_t scale_index, size_t position_index,
                              size_t channel, double correction_factor) {
    assert(scale_index < list_per_scale_.size());
    assert(position_index < list_per_scale_[scale_index].positions.size());
    assert(channel < n_frequencies_);
    float& value = list_per_scale_[scale_index]
                       .values[channel + position_index * n_frequencies_];
    value *= correction_factor;
  }

  /**
   * @brief Get vector of positions per scale index.
   */
  const aocommon::UVector<Position>& GetPositions(size_t scale_index) const {
    assert(scale_index < list_per_scale_.size());
    return list_per_scale_[scale_index].positions;
  }

  size_t NScales() const { return list_per_scale_.size(); }

  size_t NFrequencies() const { return n_frequencies_; }

  void SetNScales(size_t n_scales) { list_per_scale_.resize(n_scales); }

 private:
  struct ScaleList {
    /**
     * This list contains nFrequencies values for each
     * component, such that _positions[i] corresponds with the values
     * starting at _values[i * n_frequencies_].
     */
    aocommon::UVector<float> values;
    aocommon::UVector<Position> positions;
  };

  void Write(const std::string& filename,
             const schaapcommon::fitters::SpectralFitter& fitter,
             const aocommon::UVector<double>& scale_sizes,
             long double pixel_scale_x, long double pixel_scale_y,
             long double phase_centre_ra, long double phase_centre_dec,
             long double l_shift, long double m_shift) const;

  void LoadFromImageSet(ImageSet& image_set, size_t scale_index);

  void MergeDuplicates(size_t scale_index) {
    ScaleList& list = list_per_scale_[scale_index];
    aocommon::UVector<float> new_values;
    aocommon::UVector<Position> new_positions;

    std::vector<aocommon::Image> images(n_frequencies_);
    for (aocommon::Image& image : images)
      image = aocommon::Image(width_, height_, 0.0);
    size_t value_index = 0;
    for (size_t index = 0; index != list.positions.size(); ++index) {
      size_t position =
          list.positions[index].x + list.positions[index].y * width_;
      for (size_t frequency = 0; frequency != n_frequencies_; ++frequency) {
        images[frequency][position] += list.values[value_index];
        value_index++;
      }
    }

    list.values.clear();
    list.positions.clear();

    for (size_t image_index = 0; image_index != images.size(); ++image_index) {
      aocommon::Image& image = images[image_index];
      size_t pos_index = 0;
      for (size_t y = 0; y != height_; ++y) {
        for (size_t x = 0; x != width_; ++x) {
          if (image[pos_index] != 0.0) {
            for (size_t i = 0; i != images.size(); ++i) {
              new_values.push_back(images[i][pos_index]);
              images[i][pos_index] = 0.0;
            }
            new_positions.emplace_back(x, y);
          }
          ++pos_index;
        }
      }
    }
    std::swap(list_per_scale_[scale_index].values, new_values);
    std::swap(list_per_scale_[scale_index].positions, new_positions);
  }

  size_t width_;
  size_t height_;
  size_t n_frequencies_;
  size_t components_added_since_last_merge_;
  size_t max_components_before_merge_;
  std::vector<ScaleList> list_per_scale_;
};
}  // namespace radler
#endif
