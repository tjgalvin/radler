// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_MATH_DIJKSTRASPLITTER_H_
#define RADLER_MATH_DIJKSTRASPLITTER_H_

#include <algorithm>

namespace radler::math {

class DijkstraSplitter {
 public:
  DijkstraSplitter(std::size_t width, std::size_t height)
      : width_(width), height_(height) {}

  struct Coord {
    Coord() = default;
    Coord(std::size_t _x, std::size_t _y) : x(_x), y(_y) {}
    std::size_t x = 0;
    std::size_t y = 0;
  };

  struct Visit {
    float distance;
    Coord to;
    Coord from;
    bool operator<(const Visit& rhs) const { return distance > rhs.distance; }
  };

  void AddVerticalDivider(const float* image, float* scratch, float* output,
                          std::size_t x1, std::size_t x2) const;

  /**
   * Find the shortest vertical path through an image. The path is
   * constrained to lie between vertical lines given by x1 and x2.
   * The output is set to 1 for pixels that are part of the path, and
   * set to 0 otherwise. The reason it's a floating point is because
   * it is also used as scratch.
   */
  void DivideVertically(const float* image, float* output, std::size_t x1,
                        std::size_t x2) const;

  void AddHorizontalDivider(const float* image, float* scratch, float* output,
                            std::size_t y1, std::size_t y2) const;

  /**
   * Like DivideVertically, but for horizontal paths. The path is
   * constrained to lie between horizontal lines given by y1 and y2.
   */
  void DivideHorizontally(const float* image, float* output, std::size_t y1,
                          std::size_t y2) const;

  /**
   * Mask the space between (typically) two vertical divisions.
   * @param subdivision An image that is the result of earlier calls
   * to DivideVertically().
   * @param subimage_x An x-position that is in between the two splits.
   * @param mask A mask image for which pixels will be set to true if
   *   and only if they are part of the area specified by the
   *   two divisions.
   * @param [out] x The left side of the bounding box of the divisions.
   * @param [out] subwidth The bounding width of the two divisions.
   */
  void FloodVerticalArea(const float* subdivision, std::size_t subimage_x,
                         bool* mask, std::size_t& x,
                         std::size_t& subwidth) const;

  /**
   * Like FloodVerticalArea(), but for horizontal flooding.
   */
  void FloodHorizontalArea(const float* subdivision, std::size_t subimage_y,
                           bool* mask, std::size_t& y,
                           std::size_t& subheight) const;

  /**
   * Combines a horizontally and vertically filled area and extracts a
   * single mask where the areas overlap.
   * @param vertical_mask Mask returned by FloodHorizontalArea(), but trimmed
   * to have the specified width.
   * @param vertical_mask_x x-value returned by FloodHorizontalArea().
   * @param vertical_mask_width Width return by FloodHorizontalArea(), and width
   * of vertical_mask.
   * @param horizontal_mask Mask returned by FloodVerticalArea().
   * @param [in,out] mask Result
   * @param [out] sub_x Bounding box x-value
   * @param [out] sub_y Bounding box y-value
   * @param [out] subwidth Bounding box width
   * @param [out] subheight Bounding box height
   */
  void GetBoundingMask(const bool* vertical_mask, std::size_t vertical_mask_x,
                       std::size_t vertical_mask_width,
                       const bool* horizontal_mask, bool* mask,
                       std::size_t& sub_x, std::size_t& sub_y,
                       std::size_t& subwidth, std::size_t& subheight) const;

 private:
  /**
   * This function sets a rectangular area given by 0 <= y < height and x_start
   * <= x < x_end.
   */
  void FillColumns(float* output, std::size_t x_start, std::size_t x_end,
                   float new_value) const;

  const std::size_t width_;
  const std::size_t height_;
};
}  // namespace radler::math
#endif
