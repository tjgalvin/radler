// SPDX-License-Identifier: LGPL-3.0-only

#include "dijkstra_splitter.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <queue>

#include <aocommon/uvector.h>

namespace radler::math {

void DijkstraSplitter::AddVerticalDivider(const float* image, float* scratch,
                                          float* output, std::size_t x1,
                                          std::size_t x2) const {
  DivideVertically(image, scratch, x1, x2);
  for (std::size_t y = 0; y != height_; ++y) {
    for (std::size_t i = y * width_ + x1; i != y * width_ + x2; ++i)
      output[i] += scratch[i];
  }
}

void DijkstraSplitter::AddHorizontalDivider(const float* image, float* scratch,
                                            float* output, std::size_t y1,
                                            std::size_t y2) const {
  DivideHorizontally(image, scratch, y1, y2);
  for (std::size_t i = y1 * width_; i != y2 * width_; ++i)
    output[i] += scratch[i];
}

void DijkstraSplitter::DivideVertically(const float* image, float* output,
                                        std::size_t x1, std::size_t x2) const {
  using VisitSet = std::priority_queue<Visit>;
  VisitSet visits;

  for (std::size_t x = x1; x != x2; ++x) {
    Visit visit;
    visit.distance = 0.0;
    visit.to = Coord(x, 0);
    visit.from = Coord(x, 0);
    visits.push(visit);
  }
  aocommon::UVector<Coord> path((x2 - x1) * height_);
  FillColumns(output, x1, x2, std::numeric_limits<float>::max());
  Visit visit;
  while (!visits.empty()) {
    visit = visits.top();
    visits.pop();
    const std::size_t x = visit.to.x;
    const std::size_t y = visit.to.y;
    if (y == height_) break;
    const float current_distance = output[x + y * width_];
    const float new_distance =
        visit.distance + std::fabs(image[x + y * width_]);
    if (new_distance < current_distance) {
      output[x + y * width_] = new_distance;
      path[x - x1 + y * (x2 - x1)] = visit.from;
      visit.distance = new_distance;
      visit.from = visit.to;
      if (x > x1) {
        visit.to = Coord(x - 1, y + 1);
        visits.push(visit);
        visit.to = Coord(x - 1, y);
        visits.push(visit);
      }
      visit.to = Coord(x, y + 1);
      visits.push(visit);
      if (x < x2 - 1) {
        visit.to = Coord(x + 1, y + 1);
        visits.push(visit);
        visit.to = Coord(x + 1, y);
        visits.push(visit);
      }
    }
  }
  FillColumns(output, x1, x2, 0.0);
  Coord path_coord = visit.from;
  while (path_coord.y > 0) {
    output[path_coord.x + path_coord.y * width_] = 1.0;
    path_coord = path[path_coord.x - x1 + path_coord.y * (x2 - x1)];
  }
  output[path_coord.x] = 1.0;
}

void DijkstraSplitter::DivideHorizontally(const float* image, float* output,
                                          std::size_t y1,
                                          std::size_t y2) const {
  using VisitSet = std::priority_queue<Visit>;
  VisitSet visits;

  for (std::size_t y = y1; y != y2; ++y) {
    Visit visit;
    visit.distance = 0.0;
    visit.to = Coord(0, y);
    visit.from = Coord(0, y);
    visits.push(visit);
  }
  aocommon::UVector<Coord> path(width_ * (y2 - y1));
  std::fill(output + y1 * width_, output + y2 * width_,
            std::numeric_limits<float>::max());
  Visit visit;
  while (!visits.empty()) {
    visit = visits.top();
    visits.pop();
    const std::size_t x = visit.to.x;
    const std::size_t y = visit.to.y;
    if (x == width_) break;
    const float current_distance = output[x + y * width_];
    const float new_distance =
        visit.distance + std::fabs(image[x + y * width_]);
    if (new_distance < current_distance) {
      output[x + y * width_] = new_distance;
      path[x + (y - y1) * width_] = visit.from;
      visit.distance = new_distance;
      visit.from = visit.to;
      if (y > y1) {
        visit.to = Coord(x + 1, y - 1);
        visits.push(visit);
        visit.to = Coord(x, y - 1);
        visits.push(visit);
      }
      visit.to = Coord(x + 1, y);
      visits.push(visit);
      if (y < y2 - 1) {
        visit.to = Coord(x + 1, y + 1);
        visits.push(visit);
        visit.to = Coord(x, y + 1);
        visits.push(visit);
      }
    }
  }
  std::fill(output + y1 * width_, output + y2 * width_, 0.0);
  Coord path_coord = visit.from;
  while (path_coord.x > 0) {
    output[path_coord.x + path_coord.y * width_] = 1.0;
    path_coord = path[path_coord.x + (path_coord.y - y1) * width_];
  }
  output[path_coord.y * width_] = 1.0;
}

void DijkstraSplitter::FloodVerticalArea(const float* subdivision,
                                         std::size_t subimage_x, bool* mask,
                                         std::size_t& x,
                                         std::size_t& subwidth) const {
  std::fill(mask, mask + width_ * height_, false);
  x = width_;
  std::size_t x2 = 0;
  for (std::size_t y = 0; y != height_; ++y) {
    const float* division_row = &subdivision[y * width_];
    bool* mask_row = &mask[y * width_];
    int x_iter = subimage_x;
    // Move to the left until a border is hit
    while (x_iter >= 0 && division_row[x_iter] == 0.0) {
      mask_row[x_iter] = true;
      --x_iter;
    }
    // Continue to the left through the border
    while (x_iter >= 0 && division_row[x_iter] != 0.0) {
      mask_row[x_iter] = true;
      --x_iter;
    }
    x = std::min<std::size_t>(x, x_iter + 1);
    x_iter = subimage_x + 1;
    // Move to the right until a border is hit
    while (std::size_t(x_iter) < width_ && division_row[x_iter] == 0.0) {
      mask_row[x_iter] = true;
      ++x_iter;
    }
    x2 = std::max<std::size_t>(x2, x_iter);
  }
  if (x2 < x) {
    subwidth = 0;
  } else {
    subwidth = x2 - x;
  }
}

void DijkstraSplitter::FloodHorizontalArea(const float* subdivision,
                                           std::size_t subimage_y, bool* mask,
                                           std::size_t& y,
                                           std::size_t& subheight) const {
  std::fill(mask, mask + width_ * height_, false);
  y = height_;
  std::size_t y2 = 0;
  for (std::size_t x = 0; x != width_; ++x) {
    int y_iter = subimage_y;
    // Move up until a border is hit
    while (y_iter >= 0 && subdivision[y_iter * width_ + x] == 0.0) {
      mask[y_iter * width_ + x] = true;
      --y_iter;
    }
    // Continue to the left through the border
    while (y_iter >= 0 && subdivision[y_iter * width_ + x] != 0.0) {
      mask[y_iter * width_ + x] = true;
      --y_iter;
    }
    y = std::min<std::size_t>(y, y_iter + 1);
    y_iter = subimage_y + 1;
    // Move to the right until a border is hit
    while (std::size_t(y_iter) < height_ &&
           subdivision[y_iter * width_ + x] == 0.0) {
      mask[y_iter * width_ + x] = true;
      ++y_iter;
    }
    y2 = std::max<std::size_t>(y2, y_iter);
  }
  if (y2 < y) {
    subheight = 0;
  } else {
    subheight = y2 - y;
  }
}

void DijkstraSplitter::GetBoundingMask(const bool* vertical_mask,
                                       std::size_t vertical_mask_x,
                                       std::size_t vertical_mask_width,
                                       const bool* horizontal_mask, bool* mask,
                                       std::size_t& sub_x, std::size_t& sub_y,
                                       std::size_t& subwidth,
                                       std::size_t& subheight) const {
  sub_x = vertical_mask_width + vertical_mask_x;
  sub_y = height_;
  std::size_t sub_x2 = 0;
  std::size_t sub_y2 = 0;
  for (std::size_t y = 0; y != height_; ++y) {
    const bool* vertical_mask_row = &vertical_mask[y * vertical_mask_width];
    const bool* horizontal_mask_row = &horizontal_mask[y * width_];
    bool* mask_row = &mask[y * width_];
    for (std::size_t x = 0; x != vertical_mask_width; ++x) {
      std::size_t hx = x + vertical_mask_x;
      if (vertical_mask_row[x] && horizontal_mask_row[hx]) {
        mask_row[hx] = true;
        sub_x = std::min(hx, sub_x);
        sub_y = std::min(y, sub_y);
        sub_x2 = std::max(hx, sub_x2);
        sub_y2 = std::max(y, sub_y2);
      } else {
        mask_row[hx] = false;
      }
    }
  }
  if (sub_x2 < sub_x) {
    subwidth = 0;
    subheight = 0;
  } else {
    subwidth = sub_x2 + 1 - sub_x;
    subheight = sub_y2 + 1 - sub_y;
  }
  // If dimensions start off even, keep subimages even too
  if (width_ % 2 == 0) {
    if (subwidth % 2 != 0) {
      ++subwidth;
      if (subwidth + sub_x >= width_) --sub_x;
    }
  }
  if (height_ % 2 == 0) {
    if (subheight % 2 != 0) {
      ++subheight;
      if (subheight + sub_y >= height_) --sub_y;
    }
  }
}

void DijkstraSplitter::FillColumns(float* output, std::size_t x_start,
                                   std::size_t x_end, float new_value) const {
  assert(x_start <= x_end);
  for (std::size_t y = 0; y != height_; ++y) {
    std::fill(output + width_ * y + x_start, output + width_ * y + x_end,
              new_value);
  }
}

}  // namespace radler::math