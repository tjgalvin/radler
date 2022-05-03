// SPDX-License-Identifier: LGPL-3.0-only

#include "algorithms/iuwt/image_analysis.h"

#include <stack>

namespace radler::algorithms::iuwt::image_analysis {

static bool ExceedsThreshold(float val, float threshold) {
  if (threshold >= 0.0) {
    return val > threshold;
  } else {
    return val < threshold || val > -threshold;
  }
}

static bool ExceedsThresholdAbs(float val, float threshold) {
  return std::fabs(val) > threshold;
}

bool IsHighestOnScale0(const IuwtDecomposition& iuwt, IuwtMask& marked_mask,
                       size_t& x, size_t& y, size_t end_scale,
                       float& highest_scale0) {
  const size_t width = iuwt.Width(), height = iuwt.Height();
  Component component(x, y, 0);
  std::stack<Component> todo;
  todo.push(component);
  marked_mask[0][x + y * width] = false;
  highest_scale0 = iuwt[0][x + y * width];
  float highestOtherScales = 0.0;
  while (!todo.empty()) {
    Component c = todo.top();
    todo.pop();
    size_t index = c.x + c.y * width;
    if (c.scale == 0) {
      if (ExceedsThreshold(iuwt[0][index], highest_scale0)) {
        highest_scale0 = iuwt[0][index];
        x = c.x;
        y = c.y;
      }
    } else {
      if (ExceedsThreshold(iuwt[c.scale][index], highestOtherScales)) {
        highestOtherScales = iuwt[c.scale][index];
      }
    }
    if (c.x > 0) {
      if (marked_mask[c.scale][index - 1]) {
        marked_mask[c.scale][index - 1] = false;
        todo.push(Component(c.x - 1, c.y, c.scale));
      }
    }
    if (c.x < width - 1) {
      if (marked_mask[c.scale][index + 1]) {
        marked_mask[c.scale][index + 1] = false;
        todo.push(Component(c.x + 1, c.y, c.scale));
      }
    }
    if (c.y > 0) {
      if (marked_mask[c.scale][index - width]) {
        marked_mask[c.scale][index - width] = false;
        todo.push(Component(c.x, c.y - 1, c.scale));
      }
    }
    if (c.y < height - 1) {
      if (marked_mask[c.scale][index + width]) {
        marked_mask[c.scale][index + width] = false;
        todo.push(Component(c.x, c.y + 1, c.scale));
      }
    }
    if (c.scale > 0) {
      if (marked_mask[c.scale - 1][index]) {
        marked_mask[c.scale - 1][index] = false;
        todo.push(Component(c.x, c.y, c.scale - 1));
      }
    }
    if (c.scale < static_cast<int>(end_scale) - 1) {
      if (marked_mask[c.scale + 1][index]) {
        marked_mask[c.scale + 1][index] = false;
        todo.push(Component(c.x, c.y, c.scale + 1));
      }
    }
  }
  return std::fabs(highest_scale0) > std::fabs(highestOtherScales);
}

void Floodfill(const IuwtDecomposition& iuwt, IuwtMask& mask,
               const aocommon::UVector<float>& thresholds, size_t min_scale,
               size_t end_scale, const Component& component, float clean_border,
               size_t& area_size) {
  const size_t width = iuwt.Width(), height = iuwt.Height();
  size_t xBorder = clean_border * width;
  size_t yBorder = clean_border * height;
  size_t minX = xBorder, maxX = width - xBorder;
  size_t minY = yBorder, maxY = height - yBorder;

  area_size = 0;
  end_scale = std::min<size_t>(end_scale, iuwt.NScales());
  std::stack<Component> todo;
  todo.push(component);
  mask[component.scale][component.x + component.y * width] = true;
  while (!todo.empty()) {
    Component c = todo.top();
    ++area_size;
    todo.pop();
    size_t index = c.x + c.y * width;
    if (c.x > minX) {
      if (ExceedsThreshold(iuwt[c.scale][index - 1], thresholds[c.scale]) &&
          !mask[c.scale][index - 1]) {
        mask[c.scale][index - 1] = true;
        todo.push(Component(c.x - 1, c.y, c.scale));
      }
    }
    if (c.x < maxX - 1) {
      if (ExceedsThreshold(iuwt[c.scale][index + 1], thresholds[c.scale]) &&
          !mask[c.scale][index + 1]) {
        mask[c.scale][index + 1] = true;
        todo.push(Component(c.x + 1, c.y, c.scale));
      }
    }
    if (c.y > minY) {
      if (ExceedsThreshold(iuwt[c.scale][index - width], thresholds[c.scale]) &&
          !mask[c.scale][index - width]) {
        mask[c.scale][index - width] = true;
        todo.push(Component(c.x, c.y - 1, c.scale));
      }
    }
    if (c.y < maxY - 1) {
      if (ExceedsThreshold(iuwt[c.scale][index + width], thresholds[c.scale]) &&
          !mask[c.scale][index + width]) {
        mask[c.scale][index + width] = true;
        todo.push(Component(c.x, c.y + 1, c.scale));
      }
    }
    if (c.scale > static_cast<int>(min_scale)) {
      if (ExceedsThreshold(iuwt[c.scale - 1][index], thresholds[c.scale - 1]) &&
          !mask[c.scale - 1][index]) {
        mask[c.scale - 1][index] = true;
        todo.push(Component(c.x, c.y, c.scale - 1));
      }
    }
    if (c.scale < static_cast<int>(end_scale) - 1) {
      if (ExceedsThreshold(iuwt[c.scale + 1][index], thresholds[c.scale + 1]) &&
          !mask[c.scale + 1][index]) {
        mask[c.scale + 1][index] = true;
        todo.push(Component(c.x, c.y, c.scale + 1));
      }
    }
  }
}

void MaskedFloodfill(const IuwtDecomposition& iuwt, IuwtMask& mask,
                     const aocommon::UVector<float>& thresholds,
                     size_t min_scale, size_t end_scale,
                     const Component& component, float clean_border,
                     const bool* prior_mask, size_t& area_size) {
  const size_t width = iuwt.Width(), height = iuwt.Height();
  size_t xBorder = clean_border * width;
  size_t yBorder = clean_border * height;
  size_t minX = xBorder, maxX = width - xBorder;
  size_t minY = yBorder, maxY = height - yBorder;

  area_size = 0;
  end_scale = std::min<size_t>(end_scale, iuwt.NScales());
  std::stack<Component> todo;
  todo.push(component);
  mask[component.scale][component.x + component.y * width] = true;
  while (!todo.empty()) {
    Component c = todo.top();
    ++area_size;
    todo.pop();
    size_t index = c.x + c.y * width;
    if (c.x > minX) {
      if (ExceedsThreshold(iuwt[c.scale][index - 1], thresholds[c.scale]) &&
          !mask[c.scale][index - 1] && prior_mask[index - 1]) {
        mask[c.scale][index - 1] = true;
        todo.push(Component(c.x - 1, c.y, c.scale));
      }
    }
    if (c.x < maxX - 1) {
      if (ExceedsThreshold(iuwt[c.scale][index + 1], thresholds[c.scale]) &&
          !mask[c.scale][index + 1] && prior_mask[index + 1]) {
        mask[c.scale][index + 1] = true;
        todo.push(Component(c.x + 1, c.y, c.scale));
      }
    }
    if (c.y > minY) {
      if (ExceedsThreshold(iuwt[c.scale][index - width], thresholds[c.scale]) &&
          !mask[c.scale][index - width] && prior_mask[index - width]) {
        mask[c.scale][index - width] = true;
        todo.push(Component(c.x, c.y - 1, c.scale));
      }
    }
    if (c.y < maxY - 1) {
      if (ExceedsThreshold(iuwt[c.scale][index + width], thresholds[c.scale]) &&
          !mask[c.scale][index + width] && prior_mask[index + width]) {
        mask[c.scale][index + width] = true;
        todo.push(Component(c.x, c.y + 1, c.scale));
      }
    }
    if (c.scale > static_cast<int>(min_scale)) {
      if (ExceedsThreshold(iuwt[c.scale - 1][index], thresholds[c.scale - 1]) &&
          !mask[c.scale - 1][index] && prior_mask[index]) {
        mask[c.scale - 1][index] = true;
        todo.push(Component(c.x, c.y, c.scale - 1));
      }
    }
    if (c.scale < static_cast<int>(end_scale) - 1) {
      if (ExceedsThreshold(iuwt[c.scale + 1][index], thresholds[c.scale + 1]) &&
          !mask[c.scale + 1][index] && prior_mask[index]) {
        mask[c.scale + 1][index] = true;
        todo.push(Component(c.x, c.y, c.scale + 1));
      }
    }
  }
}

void SelectStructures(const IuwtDecomposition& iuwt, IuwtMask& mask,
                      const aocommon::UVector<float>& thresholds,
                      size_t min_scale, size_t end_scale, float clean_border,
                      const bool* prior_mask, size_t& area_size) {
  const size_t width = iuwt.Width(), height = iuwt.Height();
  const size_t xBorder = clean_border * width, yBorder = clean_border * height,
               minX = xBorder, maxX = width - xBorder, minY = yBorder,
               maxY = height - yBorder;

  area_size = 0;

  for (size_t scale = min_scale; scale != end_scale; ++scale) {
    for (size_t y = minY; y != maxY; ++y) {
      for (size_t x = minX; x != maxX; ++x) {
        size_t index = x + y * width;
        bool isInPriorMask = (prior_mask == nullptr) || prior_mask[index];
        if (ExceedsThreshold(iuwt[scale][index], thresholds[scale]) &&
            !mask[scale][index] && isInPriorMask) {
          Component component(x, y, scale);
          size_t subAreaSize = 0;
          if (prior_mask == nullptr) {
            Floodfill(iuwt, mask, thresholds, min_scale, end_scale, component,
                      clean_border, subAreaSize);
          } else {
            MaskedFloodfill(iuwt, mask, thresholds, min_scale, end_scale,
                            component, clean_border, prior_mask, subAreaSize);
          }
          area_size += subAreaSize;
        }
      }
    }
  }
}

void FloodFill2D(const float* image, bool* mask, float threshold,
                 const Component2D& component, size_t width, size_t height,
                 size_t& area_size) {
  area_size = 0;
  std::stack<Component2D> todo;
  todo.push(component);
  mask[component.x + component.y * width] = true;
  while (!todo.empty()) {
    Component2D c = todo.top();
    ++area_size;
    todo.pop();
    size_t index = c.x + c.y * width;
    if (c.x > 0) {
      if (ExceedsThreshold(image[index - 1], threshold) && !mask[index - 1]) {
        mask[index - 1] = true;
        todo.push(Component2D(c.x - 1, c.y));
      }
    }
    if (c.x < width - 1) {
      if (ExceedsThreshold(image[index + 1], threshold) && !mask[index + 1]) {
        mask[index + 1] = true;
        todo.push(Component2D(c.x + 1, c.y));
      }
    }
    if (c.y > 0) {
      if (ExceedsThreshold(image[index - width], threshold) &&
          !mask[index - width]) {
        mask[index - width] = true;
        todo.push(Component2D(c.x, c.y - 1));
      }
    }
    if (c.y < height - 1) {
      if (ExceedsThreshold(image[index + width], threshold) &&
          !mask[index + width]) {
        mask[index + width] = true;
        todo.push(Component2D(c.x, c.y + 1));
      }
    }
  }
}

void FloodFill2D(const float* image, bool* mask, float threshold,
                 const Component2D& component, size_t width, size_t height,
                 std::vector<Component2D>& area) {
  area.clear();
  std::stack<Component2D> todo;
  todo.push(component);
  mask[component.x + component.y * width] = true;
  while (!todo.empty()) {
    Component2D c = todo.top();
    area.push_back(c);
    todo.pop();
    size_t index = c.x + c.y * width;
    if (c.x > 0) {
      if (ExceedsThresholdAbs(image[index - 1], threshold) &&
          !mask[index - 1]) {
        mask[index - 1] = true;
        todo.push(Component2D(c.x - 1, c.y));
      }
    }
    if (c.x < width - 1) {
      if (ExceedsThresholdAbs(image[index + 1], threshold) &&
          !mask[index + 1]) {
        mask[index + 1] = true;
        todo.push(Component2D(c.x + 1, c.y));
      }
    }
    if (c.y > 0) {
      if (ExceedsThresholdAbs(image[index - width], threshold) &&
          !mask[index - width]) {
        mask[index - width] = true;
        todo.push(Component2D(c.x, c.y - 1));
      }
    }
    if (c.y < height - 1) {
      if (ExceedsThresholdAbs(image[index + width], threshold) &&
          !mask[index + width]) {
        mask[index + width] = true;
        todo.push(Component2D(c.x, c.y + 1));
      }
    }
  }
}
}  // namespace radler::algorithms::iuwt::image_analysis
