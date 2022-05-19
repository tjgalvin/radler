// SPDX-License-Identifier: LGPL-3.0-only

#include "work_table.h"

#include <algorithm>
#include <cassert>

namespace radler {

WorkTable::WorkTable(std::size_t n_original_groups,
                     std::size_t n_deconvolution_groups,
                     std::size_t channel_index_offset)
    : entries_(),
      channel_index_offset_(channel_index_offset),
      original_groups_(
          std::max(n_original_groups, static_cast<std::size_t>(1))),
      deconvolution_groups_(
          (n_deconvolution_groups == 0)
              ? original_groups_.size()
              // The number of deconvolution groups is always less or equal to
              // the number of original groups.
              : std::min(original_groups_.size(), n_deconvolution_groups)) {
  // Create an entry in deconvolution_groups for each original group.
  for (std::size_t i = 0; i < original_groups_.size(); ++i) {
    std::size_t deconvolution_index =
        i * deconvolution_groups_.size() / original_groups_.size();
    deconvolution_groups_[deconvolution_index].push_back(i);
  }
}

void WorkTable::AddEntry(std::unique_ptr<WorkTableEntry> entry) {
  const std::size_t original_channel_index = entry->original_channel_index;
  assert(original_channel_index < original_groups_.size());

  entry->index = entries_.size();
  entries_.push_back(std::move(entry));

  original_groups_[original_channel_index].push_back(entries_.back().get());
}
}  // namespace radler
