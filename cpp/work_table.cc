// SPDX-License-Identifier: LGPL-3.0-only

#include "work_table.h"

#include <aocommon/throwruntimeerror.h>

#include <algorithm>
#include <cassert>
#include <stdexcept>

namespace radler {

WorkTable::WorkTable(std::vector<PsfOffset> psf_offsets,
                     std::size_t n_original_groups,
                     std::size_t n_deconvolution_groups,
                     std::size_t channel_index_offset)
    : entries_(),
      psf_offsets_(std::move(psf_offsets)),
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

void WorkTable::ValidatePsfs() const {
  const size_t n_psfs = std::max<size_t>(1, psf_offsets_.size());

  if (!entries_.empty()) {
    // The Front() entry should always have PSFs.
    if (Front().psf_accessors.size() != n_psfs) {
      aocommon::ThrowRuntimeError(
          "WorkTable: Expected ", n_psfs,
          " PSF accessors in the first entry, but found ",
          Front().psf_accessors.size(), " PSF accessors.");
    }

    for (const Group& group : original_groups_) {
      for (size_t i = 0; i < group.size(); ++i) {
        const WorkTableEntry& entry = *group[i];

        // The first entry of each OriginalGroup must have PSFs.
        if (i == 0) {
          if (entry.psf_accessors.size() != n_psfs) {
            aocommon::ThrowRuntimeError(
                "WorkTable: Expected ", n_psfs,
                " PSF accessors per entry, but found an entry with ",
                entry.psf_accessors.size(), " PSF accessors.");
          }
          for (size_t psf_index = 0; psf_index < n_psfs; ++psf_index) {
            const size_t width = entry.psf_accessors[psf_index]->Width();
            const size_t height = entry.psf_accessors[psf_index]->Height();
            if (width == 0 || height == 0) {
              aocommon::ThrowRuntimeError(
                  "WorkTable: Found an entry with an empty image for PSF "
                  "accessor ",
                  psf_index, ".");
            }

            if (width != Front().psf_accessors[psf_index]->Width() ||
                height != Front().psf_accessors[psf_index]->Height()) {
              aocommon::ThrowRuntimeError(
                  "WorkTable: Found an entry with a different size for PSF "
                  "accessor ",
                  psf_index, ".");
            }
          }
        } else {
          if (!entry.psf_accessors.empty()) {
            throw std::runtime_error(
                "WorkTable: Only the first entry for a channel may have PSF "
                "accessors.");
          }
        }
      }
    }
  }
}

}  // namespace radler
