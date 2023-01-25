// SPDX-License-Identifier: LGPL-3.0-only

#include "work_table.h"

#include <array>

#include <boost/test/unit_test.hpp>

#include "test/smartptr.h"

namespace radler {

namespace {
/**
 * Helper class for testing ValidatePsfs(), which should only use the Width()
 * and Height() members of the ImageAccessors for the PSF images.
 */
class SizeOnlyAccessor : public aocommon::ImageAccessor {
 public:
  explicit SizeOnlyAccessor(const size_t width, const size_t height)
      : width_(width), height_(height) {}

  size_t Width() const override { return width_; }
  size_t Height() const override { return height_; }
  void Load(float*) const override {
    BOOST_FAIL("Unexpected SizeOnlyAccessor::Load call");
  }
  void Store(const float*) override {
    BOOST_FAIL("Unexpected SizeOnlyAccessor::Store call");
  }

 private:
  size_t width_;
  size_t height_;
};
}  // namespace

BOOST_AUTO_TEST_SUITE(work_table)

BOOST_AUTO_TEST_CASE(constructor) {
  const size_t kTableSize = 42;

  WorkTable table({}, kTableSize, kTableSize);

  BOOST_TEST(table.OriginalGroups().size() == kTableSize);
  for (const WorkTable::Group& group : table.OriginalGroups()) {
    BOOST_TEST(group.empty());
  }

  BOOST_TEST_REQUIRE(table.DeconvolutionGroups().size() == kTableSize);
  for (size_t index = 0; index < kTableSize; ++index) {
    BOOST_TEST(table.DeconvolutionGroups()[index] ==
               std::vector<size_t>(1, index));
  }

  BOOST_TEST((table.Begin() == table.End()));
  BOOST_TEST(table.Size() == 0);
}

BOOST_AUTO_TEST_CASE(single_deconvolution_group) {
  WorkTable table({}, 7, 1);
  const std::vector<std::vector<size_t>> kExpectedGroups{{0, 1, 2, 3, 4, 5, 6}};
  BOOST_TEST_REQUIRE(table.DeconvolutionGroups() == kExpectedGroups);
}

BOOST_AUTO_TEST_CASE(multiple_deconvolution_groups) {
  WorkTable table({}, 7, 3);
  const std::vector<std::vector<size_t>> kExpectedGroups{
      {0, 1, 2}, {3, 4}, {5, 6}};
  BOOST_TEST_REQUIRE(table.DeconvolutionGroups() == kExpectedGroups);
}

BOOST_AUTO_TEST_CASE(too_many_deconvolution_groups) {
  WorkTable table({}, 7, 42);
  const std::vector<std::vector<size_t>> kExpectedGroups{{0}, {1}, {2}, {3},
                                                         {4}, {5}, {6}};
  BOOST_TEST_REQUIRE(table.DeconvolutionGroups() == kExpectedGroups);
}

BOOST_AUTO_TEST_CASE(add_entries) {
  WorkTable table({}, 3, 1);

  std::array<test::UniquePtr<WorkTableEntry>, 3> entries;
  entries[0]->original_channel_index = 1;
  entries[1]->original_channel_index = 0;
  entries[2]->original_channel_index = 1;

  for (test::UniquePtr<WorkTableEntry>& entry : entries) {
    table.AddEntry(entry.take());
    // table.AddEntry(std::move(entry));
  }

  // Check if the OriginalGroups have the correct size and correct entries.
  const std::vector<WorkTable::Group>& original_groups = table.OriginalGroups();
  BOOST_TEST_REQUIRE(original_groups.size() == 3);

  BOOST_TEST_REQUIRE(original_groups[0].size() == 1);
  BOOST_TEST_REQUIRE(original_groups[1].size() == 2);
  BOOST_TEST(original_groups[2].empty());

  BOOST_TEST(original_groups[0][0] == entries[1].get());
  BOOST_TEST(original_groups[1][0] == entries[0].get());
  BOOST_TEST(original_groups[1][1] == entries[2].get());

  // Check if a range based loop, which uses begin() and end(), yields the
  // entries.
  size_t index = 0;
  for (const WorkTableEntry& entry : table) {
    BOOST_TEST(&entry == entries[index].get());
    ++index;
  }

  // Finally, check Front() and Size().
  BOOST_TEST(&table.Front() == entries.front().get());
  BOOST_TEST(table.Size() == entries.size());
}

BOOST_AUTO_TEST_CASE(print_no_entries_no_psfs) {
  WorkTable table({}, 0, 0);

  std::stringstream output;
  output << table;
  BOOST_CHECK_EQUAL(output.str(),
                    R"(=== IMAGING TABLE ===
Original groups       1
Deconvolution groups  1
Channel index         0
)");
}

BOOST_AUTO_TEST_CASE(print_entries_no_psfs) {
  WorkTable table({}, 3, 1);
  table.AddEntry(std::make_unique<WorkTableEntry>(
      WorkTableEntry{0, 5'000'000.0, 10'000'000.0,
                     aocommon::PolarizationEnum::StokesQ, 0, 4, 12, 1.234}));
  table.AddEntry(std::make_unique<WorkTableEntry>(
      WorkTableEntry{1, 1'000'000.0, 2'000'000.0,
                     aocommon::PolarizationEnum::StokesI, 1, 8, 13, 1.01}));
  table.AddEntry(std::make_unique<WorkTableEntry>(
      WorkTableEntry{1, 100'000'000.0, 200'000'000.0,
                     aocommon::PolarizationEnum::StokesU, 0, 16, 12, 1.1}));

  std::stringstream output;
  output << table;
  BOOST_CHECK_EQUAL(output.str(),
                    R"(=== IMAGING TABLE ===
Original groups       3
Deconvolution groups  1
Channel index         0
   # Pol Ch Mask Interval Weight Freq(MHz)
   0   Q  0   12        4  1.234 5-10
   1   I  1   13        8   1.01 1-2
   2   U  0   12       16    1.1 100-200
)");
}

BOOST_AUTO_TEST_CASE(print_entries_psfs) {
  WorkTable table({PsfOffset{1, 2}, PsfOffset{3, 4}}, 3, 1);
  table.AddEntry(std::make_unique<WorkTableEntry>(
      WorkTableEntry{0, 5'000'000.0, 10'000'000.0,
                     aocommon::PolarizationEnum::StokesQ, 0, 4, 12, 1.234}));
  table.AddEntry(std::make_unique<WorkTableEntry>(
      WorkTableEntry{1, 1'000'000.0, 2'000'000.0,
                     aocommon::PolarizationEnum::StokesI, 1, 8, 13, 1.01}));
  table.AddEntry(std::make_unique<WorkTableEntry>(
      WorkTableEntry{1, 100'000'000.0, 200'000'000.0,
                     aocommon::PolarizationEnum::StokesU, 0, 16, 12, 1.1}));

  std::stringstream output;
  output << table;
  BOOST_CHECK_EQUAL(output.str(),
                    R"(=== IMAGING TABLE ===
Original groups       3
Deconvolution groups  1
Channel index         0
   # Pol Ch Mask Interval Weight Freq(MHz)
   0   Q  0   12        4  1.234 5-10
   1   I  1   13        8   1.01 1-2
   2   U  0   12       16    1.1 100-200
=== PSFs ===
[x: 1, y: 2]
[x: 3, y: 4]
)");
}

BOOST_AUTO_TEST_CASE(validate_psfs_valid_no_entries) {
  const WorkTable table({}, 1, 1);
  // Valid since there are no entries.
  table.ValidatePsfs();
}

BOOST_AUTO_TEST_CASE(validate_psfs_valid_single_entry) {
  {
    WorkTable table({PsfOffset{1, 1}}, 1, 1);
    auto entry = std::make_unique<WorkTableEntry>(
        WorkTableEntry{0, 5'000'000.0, 10'000'000.0,
                       aocommon::PolarizationEnum::StokesQ, 0, 4, 12, 1.234});
    entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 6));
    table.AddEntry(std::move(entry));
    table.ValidatePsfs();
  }
  for (int n_psfs = 2; n_psfs < 10; ++n_psfs) {
    WorkTable table(std::vector<PsfOffset>(n_psfs), 1, 1);
    auto entry = std::make_unique<WorkTableEntry>(
        WorkTableEntry{0, 5'000'000.0, 10'000'000.0,
                       aocommon::PolarizationEnum::StokesQ, 0, 4, 12, 1.234});

    for (int i = 0; i < n_psfs; ++i) {
      entry->psf_accessors.emplace_back(
          std::make_unique<SizeOnlyAccessor>(4 + i, 6));
    }
    table.AddEntry(std::move(entry));

    table.ValidatePsfs();
  }
}

BOOST_AUTO_TEST_CASE(validate_psfs_valid_multiple_entries) {
  const size_t kNChannels = 3;
  const size_t kNPsfs = 2;
  const size_t kWidth = 42;
  const size_t kHeight = 43;

  WorkTable table(std::vector<PsfOffset>(kNPsfs), kNChannels, kNChannels);

  for (size_t ch = 0; ch < kNChannels; ch++) {
    auto i_entry = std::make_unique<WorkTableEntry>();
    i_entry->polarization = aocommon::PolarizationEnum::StokesI;
    i_entry->original_channel_index = ch;
    for (size_t psf_index = 0; psf_index < kNPsfs; ++psf_index) {
      i_entry->psf_accessors.emplace_back(
          std::make_unique<SizeOnlyAccessor>(kWidth, kHeight));
    }
    table.AddEntry(std::move(i_entry));

    auto q_entry = std::make_unique<WorkTableEntry>();
    q_entry->polarization = aocommon::PolarizationEnum::StokesQ;
    q_entry->original_channel_index = ch;
    // Do not add PSF accessor
    table.AddEntry(std::move(q_entry));

    auto u_entry = std::make_unique<WorkTableEntry>();
    u_entry->polarization = aocommon::PolarizationEnum::StokesU;
    u_entry->original_channel_index = ch;
    // Do not add PSF accessor
    table.AddEntry(std::move(u_entry));
  }
  BOOST_CHECK_NO_THROW(table.ValidatePsfs());
}

BOOST_AUTO_TEST_CASE(validate_psfs_too_few_accessors) {
  WorkTable table({}, 1, 1);
  table.AddEntry(std::make_unique<WorkTableEntry>(
      WorkTableEntry{0, 5'000'000.0, 10'000'000.0,
                     aocommon::PolarizationEnum::StokesQ, 0, 4, 12, 1.234}));

  BOOST_CHECK_THROW(table.ValidatePsfs(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(validate_psfs_too_many_accessors) {
  WorkTable table({PsfOffset{1, 1}}, 1, 1);
  auto entry = std::make_unique<WorkTableEntry>(
      WorkTableEntry{0, 5'000'000.0, 10'000'000.0,
                     aocommon::PolarizationEnum::StokesQ, 0, 4, 12, 1.234});
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 5));
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(6, 7));
  table.AddEntry(std::move(entry));
  BOOST_CHECK_THROW(table.ValidatePsfs(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(validate_psfs_zero_width) {
  WorkTable table({}, 1, 1);
  auto entry = std::make_unique<WorkTableEntry>();
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(0, 5));
  table.AddEntry(std::move(entry));
  BOOST_CHECK_THROW(table.ValidatePsfs(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(validate_psfs_zero_height) {
  WorkTable table({}, 1, 1);
  auto entry = std::make_unique<WorkTableEntry>();
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 0));
  table.AddEntry(std::move(entry));
  BOOST_CHECK_THROW(table.ValidatePsfs(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(validate_psfs_second_group_without_accessors) {
  WorkTable table({}, 2, 1);

  auto entry = std::make_unique<WorkTableEntry>();
  entry->original_channel_index = 0;
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 5));
  table.AddEntry(std::move(entry));

  entry = std::make_unique<WorkTableEntry>();
  entry->original_channel_index = 1;
  table.AddEntry(std::move(entry));

  BOOST_CHECK_THROW(table.ValidatePsfs(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(validate_psfs_second_group_with_too_many_accessors) {
  WorkTable table({}, 2, 1);

  auto entry = std::make_unique<WorkTableEntry>();
  entry->original_channel_index = 0;
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 5));
  table.AddEntry(std::move(entry));

  entry = std::make_unique<WorkTableEntry>();
  entry->original_channel_index = 1;
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 5));
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 5));
  table.AddEntry(std::move(entry));

  BOOST_CHECK_THROW(table.ValidatePsfs(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(validate_psfs_two_groups_incorrect_width) {
  WorkTable table({}, 2, 1);

  auto entry = std::make_unique<WorkTableEntry>();
  entry->original_channel_index = 0;
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 5));
  table.AddEntry(std::move(entry));

  entry = std::make_unique<WorkTableEntry>();
  entry->original_channel_index = 1;
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(10, 5));
  table.AddEntry(std::move(entry));

  BOOST_CHECK_THROW(table.ValidatePsfs(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(validate_psfs_two_groups_incorrect_height) {
  WorkTable table({}, 2, 1);

  auto entry = std::make_unique<WorkTableEntry>();
  entry->original_channel_index = 0;
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 5));
  table.AddEntry(std::move(entry));

  entry = std::make_unique<WorkTableEntry>();
  entry->original_channel_index = 1;
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 10));
  table.AddEntry(std::move(entry));

  BOOST_CHECK_THROW(table.ValidatePsfs(), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(validate_psfs_single_group_second_entry_has_accessors) {
  WorkTable table({}, 1, 1);

  auto entry = std::make_unique<WorkTableEntry>();
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 5));
  table.AddEntry(std::move(entry));

  entry = std::make_unique<WorkTableEntry>();
  entry->psf_accessors.emplace_back(std::make_unique<SizeOnlyAccessor>(4, 5));
  table.AddEntry(std::move(entry));

  BOOST_CHECK_THROW(table.ValidatePsfs(), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace radler
