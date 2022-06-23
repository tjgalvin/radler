// SPDX-License-Identifier: LGPL-3.0-only

#include "work_table.h"

#include <array>

#include <boost/test/unit_test.hpp>

#include "test/smartptr.h"

namespace radler {

BOOST_AUTO_TEST_SUITE(work_table)

BOOST_AUTO_TEST_CASE(constructor) {
  const size_t kTableSize = 42;

  WorkTable table(kTableSize, kTableSize);

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
  WorkTable table(7, 1);
  const std::vector<std::vector<size_t>> kExpectedGroups{{0, 1, 2, 3, 4, 5, 6}};
  BOOST_TEST_REQUIRE(table.DeconvolutionGroups() == kExpectedGroups);
}

BOOST_AUTO_TEST_CASE(multiple_deconvolution_groups) {
  WorkTable table(7, 3);
  const std::vector<std::vector<size_t>> kExpectedGroups{
      {0, 1, 2}, {3, 4}, {5, 6}};
  BOOST_TEST_REQUIRE(table.DeconvolutionGroups() == kExpectedGroups);
}

BOOST_AUTO_TEST_CASE(too_many_deconvolution_groups) {
  WorkTable table(7, 42);
  const std::vector<std::vector<size_t>> kExpectedGroups{{0}, {1}, {2}, {3},
                                                         {4}, {5}, {6}};
  BOOST_TEST_REQUIRE(table.DeconvolutionGroups() == kExpectedGroups);
}

BOOST_AUTO_TEST_CASE(add_entries) {
  WorkTable table(3, 1);

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

BOOST_AUTO_TEST_CASE(print_empty) {
  WorkTable table(0, 0);

  std::stringstream output;
  output << table;
  BOOST_CHECK_EQUAL(output.str(),
                    R"(=== IMAGING TABLE ===
Original groups       1
Deconvolution groups  1
Channel index         0
)");
}

BOOST_AUTO_TEST_CASE(print_not_empty) {
  WorkTable table(3, 1);
  table.AddEntry(std::make_unique<WorkTableEntry>(
      WorkTableEntry{0, 5'000'000.0, 10'000'000.0,
                     aocommon::PolarizationEnum::StokesQ, 0, 4, 1.234}));
  table.AddEntry(std::make_unique<WorkTableEntry>(
      WorkTableEntry{1, 1'000'000.0, 2'000'000.0,
                     aocommon::PolarizationEnum::StokesI, 1, 8, 1.01}));
  table.AddEntry(std::make_unique<WorkTableEntry>(
      WorkTableEntry{1, 100'000'000.0, 200'000'000.0,
                     aocommon::PolarizationEnum::StokesU, 0, 16, 1.1}));

  std::stringstream output;
  output << table;
  BOOST_CHECK_EQUAL(output.str(),
                    R"(=== IMAGING TABLE ===
Original groups       3
Deconvolution groups  1
Channel index         0
   # Pol Ch Interval Weight Freq(MHz)
   0   Q  0        4  1.234 5-10
   1   I  1        8   1.01 1-2
   2   U  0       16    1.1 100-200
)");
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace radler
