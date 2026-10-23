// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_UTILS_PHASE_TIMER_H_
#define RADLER_UTILS_PHASE_TIMER_H_

// Profiling-only instrumentation: accumulates wall-clock time per named phase
// over all threads and prints a summary to stderr at program exit.

#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <map>
#include <mutex>
#include <utility>

namespace radler::utils {

enum class Phase : size_t {
  kPeakOnlyPassWall,
  kMainPassWall,
  kSetupConvolvePsfs,
  kInitialScaleSearch,
  kTransformList,
  kSubMinorLoopRun,
  kCorrectResidualDirty,
  kModelUpdate,
  kScaleSearch,
  kCount
};

enum class Counter : size_t {
  kModelComponents,
  kModelUpdatesDirect,
  kModelUpdatesFft,
  kSubtractionsDirect,
  kSubtractionsFft,
  kCount
};

class PhaseTimers {
 public:
  static PhaseTimers& Instance() {
    static PhaseTimers instance;
    return instance;
  }

  void Add(Phase phase, int64_t nanoseconds) {
    const size_t index = static_cast<size_t>(phase);
    nanoseconds_[index].fetch_add(nanoseconds, std::memory_order_relaxed);
    calls_[index].fetch_add(1, std::memory_order_relaxed);
  }

  void Count(Counter counter, int64_t value) {
    counts_[static_cast<size_t>(counter)].fetch_add(value,
                                                    std::memory_order_relaxed);
  }

  /// Counts whether an algorithm object is run with the same image dimensions
  /// as on its previous run.
  void CountDimensions(const void* algorithm, size_t width, size_t height) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto iter = dimensions_.find(algorithm);
    if (iter == dimensions_.end()) {
      dimensions_.emplace(algorithm, std::make_pair(width, height));
      ++first_runs_;
    } else {
      if (iter->second == std::make_pair(width, height)) {
        ++same_dimensions_;
      } else {
        ++changed_dimensions_;
      }
      iter->second = std::make_pair(width, height);
    }
  }

  ~PhaseTimers() {
    static constexpr const char* kNames[] = {
        "peak-only pass (wall)",   "main pass (wall)",
        "setup ConvolvePsfs",      "initial scale search",
        "transformList Transform", "SubMinorLoop::Run",
        "CorrectResidualDirty",    "model update",
        "scale search (in loop)"};
    std::fprintf(stderr,
                 "\n== Radler phase timers (seconds, summed over "
                 "threads unless 'wall') ==\n");
    for (size_t i = 0; i != static_cast<size_t>(Phase::kCount); ++i) {
      std::fprintf(stderr, "%-26s %10.2f s  %8lld calls\n", kNames[i],
                   nanoseconds_[i].load() * 1e-9,
                   static_cast<long long>(calls_[i].load()));
    }
    std::fprintf(stderr,
                 "sub-image runs: %zu first, %zu same dimensions as previous, "
                 "%zu changed\n",
                 first_runs_, same_dimensions_, changed_dimensions_);
    std::fprintf(stderr,
                 "model components: %lld, model updates: %lld direct, %lld "
                 "fft\n",
                 static_cast<long long>(counts_[0].load()),
                 static_cast<long long>(counts_[1].load()),
                 static_cast<long long>(counts_[2].load()));
    std::fprintf(stderr, "residual subtractions: %lld direct, %lld fft\n",
                 static_cast<long long>(counts_[3].load()),
                 static_cast<long long>(counts_[4].load()));
  }

 private:
  PhaseTimers() = default;
  std::mutex mutex_;
  std::map<const void*, std::pair<size_t, size_t>> dimensions_;
  size_t first_runs_ = 0;
  size_t same_dimensions_ = 0;
  size_t changed_dimensions_ = 0;
  std::array<std::atomic<int64_t>, static_cast<size_t>(Phase::kCount)>
      nanoseconds_{};
  std::array<std::atomic<int64_t>, static_cast<size_t>(Phase::kCount)> calls_{};
  std::array<std::atomic<int64_t>, static_cast<size_t>(Counter::kCount)>
      counts_{};
};

class ScopedPhaseTimer {
 public:
  explicit ScopedPhaseTimer(Phase phase)
      : phase_(phase), start_(std::chrono::steady_clock::now()) {}
  ~ScopedPhaseTimer() {
    const auto end = std::chrono::steady_clock::now();
    PhaseTimers::Instance().Add(
        phase_,
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_)
            .count());
  }

 private:
  Phase phase_;
  std::chrono::steady_clock::time_point start_;
};

}  // namespace radler::utils

#endif
