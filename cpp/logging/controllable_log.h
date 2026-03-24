// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RADLER_LOGGING_CONTROLLABLE_LOG_H_
#define RADLER_LOGGING_CONTROLLABLE_LOG_H_

#include <mutex>
#include <string>
#include <vector>

#include <aocommon/logger.h>

namespace radler::logging {

class ControllableLog final : public aocommon::LogReceiver {
 public:
  ControllableLog(std::mutex* mutex)
      : LogReceiver(), _mutex(mutex), _isMuted(false), _isActive(true) {}

  ControllableLog(const ControllableLog&) = default;
  ControllableLog(ControllableLog&&) = default;
  ControllableLog& operator=(const ControllableLog&) = default;
  ControllableLog& operator=(ControllableLog&&) = default;

  void Mute(bool mute) { _isMuted = mute; }
  bool IsMuted() const { return _isMuted; }

  void Activate(bool active) { _isActive = active; }
  bool IsActive() const { return _isActive; }

  void SetTag(const std::string& tag) { _tag = tag; }
  void SetOutputOnce(const std::string& str) { _outputOnce = str; }

 protected:
  void Output(aocommon::LogLevel level, const std::string& str) override {
    if (!Skip(level) && !str.empty()) {
      std::lock_guard<std::mutex> lock(*_mutex);

      _lineBuffer += str;
      if (_lineBuffer.back() == '\n') {
        if (!_outputOnce.empty()) {
          LogReceiver::Output(level, _outputOnce);
          _outputOnce.clear();
        }
        LogReceiver::Output(level, _tag);
        LogReceiver::Output(level, _lineBuffer);
        _lineBuffer.clear();
      }
    }
  }
  void Flush(aocommon::LogLevel level) override {
    if (!Skip(level)) {
      std::lock_guard<std::mutex> lock(*_mutex);
      LogReceiver::Flush(level);
    }
  }

 private:
  bool Skip(aocommon::LogLevel level) const {
    return ((level == aocommon::LogLevel::kDebug ||
             level == aocommon::LogLevel::kInfo) &&
            _isMuted) ||
           (level == aocommon::LogLevel::kDebug &&
            !aocommon::Logger::IsVerbose());
  }
  std::mutex* _mutex;
  std::string _tag;
  bool _isMuted;
  bool _isActive;
  std::string _lineBuffer;
  std::string _outputOnce;
};

}  // namespace radler::logging

#endif
