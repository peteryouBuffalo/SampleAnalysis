#include "LumiMaskFilter.h"

#include <algorithm>

#include <yaml-cpp/yaml.h>
//#include <yaml.h>


//LumiMaskFilter::LumiMaskFilter(std::string const &path)
//    : cachedRun_{-1}, cachedLS_{-1} {
//  LoadLumiMask(path);
//}

void LumiMaskFilter::Set(std::string const &path) {
  cachedRun_ = -1;
  cachedLS_ = -1;
  LoadLumiMask(path);
}

bool LumiMaskFilter::Pass(int64_t run, int64_t ls) const {
  // Events are read consequently, which means that they will often come from
  // the same LS as in the previous call. Implement a caching to take advantage
  // of this.
  if (run != cachedRun_ or ls != cachedLS_) {
    cachedResult_ = PassImpl(run, ls);
    cachedRun_ = run;
    cachedLS_ = ls;
  }

  return cachedResult_;
}


void LumiMaskFilter::SetSourceBranches(TTreeReader &reader) {
  srcRunNumber_.emplace(reader, "run");
  srcLuminosityBlock_.emplace(reader, "luminosityBlock");
}


void LumiMaskFilter::LoadLumiMask(std::string const &path) {
  // Since JSON is a strict subset of YAML, can use a YAML parser
  auto const parsedMask = YAML::LoadFile(path);

  if (not parsedMask.IsMap()) {
    std::ostringstream message;
    message << "File \"" << path << "\" does not contain a mapping.";
    throw std::runtime_error(message.str());
  }

  for (auto const &runInfo: parsedMask) {
    auto const run = runInfo.first.as<int64_t>();
    std::vector<LumiRange> lumis;

    auto const &lumiRangesNode = runInfo.second;

    if (not lumiRangesNode.IsSequence()) {
      std::ostringstream message;
      message << "In file \"" << path << "\" run " << run
          << " is not mapped to a sequence.";
      throw std::runtime_error(message.str());
    }

    for (auto const &rangeNode: lumiRangesNode) {
      if (not rangeNode.IsSequence() or rangeNode.size() != 2) {
        std::ostringstream message;
        message << "In file \"" << path << "\" at least one of lumi ranges for "
            "run " << run << " is not a sequence or contains a wrong number of "
            "entries.";
        throw std::runtime_error(message.str());
      }

      auto const lsStart = rangeNode[0].as<int64_t>();
      auto const lsStop = rangeNode[1].as<int64_t>();

      if (lsStart > lsStop) {
        std::ostringstream message;
        message << "In file \"" << path << "\" entry for run " << run
            << " contains a lumi range with a wrong ordering: ["
            << lsStart << ", " << lsStop << "].";
        throw std::runtime_error(message.str());
      }

      lumis.emplace_back(lsStart, lsStop);
    }

    // Make sure that the lumi ranges are sorted
    std::sort(
        lumis.begin(), lumis.end(),
        [](auto const &lhs, auto const &rhs){return (lhs.first < rhs.first);});

    // Check for overlaps
    for (int i = 0; i < int(lumis.size()) - 1; ++i) {
      if (lumis[i].second > lumis[i + 1].first) {
        std::ostringstream message;
        message << "In file \"" << path << "\" entry for run " << run
            << " contains overlapping ranges of LS.";
        throw std::runtime_error(message.str());
      }
    }
    if (lumiMask_.find(run) != lumiMask_.end()) {
      std::ostringstream message;
      message << "File \"" << path << "\" contains multiple entries for run "
          << run << ".";
      throw std::runtime_error(message.str());
    }

    lumiMask_.emplace(std::make_pair(run, lumis));
  }
}


bool LumiMaskFilter::PassImpl(int64_t run, int64_t ls) const {
  auto const res = lumiMask_.find(run);

  if (res == lumiMask_.end())
    return false;

  auto const &lumis = res->second;

  // Find the first LS whose right boundary is larger than or equal to the given
  // LS. This is the only LS range that could contain the given event.
  auto const rangeIt = std::lower_bound(
      lumis.begin(), lumis.end(), ls,
      [](LumiRange const &r, int64_t l){return (r.second < l);});

  if (rangeIt == lumis.end())
    return false;

  // Only left to check the left boundary of the LS range
  return (rangeIt->first <= ls);
}

