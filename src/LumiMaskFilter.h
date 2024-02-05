#ifndef SHEARS_PRUNER_BIN_LUMIMASKFILTER_H_
#define SHEARS_PRUNER_BIN_LUMIMASKFILTER_H_

#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>


/**
 * \brief Implements filtering accoring to a luminosity mask
 *
 * The luminosity mask must be given as a JSON file in the standard PPD format.
 * The filter can be used in two ways: either relevant part of the event ID is
 * provided explicitly to \ref Pass(int64_t, int64_t) or ID of the current event
 * is read in \ref Pass(). In the latter case source branches for event ID must
 * be set up with SetSourceBranches before the first call to \ref Pass().
 */
class LumiMaskFilter {
 public:
  /// Constructor from a path to JSON file with luminosity mask
  //LumiMaskFilter(std::string const &path);
  LumiMaskFilter() {};
  void Set(std::string const &path);


  /// Check if given run and luminosity section are included in the mask
  bool Pass(int64_t run, int64_t ls) const;

  /**
   * \brief Check if the current event is included in the mask
   *
   * Before this method can be used, SetSourceBranches must be called to set up
   * branches from which event ID will be read.
   */
  bool Pass() const {
    return Pass(**srcRunNumber_, **srcLuminosityBlock_);
  }

  /**
   * \brief Set up source branches to read event ID
   *
   * Needed to use argumentless version of \ref Pass().
   */
  void SetSourceBranches(TTreeReader &reader);

 private:
  using LumiRange = std::pair<int64_t, int64_t>;

  /// Read luminosity mask from given JSON file
  void LoadLumiMask(std::string const &path);

  /// Implementation for the check against the luminosity mask
  bool PassImpl(int64_t run, int64_t ls) const;

  /**
   * \brief Luminosity mask
   *
   * Vector of LS ranges is sorted.
   */
  std::map<int64_t, std::vector<LumiRange>> lumiMask_;

  mutable std::optional<TTreeReaderValue<UInt_t>> srcRunNumber_;
  mutable std::optional<TTreeReaderValue<UInt_t>> srcLuminosityBlock_;

  /// Cached run and LS numbers from most recent call to Pass
  mutable int64_t cachedRun_, cachedLS_;

  /// Cached result from most recent call to \ref Pass(int64_t, int64_t)
  mutable bool cachedResult_;
};

#endif  // SHEARS_PRUNER_BIN_LUMIMASKFILTER_H_

