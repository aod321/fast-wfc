#include "wave.hpp"

#include <limits>

namespace {
    /**
 * Normalize a vector so the sum of its elements is equal to 1.0f
 */
    std::vector<double> normalize(std::vector<double> v) {
        double sum = 0.0;
        std::vector<double> result;
        for (auto e : v) {
            sum += e;
        }
        for (auto e : v) {
            result.push_back(e /= sum);
        }
        return result;
    }

/**
 * Return distribution * log(distribution).
 */
std::vector<double>
get_plogp(const std::vector<double> &distribution) noexcept {
  std::vector<double> plogp;
  for (unsigned i = 0; i < distribution.size(); i++) {
    plogp.push_back(distribution[i] * log(distribution[i]));
  }
  return plogp;
}

/**
 * Return min(v) / 2.
 */
double get_min_abs_half(const std::vector<double> &v) noexcept {
  double min_abs_half = std::numeric_limits<double>::infinity();
  for (unsigned i = 0; i < v.size(); i++) {
    min_abs_half = std::min(min_abs_half, std::abs(v[i] / 2.0));
  }
  return min_abs_half;
}

} // namespace


Wave::Wave(unsigned height, unsigned width,
     const std::vector<double> &patterns_frequencies) noexcept
  : patterns_frequencies(patterns_frequencies),
    plogp_patterns_frequencies(get_plogp(patterns_frequencies)),
    min_abs_half_plogp(get_min_abs_half(plogp_patterns_frequencies)),
    is_impossible(false), nb_patterns(patterns_frequencies.size()),
    data(width * height, nb_patterns, 1), width(width), height(height),
    cell_partterns_weights(width * height, patterns_frequencies.size()),
    normalized_cell_partterns_weights(width * height, patterns_frequencies.size()),
    size(height * width) {
  // Initialize the memoisation of entropy.
  double base_entropy = 0;
  double base_s = 0;
  for (unsigned i = 0; i < nb_patterns; i++) {
    base_entropy += plogp_patterns_frequencies[i];
    base_s += patterns_frequencies[i];
  }
  double log_base_s = log(base_s);
  double entropy_base = log_base_s - base_entropy / base_s;
  memoisation.plogp_sum = std::vector<double>(width * height, base_entropy);
  memoisation.sum = std::vector<double>(width * height, base_s);
  memoisation.log_sum = std::vector<double>(width * height, log_base_s);
  memoisation.nb_patterns =
    std::vector<unsigned>(width * height, static_cast<unsigned>(nb_patterns));
  memoisation.entropy = std::vector<double>(width * height, entropy_base);
  auto normalized_patterns_frequencies = normalize(patterns_frequencies);
  for (unsigned i = 0; i < width * height; i++) {
    for (unsigned j = 0; j < nb_patterns; j++) {
        cell_partterns_weights.get(i, j) = patterns_frequencies[j];
        normalized_cell_partterns_weights.get(i, j) = normalized_patterns_frequencies[j];
    }
  }
}

void Wave::set_cell_partterns_weights(Array2D<double> weights) noexcept{
    cell_partterns_weights = weights;
    for (unsigned i = 0; i < width * height; i++) {
        std::vector<double> normalized_weights;
        for (unsigned j = 0; j < nb_patterns; j++) {
            normalized_weights.push_back(cell_partterns_weights.get(i, j));
        }
        normalized_weights = normalize(normalized_weights);
        for (unsigned j = 0; j < nb_patterns; j++) {
            normalized_cell_partterns_weights.get(i, j) = normalized_weights[j];
        }
    }
}

void Wave::set(unsigned index, unsigned pattern, bool value) noexcept {
  bool old_value = data.get(index, pattern);
  // If the value isn't changed, nothing needs to be done.
  if (old_value == value) {
    return;
  }
  // Otherwise, the memoisation should be updated.
  data.get(index, pattern) = value;
  memoisation.plogp_sum[index] -= plogp_patterns_frequencies[pattern];
  memoisation.sum[index] -= patterns_frequencies[pattern];
  memoisation.log_sum[index] = log(memoisation.sum[index]);
  memoisation.nb_patterns[index]--;
  memoisation.entropy[index] =
    memoisation.log_sum[index] -
    memoisation.plogp_sum[index] / memoisation.sum[index];
  // If there is no patterns possible in the cell, then there is a
  // contradiction.
  if (memoisation.nb_patterns[index] == 0) {
    is_impossible = true;
  }
}


int Wave::get_min_entropy(std::minstd_rand &gen) const noexcept {
  if (is_impossible) {
    return -2;
  }

  std::uniform_real_distribution<> dis(0, min_abs_half_plogp);

  // The minimum entropy (plus a small noise)
  double min = std::numeric_limits<double>::infinity();
  int argmin = -1;

  for (unsigned i = 0; i < size; i++) {

    // If the cell is decided, we do not compute the entropy (which is equal
    // to 0).
    double nb_patterns_local = memoisation.nb_patterns[i];
    if (nb_patterns_local == 1) {
      continue;
    }

    // Otherwise, we take the memoised entropy.
    double entropy = memoisation.entropy[i];

    // We first check if the entropy is less than the minimum.
    // This is important to reduce noise computation (which is not
    // negligible).
    if (entropy <= min) {

      // Then, we add noise to decide randomly which will be chosen.
      // noise is smaller than the smallest p * log(p), so the minimum entropy
      // will always be chosen.
      double noise = dis(gen);
      if (entropy + noise < min) {
        min = entropy + noise;
        argmin = i;
      }
    }
  }

  return argmin;
}
