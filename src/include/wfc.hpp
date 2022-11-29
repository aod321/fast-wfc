#ifndef FAST_WFC_WFC_HPP_
#define FAST_WFC_WFC_HPP_

#include <optional>
#include <random>
#include <set>


#include "utils/array2D.hpp"
#include "propagator.hpp"
#include "wave.hpp"

/**
 * Class containing the generic WFC algorithm.
 */
class WFC {
protected:
  /**
   * The distribution of the patterns as given in input.
   */
  /**
   * The wave, indicating which patterns can be put in which cell.
   */
  Wave wave;

private:
  /**
   * The random number generator.
   */
  std::minstd_rand gen;

  /**
   * The number of distinct patterns.
   */
  const size_t nb_patterns;

  /**
   * The propagator, used to propagate the information in the wave.
   */
  Propagator propagator;

  std::set<unsigned > ramp_ids;
std::vector<double> patterns_frequencies;


    /**
   * Transform the wave to a valid output (a 2d array of patterns that aren't in
   * contradiction). This function should be used only when all cell of the wave
   * are defined.
   */
  Array2D<unsigned> wave_to_output() const noexcept;

public:
  /**
   * Basic constructor initializing the algorithm.
   */
  WFC(bool periodic_output, int seed, std::vector<double> patterns_frequencies,
      Propagator::PropagatorState propagator, unsigned wave_height,
      unsigned wave_width)
    noexcept;

    /**
   * Basic constructor initializing the algorithm.
   */
    WFC(bool periodic_output, int seed, std::vector<double> patterns_frequencies,
        Propagator::PropagatorState propagator,
        unsigned wave_height,
        unsigned wave_width, Propagator::NeghborWeights neghbor_weights,
        const std::set<unsigned>& ramp_ids)
    noexcept;

  Wave get_wave() const noexcept {
    return wave;
  }

  void set_seed(unsigned seed) noexcept {
    gen.seed(seed);
  }

  Propagator get_propagator() const noexcept {
    return propagator;
  }


  void remove_border_ramp();
  /**
   * Run the algorithm, and return a result if it succeeded.
   */
  std::optional<Array2D<unsigned>> run() noexcept;

  void mutate(Wave base_wave, double new_weight) noexcept;

  void set_propagator(Propagator base_propagator) noexcept {
      propagator.set_propagator_state(base_propagator.get_propagator_state());
      propagator.set_neghbor_weights(base_propagator.get_neghbor_weights());
  }

  /**
   * Return value of observe.
   */
  enum ObserveStatus {
    success,    // WFC has finished and has succeeded.
    failure,    // WFC has finished and failed.
    to_continue // WFC isn't finished.
  };

  /**
   * Define the value of the cell with lowest entropy.
   */
  ObserveStatus observe() noexcept;

  /**
   * Propagate the information of the wave.
   */
  void propagate() noexcept { propagator.propagate(wave); }

  /**
   * Remove pattern from cell (i,j).
   */
  void remove_wave_pattern(unsigned i, unsigned j, unsigned pattern) noexcept {
    if (wave.get(i, j, pattern)) {
      wave.set(i, j, pattern, false);
      propagator.add_to_propagator(i, j, pattern);
    }
  }
};

#endif // FAST_WFC_WFC_HPP_
