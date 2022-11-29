#include "wfc.hpp"
#include <limits>
#include <utility>


namespace {
    /**
     * Normalize a vector so the sum of its elements is equal to 1.0f
     */
    std::vector<double>& normalize(std::vector<double>& v) {
        double sum_weights = 0.0;
        for(double weight: v) {
            sum_weights += weight;
        }

        double inv_sum_weights = 1.0/sum_weights;
        for(double& weight: v) {
            weight *= inv_sum_weights;
        }

        return v;
    }

    // C++ template,  if an element of a vector exists, vector
    template<typename T>
    bool exists(const std::vector<T>& v, const T& e) {
        return std::find(v.begin(), v.end(), e) != v.end();
    }

    // C++ template,  if an element of a vector exists, set
    template<typename T>
    bool exists(const std::set<T>& v, const T& e) {
        return v.find(e) != v.end();
    }

}


Array2D<unsigned> WFC::wave_to_output() const noexcept {
  Array2D<unsigned> output_patterns(wave.height, wave.width);
  for (unsigned i = 0; i < wave.size; i++) {
    for (unsigned k = 0; k < nb_patterns; k++) {
      if (wave.get(i, k)) {
        output_patterns.data[i] = k;
      }
    }
  }
  return output_patterns;
}

/*
 *  Mutate the cell weights
 */
void WFC::mutate(Wave base_wave, double new_weight) noexcept{
    auto cell_weights = base_wave.get_cell_partterns_weights();
    for (unsigned i = 0; i < base_wave.size; i++) {
        for (unsigned k = 0; k < nb_patterns; k++) {
        if (base_wave.get(i, k)) {
            cell_weights.get(i,k) *= new_weight;
        }
        }
    }
    wave.set_cell_partterns_weights(cell_weights);
}

WFC::WFC(bool periodic_output, int seed,
         std::vector<double> patterns_frequencies,
         Propagator::PropagatorState propagator, unsigned wave_height,
         unsigned wave_width)
noexcept
        : gen(seed),
          patterns_frequencies(patterns_frequencies),
          wave(wave_height, wave_width, patterns_frequencies),
          nb_patterns(propagator.size()),
          propagator(wave.height, wave.width, periodic_output, propagator, gen) {}

WFC::WFC(bool periodic_output, int seed, std::vector<double> patterns_frequencies,
         Propagator::PropagatorState propagator,
         unsigned wave_height,
         unsigned wave_width, Propagator::NeghborWeights neghbor_weights, const std::set<unsigned>& ramp_ids)
noexcept
        : gen(seed),
          wave(wave_height, wave_width, patterns_frequencies),
          nb_patterns(propagator.size()),
          ramp_ids(ramp_ids),
          propagator(wave.height, wave.width, periodic_output, std::move(neghbor_weights), propagator,gen){
}

std::optional<Array2D<unsigned>> WFC::run() noexcept {
  while (true) {

    // Define the value of an undefined cell.
    ObserveStatus result = observe();

    // Check if the algorithm has terminated.
    if (result == failure) {
      return std::nullopt;
    } else if (result == success) {
      return wave_to_output();
    }

    // Propagate the information.
    propagator.propagate(wave);
    propagator.neghbour_propagate(wave, ramp_ids);
  }
}


WFC::ObserveStatus WFC::observe() noexcept {
    auto border_list = propagator.get_border_list();
    // Get the cell with the lowest entropy.
    int argmin = wave.get_min_entropy(gen);

    // If there is a contradiction, the algorithm has failed.
    if (argmin == -2) {
      return failure;
    }

    // If the lowest entropy is 0, then the algorithm has succeeded and
    // finished.
    if (argmin == -1) {
      wave_to_output();
      return success;
    }
    // Choose an element according to the pattern distribution
    double s = 0;
    auto normlized_cell_partterns_weights = wave.get_normlized_cell_partterns_weights();
    for (unsigned k = 0; k < nb_patterns; k++) {
        s += wave.get(argmin, k) ? normlized_cell_partterns_weights.get(argmin, k) : 0;
    }

    std::uniform_real_distribution<> dis(0, s);
    double random_value = dis(gen);
    size_t chosen_value = nb_patterns - 1;
    // Avoid of boundary ramps
    // ---------
    auto while_flag = false;
    // if the cell is in the border list, then the chosen value must not be a ramp
    if(exists<unsigned>(border_list, argmin)) while_flag = true;
    do {
        for (unsigned k = 0; k < nb_patterns; k++) {
            random_value -= wave.get(argmin, k) ? normlized_cell_partterns_weights.get(argmin, k) : 0;
            if (random_value <= 0) {
                chosen_value = k;
                break;
            }
        }
    }while(while_flag && exists<unsigned >(ramp_ids, chosen_value)); // until the chosen value is not a ramp
    propagator.add_to_neb_propagator(argmin / wave.width, argmin % wave.width,chosen_value);
    // And define the cell with the pattern.
    for (unsigned k = 0; k < nb_patterns; k++) {
      if (wave.get(argmin, k) != (k == chosen_value)) {
          propagator.add_to_propagator(argmin / wave.width, argmin % wave.width,
                                       k);
          wave.set(argmin, k, false);
      }
    }

    return to_continue;
  }
