#ifndef FAST_WFC_PROPAGATOR_HPP_
#define FAST_WFC_PROPAGATOR_HPP_

#include "direction.hpp"
#include "utils/array3D.hpp"
#include <tuple>
#include <vector>
#include <array>
#include <random>
#include <set>


class Wave;

/**
 * Propagate information about patterns in the wave.
 */
class Propagator {
public:
  using PropagatorState = std::vector<std::array<std::vector<unsigned>, 4>>;
  using NeghborWeights = std::vector<std::array<std::vector<double>, 4>>;

/**
 * True if the wave and the output is toric.
 */
private:
  /**
   * The size of the patterns.
   */
  const std::size_t patterns_size;
  const bool periodic_output;
  std::minstd_rand gen;

  /**
   * propagator[pattern1][direction] contains all the patterns that can
   * be placed in next to pattern1 in the direction direction.
   */
  PropagatorState propagator_state;
  NeghborWeights neghbor_weights;

  /**
   * The wave width and height.
   */
  const unsigned wave_width;
  const unsigned wave_height;
  std::set<unsigned > border_list;

    /**
     * All the tuples (y, x, pattern) that should be propagated.
     * The tuple should be propagated when wave.get(y, x, pattern) is set to
     * false.
     */
  std::vector<std::tuple<unsigned, unsigned, unsigned>> propagating;

  /**
   * All the tuples (y, x, pattern) that its neghbour should be propagated.
   * The tuple should be propagated when cells is collapsed to a specific pattern.
   */
  std::vector<std::tuple<unsigned, unsigned, unsigned>> neb_propagating;

  /**
   * compatible.get(y, x, pattern)[direction] contains the number of patterns
   * present in the wave that can be placed in the cell next to (y,x) in the
   * opposite direction of direction without being in contradiction with pattern
   * placed in (y,x). If wave.get(y, x, pattern) is set to false, then
   * compatible.get(y, x, pattern) has every element negative or null
   */
  Array3D<std::array<int, 4>> compatible;

  /**
   * Initialize compatible.
   */
  void init_compatible() noexcept;

public:
    void clean_propagating() noexcept{
        propagating.clear();
    }
  /**
   * Constructor building the propagator and initializing compatible.
   */
   Propagator(unsigned wave_height, unsigned wave_width, bool periodic_output,
             PropagatorState propagator_state, std::minstd_rand &gen) noexcept
      : patterns_size(propagator_state.size()),
        propagator_state(propagator_state), wave_width(wave_width),
        wave_height(wave_height), periodic_output(periodic_output),
        compatible(wave_height, wave_width, patterns_size) {
    init_compatible();
      for(int x = 0; x < wave_width; x++){
          for(int y = 0; y < wave_height; y++){
              if(x == 0 || y ==0|| x == wave_width-1 || y == wave_height-1){
                  border_list.insert(y*wave_width + x);
              }
          }
      }
  }

    /**
   * Constructor building the propagator and initializing compatible.
   */
   Propagator(unsigned wave_height, unsigned wave_width, bool periodic_output,
               NeghborWeights neghbor_weights,
               PropagatorState propagator_state, std::minstd_rand &gen) noexcept
            : patterns_size(propagator_state.size()),
              neghbor_weights(neghbor_weights),
              propagator_state(propagator_state), wave_width(wave_width),
              wave_height(wave_height), periodic_output(periodic_output),
              gen(gen),
              compatible(wave_height, wave_width, patterns_size) {
        init_compatible();
        for(int x = 0; x < wave_width; x++){
            for(int y = 0; y < wave_height; y++){
                if(x == 0 || y ==0|| x == wave_width-1 || y == wave_height-1){
                    border_list.insert(y*wave_width + x);
                }
            }
        }
    }

  /**
   * Add an element to the propagator.
   * This function is called when wave.get(y, x, pattern) is set to false.
   */
  void add_to_propagator(unsigned y, unsigned x, unsigned pattern) noexcept {
    // All the direction are set to 0, since the pattern cannot be set in (y,x).
    std::array<int, 4> temp = {};
    compatible.get(y, x, pattern) = temp;
    propagating.emplace_back(y, x, pattern);
  }

void add_to_neb_propagator(unsigned y, unsigned x, unsigned pattern) noexcept {
        neb_propagating.emplace_back(y, x, pattern);
}

  PropagatorState get_propagator_state() {
      return propagator_state;
  }

  NeghborWeights get_neghbor_weights() {
      return neghbor_weights;
  }

  void set_propagator_state(PropagatorState base_propagator_state) {
      propagator_state = base_propagator_state;
  }

  void set_neghbor_weights(NeghborWeights base_neghbor_weights) {
      neghbor_weights = base_neghbor_weights;
  }

  const auto get_border_list() {
      return border_list;
  }
  /**
   * Propagate the information given with add_to_propagator.
   */
  void propagate(Wave &wave) noexcept;
  void neghbour_propagate(Wave &wave, std::set<unsigned> ramp_ids) noexcept;
};

#endif // FAST_WFC_PROPAGATOR_HPP_
