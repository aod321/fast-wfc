#include "propagator.hpp"
#include "wave.hpp"

namespace {
    // C++ template,  if an element of a vector exists, set
    template<typename T>
    bool exists(const std::set<T>& v, const T& e) {
        return v.find(e) != v.end();
    }

    // C++ template,  if an element of a vector exists, vector
    template<typename T>
    bool exists(const std::vector<T>& v, const T& e) {
        return std::find(v.begin(), v.end(), e) != v.end();
    }
}

void Propagator::init_compatible() noexcept {
  std::array<int, 4> value;
  // We compute the number of pattern compatible in all directions.
  for (unsigned y = 0; y < wave_height; y++) {
    for (unsigned x = 0; x < wave_width; x++) {
      for (unsigned pattern = 0; pattern < patterns_size; pattern++) {
        for (int direction = 0; direction < 4; direction++) {
          value[direction] =
            static_cast<unsigned>(propagator_state[pattern][get_opposite_direction(direction)]
                                  .size());
        }
        compatible.get(y, x, pattern) = value;
      }
    }
  }
}

void Propagator::propagate(Wave &wave) noexcept {

  // We propagate every element while there is element to propagate.
  while (propagating.size() != 0) {

    // The cell and pattern that has been set to false.
    unsigned y1, x1, pattern;
    std::tie(y1, x1, pattern) = propagating.back();
    propagating.pop_back();

    // We propagate the information in all 4 directions.
    for (unsigned direction = 0; direction < 4; direction++) {

      // We get the next cell in the direction direction.
      int dx = directions_x[direction];
      int dy = directions_y[direction];
      int x2, y2;
      if (periodic_output) {
        x2 = ((int)x1 + dx + (int)wave.width) % wave.width;
        y2 = ((int)y1 + dy + (int)wave.height) % wave.height;
      } else {
        x2 = x1 + dx;
        y2 = y1 + dy;
        if (x2 < 0 || x2 >= (int)wave.width) {
          continue;
        }
        if (y2 < 0 || y2 >= (int)wave.height) {
          continue;
        }
      }

      // The index of the second cell, and the patterns compatible
      unsigned i2 = x2 + y2 * wave.width;
      std::vector<unsigned> &patterns =
        propagator_state[pattern][direction];

      // For every pattern that could be placed in that cell without being in
      // contradiction with pattern1
      for (auto it = patterns.begin(), it_end = patterns.end(); it < it_end;
           ++it) {

        // We decrease the number of compatible patterns in the opposite
        // direction If the pattern was discarded from the wave, the element
        // is still negative, which is not a problem
        std::array<int, 4> &value = compatible.get(y2, x2, *it);
        value[direction]--;

        // If the element was set to 0 with this operation, we need to remove
        // the pattern from the wave, and propagate the information
        if (value[direction] == 0) {
          add_to_propagator(y2, x2, *it);
          wave.set(i2, *it, false);
        }
      }
    }
  }
}


void Propagator::neghbour_propagate(Wave &wave, std::set<unsigned> ramp_ids) noexcept {
    auto &cell_weights = wave.get_cell_partterns_weights();
    // We propagate every element while there is element to propagate.
    while (neb_propagating.size() != 0) {

        // The cell and pattern that just collapsed
        auto [y1, x1, pattern] = neb_propagating.back();
        neb_propagating.pop_back();

        // We propagate the information in all 4 directions.
        for (unsigned direction = 0; direction < 4; direction++) {

            // We get the next cell in the direction.
            int dx = directions_x[direction];
            int dy = directions_y[direction];
            int x2, y2;
            if (periodic_output) {
                x2 = ((int)x1 + dx + (int)wave.width) % wave.width;
                y2 = ((int)y1 + dy + (int)wave.height) % wave.height;
            } else {
                x2 = x1 + dx;
                y2 = y1 + dy;
                if (x2 < 0 || x2 >= (int)wave.width) {
                    continue;
                }
                if (y2 < 0 || y2 >= (int)wave.height) {
                    continue;
                }
            }

            // The index of the second cell, and the pattern neghbour weights
            unsigned i2 = x2 + y2 * wave.width;
            // If the cell is at border, it shouldn't be a ramp
            if(exists(border_list, i2)){
                for(auto ramp_id : ramp_ids){
                    cell_weights.get(i2, ramp_id) = 0.0;
                }
            }

            std::vector<double> &weights =
                    neghbor_weights[pattern][direction];
            // For every parten that related to the collapsed pattern, set its weights
            for (auto it = weights.begin(), it_end = weights.end(); it < it_end;
                 ++it) {
                    if(wave.get(i2, it-weights.begin())){
                        cell_weights.get(i2, it - weights.begin()) *= *it;
                    }
                }
        }
    }
    wave.set_cell_partterns_weights(cell_weights);
}