#ifndef FAST_WFC_XML_WFC_HPP_
#define FAST_WFC_XML_WFC_HPP_

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include "time.h"

#include "tiling_wfc.hpp"
#include "utils/array3D.hpp"
#include "wfc.hpp"
#include "external/rapidxml.hpp"
#include "image.hpp"
#include "rapidxml_utils.hpp"
#include "utils.hpp"
#include <unordered_set>


using namespace rapidxml;
using namespace std;
using namespace tilingWFC;

namespace XMLWFC
{
/**
 * Get a random seed.
 * This function use random_device on linux, but use the C rand function for
 * other targets. This is because, for instance, windows don't implement
 * random_device non-deterministically.
 */
int get_random_seed() {
#ifdef __linux__
    return random_device()();
#else
    return rand();
#endif
}

/**
 * Transform a symmetry name into its Symmetry enum
 */
Symmetry to_symmetry(const string &symmetry_name) {
    if (symmetry_name == "X") {
        return Symmetry::X;
    }
    if (symmetry_name == "T") {
        return Symmetry::T;
    }
    if (symmetry_name == "I") {
        return Symmetry::I;
    }
    if (symmetry_name == "L") {
        return Symmetry::L;
    }
    if (symmetry_name == "\\") {
        return Symmetry::backslash;
    }
    if (symmetry_name == "P") {
        return Symmetry::P;
    }
    throw symmetry_name + "is an invalid Symmetry";
}

/**
 * Read the names of the tiles in the subset in a tiling WFC problem
 */
std::optional<unordered_set<string>> read_subset_names(xml_node<> *root_node,
                                                       const string &subset) {
    unordered_set<string> subset_names;
    xml_node<> *subsets_node = root_node->first_node("subsets");
    if (!subsets_node) {
        return std::nullopt;
    }
    xml_node<> *subset_node = subsets_node->first_node("subset");
    while (subset_node &&
           rapidxml::get_attribute(subset_node, "name") != subset) {
        subset_node = subset_node->next_sibling("subset");
    }
    if (!subset_node) {
        return std::nullopt;
    }
    for (xml_node<> *node = subset_node->first_node("tile"); node;
         node = node->next_sibling("tile")) {
        subset_names.insert(rapidxml::get_attribute(node, "name"));
    }
    return subset_names;
}

/**
 * Read all tiles for a tiling problem
 */
//unordered_map<string, Tile<Color>> read_tiles(xml_node<> *root_node,
std::vector<std::pair<string, Tile<Color>>> read_tiles(xml_node<> *root_node,
                                              const string &current_dir,
                                              const string &subset,
                                              unsigned size) {
    std::optional<unordered_set<string>> subset_names =
            read_subset_names(root_node, subset);
//    unordered_map<string, Tile<Color>> tiles;
    std::vector<std::pair<string, Tile<Color>>> tiles;
    xml_node<> *tiles_node = root_node->first_node("tiles");
    for (xml_node<> *node = tiles_node->first_node("tile"); node;
         node = node->next_sibling("tile")) {
        string name = rapidxml::get_attribute(node, "name");
        if (subset_names != nullopt &&
            subset_names->find(name) == subset_names->end()) {
            continue;
        }
        Symmetry symmetry =
                to_symmetry(rapidxml::get_attribute(node, "symmetry", "X"));
        double weight = stod(rapidxml::get_attribute(node, "weight", "1.0"));
        const std::string image_path = current_dir + "/" + name + ".png";
        optional<Array2D<Color>> image = read_image(image_path);

        if (image == nullopt) {
            vector<Array2D<Color>> images;
            for (unsigned i = 0; i < nb_of_possible_orientations(symmetry); i++) {
                const std::string image_path =
                        current_dir + "/" + name + " " + to_string(i) + ".png";
                optional<Array2D<Color>> image = read_image(image_path);
                if (image == nullopt) {
                    throw "Error while loading " + image_path;
                }
                if ((image->width != size) || (image->height != size)) {
                    throw "Image " + image_path + " has wrond size";
                }
                images.push_back(*image);
            }
            Tile<Color> tile = {images, symmetry, weight};
//            tiles.insert({name, tile});
            tiles.push_back({name, tile});
        } else {
            if ((image->width != size) || (image->height != size)) {
                throw "Image " + image_path + " has wrong size";
            }

            Tile<Color> tile(*image, symmetry, weight);
//            tiles.insert({name, tile});
            tiles.push_back({name, tile});
        }
    }

    return tiles;
}


/**A
 * Read the neighbors constraints for a tiling problem.
 * A value {t1,o1,t2,o2} means that the tile t1 with orientation o1 can be
 * placed at the right of the tile t2 with orientation o2.
 */
//vector<tuple<string, unsigned, string, unsigned, double>>
auto
read_neighbors_weights(xml_node<> *root_node) {
    vector<tuple<string, unsigned, string, unsigned, double>> neighbors;
    xml_node<> *neighbor_node = root_node->first_node("neighbors");
    for (xml_node<> *node = neighbor_node->first_node("neighbor"); node;
         node = node->next_sibling("neighbor")) {
        string left = rapidxml::get_attribute(node, "left");
        string::size_type left_delimiter = left.find(" ");
        string left_tile = left.substr(0, left_delimiter);
        unsigned left_orientation = 0;
        if (left_delimiter != string::npos) {
            left_orientation = stoi(left.substr(left_delimiter, string::npos));
        }

        string right = rapidxml::get_attribute(node, "right");
        string::size_type right_delimiter = right.find(" ");
        string right_tile = right.substr(0, right_delimiter);
        unsigned right_orientation = 0;
        if (right_delimiter != string::npos) {
            right_orientation = stoi(right.substr(right_delimiter, string::npos));
        }
        string weight_str = rapidxml::get_attribute(node, "weight", "1.0");
        double weight = stod(weight_str);
        neighbors.push_back(
                {left_tile, left_orientation, right_tile, right_orientation, weight});
    }
    return neighbors;
}


/**
 * Read an instance of a tiling WFC problem.
 */
//std::optional<vector<unsigned >> read_simpletiled_instance(xml_node<> *node,
auto read_simpletiled_instance(xml_node<> *node,
                               const string &current_dir) noexcept {
    string name = rapidxml::get_attribute(node, "name");
    string subset = rapidxml::get_attribute(node, "subset", "tiles");
    bool periodic_output =
            (rapidxml::get_attribute(node, "periodic", "False") == "True");
    unsigned width = stoi(rapidxml::get_attribute(node, "width", "48"));
    unsigned height = stoi(rapidxml::get_attribute(node, "height", "48"));

//    cout << name << " " << subset << " started!" << endl;

    ifstream config_file("samples/" + name + "/data.xml");
    vector<char> buffer((istreambuf_iterator<char>(config_file)),
                        istreambuf_iterator<char>());
    buffer.push_back('\0');
    xml_document<> data_document;
    data_document.parse<0>(&buffer[0]);
    xml_node<> *data_root_node = data_document.first_node("set");
    unsigned size = stoi(rapidxml::get_attribute(data_root_node, "size"));

    vector<pair<string, Tile<Color>>> tiles_map =
            read_tiles(data_root_node, current_dir + "/" + name, subset, size);
    unordered_map<string, unsigned> tiles_id;
    vector<Tile<Color>> tiles;
    unsigned id = 0;
    for (pair<string, Tile<Color>> tile : tiles_map) {
        tiles_id.insert({tile.first, id});
        tiles.push_back(tile.second);
        id++;
    }

    vector<tuple<string, unsigned, string, unsigned, double>> neighbors =
            read_neighbors_weights(data_root_node);
    vector<tuple<unsigned, unsigned, unsigned, unsigned, double>> neighbors_ids;
    for (auto neighbor : neighbors) {
        const string &neighbor1 = get<0>(neighbor);
        const int &orientation1 = get<1>(neighbor);
        const string &neighbor2 = get<2>(neighbor);
        const int &orientation2 = get<3>(neighbor);
        const double &weight = get<4>(neighbor);
        if (tiles_id.find(neighbor1) == tiles_id.end()) {
            continue;
        }
        if (tiles_id.find(neighbor2) == tiles_id.end()) {
            continue;
        }
        neighbors_ids.push_back(make_tuple(tiles_id[neighbor1], orientation1,
                                           tiles_id[neighbor2], orientation2, weight));
    }

    int seed = get_random_seed();
    return make_tuple(tiles, neighbors_ids, width, height, periodic_output, seed);
//    TilingWFC<Color> wfc(tiles, neighbors_ids, height, width, {periodic_output},
//                          seed);

}

vector<Tile<Color>> get_tiles(const string &config_path) noexcept{
    ifstream config_file(config_path);
    // if file not exists
    if(!config_file) {
        cout << "File " << config_path << " does not exist" << endl;
        return {};
    }
    vector<char> buffer((istreambuf_iterator<char>(config_file)),
                        istreambuf_iterator<char>());
    buffer.push_back('\0');
    xml_document<> document;
    document.parse<0>(&buffer[0]);

    xml_node<> *root_node = document.first_node("samples");
    string dir_path = get_dir(config_path) + "/" + "samples";
    auto [tiles,
            neighbors_ids,
            width, height,
            periodic_output, seed
    ] = read_simpletiled_instance(root_node->first_node("sample"), dir_path);
    return tiles;
}

/**
 * Build a inital Wave
 */
std::optional<Wave> build_wave(const string &config_path) noexcept {
        ifstream config_file(config_path);
        // if file not exists
        if(!config_file) {
            cout << "File " << config_path << " does not exist" << endl;
            return {};
        }
        vector<char> buffer((istreambuf_iterator<char>(config_file)),
                            istreambuf_iterator<char>());
        buffer.push_back('\0');
        xml_document<> document;
        document.parse<0>(&buffer[0]);

        xml_node<> *root_node = document.first_node("samples");
        string dir_path = get_dir(config_path) + "/" + "samples";
        auto [tiles,
                neighbors_ids,
                width, height,
                periodic_output, seed
        ] = read_simpletiled_instance(root_node->first_node("simpletiled"), dir_path);
        TilingWFC<Color> wfc(tiles, neighbors_ids, height, width, {periodic_output},
                             seed);
        return wfc.get_wave();
}

/**
 * Build wave from id
 */
std::optional<Wave> wave_from_id(const string &config_path, vector<pair<unsigned, unsigned>> &tiles_id) noexcept {
        ifstream config_file(config_path);
        // if file not exists
        if(!config_file) {
            cout << "File " << config_path << " does not exist" << endl;
            return {};
        }
        vector<char> buffer((istreambuf_iterator<char>(config_file)),
                            istreambuf_iterator<char>());
        buffer.push_back('\0');
        xml_document<> document;
        document.parse<0>(&buffer[0]);

        xml_node<> *root_node = document.first_node("samples");
        string dir_path = get_dir(config_path) + "/" + "samples";
        auto [tiles,
                neighbors_ids,
                width, height,
                periodic_output, seed
        ] = read_simpletiled_instance(root_node->first_node("simpletiled"), dir_path);
        TilingWFC<Color> wfc(tiles, neighbors_ids, height, width, {periodic_output},
                             seed);
        auto wave = wfc.get_wave();
        // set all parterns in wave false
        for (unsigned i = 0; i < wave.size; i++) {
            for(unsigned j=0;j<wave.nb_patterns;j++){
                wave.set(i, j, false);
            }
        }
        // set input parterns true
        for (int i = 0; i < wave.height; ++i) {
            for (int j = 0; j < wave.width; ++j) {
                wave.set(i,j,tiles_id[i*wave.width+j].first, true);
            }
        }
        return wave;
}


auto read_config(const string &config_path, bool out_img=false) noexcept {
    ifstream config_file(config_path);
    // if file not exists
//    if(!config_file) {
//        cout << "File " << config_path << " does not exist" << endl;
//        return {};
//    }
    vector<char> buffer((istreambuf_iterator<char>(config_file)),
                        istreambuf_iterator<char>());
    buffer.push_back('\0');
    xml_document<> document;
    document.parse<0>(&buffer[0]);

    xml_node<> *root_node = document.first_node("samples");
    string dir_path = get_dir(config_path) + "/" + "samples";
    auto [tiles,
            neighbors_ids,
            width, height,
            periodic_output, seed
    ] = read_simpletiled_instance(root_node->first_node("simpletiled"), dir_path);
    return make_tuple(tiles, neighbors_ids, height,width, TilingWFCOptions{periodic_output}, seed);
}


class XLandWFC {
private:
    const std::string config_path;
    tilingWFC::TilingWFC<Color> wfc;
    tilingWFC::TilingWFC<Color>* current_wfc;

public:
    XLandWFC(const std::string& config_path):
            config_path(config_path),
            wfc(std::get<0>(XMLWFC::read_config(config_path)),
                std::get<1>(XMLWFC::read_config(config_path)),
                std::get<2>(XMLWFC::read_config(config_path)),
                std::get<3>(XMLWFC::read_config(config_path)),
                std::get<4>(XMLWFC::read_config(config_path)),
                std::get<5>(XMLWFC::read_config(config_path)))
    {
    };
        auto build_a_open_area_wave(){
            auto base_wave = wfc.get_wave();
            auto out = new Wave(base_wave.height, base_wave.width, base_wave.patterns_frequencies);
            for (int i = 0; i < base_wave.size; i++) {
                for (int k = 0; k < base_wave.nb_patterns; k++) {
                    out->set(i, k, false);
                }
            }
            for(int i = 0; i < base_wave.size; i++) {
                out->set(i, 0, true);
            }
            return out;
        }

        static auto build_a_open_area_wave(const std::string _conifg_path){
            auto temp_wfc = new TilingWFC<Color> (std::get<0>(XMLWFC::read_config(_conifg_path)),
                                                    std::get<1>(XMLWFC::read_config(_conifg_path)),
                                                    std::get<2>(XMLWFC::read_config(_conifg_path)),
                                                    std::get<3>(XMLWFC::read_config(_conifg_path)),
                                                    std::get<4>(XMLWFC::read_config(_conifg_path)),
                                                    std::get<5>(XMLWFC::read_config(_conifg_path)));
            auto base_wave = temp_wfc->get_wave();
            auto out = new Wave(base_wave.height, base_wave.width, base_wave.patterns_frequencies);
            for (int i = 0; i < base_wave.size; i++) {
                for (int k = 0; k < base_wave.nb_patterns; k++) {
                    out->set(i, k, false);
                }
            }
            for(int i = 0; i < base_wave.size; i++) {
                out->set(i, 0, true);
            }
            delete temp_wfc;
            return out;
        }

        auto build_wave_from_ids(vector<pair<unsigned, unsigned>> &tiles_id){
             auto base_wave = wfc.get_wave();
            auto out = new Wave(base_wave.height, base_wave.width, base_wave.patterns_frequencies);
            // set all parterns in wave false
            for (unsigned i = 0; i < out->size; i++) {
                for(unsigned j=0;j<out->nb_patterns;j++){
                    out->set(i, j, false);
                }
            }
            auto ids = wfc.oriented_tile_ids_to_id(tiles_id);
            // set input parterns true
            for (int i = 0; i < out->size; ++i) {
                for(int k=0; k< out->nb_patterns; k++){
                    out->set(i, ids.data[i], true);
                }
            }
            return out;
        }

        static auto build_wave_from_ids(vector<pair<unsigned, unsigned>> &tiles_id, const std::string _conifg_path){
            // auto base_wave = wfc.get_wave();
            auto temp_wfc = new TilingWFC<Color> (std::get<0>(XMLWFC::read_config(_conifg_path)),
                                                  std::get<1>(XMLWFC::read_config(_conifg_path)),
                                                  std::get<2>(XMLWFC::read_config(_conifg_path)),
                                                  std::get<3>(XMLWFC::read_config(_conifg_path)),
                                                  std::get<4>(XMLWFC::read_config(_conifg_path)),
                                                  std::get<5>(XMLWFC::read_config(_conifg_path)));
            auto base_wave = temp_wfc->get_wave();
            auto out = new Wave(base_wave.height, base_wave.width, base_wave.patterns_frequencies);
            // set all parterns in wave false
            for (unsigned i = 0; i < out->size; i++) {
                for(unsigned j=0;j<out->nb_patterns;j++){
                    out->set(i, j, false);
                }
            }
            auto ids =  temp_wfc->oriented_tile_ids_to_id(tiles_id);
            // set input parterns true
            for (int i = 0; i < out->size; ++i) {
                for(int k=0; k< out->nb_patterns; k++){
                    out->set(i, ids.data[i], true);
                }
            }
            delete temp_wfc;
            return out;
        }

        auto get_ids_from_wave(Wave &wave){
            Array2D<unsigned> ids(wave.height, wave.width);
            for (int i = 0; i < wave.size; ++i) {
                for(int k=0; k< wave.nb_patterns; k++){
                    if(wave.get(i, k)){
                        ids.data[i] = k;
                    }
                }
            }
            auto result = wfc.id_to_oriented_tiles(ids);
            return result;
        }

        static auto get_ids_from_wave(Wave &wave, const std::string _conifg_path){
            Array2D<unsigned> ids(wave.height, wave.width);
            for (int i = 0; i < wave.size; ++i) {
                for(int k=0; k< wave.nb_patterns; k++){
                    if(wave.get(i, k)){
                        ids.data[i] = k;
                    }
                }
            }
            auto temp_wfc = new TilingWFC<Color> (std::get<0>(XMLWFC::read_config(_conifg_path)),
                                                  std::get<1>(XMLWFC::read_config(_conifg_path)),
                                                  std::get<2>(XMLWFC::read_config(_conifg_path)),
                                                  std::get<3>(XMLWFC::read_config(_conifg_path)),
                                                  std::get<4>(XMLWFC::read_config(_conifg_path)),
                                                  std::get<5>(XMLWFC::read_config(_conifg_path)));
            auto result = temp_wfc->id_to_oriented_tiles(ids);
            delete temp_wfc;
            return result;
        }

        static std::tuple<std::vector<std::pair<unsigned, unsigned>>,
                std::optional<Array2D<Color>>>
        run(const std::string _config_path, unsigned seed, bool out_img=false) {
            auto _current_wfc = new TilingWFC<Color>(std::get<0>(XMLWFC::read_config(_config_path)),
                                                     std::get<1>(XMLWFC::read_config(_config_path)),
                                                     std::get<2>(XMLWFC::read_config(_config_path)),
                                                     std::get<3>(XMLWFC::read_config(_config_path)),
                                                     std::get<4>(XMLWFC::read_config(_config_path)),
                                                     std::get<5>(XMLWFC::read_config(_config_path)));
            _current_wfc->set_seed(seed);
            for (unsigned test = 0; test < 30; test++) {
                auto success = _current_wfc->run();
                if (success.has_value()) {
                    auto output = success.value();
                    auto oriented_tiles = _current_wfc->id_to_oriented_tiles(output);
                    if (out_img) {
    //                 cout << " finished!" << endl;
                        auto img = _current_wfc->id_to_tiling(output);
                        delete _current_wfc;
                        return make_tuple(oriented_tiles, img);
                    } else {
                        delete _current_wfc;
    //                 cout << " finished!" << endl;
                        return make_tuple(oriented_tiles, std::nullopt);
                    }
                    break;
                } else {
                    if (test == 29) {
                        delete _current_wfc;
                        _current_wfc = new TilingWFC<Color>(std::get<0>(XMLWFC::read_config(_config_path)),
                                                            std::get<1>(XMLWFC::read_config(_config_path)),
                                                            std::get<2>(XMLWFC::read_config(_config_path)),
                                                            std::get<3>(XMLWFC::read_config(_config_path)),
                                                            std::get<4>(XMLWFC::read_config(_config_path)),
                                                            std::get<5>(XMLWFC::read_config(_config_path)));
                        _current_wfc->set_seed(seed);
                        test = 0;
                    }
                }
            }
    }

    std::tuple<std::vector<std::pair<unsigned, unsigned>>,
                std::optional<Array2D<Color>>>
        mutate(Wave base_wave, double new_weight, int iter_count=1, bool out_img= false) {
            for (unsigned i = 0; i < iter_count; i++) {
                current_wfc = new TilingWFC<Color> (std::get<0>(XMLWFC::read_config(config_path)),
                                                    std::get<1>(XMLWFC::read_config(config_path)),
                                                    std::get<2>(XMLWFC::read_config(config_path)),
                                                    std::get<3>(XMLWFC::read_config(config_path)),
                                                    std::get<4>(XMLWFC::read_config(config_path)),
                                                    std::get<5>(XMLWFC::read_config(config_path)));
                for (unsigned test = 0; test < 30; test++) {
                    auto success = current_wfc->mutate(base_wave, new_weight);
                    if (success.has_value()) {
                        auto output = success.value();
                        auto oriented_tiles = current_wfc->id_to_oriented_tiles(output);
                        if (out_img) {
//                 cout << " finished!" << endl;
                            auto img = current_wfc->id_to_tiling(output);
                            delete current_wfc;
                            return make_tuple(oriented_tiles, img);
                        } else {
//                 cout << " finished!" << endl;
                            delete current_wfc;
                            return make_tuple(oriented_tiles, std::nullopt);
                        }
                        break;
                    } else {
                    if(test==29){delete current_wfc; 
                        current_wfc = new TilingWFC<Color> (std::get<0>(XMLWFC::read_config(config_path)),
                                            std::get<1>(XMLWFC::read_config(config_path)),
                                            std::get<2>(XMLWFC::read_config(config_path)),
                                            std::get<3>(XMLWFC::read_config(config_path)),
                                            std::get<4>(XMLWFC::read_config(config_path)),
                                            std::get<5>(XMLWFC::read_config(config_path)));
                        test = 0;
                    }
                    }
                }
            }
        }

        static  std::tuple<std::vector<std::pair<unsigned, unsigned>>,
                std::optional<Array2D<Color>>>
        mutate(Wave base_wave, double new_weight, unsigned seed,
                      const std::string _config_path, int iter_count=1, bool out_img= false) {
            for (unsigned i = 0; i < iter_count; i++) {
                auto _current_wfc = new TilingWFC<Color> (std::get<0>(XMLWFC::read_config(_config_path)),
                                                    std::get<1>(XMLWFC::read_config(_config_path)),
                                                    std::get<2>(XMLWFC::read_config(_config_path)),
                                                    std::get<3>(XMLWFC::read_config(_config_path)),
                                                    std::get<4>(XMLWFC::read_config(_config_path)),
                                                    std::get<5>(XMLWFC::read_config(_config_path)));
                _current_wfc->set_seed(seed);
                for (unsigned test = 0; test < 30; test++) {
                    auto success = _current_wfc->mutate(base_wave, new_weight);
                    if (success.has_value()) {
                        auto output = success.value();
                        auto oriented_tiles = _current_wfc->id_to_oriented_tiles(output);
                        if (out_img) {
    //                 cout << " finished!" << endl;
                            auto img = _current_wfc->id_to_tiling(output);
                            delete _current_wfc;
                            return make_tuple(oriented_tiles, img);
                        } else {
    //                 cout << " finished!" << endl;
                            delete _current_wfc;
                            return make_tuple(oriented_tiles, std::nullopt);
                        }
                        break;
                    } else {
                        if(test==29){delete _current_wfc;
                            _current_wfc = new TilingWFC<Color> (std::get<0>(XMLWFC::read_config(_config_path)),
                                                                std::get<1>(XMLWFC::read_config(_config_path)),
                                                                std::get<2>(XMLWFC::read_config(_config_path)),
                                                                std::get<3>(XMLWFC::read_config(_config_path)),
                                                                std::get<4>(XMLWFC::read_config(_config_path)),
                                                                std::get<5>(XMLWFC::read_config(_config_path)));
                            _current_wfc->set_seed(seed);
                            test = 0;
                        }
                    }
                }
            }
        }

    std::tuple<std::vector<std::pair<unsigned, unsigned>>,
            std::optional<Array2D<Color>>>
    run(bool out_img=false) {
        current_wfc = new TilingWFC<Color>(std::get<0>(XMLWFC::read_config(config_path)),
                                           std::get<1>(XMLWFC::read_config(config_path)),
                                           std::get<2>(XMLWFC::read_config(config_path)),
                                           std::get<3>(XMLWFC::read_config(config_path)),
                                           std::get<4>(XMLWFC::read_config(config_path)),
                                           std::get<5>(XMLWFC::read_config(config_path)));
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        auto random_seed = (time_t)ts.tv_nsec;
        current_wfc->set_seed(random_seed);
        for (unsigned test = 0; test < 30; test++) {
            auto success = current_wfc->run();
            if (success.has_value()) {
                auto output = success.value();
                auto oriented_tiles = current_wfc->id_to_oriented_tiles(output);
                if (out_img) {
//                 cout << " finished!" << endl;
                    auto img = current_wfc->id_to_tiling(output);
                    delete current_wfc;
                    return make_tuple(oriented_tiles, img);
                } else {
                    delete current_wfc;
//                 cout << " finished!" << endl;
                    return make_tuple(oriented_tiles, std::nullopt);
                }
                break;
            } else {
                if (test == 29) {
                    delete current_wfc;
                    current_wfc = new TilingWFC<Color>(std::get<0>(XMLWFC::read_config(config_path)),
                                                       std::get<1>(XMLWFC::read_config(config_path)),
                                                       std::get<2>(XMLWFC::read_config(config_path)),
                                                       std::get<3>(XMLWFC::read_config(config_path)),
                                                       std::get<4>(XMLWFC::read_config(config_path)),
                                                       std::get<5>(XMLWFC::read_config(config_path)));
                    test = 0;
                }
            }
        }
    }
};

} // namespace xmlwfc
#endif //XMLWFC_XMLWFC_H