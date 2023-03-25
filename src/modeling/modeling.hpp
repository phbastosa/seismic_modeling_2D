# ifndef MODELING_HPP
# define MODELING_HPP

# include <chrono>
# include <iostream>

# include "../utils/io/io.hpp"
# include "../model/model.hpp"
# include "../wavelet/wavelet.hpp"
# include "../geometry/geometry.hpp"

# include "../geometry/regular/regular.hpp"
# include "../geometry/streamer/streamer.hpp"

class Modeling
{
protected:

    std::chrono::_V2::system_clock::time_point ti;

    bool export_receiver_output;
    bool export_wavefield_output;

    float * receiver_output = nullptr;
    float * wavefield_output = nullptr;

    std::string receiver_output_folder;
    std::string wavefield_output_folder;

    Model * model = nullptr;
    Wavelet * wavelet = nullptr;
    Geometry * geometry = nullptr;

    void set_geometry();

    virtual void info_message() = 0;
    virtual void build_outputs() = 0;

public:

    int shot_id;
    int total_shots;
    int total_nodes;

    std::string file;

    void get_execution_time();
    void show_execution_time();

    virtual void propagation() = 0;
    virtual void set_components() = 0;
    virtual void set_wavefields() = 0;
    virtual void export_outputs() = 0;

    virtual void set_parameters() = 0;
};

# endif