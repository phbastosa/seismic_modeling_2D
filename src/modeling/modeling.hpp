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
private:

    std::chrono::_V2::system_clock::time_point ti;

protected:

    std::string title;

    int receiver_output_samples;
    int wavefield_output_samples;

    bool export_receiver_output;
    bool export_wavefield_output;    

    std::string receiver_output_file;
    std::string wavefield_output_file;
    std::string receiver_output_folder;
    std::string wavefield_output_folder;

    float * receiver_output = nullptr;
    float * wavefield_output = nullptr;

    Model * model = nullptr;
    Wavelet * wavelet = nullptr;
    Geometry * geometry = nullptr;

    void set_geometry();
    void info_message();

    virtual void build_outputs() = 0;

public:

    int shot_id;
    int total_shots;
    int total_nodes;

    std::string file;

    void export_outputs();
    void get_execution_time();
    void show_execution_time();

    virtual void propagation() = 0;
    virtual void set_components() = 0;
    virtual void set_wavefields() = 0;

    virtual void set_parameters() = 0;
};

# endif