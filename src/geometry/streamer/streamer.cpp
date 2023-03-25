# include <iostream>

# include "streamer.hpp"

void Streamer::set_parameters(std::string file)
{
    set_shot_geometry(file);
    set_node_geometry(file);
    set_relation();
}

void Streamer::set_shot_geometry(std::string file)
{
    float selev = std::stof(catch_parameter("s_elev", file));
    shots.total = std::stoi(catch_parameter("n_shot", file));

    float i_shot = std::stof(catch_parameter("i_shot", file));
    float f_shot = std::stof(catch_parameter("f_shot", file));

    auto tmp_x = linspace(i_shot, f_shot, shots.total);

    shots.x = new float[shots.total]();
    shots.z = new float[shots.total]();

    for (int i = 0; i < shots.total; i++)
    {
        shots.x[i] = tmp_x[i];
        shots.z[i] = selev;    
    }

    std::vector<float>().swap(tmp_x);
}

void Streamer::set_node_geometry(std::string file)
{
    int spread = std::stoi(catch_parameter("spread", file));

    float gelev = std::stof(catch_parameter("g_elev", file));
    nodes.total = spread * shots.total;

    float offset_min = std::stof(catch_parameter("offset_min", file));
    float offset_max = std::stof(catch_parameter("offset_max", file));

    nodes.x = new float[nodes.total]();
    nodes.z = new float[nodes.total]();

    auto tmp_x = linspace(offset_min, offset_max, spread);

    for (int j = 0; j < shots.total; j++)
    {    
        for (int i = 0; i < spread; i++)
        {
            nodes.x[i + j*spread] = shots.x[j] + tmp_x[i];
            nodes.z[i + j*spread] = gelev;
        }
    }
}

void Streamer::set_relation()
{
    iRel = new int[shots.total]();
    fRel = new int[shots.total]();    

    int spread = nodes.total / shots.total;

    for (int i = 0; i < shots.total; i++)
    {
        iRel[i] = i*spread;
        fRel[i] = i*spread + spread;
    }
}
