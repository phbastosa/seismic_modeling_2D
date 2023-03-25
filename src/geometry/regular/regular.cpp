# include <iostream>

# include "regular.hpp"

void Regular::set_parameters(std::string file)
{
    shots.total = std::stoi(catch_parameter("total_shots", file));
    nodes.total = std::stoi(catch_parameter("total_nodes", file));

    reciprocity = str2bool(catch_parameter("reciprocity", file));

    set_shot_geometry(file);
    set_node_geometry(file);

    if (reciprocity)
    {
        std::swap(shots.x, nodes.x);
        std::swap(shots.z, nodes.z);
        std::swap(shots.total, nodes.total);
    }

    set_relation();
}

void Regular::set_shot_geometry(std::string file)
{
    shots.x = new float[shots.total]();
    shots.z = new float[shots.total]();

    std::vector<std::string> shot_beg_position, shot_end_position;

    shot_beg_position = split(catch_parameter("shot_beg_position", file),',');
    shot_end_position = split(catch_parameter("shot_end_position", file),',');

    float z_shot_beg = std::stof(shot_beg_position[0]); 
    float z_shot_end = std::stof(shot_end_position[0]); 
    float x_shot_beg = std::stof(shot_beg_position[1]); 
    float x_shot_end = std::stof(shot_end_position[1]); 
    
    if (z_shot_beg == z_shot_end)
    { 
        if (x_shot_beg <= x_shot_end)
        {
            auto tmp_x = linspace(x_shot_beg, x_shot_end, shots.total);       

            for (int i = 0; i < shots.total; i++)
            {
                shots.x[i] = tmp_x[i];
                shots.z[i] = z_shot_beg;
            }

            std::vector<float>().swap(tmp_x);
        }
        else
        {
            throw std::invalid_argument("Geometry error: Parameter \033[31mshot_beg_x\033[0;0m must be greather than \033[31mshot_end_x\033[0;0m in horizontal geometry!");
        }
    }

    if (x_shot_beg == x_shot_end)
    {
        if (z_shot_beg <= z_shot_end)
        {
            auto tmp_z = linspace(z_shot_beg, z_shot_end, shots.total);       

            for (int i = 0; i < shots.total; i++)
            {
                shots.x[i] = x_shot_beg;
                shots.z[i] = tmp_z[i];
            }

            std::vector<float>().swap(tmp_z);
        }
        else
        {
            throw std::invalid_argument("Geometry error: Parameter \033[31mshot_beg_z\033[0;0m must be greather than \033[31mshot_end_z\033[0;0m in vertical geometry!");            
        }
    }
}

void Regular::set_node_geometry(std::string file)
{
    nodes.x = new float[nodes.total]();
    nodes.z = new float[nodes.total]();

    std::vector<std::string> node_beg_position, node_end_position;

    node_beg_position = split(catch_parameter("node_beg_position", file),',');
    node_end_position = split(catch_parameter("node_end_position", file),',');

    float z_node_beg = std::stof(node_beg_position[0]); 
    float z_node_end = std::stof(node_end_position[0]); 
    float x_node_beg = std::stof(node_beg_position[1]); 
    float x_node_end = std::stof(node_end_position[1]); 
    
    if (z_node_beg == z_node_end)
    { 
        if (x_node_beg <= x_node_end)
        {
            auto tmp_x = linspace(x_node_beg, x_node_end, nodes.total);       

            for (int i = 0; i < nodes.total; i++)
            {
                nodes.x[i] = tmp_x[i];
                nodes.z[i] = z_node_beg;
            }

            std::vector<float>().swap(tmp_x);
        }
        else
        {
            throw std::invalid_argument("Geometry error: Parameter \033[31mnode_beg_x\033[0;0m must be greather than \033[31mnode_end_x\033[0;0m in horizontal geometry!");
        }
    }

    if (x_node_beg == x_node_end)
    {
        if (z_node_beg <= z_node_end)
        {
            auto tmp_z = linspace(z_node_beg, z_node_end, nodes.total);       

            for (int i = 0; i < nodes.total; i++)
            {
                nodes.x[i] = x_node_beg;
                nodes.z[i] = tmp_z[i];
            }

            std::vector<float>().swap(tmp_z);
        }
        else
        {
            throw std::invalid_argument("Geometry error: Parameter \033[31mnode_beg_z\033[0;0m must be greather than \033[31mnode_end_z\033[0;0m in vertical geometry!");            
        }
    }
}

void Regular::set_relation()
{
    iRel = new int[shots.total]();
    fRel = new int[shots.total]();    

    for (int i = 0; i < shots.total; i++)
    {
        iRel[i] = 0;
        fRel[i] = nodes.total;
    }
}
