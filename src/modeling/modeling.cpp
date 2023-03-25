# include "modeling.hpp"

void Modeling::get_execution_time()
{
    ti = std::chrono::system_clock::now();
}

void Modeling::show_execution_time()
{
    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::cout<<"\nRun time: "<<elapsed_seconds.count()<<" s."<<std::endl;
}

void Modeling::set_geometry()
{
    Geometry * gtypes[] = 
    {
        new Regular(),
        new Streamer()
    };

    int geometry_type = std::stoi(catch_parameter("geometry_type", file));
    
    geometry = gtypes[geometry_type];
    
    geometry->set_parameters(file);

    total_nodes = geometry->nodes.total;
    total_shots = geometry->shots.total;

    geometry->shots.idx = new int[total_shots]();
    geometry->shots.idz = new int[total_shots]();

    for (int i = 0; i < total_shots; i++)
    {
        geometry->shots.idx[i] = (int)(geometry->shots.x[i] / model->dx) + model->nb;
        geometry->shots.idz[i] = (int)(geometry->shots.z[i] / model->dz) + model->nb;
    }

    geometry->nodes.idx = new int[total_nodes]();
    geometry->nodes.idz = new int[total_nodes]();

    for (int i = 0; i < total_nodes; i++)
    {
        geometry->nodes.idx[i] = (int)(geometry->nodes.x[i] / model->dx) + model->nb;
        geometry->nodes.idz[i] = (int)(geometry->nodes.z[i] / model->dz) + model->nb;
    }
}