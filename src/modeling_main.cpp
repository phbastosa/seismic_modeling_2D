# include <chrono>
# include <iostream>

# include "modeling/eikonal/eikonal_modeling.hpp"
# include "modeling/scalar/scalar_modeling.hpp"
# include "modeling/acoustic/acoustic_modeling.hpp"
# include "modeling/elastic/elastic_modeling.hpp"

int main(int argc, char **argv)
{    
    auto ti = std::chrono::system_clock::now();
    
    Modeling * modeling[] = 
    {
        new Eikonal_modeling(),
        new Scalar_modeling(),
        new Acoustic_modeling(),
        new Elastic_modeling()    
    };
    
    std::string file = std::string(argv[1]);

    int type = std::stoi(catch_parameter("modeling_type", file));

    modeling[type]->set_parameters(file);

    modeling[type]->set_components();

    for (int shot = 0; shot < modeling[type]->total_shots; shot++)
    {
        modeling[type]->shot_id = shot;

        modeling[type]->set_wavefields();

        modeling[type]->propagation();
        
        modeling[type]->export_outputs();
    }

    auto tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::cout<<"\nModeling run time: "<<elapsed_seconds.count()<<" s."<<std::endl;
    
    return 0;
}