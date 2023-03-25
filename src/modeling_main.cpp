# include "modeling/eikonal/eikonal_modeling.hpp"
# include "modeling/scalar/scalar_modeling.hpp"
# include "modeling/acoustic/acoustic_modeling.hpp"
# include "modeling/elastic/elastic_modeling.hpp"

int main(int argc, char **argv)
{        
    Modeling * modeling[] = 
    {
        new Eikonal_modeling(),
        new Scalar_modeling(),
        new Acoustic_modeling(),
        new Elastic_modeling()    
    };
    
    int type = std::stoi(catch_parameter("modeling_type", std::string(argv[1])));

    modeling[type]->file = std::string(argv[1]);

    modeling[type]->set_parameters();
    modeling[type]->set_components();

    modeling[type]->get_execution_time();

    for (int shot = 0; shot < modeling[type]->total_shots; shot++)
    {
        modeling[type]->shot_id = shot;

        modeling[type]->set_wavefields();

        modeling[type]->propagation();
        
        modeling[type]->export_outputs();
    }

    modeling[type]->show_execution_time();

    return 0;
}