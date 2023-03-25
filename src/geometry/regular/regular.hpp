# ifndef REGULAR_HPP
# define REGULAR_HPP

# include "../geometry.hpp"

class Regular : public Geometry
{
private:

    void set_relation();
    void set_shot_geometry(std::string file);
    void set_node_geometry(std::string file);

public:

    void set_parameters(std::string file); 
};

# endif