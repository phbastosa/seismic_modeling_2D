# ifndef STREAMER_HPP
# define STREAMER_HPP

# include "../geometry.hpp"

class Streamer : public Geometry
{
private:

    void set_relation();
    void set_shot_geometry(std::string file);
    void set_node_geometry(std::string file);

public:

    void set_parameters(std::string file); 
};

# endif