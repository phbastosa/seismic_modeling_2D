# ifndef EIKONAL_MODEL_HPP
# define EIKONAL_MODEL_HPP

# include "../model.hpp"

class Eikonal_model : public Model
{
private:

    float * vp = nullptr;

    void expand_boundaries();

public:

    void set_parameters(std::string file);
};

# endif