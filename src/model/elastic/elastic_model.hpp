# ifndef ELASTIC_MODEL_HPP
# define ELASTIC_MODEL_HPP

# include "../model.hpp"

class Elastic_model : public Model
{
private:

    float * vp = nullptr;
    float * vs = nullptr;
    float * rho = nullptr;

    void expand_boundaries();

public:

    void set_parameters(std::string file);
};

# endif
