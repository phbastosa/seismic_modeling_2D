# ifndef SCALAR_MODEL_HPP
# define SCALAR_MODEL_HPP

# include "../model.hpp"

class Scalar_model : public Model
{
private:

    float * vp = nullptr;

    void expand_boundaries();

public:

    void set_parameters(std::string file);
};

# endif
