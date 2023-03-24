# ifndef ACOUSTIC_MODEL_HPP
# define ACOUSTIC_MODEL_HPP

# include "../model.hpp"

class Acoustic_model : public Model
{
private:

    float * vp = nullptr;
    float * rho = nullptr;

    void expand_boundaries();

public:

    void set_parameters(std::string file);
};

# endif
