# ifndef MODEL_HPP
# define MODEL_HPP

# include "../utils/io/io.hpp"

class Model
{
protected:

    virtual void expand_boundaries() = 0;

public:

    int nb;
    int nx, nz;
    int nxx, nzz;
    float dx, dz;

    int nPoints;
    int nPointsB;

    float * V = nullptr;
    float * B = nullptr;
    float * K = nullptr;
    float * M = nullptr;
    float * L = nullptr;

    virtual void set_parameters(std::string file) = 0;
};

# endif