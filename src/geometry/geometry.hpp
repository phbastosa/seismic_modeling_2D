# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../utils/io/io.hpp"

class Position
{
public:

    int total;

    float * x;
    float * z;

    int * idx;
    int * idz;
};

class Geometry
{
public:

    int * iRel;
    int * fRel;

    Position shots;
    Position nodes;

    bool reciprocity;

    virtual void set_parameters(std::string file) = 0; 
};

# endif