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
protected:

    void import_geometry();
    void export_geometry();

    std::vector<float> linspace(float xi, float xf, int n);

public:

    int * iRel;
    int * fRel;

    Position shots;
    Position nodes;

    bool reciprocity;

    virtual void set_parameters(std::string file) = 0; 
};

# endif