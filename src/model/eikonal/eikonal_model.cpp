# include "eikonal_model.hpp"

void Eikonal_model::set_parameters(std::string file)
{
    nx = std::stoi(catch_parameter("nx", file));    
    nz = std::stoi(catch_parameter("nz", file));

    nb = 1;

    dx = std::stof(catch_parameter("dx", file));
    dz = std::stof(catch_parameter("dz", file));

    nxx = nx + 2*nb;
    nzz = nz + 2*nb;

    nPoints = nx * nz;
    nPointsB = nxx * nzz;

    std::string vp_file = catch_parameter("vp_file", file);

    V = new float[nPointsB]();
    vp = new float[nPoints]();

    read_binary_float(vp_file, vp, nPoints);
    
    expand_boundaries();

    delete[] vp;
}

void Eikonal_model::expand_boundaries()
{
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            V[(i + nb) + (j + nb)*nzz] = 1.0f / vp[i + j*nz];
        }
    }

    for (int i = nb; i < nzz-nb; i++)
    {
        for (int j = 0; j < nb; j++)
        {
            V[i + j*nzz] = 1.0f / vp[(i - nb) + 0*nz];
            V[i + (nxx - j - 1)*nzz] = 1.0f / vp[(i - nb) + (nx - 1)*nz];
        }
    }
    
    for (int i = 0; i < nb; i++)
    {
        for (int j = 0; j < nxx; j++)
        {
            V[i + j*nzz] = V[nb + j*nzz];
            V[(nzz - i - 1) + j*nzz] = V[(nzz - nb - 1) + j*nzz];
        }
    }
}
