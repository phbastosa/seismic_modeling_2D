# include <cmath>

# include "acoustic_model.hpp"

void Acoustic_model::set_parameters(std::string file)
{
    nx = std::stoi(catch_parameter("nx", file));    
    nz = std::stoi(catch_parameter("nz", file));

    nb = std::stoi(catch_parameter("nb", file));

    dx = std::stof(catch_parameter("dx", file));
    dz = std::stof(catch_parameter("dz", file));

    nxx = nx + 2*nb;
    nzz = nz + 2*nb;

    nPoints = nx * nz;
    nPointsB = nxx * nzz;

    std::string vp_file = catch_parameter("vp_file", file);
    std::string rho_file = catch_parameter("rho_file", file);

    K = new float[nPointsB]();
    B = new float[nPointsB]();

    vp = new float[nPoints]();
    rho = new float[nPoints]();

    read_binary_float(vp_file, vp, nPoints);
    read_binary_float(rho_file, rho, nPoints);
    
    expand_boundaries();

    delete[] vp;
    delete[] rho;
}

void Acoustic_model::expand_boundaries()
{
    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            K[(i + nb) + (j + nb)*nzz] = rho[i + j*nz] * powf(vp[i + j*nz], 2.0f);
            B[(i + nb) + (j + nb)*nzz] = 1.0f / rho[i + j*nz];
        }
    }

    for (int i = nb; i < nzz-nb; i++)
    {
        for (int j = 0; j < nb; j++)
        {
            K[i + j*nzz] = rho[(i - nb) + 0*nz] * powf(vp[(i - nb) + 0*nz], 2.0f);
            K[i + (nxx - j - 1)*nzz] = rho[(i - nb) + (nx - 1)*nz] * powf(vp[(i - nb) + (nx - 1)*nz], 2.0f);

            B[i + j*nzz] = 1.0f / rho[(i - nb) + 0*nz];
            B[i + (nxx - j - 1)*nzz] = 1.0f / rho[(i - nb) + (nx - 1)*nz];
        }
    }
    
    for (int i = 0; i < nb; i++)
    {
        for (int j = 0; j < nxx; j++)
        {
            K[i + j*nzz] = K[nb + j*nzz];
            K[(nzz - i - 1) + j*nzz] = K[(nzz - nb - 1) + j*nzz];

            B[i + j*nzz] = B[nb + j*nzz];
            B[(nzz - i - 1) + j*nzz] = B[(nzz - nb - 1) + j*nzz];
        }
    }
}
