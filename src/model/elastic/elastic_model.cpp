# include <cmath> 

# include "elastic_model.hpp"

void Elastic_model::set_parameters(std::string file)
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
    std::string vs_file = catch_parameter("vs_file", file);
    std::string rho_file = catch_parameter("rho_file", file);

    B = new float[nPointsB]();
    M = new float[nPointsB]();
    L = new float[nPointsB]();

    vp = new float[nPoints]();
    vs = new float[nPoints]();
    rho = new float[nPoints]();

    read_binary_float(vp_file, vp, nPoints);
    read_binary_float(vs_file, vs, nPoints);
    read_binary_float(rho_file, rho, nPoints);
    
    expand_boundaries();

    delete[] vp;
    delete[] vs;
    delete[] rho;
}

void Elastic_model::expand_boundaries()
{
    for(int index = 0; index < nxx*nzz; index++) 
    {
        M[index] = rho[index]*pow(vs[index],2.0f);
        L[index] = rho[index]*pow(vp[index],2.0f) - 2.0f*M[index];
    }


    for (int i = 0; i < nz; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            B[(i + nb) + (j + nb)*nzz] = 1.0f / rho[i + j*nz];

            M[(i + nb) + (j + nb)*nzz] = rho[i + j*nz]*powf(vs[i + j*nz], 2.0f);         
            
            L[(i + nb) + (j + nb)*nzz] = rho[i + j*nz]*powf(vp[i + j*nz], 2.0f) - 2.0f*M[(i + nb) + (j + nb)*nzz];
        }
    }

    for (int i = nb; i < nzz-nb; i++)
    {
        for (int j = 0; j < nb; j++)
        {
            B[i + j*nzz] = 1.0f / rho[(i - nb) + 0*nz];
            B[i + (nxx - j - 1)*nzz] = 1.0f / rho[(i - nb) + (nx - 1)*nz];

            M[i + j*nzz] = rho[i + 0*nz]*powf(vs[i + 0*nz], 2.0f);
            M[i + (nxx - j - 1)*nzz] = rho[(i - nb) + (nx - 1)*nz]*powf(vs[(i - nb) + (nx - 1)*nz], 2.0f);

            L[i + j*nzz] = rho[i + 0*nz]*powf(vp[i + 0*nz], 2.0f) - 2.0f*M[i + j*nzz];
            L[i + (nxx - j - 1)*nzz] = rho[(i - nb) + (nx - 1)*nz]*powf(vp[(i - nb) + (nx - 1)*nz], 2.0f) - 2.0f*M[i + j*nzz];
        }
    }
    
    for (int i = 0; i < nb; i++)
    {
        for (int j = 0; j < nxx; j++)
        {
            B[i + j*nzz] = B[nb + j*nzz];
            B[(nzz - i - 1) + j*nzz] = B[(nzz - nb - 1) + j*nzz];

            M[i + j*nzz] = M[nb + j*nzz];
            M[(nzz - i - 1) + j*nzz] = M[(nzz - nb - 1) + j*nzz];

            L[i + j*nzz] = L[nb + j*nzz];
            L[(nzz - i - 1) + j*nzz] = L[(nzz - nb - 1) + j*nzz];
        }
    }
}
