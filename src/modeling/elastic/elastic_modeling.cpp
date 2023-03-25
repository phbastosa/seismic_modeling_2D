# include <cmath>
# include <iostream>

# include "elastic_modeling.hpp"

void Elastic_modeling::set_parameters()
{
    model = new Elastic_model();
    wavelet = new Gaussian_1st();

    model->set_parameters(file);
    wavelet->set_parameters(file);

    factor = std::stof(catch_parameter("damping_factor", file));

    export_receiver_output = str2bool(catch_parameter("export_seismogram", file));
    export_wavefield_output = str2bool(catch_parameter("export_snapshots", file));

    i_snap = std::stoi(catch_parameter("i_snap", file));
    f_snap = std::stoi(catch_parameter("f_snap", file));
    d_snap = std::stoi(catch_parameter("d_snap", file));

    receiver_output_folder = catch_parameter("seismogram_folder", file);
    wavefield_output_folder = catch_parameter("snapshots_folder", file);

    n_snap = (f_snap - i_snap) / d_snap + 1;

    total_times = wavelet->nt;
}

void Elastic_modeling::set_components()
{
    set_geometry();
    set_abc_dampers();

    Vx = new float[model->nPointsB]();
    Vz = new float[model->nPointsB]();
    Txx = new float[model->nPointsB]();
    Tzz = new float[model->nPointsB]();
    Txz = new float[model->nPointsB]();

    receiver_output = new float[total_times * total_nodes]();
    wavefield_output = new float[n_snap * model->nPoints]();
}

void Elastic_modeling::set_abc_dampers()
{
    damp1D = new float[model->nb]();
    damp2D = new float[model->nb * model->nb]();

    for (int i = 0; i < model->nb; i++) 
    {
        damp1D[i] = expf(-powf(factor * (model->nb - i), 2.0f));
    }

    for(int i = 0; i < model->nb; i++) 
    {
        for (int j = 0; j < model->nb; j++)
        {   
            damp2D[j + i*model->nb] += damp1D[i]; // up to bottom
            damp2D[i + j*model->nb] += damp1D[i]; // left to right
        }
    }

    for (int index = 0; index < model->nb * model->nb; index++)
        damp2D[index] -= 1.0f;
}

void Elastic_modeling::set_wavefields()
{
    for (int index = 0; index < model->nPointsB; index++)
    {
        Vx[index] = 0.0f;
        Vz[index] = 0.0f;
        Txx[index] = 0.0f;
        Tzz[index] = 0.0f;
        Txz[index] = 0.0f;
    }
}

void Elastic_modeling::info_message()
{
    if (time_id % (total_times / 100) == 0)
    {
        int result = system("clear");
        
        std::cout<<"2D wave propagation in elastic media\n\n";
        
        std::cout<<"Total x model length = "<<(model->nx-1)*model->dx<<" m\n";
        std::cout<<"Total Z model length = "<<(model->nz-1)*model->dz<<" m\n\n";
        
        std::cout<<"Shot "<<shot_id+1<<" of "<<total_shots<<"\n\n";

        std::cout<<"Position (z,x) = ("<<geometry->shots.z[shot_id]<<", "<<geometry->shots.x[shot_id]<<") m\n\n";

        std::cout<<"Modeling progression: "<< 100.0f *((float)(time_id+1) / (float)(total_times)) <<" %\n\n";    
    }
}

void Elastic_modeling::propagation()
{
    snap_id = 0;

    for (int time = 0; time < total_times; time++)
    {
        time_id = time;

        info_message();
        apply_wavelet();

        kernel_propagation();        
        
        build_outputs();
    }
}

void Elastic_modeling::apply_wavelet()
{
    int sid = geometry->shots.idz[shot_id] + geometry->shots.idx[shot_id] * model->nzz;    

    Txx[sid] += wavelet->amp[time_id] / (model->dx * model->dz);
    Tzz[sid] += wavelet->amp[time_id] / (model->dx * model->dz);
}

void Elastic_modeling::kernel_propagation()
{
    float damper, Mxz, Bx, Bz; 
    float dVx_dx, dVz_dz, dVz_dx, dVx_dz;    
    float dTxx_dx, dTxz_dx, dTzz_dz, dTxz_dz;

    # pragma omp parallel for
    for (int index = 0; index < model->nPointsB; index++)
    {
        int i = (int)(index % model->nzz);
        int j = (int)(index / model->nzz);

        if ((i >= 3) && (i < model->nzz - 4) && (j >= 4) && (j < model->nxx - 3)) 
        {
            dTxx_dx = (75.0f*(Txx[i + (j-4)*model->nzz] - Txx[i + (j+3)*model->nzz]) +
                     1029.0f*(Txx[i + (j+2)*model->nzz] - Txx[i + (j-3)*model->nzz]) +
                     8575.0f*(Txx[i + (j-2)*model->nzz] - Txx[i + (j+1)*model->nzz]) +
                   128625.0f*(Txx[i + j*model->nzz]     - Txx[i + (j-1)*model->nzz])) / (model->dx*107520.0f);

            dTxz_dz = (75.0f*(Txz[(i-3) + j*model->nzz] - Txz[(i+4) + j*model->nzz]) +
                     1029.0f*(Txz[(i+3) + j*model->nzz] - Txz[(i-2) + j*model->nzz]) + 
                     8575.0f*(Txz[(i-1) + j*model->nzz] - Txz[(i+2) + j*model->nzz]) +
                   128625.0f*(Txz[(i+1) + j*model->nzz] - Txz[i + j*model->nzz])) / (model->dz*107520.0f);

            Bx = 0.5f*(model->B[i + (j+1)*model->nzz] + model->B[i + j*model->nzz]);

            Vx[index] += wavelet->dt*Bx*(dTxx_dx + dTxz_dz);  
        }
      
        if((i >= 4) && (i < model->nzz - 3) && (j >= 3) && (j < model->nxx - 4)) 
        {
            dTxz_dx = (75.0f*(Txz[i + (j-3)*model->nzz] - Txz[i + (j+4)*model->nzz]) +
                     1029.0f*(Txz[i + (j+3)*model->nzz] - Txz[i + (j-2)*model->nzz]) +
                     8575.0f*(Txz[i + (j-1)*model->nzz] - Txz[i + (j+2)*model->nzz]) +
                   128625.0f*(Txz[i + (j+1)*model->nzz] - Txz[i + j*model->nzz])) / (model->dx*107520.0f);

            dTzz_dz = (75.0f*(Tzz[(i-4) + j*model->nzz] - Tzz[(i+3) + j*model->nzz]) + 
                     1029.0f*(Tzz[(i+2) + j*model->nzz] - Tzz[(i-3) + j*model->nzz]) +
                     8575.0f*(Tzz[(i-2) + j*model->nzz] - Tzz[(i+1) + j*model->nzz]) +
                   128625.0f*(Tzz[i + j*model->nzz]     - Tzz[(i-1) + j*model->nzz])) / (model->dz*107520.0f);
        
            Bz = 0.5f*(model->B[(i+1) + j*model->nzz] + model->B[i + j*model->nzz]);

            Vz[index] += wavelet->dt*Bz*(dTxz_dx + dTzz_dz); 
        }
    }

    # pragma omp parallel for
    for (int index = 0; index < model->nPointsB; index++)
    {
        int i = (int)(index % model->nzz);
        int j = (int)(index / model->nzz);

        if ((i >= 3) && (i < model->nzz - 4) && (j >= 3) && (j < model->nxx - 4)) 
        {           
            dVx_dx = (75.0f*(Vx[i + (j-3)*model->nzz] - Vx[i + (j+4)*model->nzz]) + 
                    1029.0f*(Vx[i + (j+3)*model->nzz] - Vx[i + (j-2)*model->nzz]) +
                    8575.0f*(Vx[i + (j-1)*model->nzz] - Vx[i + (j+2)*model->nzz]) + 
                  128625.0f*(Vx[i + (j+1)*model->nzz] - Vx[i + j*model->nzz])) / (model->dx * 107520.0f);

            dVz_dz = (75.0f*(Vz[(i-3) + j*model->nzz] - Vz[(i+4) + j*model->nzz]) +   
                    1029.0f*(Vz[(i+3) + j*model->nzz] - Vz[(i-2) + j*model->nzz]) +
                    8575.0f*(Vz[(i-1) + j*model->nzz] - Vz[(i+2) + j*model->nzz]) +
                  128625.0f*(Vz[(i+1) + j*model->nzz] - Vz[i + j*model->nzz])) / (model->dz * 107520.0f);     

            Txx[index] += wavelet->dt*((model->L[index] + 2.0f*model->M[index])*dVx_dx + model->L[index]*dVz_dz);   
            Tzz[index] += wavelet->dt*((model->L[index] + 2.0f*model->M[index])*dVz_dz + model->L[index]*dVx_dx);
        }
    
        if ((i >= 4) && (i < model->nzz - 3) && (j >= 4) && (j < model->nxx - 3)) 
        {
            dVz_dx = (75.0f*(Vz[i + (j-4)*model->nzz] - Vz[i + (j+3)*model->nzz]) +
                    1029.0f*(Vz[i + (j+2)*model->nzz] - Vz[i + (j-3)*model->nzz]) +
                    8575.0f*(Vz[i + (j-2)*model->nzz] - Vz[i + (j+1)*model->nzz]) +
                  128625.0f*(Vz[i + j*model->nzz]     - Vz[i + (j-1)*model->nzz])) / (model->dx*107520.0f);

            dVx_dz = (75.0f*(Vx[(i-4) + j*model->nzz] - Vx[(i+3) + j*model->nzz]) +
                    1029.0f*(Vx[(i+2) + j*model->nzz] - Vx[(i-3) + j*model->nzz]) +
                    8575.0f*(Vx[(i-2) + j*model->nzz] - Vx[(i+1) + j*model->nzz]) +
                  128625.0f*(Vx[i + j*model->nzz]     - Vx[(i-1) + j*model->nzz])) / (model->dz*107520.0f);

            Mxz = powf(0.25f*(1.0f / model->M[(i+1) + j*model->nzz] + 
                              1.0f / model->M[i + (j+1)*model->nzz] + 
                              1.0f / model->M[(i+1) + (j+1)*model->nzz] + 
                              1.0f / model->M[i + j*model->nzz]), -1.0f); 

            Txz[index] += wavelet->dt*Mxz*(dVx_dz + dVz_dx);            
        }          

        if ((i >= model->nb) && (i < model->nzz - model->nb) && (j >= model->nb) && (j < model->nxx - model->nb))
        {
            damper = 1.0f;
        }
        else
        {
            if ((i < model->nb) && (j >= model->nb) && (j < model->nxx - model->nb)) 
            {
                damper = damp1D[i];
            }         
            else if ((i >= model->nzz - model->nb) && (j >= model->nb) && (j < model->nxx - model->nb)) 
            {
                damper = damp1D[model->nb - (i - (model->nzz - model->nb)) - 1];
            }         
            else if ((i >= model->nb) && (i < model->nzz - model->nb) && (j < model->nb))
            {
                damper = damp1D[j];
            }
            else if ((i >= model->nb) && (i < model->nzz - model->nb) && (j >= model->nxx - model->nb))
            {
                damper = damp1D[model->nb - (j - (model->nxx - model->nb)) - 1];
            }
            else if ((i < model->nb) && (j < model->nb))
            {
                damper = damp2D[i + j*model->nb];
            }
            else if ((i < model->nb) && (j >= model->nxx - model->nb))
            {
                damper = damp2D[i + (model->nb - (j - (model->nxx - model->nb)) - 1)*model->nb];
            }            
            else if ((i >= model->nzz - model->nb) && (j < model->nb))
            {
                damper = damp2D[(model->nb - (i - (model->nzz - model->nb)) - 1) + j*model->nb];
            }    
            else if ((i >= model->nzz - model->nb) && (j >= model->nxx - model->nb))
            {
                damper = damp2D[(model->nb - (i - (model->nzz - model->nb)) - 1) + (model->nb - (j - (model->nxx - model->nb)) - 1)*model->nb];
            }
        }

        Vx[index] *= damper;
        Vz[index] *= damper;
        Txx[index] *= damper;
        Tzz[index] *= damper;
        Txz[index] *= damper;
    }
}

void Elastic_modeling::build_outputs()
{
    if (export_receiver_output)
    {
        int current_node = 0;
        int iNode = geometry->iRel[shot_id];
        int fNode = geometry->fRel[shot_id];
        
        for (int node = iNode; node < fNode; node++)
        {    
            int index = geometry->nodes.idz[node] + geometry->nodes.idx[node] * model->nzz;
                
            receiver_output[time_id +  current_node * total_times] = 0.5f*(Txx[index] + Tzz[index]);

            current_node += 1;
        }
    }
        
    if (export_wavefield_output)
    {
        if (time_id % d_snap == 0)
        {
            if ((time_id >= i_snap) && (time_id <= f_snap))
            {            
                # pragma omp parallel for
                for (int index = 0; index < model->nPoints; index++)
                {
                    int i = (int)(index % model->nz);
                    int j = (int)(index / model->nz);
                
                    wavefield_output[index + snap_id * model->nPoints] = 0.5f*(Txx[(i + model->nb) + (j + model->nb) * model->nzz] + Tzz[(i + model->nb) + (j + model->nb) * model->nzz]);
                }

                snap_id += 1;
            }
        }
    }
}

void Elastic_modeling::export_outputs()
{
    std::string receiver_output_name = receiver_output_folder + "seismogram_elastic_" + std::to_string(total_times) + "x" + std::to_string(total_nodes) + "_shot_" + std::to_string(shot_id+1) + ".bin";
    std::string wavefield_output_name = wavefield_output_folder + "snapshots_elastic_" + std::to_string(n_snap) + "x" + std::to_string(model->nz) + "x" + std::to_string(model->nx) + "_shot_" + std::to_string(shot_id+1) + ".bin";
    
    if (export_receiver_output) 
        write_binary_float(receiver_output_name, receiver_output, total_times * total_nodes);
    
    if (export_wavefield_output) 
        write_binary_float(wavefield_output_name, wavefield_output, n_snap * model->nPoints);
}


