# include <cmath>
# include <iostream>

# include "eikonal_modeling.hpp"

void Eikonal_modeling::set_parameters()
{
    model = new Eikonal_model();

    model->set_parameters(file);

    export_receiver_output = str2bool(catch_parameter("export_first_arrivals", file));
    export_wavefield_output = str2bool(catch_parameter("export_travel_times", file));

    receiver_output_folder = catch_parameter("first_arrivals_folder", file);
    wavefield_output_folder = catch_parameter("travel_times_folder", file);
}

void Eikonal_modeling::set_components()
{
    set_geometry();

    T = new float[model->nPointsB]();

    receiver_output = new float[total_nodes]();
    wavefield_output = new float[model->nPoints]();
}

void Eikonal_modeling::set_wavefields()
{
    float sid = geometry->shots.idz[shot_id] + geometry->shots.idx[shot_id] * model->nzz;
    
    for (int index = 0; index < model->nPointsB; index++)
    {
        T[index] = 1e6f;

        if (index == sid) T[index] = 0.0f;      
    }
}

void Eikonal_modeling::info_message()
{
    int result = system("clear");
    
    std::cout<<"2D eikonal equation solver in scalar acoustic media\n\n";
    
    std::cout<<"Total x model length = "<<(model->nx-1)*model->dx<<" m\n";
    std::cout<<"Total Z model length = "<<(model->nz-1)*model->dz<<" m\n\n";
    
    std::cout<<"Shot "<<shot_id+1<<" of "<<total_shots<<"\n\n";

    std::cout<<"Position (z,x) = ("<<geometry->shots.z[shot_id]<<", "<<geometry->shots.x[shot_id]<<") m\n\n";
}

void Eikonal_modeling::propagation()
{
    info_message();

    kernel_propagation();

    build_outputs();
}

void Eikonal_modeling::kernel_propagation()
{
    float sx = geometry->shots.x[shot_id];
    float sz = geometry->shots.z[shot_id];

    int sidz = geometry->shots.idz[shot_id];
    int sidx = geometry->shots.idx[shot_id];

    int id = sidz + sidx * model->nzz;

    T[id + 1] = model->V[id] * sqrtf(powf((float)(sidz + 1)*model->dz - sz, 2.0f) + powf((float)(sidx)*model->dx - sx, 2.0f));
    T[id - 1] = model->V[id] * sqrtf(powf((float)(sidz - 1)*model->dz - sz, 2.0f) + powf((float)(sidx)*model->dx - sx, 2.0f));
    T[id + model->nzz] = model->V[id] * sqrtf(powf((float)(sidz)*model->dz - sz, 2.0f) + powf((float)(sidx + 1)*model->dx - sx, 2.0f));
    T[id - model->nzz] = model->V[id] * sqrtf(powf((float)(sidz)*model->dz - sz, 2.0f) + powf((float)(sidx - 1)*model->dx - sx, 2.0f));
    T[id + 1 + model->nzz] = model->V[id] * sqrtf(powf((float)(sidz + 1)*model->dz - sz, 2.0f) + powf((float)(sidx + 1)*model->dx - sx, 2.0f));
    T[id + 1 - model->nzz] = model->V[id] * sqrtf(powf((float)(sidz + 1)*model->dz - sz, 2.0f) + powf((float)(sidx - 1)*model->dx - sx, 2.0f));
    T[id - 1 + model->nzz] = model->V[id] * sqrtf(powf((float)(sidz - 1)*model->dz - sz, 2.0f) + powf((float)(sidx + 1)*model->dx - sx, 2.0f));
    T[id - 1 - model->nzz] = model->V[id] * sqrtf(powf((float)(sidz - 1)*model->dz - sz, 2.0f) + powf((float)(sidx - 1)*model->dx - sx, 2.0f));

    dz2i = 1.0 / (model->dz * model->dz);
    dx2i = 1.0 / (model->dx * model->dx);

    sgntz = 1; sgnvz = 1;
    for (i = 1; i < model->nzz; i++)
    {
        sgntx = 1; sgnvx = 1;
        for (j = 1; j < model->nxx; j++)
            inner_sweep();
        
        sgntx = -1; sgnvx = 0;
        for (j = model->nxx-2; j >= 0; j--)
            inner_sweep();
    }

    sgntz = -1; sgnvz = 0;
    for (i = model->nzz-2; i >= 0; i--)
    {
        sgntx = 1; sgnvx = 1;
        for (j = 1; j < model->nxx; j++)
            inner_sweep();

        sgntx = -1; sgnvx = 0;
        for (j = model->nxx-2; j >= 0; j--)
            inner_sweep();
    }        

    sgntx = 1; sgnvx = 1;
    for (j = 1; j < model->nxx; j++) 
    {
        sgntz = 1; sgnvz = 1;
        for (i = 1; i < model->nzz; i++)
            inner_sweep();

        sgntz = -1; sgnvz = 0;
        for (i = model->nzz-2; i >= 0; i--)
            inner_sweep();
    }
        
    sgntx = -1; sgnvx = 0;
    for (j = model->nxx-2; j >= 0; j--)
    {
        sgntz = 1; sgnvz = 1;
        for (i = 1; i < model->nzz; i++) 
            inner_sweep();

        sgntz = -1; sgnvz = 0;
        for (i = model->nzz-2; i >= 0; i--) 
            inner_sweep();
    }
}

void Eikonal_modeling::inner_sweep()
{
    int i1 = i - sgnvz;
    int j1 = j - sgnvx;

    float tv = T[(i - sgntz) + j*model->nzz];
    float te = T[i + (j - sgntx)*model->nzz];
    float tev = T[(i - sgntz) + (j - sgntx)*model->nzz];

    float Sref;
    float t1D, t2D;
    float t1d1, t1d2;

    Sref = std::min(model->V[i1 + std::max(j - 1, 1)*model->nzz], model->V[i1 + std::min(j, model->nxx - 1)*model->nzz]);
    t1d1 = tv + model->dz * Sref; 
    
    Sref = std::min(model->V[std::max(i - 1, 1) + j1*model->nzz], model->V[std::min(i, model->nzz - 1) + j1*model->nzz]);
    t1d2 = te + model->dx * Sref; 

    t1D = std::min(t1d1, t1d2);

    float t1 = 1e6f;
    float t2 = 1e6f; 
    float t3 = 1e6f;
    
    Sref = model->V[i1 + j1*model->nzz];

    if ((tv <= te + model->dx * Sref) && (te <= tv + model->dz * Sref))
    {
        if ((te - tev >= 0.0) && (tv - tev >= 0.0))
        {
            float ta = tev + te - tv;
            float tb = tev - te + tv;
        
            t1 = ((tb*dz2i + ta*dx2i) + sqrtf(4.0*Sref*Sref*(dz2i + dx2i) - dz2i*dx2i*(ta - tb)*(ta - tb))) / (dz2i + dx2i);
        }
        else if ((te - tev <= Sref*model->dz*model->dz / sqrtf(model->dx*model->dx + model->dz*model->dz)) && (te - tev >= 0.0))
        {
            t2 = te + model->dx * sqrtf(Sref*Sref - ((te - tev) / model->dz)*((te - tev) / model->dz));
        }
        else if ((tv - tev <= Sref*model->dx*model->dx / sqrtf(model->dx*model->dx + model->dz*model->dz)) && (tv - tev >= 0.0))
        {
            t3 = tv + model->dz * sqrtf(Sref*Sref - ((tv - tev) / model->dx)*((tv - tev) / model->dx));
        }
    }    
    
    t2D = std::min(t1, std::min(t2, t3));

    T[i + j*model->nzz] = std::min(T[i + j*model->nzz], std::min(t1D, t2D));    
}

void Eikonal_modeling::build_outputs()
{
    if (export_receiver_output)
    {
        int current_node = 0;
        int iNode = geometry->iRel[shot_id];
        int fNode = geometry->fRel[shot_id];
        
        for (int node = iNode; node < fNode; node++)
        {    
            int index = geometry->nodes.idz[node] + geometry->nodes.idx[node] * model->nzz;

            receiver_output[current_node] = T[index];

            current_node += 1;
        }
    }

    if (export_wavefield_output)
    {
        for (int index = 0; index < model->nPoints; index++)
        {
            int i = (int)(index % model->nz);
            int j = (int)(index / model->nz);

            wavefield_output[index] = T[(i + model->nb) + (j + model->nb)*model->nzz]; 
        }
    }
}

void Eikonal_modeling::export_outputs()
{
    std::string receiver_output_name = receiver_output_folder + "first_arrivals_" + std::to_string(geometry->fRel[0]) + "_shot_" + std::to_string(shot_id+1) + ".bin";
    std::string wavefield_output_name = wavefield_output_folder + "travel_time_" + std::to_string(model->nz) + "x" + std::to_string(model->nx) + "_shot_" + std::to_string(shot_id+1) + ".bin";

    if (export_receiver_output) 
        write_binary_float(receiver_output_name, receiver_output, total_nodes);
    
    if (export_wavefield_output) 
        write_binary_float(wavefield_output_name, wavefield_output, model->nPoints);
}











