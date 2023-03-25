# include <cmath>

# include "gaussian_2nd.hpp"

void Gaussian_2nd::set_parameters(std::string file)
{
    nt = std::stoi(catch_parameter("nt", file));
    dt = std::stof(catch_parameter("dt", file));

    fmax = std::stof(catch_parameter("fmax", file));
    tlag = std::stof(catch_parameter("tlag", file));
    gain = std::stof(catch_parameter("gain", file));

    wavelet_file = catch_parameter("wavelet_file", file);
    import_wavelet = str2bool(catch_parameter("import_wavelet", file));

    amp = new float[nt]();

    if (import_wavelet) 
        read_binary_float(wavelet_file, amp, nt);
    else
        build_amplitudes();
}

void Gaussian_2nd::build_amplitudes()
{
    int n_freq = 10000;

    float factor_real;
    float factor_imag;

    float pi = 4.0f * atanf(1.0f);  
    float fc = fmax / (3.0f * sqrtf(pi));

    float * aux_amp = new float[nt]();

    for (int n = 0; n < nt; n++)
    {        
        float aux1 = 1.0f - 2.0f*pi*powf(n*dt - tlag, 2.0f) * powf(fc, 2.0f) * powf(pi, 2.0f);
        float aux2 = expf(-pi * powf(n*dt - tlag, 2.0f) * powf(fc, 2.0f) * powf(pi, 2.0f));    
        
        aux_amp[n] = aux1 * aux2;  
    }

    float * omega = new float[n_freq]();

    for (int w = 0; w < n_freq; w++)
        omega[w] = -pi + 2*pi*w/n_freq;

    float dw = fabsf(fabsf(omega[1]) - fabsf(omega[0])); 

    float * input_real = new float[n_freq]();
    float * input_imag = new float[n_freq]();
    
    for (int n = 0; n < nt; n++)
    {  
        for (int w = 0; w < n_freq; w++)
        {
            input_real[w] += aux_amp[n] * cosf(omega[w] * n);
            input_imag[w] -= aux_amp[n] * sinf(omega[w] * n);
        }
    }

    for (int w = 0; w < n_freq; w++)
    {
        float factor = 0.5f * sqrtf(2.0f * fabsf(omega[w]));

        if (omega[w] < 0)
        {
            input_real[w] = factor * (input_real[w] + input_imag[w]);
            input_imag[w] = factor * (input_real[w] - input_imag[w]);
        }
        else
        {
            input_real[w] = factor * (input_real[w] - input_imag[w]);
            input_imag[w] = factor * (input_real[w] + input_imag[w]);
        }
    }
    
    float max = 0.0f;
    for (int n = 0; n < nt; n++)
    {  
        for (int w = 0; w < n_freq; w++)
            amp[n] += (input_real[w] * cosf(omega[w] * n) - input_imag[w] * sinf(omega[w] * n)) * dw;
        
        if (amp[n] > max) max = amp[n];
    }    
    
    for (int n = 0; n < nt; n++)
        amp[n] *= gain / max;

    delete[] omega;
    delete[] aux_amp;
    delete[] input_real;
    delete[] input_imag;
}