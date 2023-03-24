# ifndef WAVELET_HPP
# define WAVELET_HPP

# include "../utils/io/io.hpp"

class Wavelet
{
protected:

    float fmax;
    float tlag;
    float gain;

    bool import_wavelet;
    std::string wavelet_file;

    virtual void build_amplitudes() = 0;

public:

    int nt;
    float dt;

    float * amp = nullptr;

    virtual void set_parameters(std::string file) = 0;
};

# endif