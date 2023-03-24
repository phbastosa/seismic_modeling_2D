# ifndef GAUSSIAN_1ST_HPP
# define GAUSSIAN_1ST_HPP

# include "../wavelet.hpp"

class Gaussian_1st : public Wavelet
{
public:

    void build_amplitudes();

    void set_parameters(std::string file);
};

# endif