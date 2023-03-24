# ifndef GAUSSIAN_2ND_HPP
# define GAUSSIAN_2ND_HPP

# include "../wavelet.hpp"

class Gaussian_2nd : public Wavelet
{
public:
    
    void set_parameters(std::string file);
    void build_amplitudes();
};

# endif