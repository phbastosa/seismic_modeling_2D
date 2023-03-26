# ifndef EIKONAL_MODELING_HPP
# define EIKONAL_MODELING_HPP

# include "../modeling.hpp"

# include "../../model/eikonal/eikonal_model.hpp"

# include "../../wavelet/gaussian_2nd/gaussian_2nd.hpp"

# include "../../geometry/regular/regular.hpp"
# include "../../geometry/streamer/streamer.hpp"

class Eikonal_modeling : public Modeling
{
private:    

    int i, j;
    int sgnvx, sgnvz;
    int sgntx, sgntz;
    float dx2i, dz2i;

    float * T = nullptr;

    void inner_sweep();
    void build_outputs();
    void kernel_propagation();

public:

    void propagation();
    void set_components();
    void set_wavefields();
    void set_parameters();
};

# endif