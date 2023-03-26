# ifndef SCALAR_MODELING_HPP
# define SCALAR_MODELING_HPP

# include "../modeling.hpp"

# include "../../model/scalar/scalar_model.hpp"

# include "../../wavelet/gaussian_2nd/gaussian_2nd.hpp"

# include "../../geometry/regular/regular.hpp"
# include "../../geometry/streamer/streamer.hpp"

class Scalar_modeling : public Modeling
{
private:    

    int time_id;
    int total_times;

    int snap_id;
    int i_snap, f_snap;
    int d_snap, n_snap;

    float * U_pre = nullptr;
    float * U_pas = nullptr;
    float * U_fut = nullptr;

    float factor;
    float * damp1D = nullptr;
    float * damp2D = nullptr;

    void build_outputs();
    void apply_wavelet();
    void set_abc_dampers();
    void update_wavefield();
    void kernel_propagation();

public:

    void propagation();
    void set_components();
    void set_wavefields();
    void set_parameters();
};

# endif