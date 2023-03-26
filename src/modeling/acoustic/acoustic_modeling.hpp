# ifndef ACOUSTIC_MODELING_HPP
# define ACOUSTIC_MODELING_HPP

# include "../modeling.hpp"

# include "../../model/acoustic/acoustic_model.hpp"

# include "../../wavelet/gaussian_1st/gaussian_1st.hpp"

# include "../../geometry/regular/regular.hpp"
# include "../../geometry/streamer/streamer.hpp"

class Acoustic_modeling : public Modeling
{
private:    

    int time_id;
    int total_times;

    int snap_id;
    int i_snap, f_snap;
    int d_snap, n_snap;

    float * P = nullptr;
    float * Vx = nullptr;
    float * Vz = nullptr;

    float factor;
    float * damp1D = nullptr;
    float * damp2D = nullptr;

    void set_abc_dampers();

    void build_outputs();
    void apply_wavelet();
    void kernel_propagation();

public:

    void propagation();
    void set_components();
    void set_wavefields();

    void set_parameters();
};

# endif