# ifndef ELASTIC_MODELING_HPP
# define ELASTIC_MODELING_HPP

# include "../modeling.hpp"

# include "../../model/elastic/elastic_model.hpp"

# include "../../wavelet/gaussian_1st/gaussian_1st.hpp"

class Elastic_modeling : public Modeling
{
private:    

    int time_id;
    int total_times;

    int snap_id;
    int i_snap, f_snap;
    int d_snap, n_snap;

    float * Vx = nullptr;
    float * Vz = nullptr;
    float * Txx = nullptr;
    float * Tzz = nullptr;
    float * Txz = nullptr;

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