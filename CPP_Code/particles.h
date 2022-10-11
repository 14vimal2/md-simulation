#include "para.h"
/// @brief Particles class contains information of positions and previous positions at all instant of time, it also provides all methods to calculate necessary quantities
template <class _int_type_, class _dec_type_>
class Particles
{
private:
    /**
     * @brief Head pointer of array containing positions of the particles
     * For N particles and d dimensions it is stored like P_0_0, P_0_1, ..., P_i_j, ..., P_(N-1)_(d-1), i refers to ith particle and j refers to jth dimension
     * i.e. - 3 particles in 2 dimensiont P_0_0, P_0_1, P_1_0, P_1_1, P_2_0, P_2_1
     */
    _dec_type_ *r;

    /**
     * @brief Head pointer of array containing previous positions of the particles
     * For N particles and d dimensions it is stored like P_0_0, P_0_1, ..., P_i_j, ..., P_(N-1)_(d-1), i refers to ith particle and j refers to jth dimension
     * i.e. - 3 particles in 2 dimensiont P_0_0, P_0_1, P_1_0, P_1_1, P_2_0, P_2_1
     */
    _dec_type_ *rp;
    //
    /**
     * @brief Head pointer of array containing force on the particles
     * For N particles and d dimensions it is stored like P_0_0, P_0_1, ..., P_i_j, ..., P_(N-1)_(d-1), i refers to ith particle and j refers to jth dimension
     * i.e. - 3 particles in 2 dimensiont P_0_0, P_0_1, P_1_0, P_1_1, P_2_0, P_2_1
     */
    _dec_type_ *F;
    /**
     * @brief Head pointer of array containing temporary variables of the particles
     * It is used as reserve memory for temporary variables
     * For N particles and d dimensions it is stored like P_0_0, P_0_1, ..., P_i_j, ..., P_(N-1)_(d-1), i refers to ith particle and j refers to jth dimension
     * i.e. - 3 particles in 2 dimensiont P_0_0, P_0_1, P_1_0, P_1_1, P_2_0, P_2_1
     */
    _dec_type_ *temp_var;
public:
    Particles(){
        r = new _dec_type_[N_times_d];
        rp = new _dec_type_[N_times_d];
        F = new _dec_type_[N_times_d];
        temp_var = new _dec_type_[N_times_d];
    };
    /** 
     * Calculates current position at each time-step
     * Use equation 𝒓(𝑡+Δ𝑡)=2 𝒓(𝑡)−𝒓(𝑡−Δ𝑡)+ f ((𝑟_𝑖𝑗 ) ⃗ )Δt^2
    */
    void time_advance_single_step()
    {
        for (_int_type_ i = 0; i < N_times_d; i++)
        {
            *(temp_var+i) = *(r+i);
            *(r + i) = 2 * (r + i) - *(rp + i) + 2 * L * ((*(r + i) - *(rp + i)) < L / 2) - 2 * L * ((*(r + i) - *(rp + i)) > L / 2) + *(F + i) * dt * dt;
            *(r + i) = *(r + i) + (*(r + i) < 0) * L - (*(r + i) > L) * L;
            *(rp+i) = *(temp_var+i);
        }
    }

    ~Particles(){
        delete[] r;
        delete[] rp;
        delete[] F;
        delete[] temp_var;
    };
};