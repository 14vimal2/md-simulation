#include "para.h"
/// @brief Particles class contains information of positions and previous positions at all instant of time, it also provides all methods to calculate necessary quantities
template <class _int_type_, class _dec_type_>
class Particles
{
private:
    /**
     * @brief Head pointer of array containing positions of the particles
     * For N particles and d dimensions it is stored like P_0_0, P_0_1, ..., P_i_j, ..., P_(N-1)_(d-1), i refers to ith particle and j refers to jth dimension
     * i.e. - 3 particles in 2 dimensions P_0_0, P_0_1, P_1_0, P_1_1, P_2_0, P_2_1
     */
    _dec_type_ *r;

    /**
     * @brief Head pointer of array containing previous positions of the particles
     * For N particles and d dimensions it is stored like P_0_0, P_0_1, ..., P_i_j, ..., P_(N-1)_(d-1), i refers to ith particle and j refers to jth dimension
     * i.e. - 3 particles in 2 dimensions P_0_0, P_0_1, P_1_0, P_1_1, P_2_0, P_2_1
     */
    _dec_type_ *rp;
    //
    /**
     * @brief Head pointer of array containing force on the particles
     * For N particles and d dimensions it is stored like P_0_0, P_0_1, ..., P_i_j, ..., P_(N-1)_(d-1), i refers to ith particle and j refers to jth dimension
     * i.e. - 3 particles in 2 dimensions P_0_0, P_0_1, P_1_0, P_1_1, P_2_0, P_2_1
     */
    _dec_type_ *F;
    /** 
     * @brief Head pointer of array containing temporary variables of the particles
     * It is used as reserve memory for temporary variables
     * For N particles and d dimensions it is stored like P_0_0, P_0_1, ..., P_i_j, ..., P_(N-1)_(d-1), i refers to ith particle and j refers to jth dimension
     * i.e. - 3 particles in 2 dimensions P_0_0, P_0_1, P_1_0, P_1_1, P_2_0, P_2_1
     */
    _dec_type_ *temp_var;

public:
    Particles()
    {
        r = new _dec_type_[N_times_d];
        rp = new _dec_type_[N_times_d];
        F = new _dec_type_[N_times_d];
        temp_var = new _dec_type_[N_times_d];
    };
    /**
     * Calculates current position at each time-step
     * Uses equations 𝒓(𝑡+Δ𝑡)=2 𝒓(𝑡)−𝒓(𝑡−Δ𝑡)+𝒂(𝒓_𝑖𝑗 )Δt^2, where 𝒂(𝒓_𝑖𝑗 ) = 𝒇/𝑚, with m = 1.
     * For performance reasons we are using the concept of branchless programming.
     * So, + 2 * L * ((r[i] - rp[i]) < L / 2), - 2 * L * ((r[i] - rp[i]) > L / 2), + (r[i] < 0) * L, and - (r[i] > L) * L  terms are used for dealing with periodic boundary conditions
     */
    void time_advance_single_step()
    {
        for (_int_type_ i = 0; i < N_times_d; i++)
        {
            temp_var[i] = r[i];
            r[i] = 2 * r[i] - rp[i] + 2 * L * ((r[i] - rp[i]) < L / 2) - 2 * L * ((r[i] - rp[i]) > L / 2) + F[i] * dt * dt;
            r[i] = r[i] + (r[i] < 0) * L - (r[i] > L) * L;
            rp[i] = temp_var[i];
            // *(temp_var + i) = *(r + i);
            // *(r + i) = 2 * (r + i) - *(rp + i) + 2 * L * ((*(r + i) - *(rp + i)) < L / 2) - 2 * L * ((*(r + i) - *(rp + i)) > L / 2) + *(F + i) * dt * dt;
            // *(r + i) = *(r + i) + (*(r + i) < 0) * L - (*(r + i) > L) * L;
            // *(rp + i) = *(temp_var + i);
        }
    }
    /**
     * @brief Get the velocities of the particles
     * Calculates velocities using equation 𝒗(𝑡)= (𝒓(𝑡+Δ𝑡)−𝒓(𝑡−Δ𝑡))/(2 Δ𝑡)
     * @return _dec_type_* Head pointer to the array containing velocities
     */
    _dec_type_ *get_velocities()
    {
        for (_int_type_ i = 0; i < N_times_d; i++)
        {
            temp_var[i] = (r[i] - rp[i] + L * (r[i] - rp[i] < -L / 2) - L * (r[i] - rp[i] > L / 2)) / (2 * dt);
        }
        return temp_var;
    }
    /**
     * @brief Get Kinetic Energy of all particles
     * Calculated using KE=(∑_(𝑖=0)^𝑁▒𝑣_𝑖^2 )/2
     * @return _dec_type_ Kinetic Energy
     */
    _dec_type_ get_KE()
    {
        get_velocities();
        _dec_type_ KE = 0;
        for (int i = 0; i < N_times_d; i++)
        {
            KE += temp_var[i] * temp_var[i];
        }
        return KE / 2;
    }
    /**
     * @brief Get Temperature of the particles
     * Calculates using equation 𝑇=(∑_(𝑖=0)^𝑁▒𝑣_𝑖^2 )/(𝑁 × 𝑑)
     * @return _dec_type_ 
     */
    _dec_type_ get_T(){
        _dec_type_ T = get_KE();
        return 2*T/N_times_d;
    }
    ~Particles()
    {
        delete[] r;
        delete[] rp;
        delete[] F;
        delete[] temp_var;
    };
};