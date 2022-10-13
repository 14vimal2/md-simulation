#include "para.h"
#include "bindata.h"
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
    /**
     * @brief Head pointer of array containing velocities of the particles
     * For N particles and d dimensions it is stored like P_0_0, P_0_1, ..., P_i_j, ..., P_(N-1)_(d-1), i refers to ith particle and j refers to jth dimension
     * i.e. - 3 particles in 2 dimensions P_0_0, P_0_1, P_1_0, P_1_1, P_2_0, P_2_1
     */
    _dec_type_ *v;
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

    bool velocities_updated;

public:
    Particles()
    {
        r = new _dec_type_[N_times_d];
        rp = new _dec_type_[N_times_d];
        v = new _dec_type_[N_times_d];
        F = new _dec_type_[N_times_d];
        temp_var = new _dec_type_[N_times_d];
        velocities_updated = true;
        init();
    };
    void init()
    {
        readdata("CPP_Code/input_files/" + positions_input_file, r, variables_datasize);
        readdata("CPP_Code/input_files/" + velocities_input_file, v, variables_datasize);
        for (_int_type_ i = 0; i < N_times_d; i++)
        {
            rp[i] = r[i] - v[i] * dt;
        }
    }

    void computeForce()
    {
        fill_n(F, N_times_d, 0);
    }

    void DisplayPositions()
    {
        cout << "Positions" << endl;
        for (int i = 0; i < N_times_d; i++)
        {
            cout << r[i] << " ";
        }
        cout << endl;
    }
    void DisplayTempVars()
    {
        cout << "TempVars" << endl;
        for (int i = 0; i < N_times_d; i++)
        {
            cout << temp_var[i] << " ";
        }
        cout << endl;
    }
    void DisplayVelocities()
    {
        cout << "velocities" << endl;
        for (int i = 0; i < N_times_d; i++)
        {
            cout << v[i] << " ";
        }
        cout << endl;
    }
    /**
     * Calculates current position at each time-step
     * Uses equations 𝒓(𝑡+Δ𝑡)=2 𝒓(𝑡)−𝒓(𝑡−Δ𝑡)+𝒂(𝒓_𝑖𝑗 )Δt^2, where 𝒂(𝒓_𝑖𝑗 ) = 𝒇/𝑚, with m = 1.
     * For performance reasons we are using the concept of branchless programming.
     * So, + 2 * L * ((r[i] - rp[i]) < L / 2), - 2 * L * ((r[i] - rp[i]) > L / 2), + (r[i] < 0) * L, and - (r[i] > L) * L  terms are used for dealing with periodic boundary conditions
     */
    void time_advance_single_step(bool update_velocity = false)
    {
        velocities_updated = update_velocity;
        if (update_velocity)
        {
            for (_int_type_ i = 0; i < N_times_d; i++)
            {
                temp_var[i] = r[i];                                                                                             // saves current position
                r[i] = 2 * r[i] - rp[i] + 2 * L * ((r[i] - rp[i]) < L / 2) - 2 * L * ((r[i] - rp[i]) > L / 2) + F[i] * dt * dt; // current position gets  new positions
                r[i] = r[i] + (r[i] < 0) * L - (r[i] > L) * L;                                                                  // fixed boundary conditions of new positions
                v[i] = r[i] - rp[i];                                                                                            // new position - previous position
                v[i] = (v[i] + L * (v[i] < -L / 2) - L * (v[i] > L / 2)) / (2 * dt);
                rp[i] = temp_var[i]; // previous position gets current position
            }
        }
        else
        {
            for (_int_type_ i = 0; i < N_times_d; i++)
            {
                temp_var[i] = r[i];                                                                                             // saves current position
                r[i] = 2 * r[i] - rp[i] + 2 * L * ((r[i] - rp[i]) < L / 2) - 2 * L * ((r[i] - rp[i]) > L / 2) + F[i] * dt * dt; // current position gets  new positions
                r[i] = r[i] + (r[i] < 0) * L - (r[i] > L) * L;                                                                  // fixed boundary conditions of new positions
                rp[i] = temp_var[i];                                                                                            // previous position gets current position                                                                                           // previous position gets current position
            }
        }
    }
    /**
     * @brief Get the velocities of the particles
     * Calculates velocities using equation 𝒗(𝑡)= (𝒓(𝑡+Δ𝑡)−𝒓(𝑡−Δ𝑡))/(2 Δ𝑡)
     * @return _dec_type_* Head pointer to the array containing velocities
     */
    _dec_type_ *get_velocities()
    {
        if (velocities_updated)
        {
            return v;
        }
        cout << "calculating velocities" << endl;
        for (_int_type_ i = 0; i < N_times_d; i++)
        {
            temp_var[i] = 2 * r[i] - rp[i] + 2 * L * ((r[i] - rp[i]) < L / 2) - 2 * L * ((r[i] - rp[i]) > L / 2) + F[i] * dt * dt;
            temp_var[i] = temp_var[i] + (temp_var[i] < 0) * L - (temp_var[i] > L) * L;
            v[i] = (temp_var[i] - rp[i] + L * (temp_var[i] - rp[i] < -L / 2) - L * (temp_var[i] - rp[i] > L / 2)) / (2 * dt);
            cout << v[i] << " ";
        }
        return v;
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
            KE += v[i] * v[i];
        }
        return KE / 2;
    }
    /**
     * @brief Get Temperature of the particles
     * Calculates using equation 𝑇=(∑_(𝑖=0)^𝑁▒𝑣_𝑖^2 )/(𝑁 × 𝑑)
     * @return _dec_type_
     */
    _dec_type_ get_T()
    {
        _dec_type_ T = get_KE();
        return 2 * T / N_times_d;
    }
    ~Particles()
    {
        delete[] r;
        delete[] rp;
        delete[] F;
        delete[] temp_var;
    };
};