#include <iostream>
#include <fstream>
#include <ctime>
#include <random>
#include "para.h"
#include "bindata.h"


///@brief Create position data
///@param N number of the particles only required for getting type
///@param boxsize size of the box only required for getting type
template <class _int_type_, class _dec_type_>
void create_position_data(_int_type_ N, _dec_type_ boxsize)
{
    _int_type_ sn = ceil(pow(N, 1.0 / d));
    _dec_type_ *pos = new _dec_type_[N_times_d];
    for (_int_type_ i = 0; i < N_times_d; i++)
    {
        *(pos + i) = (( (i / d) / (_int_type_)(pow(sn, i % d))) % sn + (1.0 / 2.0)) * L / sn;
    }

    string pif = "input_files/postions_N" + to_string(N) + "_d" + to_string(d) + "_rho" + to_string(rho) + "_squaregrid.dat";
    savedata(pif, pos, variables_datasize);
    for (_int_type_ i = 0; i < N_times_d; i++)
    {
        *(pos + i) = 0;
    }
    readdata(pif, pos, variables_datasize);
    cout << "postions" << endl;
    for (_int_type_ i = 0; i < N_times_d; i++)
    {
        cout << *(pos + i) << " ";
    }
    cout << endl
         << pif << " created";
    delete[] pos;
}

///@brief Create velocity data
///@param N number of the particles only required for getting type
///@param temperature temperature of the particles only required for getting type
template <class _int_type_, class _dec_type_>
void create_velocity_data(_int_type_ N, _dec_type_ temperature)
{
    string vif = "input_files/velocities_N" + to_string(N) + "_d" + to_string(d) + "_T" + to_string(T) + "_random.dat";
    srand(time(NULL));
    _dec_type_ *vel = new _dec_type_[N_times_d];
    _dec_type_ *vcm = new _dec_type_[d];
    fill_n(vcm, d, 0);
    for (_int_type_ i = 0; i < N_times_d; i++)
    {
        *(vel + i) = (rand() / (_dec_type_)(RAND_MAX)) - 1 / 2.0;
        *(vcm + (i % d)) += *(vel + i);
    }

    for (_int_type_ j = 0; j < d; j++)
    {
        *(vcm + j) /= N;
    }

    _dec_type_ KE = 0;

    for (_int_type_ i = 0; i < N_times_d; i++)
    {
        *(vel + i) -= vcm[i % d];
        KE += *(vel + i) * *(vel + i);
    }

    _dec_type_ scalefactor = sqrt(d * T * N / KE);
    KE = 0;

    for (_int_type_ i = 0; i < N_times_d; i++)
    {
        *(vel + i) *= scalefactor;
        KE += *(vel + i) * *(vel + i);
    }

    cout << "temperature " << KE / N / d << endl;

    savedata(vif, vel, variables_datasize);

    fill_n(vel, N_times_d, 0);
    readdata(vif, vel, variables_datasize);
    cout << "velocities" << endl;
    for (_int_type_ i = 0; i < N_times_d; i++)
    {
        cout << *(vel + i) << " ";
    }
    delete[] vel;
    delete[] vcm;
}

int main()
{
    create_position_data(N, L);
    cout << endl
         << "boxsize = " << L << endl;
    create_velocity_data(N, T);
}