#include <iostream>
#include <fstream>
#include "CPP_Code/para.h"

int main()
{
    int sn = ceil(pow(N, 1.0 / d));
    cout << positions_input_file << endl;
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < d; j++)
    //     {
    //         cout << ((i / (int)(pow(sn, j))) % sn + 1 / 2.0) * L / sn << " ";
    //     }
    // }
    // cout << endl;

    float *arr = new float[N_times_d];
    for (int i = 0; i < N_times_d; i++)
    {
        *(arr+i) = (((i / d) / (int)(pow(sn, i % d))) % sn + 0.5) * L / sn;
    }

    ofstream os("CPP_Code/inputpositions.dat", ios::binary);
    os.write(reinterpret_cast<char*> (arr), N_times_d*sizeof(float));
    os.close();

    for (int i = 0; i < N_times_d; i++)
    {
        *(arr+i) = 0;
    }

    ifstream is("CPP_Code/inputpositions.dat", ios::binary);
    is.read(reinterpret_cast<char*>(arr), N_times_d*sizeof(float));

    for (int i = 0; i < N_times_d; i++)
    {
        cout << *(arr+i) << " ";
    }
    
    

    cout << endl;
    printf("%f", L);
}