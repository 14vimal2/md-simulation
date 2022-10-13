#ifndef _INCL_BIN_DATA
#define _INCL_BIN_DATA

#include <string>
#include <fstream>

using namespace std;

///@brief Save data in binary format
///@param filename filename e.g "positions.dat"
///@param ptr pointer to the first value of data to be saved
///@param size total size of data i.e. if 5 int then 5 * sizeof(5)
template <typename _type_>
void savedata(string filename, _type_ *ptr, unsigned long long size){
    ofstream f(filename, ios::binary);
    f.write(reinterpret_cast<char*> (ptr), size);
    f.close();
}

///@brief Read data in binary format
///@param filename filename e.g "positions.dat"
///@param ptr pointer to the first value of data to be filled
///@param size total size of data i.e. if 5 int then 5 * sizeof(5)
template <typename _type_>
void readdata(string filename, _type_ *ptr, unsigned long long size){
    ifstream f(filename, ios::binary);
    f.read(reinterpret_cast<char*>(ptr), size);
    f.close();
}

#endif

