/// @brief Particles class contains information of positions and previous positions at all instant of time, it also provides all methods to calculate necessary quantities
class Particles
{
private:
    //Positions of particles
    int *r;
    //Previous positions of particles
    int *rp;
    //

    /* data */
public:
    Particles(/* args */);
    ~Particles();
};

Particles::Particles(/* args */)
{
}

Particles::~Particles()
{
}
