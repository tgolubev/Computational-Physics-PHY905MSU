#ifndef LENNARDJONES_H
#define LENNARDJONES_H

class LennardJones
{
private:
    double m_sigma = 1.0;  //use Angstroms. 3.405A is for Ar
    double m_epsilon = 1.0;
    double m_potentialEnergy = 0;
    double m_four_epsilon = 4.0*m_epsilon;
    double m_twntyfour_epsilon = 24.0*m_epsilon;

public:
    LennardJones() { }
    void calculateForces(class System &system);
    double potentialEnergy() const;
    double sigma() const;
    void setSigma(double sigma);
    double epsilon() const;
    //double twntyfour_epsilon() const;
    void setEpsilon(double epsilon);
};
#endif
