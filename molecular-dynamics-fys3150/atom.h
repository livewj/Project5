#ifndef ATOM_H
#define ATOM_H
#include "math/vec3.h"
#include <vector>
using std::vector;

class Atom
{
private:
    float m_mass;
public:
    vec3 position;
    vec3 velocity;
    vec3 force;
    vec3 initialPosition; //FFC

    Atom(double mass);
    void resetForce();
    void resetVelocityMaxwellian(double temperature);

    double mass() { return m_mass; }
    void setMass(double mass) { m_mass = mass; }
};
#endif
