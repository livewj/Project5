#include "lennardjones.h"
#include "system.h"
#include <math.h>

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System &system)
//computes the force betweet each set of atoms
{
    m_potentialEnergy = 0; // Remember to compute this in the loop

    double rCut = 2.5*m_sigma; //????
    double L = system.systemSize().x();
    double sigma6 = pow(m_sigma, 6); //sigma^6
    double energyAtRCut = 4*m_epsilon*pow(rCut, -6)*sigma6*(pow(rCut, -6)*sigma6 - 1.0);

    for (Atom *atom1 : system.atoms())
    {
        //By setting variables for atom1 we get a more tidy code
        double x1 = atom1->position.x();
        double y1 = atom1->position.y();
        double z1 = atom1->position.z();

        for (Atom *atom2 : atom1->neighbors)
        {
            double x2 = atom2->position.x();
            double deltax2x1 = x2 - x1;
            double y2 = atom2->position.y();
            double deltay2y1 = y2 - y1;
            double z2 = atom2->position.z();
            double deltaz2z1 = z2 - z1;
            //Minimum image convention:
            if (deltax2x1 < -0.5*L) deltax2x1 += L;  //x
            else if (deltax2x1 > 0.5*L) deltax2x1 -= L;
            if (deltay2y1 < -0.5*L) deltay2y1 += L;  //y
            else if (deltay2y1 > 0.5*L) deltay2y1 -= L;
            if (deltaz2z1 < -0.5*L) deltaz2z1 += L;  //z
            else if (deltaz2z1 > 0.5*L) deltaz2z1 -= L;

            double dr2 = deltax2x1*deltax2x1 + deltay2y1*deltay2y1 + deltaz2z1*deltaz2z1; // r^2 = x^2 + y^2 + z^2
            if (dr2 > rCut*rCut) continue;  //compute forces

            double dr = sqrt(dr2);
            double dr6 = dr2*dr2*dr2;
            double F = -4*m_epsilon*(-12.0*(sigma6*sigma6/(dr6*dr6*dr)) + 6.0*sigma6/(dr6*dr))*(1/dr); //Force KJHDHFS

            //Calculate forces between atom1 and atom2
            atom1->force[0] += F*deltax2x1;
            atom1->force[1] += F*deltay2y1;
            atom1->force[2] += F*deltaz2z1;

            atom2->force[0] -= F*deltax2x1;
            atom2->force[1] -= F*deltay2y1;
            atom2->force[2] -= F*deltaz2z1;

            m_potentialEnergy += 4*m_epsilon*(sigma6*sigma6/(dr6*dr6) - sigma6/dr6) - energyAtRCut;
            //PRESSURE?
        }
    }

}


























