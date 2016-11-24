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
/*
//finding the nearest neighbors of the atom
void LennardJones::CreateListOfNeighbors(System &system)
{
    for(Atom *atom : system.atoms()) {
        atom->neighbors.clear();
        atom->ListOfNeighborsVec = atom->position;
    }
    system.cellList().cutoffDistance = rShell + rCut;
    system.cellList().build(&system);
    CellList &cellList = system.cellList();
    double L = system.systemSize().x();
    double neighborListRadiusSquared = (rShell+rCut)*(rShell+rCut);
    for(int cx = 0; cx<cellList.numberOfCellsPerDimension; cx++) {
        for(int cy = 0; cy<cellList.numberOfCellsPerDimension; cy++) {
            for(int cz = 0; cz<cellList.numberOfCellsPerDimension; cz++) {
                vector<Atom*> &cell1 = cellList.cell(cx,cy,cz);
                for(int dx = -1; dx<=1; dx++) {
                    for(int dy = -1; dy<=1; dy++) {
                        for(int dz = -1; dz<=1; dz++) {
                            int i = (cx + dx + cellList.numberOfCellsPerDimension) % cellList.numberOfCellsPerDimension;
                            int j = (cy + dy + cellList.numberOfCellsPerDimension) % cellList.numberOfCellsPerDimension;
                            int k = (cz + dz + cellList.numberOfCellsPerDimension) % cellList.numberOfCellsPerDimension;
                            vector<Atom*> &cell2 = cellList.cell(i,j,k);

                            for(Atom *atomi : cell1) {
                                for(Atom *atomj : cell2) {
                                    if(atomi<=atomj) continue;
                                    double dx = atomi->position.x() - atomj->position.x();
                                    double dy = atomi->position.y() - atomj->position.y();
                                    double dz = atomi->position.z() - atomj->position.z();
                                    if(dx < -0.5*L) dx += L;
                                    else if(dx > 0.5*L) dx -= L;

                                    if(dy < -0.5*L) dy += L;
                                    else if(dy > 0.5*L) dy -= L;

                                    if(dz < -0.5*L) dz += L;
                                    else if(dz > 0.5*L) dz -= L;

                                    double dr2 = dx*dx + dy*dy + dz*dz;
                                    if(dr2<neighborListRadiusSquared) {
                                        atomi->neighbors.push_back(atomj);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

*/

void LennardJones::calculateForces(System &system)
//computes the force betweet each set of atoms
{
    m_potentialEnergy = 0; // Remember to compute this in the loop

    double L = system.systemSize().x();
    double sigma6 = pow(m_sigma, 6); //sigma^6

    for (Atom *atomi : system.atoms())
    {
        //By setting variables for atomi=atom1 and atomj=atom2 we get a more tidy code
        double x1 = atomi->position.x();
        double y1 = atomi->position.y();
        double z1 = atomi->position.z();

        //for (Atom *atomj : atomi->neighbors)
        for (Atom *atomj : system.atoms())
        {
            if(atomi<=atomj) continue;

            double x2 = atomj->position.x();
            double y2 = atomj->position.y();
            double z2 = atomj->position.z();

            //Distance between atoms
            double deltax2x1 = x2 - x1;
            double deltay2y1 = y2 - y1;
            double deltaz2z1 = z2 - z1;

            //Minimum image convention: Finding the distance between atom i and j
            if (deltax2x1 < -0.5*L) deltax2x1 += L;  //x
            else if (deltax2x1 > 0.5*L) deltax2x1 -= L;
            if (deltay2y1 < -0.5*L) deltay2y1 += L;  //y
            else if (deltay2y1 > 0.5*L) deltay2y1 -= L;
            if (deltaz2z1 < -0.5*L) deltaz2z1 += L;  //z
            else if (deltaz2z1 > 0.5*L) deltaz2z1 -= L;

            double dr2 = deltax2x1*deltax2x1 + deltay2y1*deltay2y1 + deltaz2z1*deltaz2z1; // dr^2 = x^2 + y^2 + z^2
            double dr = sqrt(dr2);
            double dr6 = pow(dr, 6); // dr2*dr2*dr2;
            double F = -4*m_epsilon*(-12.0*(sigma6*sigma6/(dr6*dr6*dr)) + 6.0*sigma6/(dr6*dr))*(1/dr); //Force

            //Calculate forces between atomi and atomj
            atomi->force[0] += F*deltax2x1;
            atomi->force[1] += F*deltay2y1;
            atomi->force[2] += F*deltaz2z1;

            atomj->force[0] -= F*deltax2x1;
            atomj->force[1] -= F*deltay2y1;
            atomj->force[2] -= F*deltaz2z1;

            m_potentialEnergy += 4*m_epsilon*(sigma6*sigma6/(dr6*dr6) - sigma6/dr6);
        }
    }
}

