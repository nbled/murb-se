#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "SimulationNBodyBarnesHut.hpp"

SimulationNBodyBarnesHut::SimulationNBodyBarnesHut(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit)
{
    this->flopsPerIte = 30.f * ((float)this->getBodies().getN() * (float)this->getBodies().getN() - (float)this->getBodies().getN())/2;
    this->accelerations.resize(this->getBodies().getN());
    this->softSquared = this->soft * this->soft;
}

void SimulationNBodyBarnesHut::freeTree(Octree *tree) {
    if(tree) {
        if (tree->internal == false) {
            free(tree);
            return;
        }

        for (int i=0;i<8;i++) {
            
            this->freeTree(tree->data.internal.children[i]);
        }

        free(tree);
    }
}

void SimulationNBodyBarnesHut::initIteration()
{
    for (unsigned long iBody = 0; iBody < this->getBodies().getN(); iBody++) {
        this->accelerations[iBody].ax = 0.f;
        this->accelerations[iBody].ay = 0.f;
        this->accelerations[iBody].az = 0.f;
    }
    this->computeOctree();
}
void SimulationNBodyBarnesHut::getBoundingBox(float *min_x_,float *max_x_,float *min_y_,float *max_y_,float *min_z_,float *max_z_) {
    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();
    float min_x = d[0].qx, max_x = d[0].qx, min_y = d[0].qy, max_y = d[0].qy, min_z = d[0].qz, max_z = d[0].qz;

    unsigned long n_bodies = this->getBodies().getN();

    for (unsigned long i=0;i<n_bodies;i++) {
        if (d[i].qx > max_x)
            max_x = d[i].qx;
        if (d[i].qx < min_x)
            min_x = d[i].qx;
        if (d[i].qy > max_y)
            max_y = d[i].qy;
        if (d[i].qy < min_y)
            min_y = d[i].qy;
        if (d[i].qz > max_z)
            max_z = d[i].qz;
        if (d[i].qz < min_z)
            min_z = d[i].qz;
    }

    *min_x_ = min_x;
    *max_x_ = max_x;
    *min_y_ = min_y;
    *max_y_ = max_y;
    *min_z_ = min_z;
    *max_z_ = max_z;


}

void SimulationNBodyBarnesHut::insertBody(Octree* tree, const dataAoS_t<float> *body) {
    if (tree->internal == false) { //End case : external node
        if (tree->data.external.body == NULL) {
            tree->data.external.body = body;
            tree->CoMx = body->qx;   //update right now CoM and mass
            tree->CoMy = body->qy;
            tree->CoMz = body->qz;
            tree->mass = body->m;
        } else {
            const dataAoS_t<float> *other_body = tree->data.external.body;
            // printf("body 1 : %e %e %e\n",body->qx,body->qy,body->qz);
            // printf("body 2 : %e %e %e\n",other_body->qx,other_body->qy,other_body->qz);

            tree->internal = true;
            // printf("constructing child %e\n",tree->size);
            for (int i=0;i<8;i++) {
                Octree *child = (Octree*) malloc(sizeof(Octree));
                child->internal = false; //Newly created child are external
                child->size = tree->size / 2;
                const float size_increment = child->size / 2;
                child->mass = 0; //we put the mass at 0 to simplify all later computation
                if (i % 2) { //new_x > x
                    child->qx = tree->qx + size_increment;
                } else {
                    child->qx = tree->qx - size_increment;
                }

                if (i % 4 > 1) {
                    child->qy = tree->qy + size_increment;
                } else {
                    child->qy = tree->qy - size_increment;
                }

                if (i > 3) {
                    child->qz = tree->qz + size_increment;
                } else {
                    child->qz = tree->qz - size_increment;
                }

                child->data.external.body = NULL;
                tree->data.internal.children[i] = child;
            }

            this->insertBody(tree,other_body); //Now that we have changed the type of the node we reinsert the original body, and the new one
            this->insertBody(tree,body);       //They will go on to the good external node.
            // this->updateTree(tree);

        }
        return;
    }
    //we go to the next group.
    int index = 0;
    index += (body->qx > tree->qx) * 1; //Permet de donner un indice différent selon tout les cas, de 0 à 7
    index += (body->qy > tree->qy) * 2;
    index += (body->qz > tree->qz) * 4;

    // printf("center : %e %e %e\n",tree->qx,tree->qy,tree->qz);
    // printf("%e %e %e going to %d\n",body->qx,body->qy,body->qz,index);

    this->insertBody(tree->data.internal.children[index], body);
    // this->updateTree(tree);
}

void SimulationNBodyBarnesHut::updateTree(Octree* tree) { //Could probably be fused with `insertBody`
    if (tree->internal == false) 
        return;


    float CoMx = 0,CoMy=0, CoMz=0, mass = 0;
    for (int i=0;i<8;i++) {
        this->updateTree(tree->data.internal.children[i]);
        float child_mass = tree->data.internal.children[i]->mass;
        CoMx += tree->data.internal.children[i]->CoMx * child_mass;
        CoMy += tree->data.internal.children[i]->CoMy * child_mass;
        CoMz += tree->data.internal.children[i]->CoMz * child_mass;
        mass += child_mass;

    }

    CoMx /= mass;
    CoMy /= mass;
    CoMz /= mass;
    tree->mass = mass;
    tree->CoMx = CoMx;
    tree->CoMy = CoMy;
    tree->CoMz = CoMz;
}

void SimulationNBodyBarnesHut::computeOctree() {
    // clock_t begin = clock();
    float min_x,max_x,min_y,max_y,min_z,max_z;
    this->getBoundingBox(&min_x,&max_x,&min_y,&max_y,&min_z,&max_z);

    float size_x = max_x - min_x;
    float size_y = max_y - min_y;
    float size_z = max_z - min_z;

    float size;
    if (size_x > size_y)
        size = size_x;
    else
        size = size_y;

    if (size_z > size)
        size = size_z;


    this->tree = (Octree*) malloc(sizeof(Octree));
    this->tree->size = size;
    this->tree->qx = (max_x - min_x) / 2 + min_x;
    this->tree->qy = (max_y - min_y) / 2 + min_y;
    this->tree->qz = (max_z - min_z) / 2 + min_z;
    this->tree->internal = false;
    this->tree->data.external.body = NULL;
    // printf("bouding box : %e %e,  %e %e,  %e %e\n",min_x,max_x,min_y,max_y,min_z,max_z);
    // printf("Created tree of size %e, at %e %e %e\n",size,this->tree->qx,this->tree->qy,this->tree->qz );

    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();

    unsigned long n_bodies = this->getBodies().getN();
    for (unsigned long i=0;i<n_bodies;i++) {
        // printf("inserting body %e %e %e\n",d[i].qx,d[i].qy,d[i].qz);
        this->insertBody(this->tree,&d[i]);
    }
    // clock_t end_1 = clock();
    this->updateTree(this->tree);
    // clock_t end = clock();
    // double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    // double time_spent_1 = (double)(end_1 - begin) / CLOCKS_PER_SEC;
    // double time_spent_2 = (double)(end - end_1) / CLOCKS_PER_SEC;

    // printf("boundingbox and insert : %f\n",time_spent_1);
    // printf("updateTree : %f\n",time_spent_2);
    // printf("compute Octree total: %f\n",time_spent);


}

void SimulationNBodyBarnesHut::computeBodyAcceleration(Octree *tree,const dataAoS_t<float> *body,float *ax, float *ay, float *az, int depth) {
    // printf("    entering node %e %e %e , size=%e, depth=%d\n",tree->qx,tree->qy,tree->qz,tree->size,depth);
    if (tree->internal == false) {
        if (tree->mass != 0) { 
            // printf("        with body %e %e %e of mass %e\n",tree->CoMx,tree->CoMy,tree->CoMz,tree->mass);
            const float rijx = tree->CoMx - body->qx; // 1 flop
            const float rijy = tree->CoMy - body->qy; // 1 flop
            const float rijz = tree->CoMz - body->qz; // 1 flop

            const float rijSquared = rijx * rijx + rijy * rijy + rijz * rijz; // 5 flops

            // compute the acceleration value between body i and body j: || ai || = G.mj / (|| rij ||² + e²)^{3/2}
            const float x = this->G / ((rijSquared + softSquared) * std::sqrt(rijSquared + softSquared));
            const float ai = x * tree->mass; // 1 flops

            *ax += ai * rijx;
            *ay += ai * rijy;
            *az += ai * rijz;
            // printf("    ax += %e  = %e\n",ai * rijx, *ax);
        } else {
            // printf("        empty\n");
        }
    } else {
        const float rijx = tree->CoMx - body->qx; // 1 flop
        const float rijy = tree->CoMy - body->qy; // 1 flop
        const float rijz = tree->CoMz - body->qz; // 1 flop

        const float rijSquared = rijx * rijx + rijy * rijy + rijz * rijz; // 5 flops

        const float s = tree->size * tree->size;

        if (s / rijSquared <= THETA * THETA) { // s/d < theta <=> s^2 / d^2 < theta^2, because everything is >= 0
            // printf("        (d=%d) with group %e %e %e of mass %e\n",depth,tree->CoMx,tree->CoMy,tree->CoMz,tree->mass);
            const float x = this->G / ((rijSquared + softSquared) * std::sqrt(rijSquared + softSquared));
            const float ai = x * tree->mass; // 1 flops

            *ax += ai * rijx;
            *ay += ai * rijy;
            *az += ai * rijz;

        } else { //We don't use this group
            // printf("        (d=%d) we go see the childs\n",depth);        
            for (int i=0;i<8;i++) {
                this->computeBodyAcceleration(tree->data.internal.children[i],body,ax,ay,az,depth + 1);
            }
        }
    }
}

void SimulationNBodyBarnesHut::computeBodiesAcceleration()
{
    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();

    unsigned long n_bodies = this->getBodies().getN();

    for (unsigned long iBody = 0; iBody < n_bodies; iBody++) {
        // printf("Computing for %e %e %e\n",d[iBody].qx,d[iBody].qy,d[iBody].qz);
        this->computeBodyAcceleration(this->tree,&d[iBody],&this->accelerations[iBody].ax,&this->accelerations[iBody].ay,&this->accelerations[iBody].az,0);
    }
}

void SimulationNBodyBarnesHut::computeOneIteration()
{
    // printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    this->initIteration();
    // clock_t begin = clock();
    this->computeBodiesAcceleration();
    // clock_t end = clock();
    // double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    // printf("compute Accelerations : %f\n",time_spent);

    this->freeTree(this->tree);
    // time integration
    this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
}
