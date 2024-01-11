#ifndef SIMULATION_N_BODY_BARNES_HUT_HPP_
#define SIMULATION_N_BODY_BARNES_HUT_HPP_

#include <string>

#include "core/SimulationNBodyInterface.hpp"
#include "core/Bodies.hpp"

#define THETA 1.0f //0.3 doesn't pass the tests

struct Octree;

struct internal_node {
  Octree* children[8];
};

struct external_node {
  const dataAoS_t<float>* body;
};


struct Octree {
  bool internal; //1 for internal node, 0 for external
  float size;
  float qx,qy,qz; //Position of the middle of the group.
  float CoMx,CoMy,CoMz; //Center of mass of the group
  float mass; //Total mass of the group  
  union {
    struct internal_node internal;
    struct external_node external;
  } data;
};


class SimulationNBodyBarnesHut : public SimulationNBodyInterface {
  protected:
    std::vector<accAoS_t<float>> accelerations; /*!< Array of body acceleration structures. */
    Octree *tree;
    float softSquared;

  public:
    SimulationNBodyBarnesHut(const unsigned long nBodies, const std::string &scheme = "galaxy", const float soft = 0.035f,
                         const unsigned long randInit = 0);
    virtual ~SimulationNBodyBarnesHut() = default;
    virtual void computeOneIteration();

  protected:
    void initIteration();
    virtual void computeBodiesAcceleration();
    void getBoundingBox(float*,float*,float*,float*,float*,float*);
    void insertBody(Octree* tree, const dataAoS_t<float> *body);
    void updateTree(Octree* tree);
    void computeOctree();
    void freeTree(Octree *tree);
    void computeBodyAcceleration(Octree *tree,const dataAoS_t<float> *body,float *ax, float *ay, float *az, int depth);

};

#endif /* SIMULATION_N_BODY_OPTIM_HPP_ */
