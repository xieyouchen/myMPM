#pragma once
#include <iostream>
#include <Eigen/Core>
#include <istream>
#include <string>
#include <fstream>
#include <vector>

#define MPM_ASSERT(condition, statement)                                       \
    if (!(condition)) {                                                        \
      cout << statement << endl;                                               \
      assert(condition);                                                       \
    }                                                                          \


using namespace std;
using namespace Eigen;

struct SimulatorConfiguration
{
    int W, H, L, grid_size;
};

class Water
{
public:
    float E = 500.0f; // 杨氏模量
    float nu = 0.4f; // 泊松比
    float mass = 0.001f; // 质量
    float density = 1.0f; //密度
    Water() {}
};