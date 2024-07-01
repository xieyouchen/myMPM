#pragma once
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <istream>
#include <string>
#include <fstream>
#include <vector>

#define MPM_ASSERT(condition, statement)                                       \
    if (!(condition)) {                                                        \
      cout << statement << endl;                                               \
      assert(condition);                                                       \
    }                                                                          \

#define MPM_INFO(statement) \
    cout << "- " statement << endl;

using namespace std;
using namespace Eigen;

struct SimulatorConfiguration
{
    Vector3f gravity{0.0f, -9.8f, 0.0f};
    Vector3f area{1.0f, 1.0f, 1.0f};
    int W, H, L, grid_size;
    float grid_interval = 0.02f;
};

struct Grid
{
    Vector3f velocity = Vector3f::Zero();
    float mass = 0;
    Vector3f force = Vector3f::Zero();
};

class Water
{
public:
    float E = 500.0f; // 杨氏模量
    float nu = 0.4f; // 泊松比
    float mass = 0.001f; // 质量
    float density = 1.0f; //密度

    // 通过上面四个属性计算以下四个属性值
    float mu = 0.5f * E / (1 + nu);
    float lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    float volume = mass / density; //体积
    float K = E / (1 - 2 * nu) / 3.0f;


};

struct Particle
{
    float mass = 0;
    Vector3f velocity = Vector3f::Zero();
    Matrix3f F = Matrix3f::Identity(); //单位矩阵，形变梯度
    Water* water = new Water();
    Matrix3f piola;

    Particle() {
        auto m = water;
        auto J = F.determinant(); // 行列式
        piola = m->mu * (F - F.transpose().inverse()) +
                m->lambda * log(J) * F.transpose().inverse();
    }
};