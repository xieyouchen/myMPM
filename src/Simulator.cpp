//
// Created by 谢宥辰 on 2024/6/23.
//

#include "Simulator.h"
#include "defines.h"
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

Simulator::Simulator(Vector3f area, float h)
{
    // 网格初始化
    grid_interval = h;
    simulator_configuration.W = area[0] / grid_interval + 1;
    simulator_configuration.H = area[1] / grid_interval + 1;
    simulator_configuration.L = area[2] / grid_interval + 1;
    simulator_configuration.grid_size = simulator_configuration.W * simulator_configuration.H * simulator_configuration.L;

    grid = new Grid[simulator_configuration.grid_size];
    particles = new Particle[simulator_configuration.grid_size];
}

tuple<Vector3i, Matrix3f> Simulator::quadratic_interpolation(const Vector3f &particle_pos) {
    Vector3i base_node = floor(particle_pos.array() - 0.5f).cast<int>();
    Matrix3f wp;
    wp << calc_quadratic(base_node(0), particle_pos(0)),
            calc_quadratic(base_node(1), particle_pos(1)),
            calc_quadratic(base_node(2), particle_pos(2));

    return {base_node, wp};
}

inline Vector3f Simulator::calc_quadratic(float o, float x) {
    // +-(o)------(o+1)--(x)--(o+2)-+
    float d0 = x - o;
    float d1 = d0 - 1;
    float d2 = 1 - d1;

    return {0.5f * (1.5f - d0) * (1.5f - d0), 0.75f - d1 * d1,
            0.5f * (1.5f - d2) * (1.5f - d2)};
}

void Simulator::add_object(const vector<Vector3f> &pos, const vector<Vector3f> v, const Water *water) {
    positions.insert(positions.end(), pos.begin(), pos.end());
    velocity.insert(velocity.end(), v.begin(), v.end());
    particles_size += positions.size();

    for(auto i = 0 ; i < particles_size ; i++) {
        particles[i].mass = water->mass;
        particles[i].velocity = v[i];
    }
}

void Simulator::transfer_P2G() {
    for(int iter = 0 ; iter < this->particles_size ; iter++) {
        // 得到粒子在网格中的位置
        Vector3f particle_pos_grid = positions[iter] / grid_interval;
        // 得到3×3网格的左下角网格下标
        auto [base_node, wp] =
                quadratic_interpolation(particle_pos_grid);

        auto p = particles[iter];
        // 网格节点循环
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    // 遍历每个网格节点进行计算
                    Vector3i node = base_node + Vector3i(i, j, k);
                    // 网格节点转化为一维数组的粒子下标
                    int index = node(0) * simulator_configuration.H * simulator_configuration.L +
                                node(1) * simulator_configuration.L +
                                node(2);

                    MPM_ASSERT(index >= 0 && index < simulator_configuration.grid_size,
                               "PARTICLE OUT OF GRID at Transfer_P2G");

                    // 权重系数累积汇总
                    float wijk = wp(i, 0) * wp(j, 1) * wp(k, 2);
                    grid[index].velocity += wijk * p.mass * p.velocity;
                    grid[index].mass += wijk * p.mass;
                }
            }
        }
    }
    // 如果网格节点的质量不够，那么直接无视
    for(int i = 0 ; i < simulator_configuration.grid_size ; i++) {
        if(grid[i].mass > 1e-15) {
            active_nodes.push_back(i);
            grid[i].velocity /= grid[i].mass;
        }
        else grid[i].velocity = Vector3f::Zero();
    }
}

void Simulator::substeps(float dt) {
    transfer_P2G();
    // 需要检测一个函数的输出结果，看看是不是正常的
}