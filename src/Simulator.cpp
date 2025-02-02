//
// Created by 谢宥辰 on 2024/6/23.
//

#include "Simulator.h"
#include "defines.h"
#include <Eigen/Core>
#include "math.h"

using namespace std;
using namespace Eigen;

Simulator::Simulator() {
    // 网格初始化
    simulator_configuration.W = simulator_configuration.area[0] / simulator_configuration.grid_interval + 1;
    simulator_configuration.H = simulator_configuration.area[1] / simulator_configuration.grid_interval + 1;
    simulator_configuration.L = simulator_configuration.area[2] / simulator_configuration.grid_interval + 1;
    simulator_configuration.grid_size = simulator_configuration.W * simulator_configuration.H * simulator_configuration.L;

    grid = new Grid[simulator_configuration.grid_size];
}

Simulator::Simulator(Vector3f area, float h)
{
    // 网格初始化
    simulator_configuration.grid_interval = h;
    simulator_configuration.W = area[0] / h + 1;
    simulator_configuration.H = area[1] / h + 1;
    simulator_configuration.L = area[2] / h + 1;
    simulator_configuration.grid_size = simulator_configuration.W * simulator_configuration.H * simulator_configuration.L;

    grid = new Grid[simulator_configuration.grid_size];
    particles = new Particle[simulator_configuration.grid_size];
}

tuple<Vector3i, Matrix3f, Matrix3f> Simulator::quadratic_interpolation(const Vector3f &particle_posistion) {
    // 得到粒子在网格中的位置
    Vector3f particle_pos = particle_posistion / simulator_configuration.grid_interval;

    Vector3i base_node = floor(particle_pos.array() - 0.5f).cast<int>();
    Matrix3f wp;
    wp << calc_quadratic(base_node(0), particle_pos(0)),
            calc_quadratic(base_node(1), particle_pos(1)),
            calc_quadratic(base_node(2), particle_pos(2));
    Matrix3f dwp;
    dwp << calc_quadratic_grad(base_node(0), particle_pos(0)),
            calc_quadratic_grad(base_node(1), particle_pos(1)),
            calc_quadratic_grad(base_node(2), particle_pos(2));

    return {base_node, wp, dwp};
}

inline Vector3f Simulator::calc_quadratic(float o, float x) {
    // +-(o)------(o+1)--(x)--(o+2)-+
    float d0 = x - o;
    float d1 = d0 - 1;
    float d2 = 1 - d1;

    return {0.5f * (1.5f - d0) * (1.5f - d0), 0.75f - d1 * d1,
            0.5f * (1.5f - d2) * (1.5f - d2)};
}

inline Vector3f Simulator::calc_quadratic_grad(float o, float x) {
    float d0 = x - o;
    float d1 = d0 - 1;
    float d2 = 1 - d1;

    return {d0 - 1.5f, -2 * d1, 1.5f - d2};
}

void Simulator::add_object(const vector<Vector3f> &pos, const vector<Vector3f> v, const Water *water) {
    positions.insert(positions.end(), pos.begin(), pos.end());
    velocity.insert(velocity.end(), v.begin(), v.end());

    int new_size = particles_size + positions.size();
    if(particles) {
        Particle* new_part = new Particle[new_size];
        std::memcpy(new_part, particles, sizeof(Particle) * particles_size); // 清空原particles指针指向的所有内容
        delete[] particles; // 释放内存

        particles = new_part;
    }
    else {
        particles = new Particle[new_size];
    }

    for(auto i = particles_size ; i < new_size ; i++) {
        particles[i].mass = water->mass;
        particles[i].velocity = v[i];
    }
    particles_size = new_size;
}

void Simulator::transfer_P2G() {
    auto pos = positions;
    for(int iter = 0 ; iter < this->particles_size ; iter++) {
        // 得到3×3网格的左下角网格下标
        auto [base_node, wp, dwp] =
                quadratic_interpolation(positions[iter]);

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

                    MPM_ASSERT(0 <= index && index < simulator_configuration.grid_size,
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

//    printf("size: %d ---> %d %f\n", active_nodes.size(), grid[active_nodes[0]].mass > 1e-6, grid[active_nodes[0]].mass);
//    for(int i = 0 ; i < active_nodes.size() ; i++) printf("%d ", active_nodes[i]);

}

void Simulator::add_gravity() {
    // 给每个网格节点的外力添加重力效果
    for(int i = 0 ; i < active_nodes.size() ; i++) {
        int idx = active_nodes[i];
        grid[idx].force += simulator_configuration.gravity * grid[idx].mass;
    }
}

void Simulator::update_grid_force() {
    for(int iter = 0 ; iter < particles_size ; iter++) {
        // 变形梯度 F 开始计算
        auto F = particles[iter].F;
        auto volume = particles[iter].water->volume;

        // 多线程中可能会出现问题，计算 piola 的时候
        particles[iter].cal_stress_tensor();
        Matrix3f piola = particles[iter].piola;
        auto [base_node, wp, dwp] = quadratic_interpolation(positions[iter]);

        for(int i = 0 ; i < 3 ; i++) {
            for(int j = 0 ; j < 3 ; j++) {
                for(int k = 0 ; k < 3 ; k++) {
                    Vector3i node = base_node + Vector3i(i, j, k);
                    // 梯度权重系数
                    auto h = simulator_configuration.grid_interval;
                    Vector3f grad_wip{dwp(i, 0) * wp(j, 1) * wp(k, 2) / h,
                                      wp(i, 0) * dwp(j, 1) * wp(k, 2) / h,
                                      wp(i, 0) * wp(j, 1) * dwp(k, 2) / h};

                    // 网格节点转化为一维数组的粒子下标
                    int index = node(0) * simulator_configuration.H * simulator_configuration.L +
                                node(1) * simulator_configuration.L +
                                node(2);

                    MPM_ASSERT(0 <= index && index < simulator_configuration.grid_size,
                               "PARTICLE OUT OF GRID");

                    grid[index].force -= volume * (piola * F.transpose()) * grad_wip;
//                    std::cout << "index: " << index << " --> " <<  grid[index].force.transpose() << std::endl;
//                    IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
//                    std::cout << "piola:" << piola.format(CleanFmt) << std::endl;
                }
            }
        }
    }
}

void Simulator::update_grid_velocity() {
    // 只需要处理有效节点
    for(int iter = 0 ; iter < active_nodes.size() ; iter++) {
        int idx = active_nodes[iter];

//        cout << "idx: " << idx << "---> v:" << grid[idx].velocity(0) << endl;
//        assert(isnormal(grid[idx].velocity(0)) == 1);
        grid[idx].velocity += simulator_configuration.dt * grid[idx].force / grid[idx].mass;
//        cout << "idx: " << idx << "---> v:" << grid[idx].velocity.transpose() << "   --> mass:" << grid[idx].mass
//        << "  ---> f:" << grid[idx].force.transpose() << endl;
    }
}

void Simulator::solve_grid_boundary() {
  // MPM_PROFILE_FUNCTION();
  // Sticky boundary
  int thickness = simulator_configuration.thickness;
  auto [W, H, L] = std::tie(simulator_configuration.W, simulator_configuration.H, simulator_configuration.L);
  // check x-axis bound
  for (int i = 0; i < thickness; i++) {
    for (int j = 0; j < H; j++) {
      for (int k = 0; k < L; k++) {
        int index1 = i * H * L + j * L + k;
        int index2 = (W - i - 1) * H * L + j * L + k;
        if (grid[index1].velocity[0] < 0) {
          grid[index1].velocity[0] = 0.0f;
        }
        if (grid[index2].velocity[0] > 0) {
          grid[index2].velocity[0] = 0.0f;
        }
      }
    }
  }
  // check y-axis bound
  for (int i = 0; i < W; i++) {
    for (int j = 0; j < thickness; j++) {
      for (int k = 0; k < L; k++) {
        int index1 = i * H * L + j * L + k;
        int index2 = i * H * L + (H - j - 1) * L + k;
        if (grid[index1].velocity[1] < 0) {
          grid[index1].velocity[1] = 0.0f;
        }
        if (grid[index2].velocity[1] > 0) {
          grid[index2].velocity[1] = 0.0f;
        }
      }
    }
  }
  // check z-axis bound
  for (int i = 0; i < W; i++) {
    for (int j = 0; j < H; j++) {
      for (int k = 0; k < thickness; k++) {
        int index1 = i * H * L + j * L + k;
        int index2 = i * H * L + j * L + (L - k - 1);
        if (grid[index1].velocity[2] < 0) {
          grid[index1].velocity[2] = 0.0f;
        }
        if (grid[index2].velocity[2] > 0) {
          grid[index2].velocity[2] = 0.0f;
        }
      }
    }
  }
}

void Simulator::update_F() {
    for(int iter = 0 ; iter < particles_size ; iter++) {
        auto F = particles[iter].F;
        auto [base_node, wp, dwp] =
                quadratic_interpolation(positions[iter]);

        // weight 就是仿射速度场，用于计算形变梯度 F 最关键的一个指标 Cp
        Matrix3f weight = Matrix3f::Zero();
        for(int i = 0 ; i < 3 ; i++) {
            for(int j = 0 ; j < 3 ; j++) {
                for(int k = 0 ; k < 3 ; k++) {
                    Vector3i node = base_node + Vector3i(i, j, k);
                    auto h = simulator_configuration.grid_interval;
                    Vector3f grad_wip{dwp(i, 0) * wp(j, 1) * wp(k, 2) / h,
                                      wp(i, 0) * dwp(j, 1) * wp(k, 2) / h,
                                      wp(i, 0) * wp(j, 1) * dwp(k, 2) / h};
                }
            }
        }
        particles[iter].F += simulator_configuration.dt * weight * F;
        if (particles[iter].F.determinant() < 0) {
            spdlog::error("particles[{}]'s determinat(F) is negative!\n{}", iter);
            assert(false);
        }
    }
}

void Simulator::transfer_G2P() {
    for(int iter = 0 ; iter < particles_size ; iter++) {
        auto [base_node, wp, dwp] = quadratic_interpolation(positions[iter]);

        Vector3f velocity_PIC = Vector3f::Zero();
        for(int i = 0 ; i < 3 ; i++) {
            for(int j = 0 ; j < 3 ; j++) {
                for(int k = 0 ; k < 3 ; k++) {
                  Vector3i node = base_node + Vector3i(i, j, k);
                  auto wijk = wp(i,0) * wp(j,1) * wp(k,2);
                  auto index = node(0) * simulator_configuration.H * simulator_configuration.L +
                          node(1) * simulator_configuration.L + node(2);
                  MPM_ASSERT(0 <= index && index < simulator_configuration.grid_size, "PARTICLE OUT OF GRID");

                  velocity_PIC += wijk * grid[index].velocity; // APIC 就是最简单的系数 * 速度
                }
            }
        }
        switch(transfer_scheme) {
          case TransferScheme::APIC:
            particles[iter].velocity = velocity_PIC;
            break;
        }
    }
}

void Simulator::advection() {
  for(int iter = 0 ; iter < particles_size ; iter++) {
    positions[iter] += simulator_configuration.dt * particles[iter].velocity;
  }
}

void Simulator::substeps() {
    transfer_P2G();
    add_gravity();
    update_grid_force();
    update_grid_velocity();
    // 处理网格的边界问题，后面再写，看看就这样出来的效果
    update_F(); // 更新形变梯度
    transfer_G2P();
    advection(); //平流，就是更新粒子的位置
    simulator_configuration.current_step++;
}