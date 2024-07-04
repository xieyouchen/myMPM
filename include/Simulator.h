//
// Created by 谢宥辰 on 2024/6/23.
//
#pragma once
#include "defines.h"

class Simulator
{
public:
    enum TransferScheme { FLIP99, FLIP95, APIC };
    SimulatorConfiguration simulator_configuration;
    Grid *grid;
    Particle *particles = nullptr;
    vector<int> active_nodes;
    int particles_size = 0;

    std::vector<Vector3f> positions;
    vector<Vector3f> velocity;
    TransferScheme transfer_scheme = TransferScheme::APIC;

    Simulator();
    Simulator(Vector3f area, float h);

    tuple<Vector3i, Matrix3f, Matrix3f> quadratic_interpolation(const Vector3f &particle_pos);

    inline Vector3f calc_quadratic(float o, float x);
    inline Vector3f calc_quadratic_grad(float o, float x);

    void add_object(const vector<Vector3f> &pos, const vector<Vector3f> v, const Water *water);

    void transfer_P2G();
    void add_gravity();
    void update_grid_force();
    void solve_grid_boundary();
    void update_grid_velocity();
    void update_F();
    void transfer_G2P();
    void advection();


    void substeps();
};
