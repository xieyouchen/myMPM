//
// Created by 谢宥辰 on 2024/6/23.
//
#pragma once
#include "defines.h"

class Simulator
{
public:
    SimulatorConfiguration sim_config;
    int particles_size = 0;

    std::vector<Vector3f> positions;
    vector<Vector3f> velocity;
    float grid_interval = 0;

    Simulator(Vector3f area, float h);

    tuple<Vector3i, Matrix3f> quadratic_interpolation(const Vector3f &particle_pos);

    inline Vector3f calc_quadratic(float o, float x);

    void add_object(const vector<Vector3f> &pos, const vector<Vector3f> v, const Water *water);

    void transfer_P2G();

    void substeps(float dt);
};
