//
// Created by 谢宥辰 on 2024/6/22.
//

#include "utils.h"
#include "Partio.h"
#include "defines.h"

bool write_particles_position(const std::string &path, const vector<Vector3f> &positions) {
    Partio::ParticlesDataMutable *particles = Partio::create();
    // 添加粒子属性
    Partio::ParticleAttribute position_atr = particles->addAttribute("position", Partio::VECTOR, 3);
    Partio::ParticleAttribute index_atr = particles->addAttribute("index", Partio::INT, 1);

    // 逐个添加粒子
    for (int i = 0; i < positions.size(); i++) {
        // 添加粒子，并获得该粒子的下标值
        Partio::ParticleIndex index_particle = particles->addParticle();
        // 获得初始化后粒子的属性指针
        auto *position = particles->dataWrite<Vector3f>(position_atr, index_particle);
        auto *index = particles->dataWrite<int>(index_atr, index_particle);

        *position = positions[i];
        *index = i;
    }
    write(path.c_str(), *particles);
    particles->release();
    cout << "- Write Particles Positions SUCCESS \n  To" << path << endl;
    return true;
}

bool read_model_positions(const string &model_path, vector<Vector3f> &positions) {
    ifstream input(model_path);
    string line;
    Vector3f pos;
    if(input) {
        while(getline(input, line)) {
            if(line[0] == 'v') {
                sscanf(line.c_str(), "v %f %f %f", &pos[0], &pos[1], &pos[2]);
                positions.push_back(pos);
            }
        }
        cout << "- Read Particles Positions SUCCESS \n  The Particles Num is " << positions.size() << endl;
        return true;
    }
    else {
        cout << "ERROR" << endl;
        return false;
    }

}




