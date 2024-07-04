#include <iostream>
#include "utils.h"
#include "Simulator.h"

using namespace std;
using namespace Eigen;


int main() {
    // y 方向作为竖直方向，初始速度设置
    Vector3f velocity{-2.5f, 0.5f, -0.3f};
    //  选择本构模型，流体2

    // 选择平流模型，FLIP99
    float alpha = 0.99f; // flip/pic

    // 选择材质为水
    Water* water = new Water();
    cout << "材质初始化完成，所选材质属性E：" << water->E << endl;

    // 初始化一个模拟器
    Simulator* simulator = new Simulator();
    // 得到物体的初始位置
    vector<Vector3f> positions;
    string model_path = "/Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/small_cube.obj";
    read_model_positions(model_path, positions);

    // 尝试将初始位置写入 .bgeo 文件
    string output_dir("/Users/xieyouchen/usually/Projects/LiuLong/Project/myMPM/output/test/");
    write_particles_position(output_dir+"0.bgeo", positions);

    // 帧渲染
    int total_frame = 10;
    for(int frame = 0 ; frame < total_frame ; frame++) {
      simulator->simulator_configuration.current_step = 0;
//        MPM_INFO("Start to Frame Compute");
      if(frame % 50 == 0) {
          // 计算到第 50 帧，粒子流动趋于稳定，再次增加一个物体进行渲染展示
          vector<Vector3f> v =  vector<Vector3f>(positions.size(), velocity);
          simulator->add_object(positions, v, water);
//            MPM_INFO("Add_Object() SUCCESS");
      }

      // 帧率60
      int frame_rate = 60;
      float dt = 1e-4f;
      int steps_per_frame = (int)ceil(1.0f / frame_rate / dt);
      for(int i = 0 ; i < steps_per_frame ; i++) {
          // 一帧要计算 160 次的，每一帧就是一个 bgeo 文件
          simulator->substeps();
      }

      string output = output_dir + to_string(frame+1) + ".bgeo";
      write_particles_position(output, positions);
    }

    return 0;
}
