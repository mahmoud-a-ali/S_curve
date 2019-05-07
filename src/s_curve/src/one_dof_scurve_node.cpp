/**
\file   one_dof_scurve_node.cpp
\brief  to generate scurve trajectory for one dof
 *
\author  Mahmoud Ali
\date   25/4/2019
*/



#include "std_msgs/String.h"
#include<vector>
#include<string>
#include "algorithm"
#include<math.h>

#include<iostream>
#include "ros/ros.h"
#include<yaml-cpp/yaml.h>
#include<trajectory_msgs/JointTrajectory.h>
#include<trajectory_msgs/JointTrajectoryPoint.h>
#include<std_msgs/Header.h>
 // different trajectory to test different case
//#include<normal_toppra_traj_instant_1.h>
#include<normal_toppra_traj_instant_2.h>

#include<scurve_one_dof_coefs_and_sampling.cpp>


#include <python2.7/Python.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
const double jm = 985,   vm = 130.00,  am= 250.0/4, sm =185.00;
const double eps= 0.0001;

int main(int argc, char **argv)
{
    ros::init(argc, argv, "one_dof_scurve_node");
    ros::NodeHandle nh;
    ROS_INFO(" test: one_dof_scurve_node ...  ");

    //============read the trajectory points generated by python script ==========
    trajectory_msgs::JointTrajectory  traj;
    traj = generate_traj();
    int n_jts = traj.joint_names.size(); //number of dof
    int n_pts = traj.points.size();      // number of waypoints
    std::cout<< "number of joints, and points: "<< n_jts << ", "<<n_pts << std::endl;


    // waypoints for specific joint
    int jt=4; // here we chose a path of one joint as it is still one_dof
    std::vector<double> P_wpt;    //vector of waypoints
    for(int pt=0; pt<n_pts; pt++){
        P_wpt.push_back( traj.points[pt].positions[jt] );
    }

    std::vector< std::vector<double> >  traj_T;
    traj_T.resize(4); // Tj, Ta, Tv, T
    std::vector<double> jrk, seg_idx;

    double ref_T = 2;
    //  compute the coef of scurve
    // if ref_T=0, it will continue with minimum optimal time ( max jerk)
    compute_1dof_scurve_coef(P_wpt, ref_T, traj_T, jrk, seg_idx, sm, vm, am, jm );
    double T_tot =0;
    for (auto& t : traj_T[3]){
//        ROS_INFO_STREAM( "t: " << t );
        T_tot += t;
    }
    ROS_INFO_STREAM( "time: " << T_tot );
    ROS_INFO_STREAM( "jrk: " << jrk[0] );

    // use sampling function
    std::vector<double>  TPVA, Ts, Ps, Vs, As;
    TPVA.resize(4);

    for (auto& p : seg_idx) // to check at which point the direction has been changed
        ROS_INFO_STREAM( "change of direction occured in: " << p );

    double tg = 0, ts =0.008;
    while(tg < T_tot){
        // second function to sample the trajectory
        bool check_sample =  sample_1dof_scurve( tg,  traj_T, TPVA,  jrk, seg_idx);
        Ts.push_back(tg);
        Ps.push_back(TPVA[1]);
        Vs.push_back(TPVA[2]);
        As.push_back(TPVA[3]);
//        ROS_INFO_STREAM(  "  t: " << TPVA[0] << "    P: "<< TPVA[1] << "   V: "<< TPVA[2] << "   A: "<< TPVA[3] );
        tg+=ts;
    }

    plt::plot(Ts, Ps);
    plt::xlabel("time");
    plt::ylabel("pos");
    plt::title("position_one_dof");
    plt::legend();
    plt::show();

    plt::plot(Ts, Vs);
    plt::xlabel("time");
    plt::ylabel("vel");
    plt::title("velocity_one_dof");
    plt::legend();
    plt::show();

    plt::plot(Ts, As);
    plt::xlabel("time");
    plt::ylabel("acc");
    plt::title("acceleration_one_dof");
    plt::legend();
    plt::show();




return 0;


}
