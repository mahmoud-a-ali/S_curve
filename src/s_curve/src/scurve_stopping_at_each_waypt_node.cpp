/**
\file   scurve_stopping_at_each_waypt_node.cpp
\brief  to scurve trajectory for ndof.
 * this implementation of scurve considers that the velocity will be null at each waypoint.
 * all the n dof are synchronized in motion.
 * timing for each points could be a ref_time passed to the compute_ndof_scurve_Coefs_per_waypt function, 
 * otherwise the optimal minimum time -for indicidual waypoints- will be computed.
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
// different trajectory to test different cases
//#include<traj_test_instant.h>
//#include<normal_toppra_traj_instant_1.h>
//#include<normal_toppra_traj_instant_2.h>
#include<normal_toppra_traj_instant_3.h>

#include<scurve_ndof_stopping_waypts_coefs_and_sampling.cpp>


#include <python2.7/Python.h>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
const double jm = 985,   vm = 130.00,  am= 250.0/4, sm =185.00;
const double eps= 0.0001;



int main(int argc, char **argv)
{
    ros::init(argc, argv, "scurve_stopping_at_each_waypt_node");
    ros::NodeHandle nh;
    ROS_INFO(" test: scurve_stopping_at_each_waypt_node ...  ");


    //============define traj and read the one generated by python script==========
    trajectory_msgs::JointTrajectory  traj;
    traj = generate_traj();
    int n_jts = traj.joint_names.size();
    int n_pts = traj.points.size();
    ROS_INFO_STREAM( "number of joints: "<< n_jts << "  ,   number of waypoints: "<<n_pts );

    // waypoints of trajectory
    std::vector< std::vector<double> > P_jt_wpt;
     P_jt_wpt.resize(n_jts);
     for(int jt=0; jt<n_jts; jt++){
         for(int pt=0; pt<n_pts; pt++){
             P_jt_wpt[jt].push_back( traj.points[pt].positions[jt] );
         }
     }



    // ref_T: reference time for each waypoint (if we prefere a ref_T higher than optimal time)
    // t_start: start_time of trajectory, tg: global time variable, ts: sampling time
     double ref_T =0, t_start =0, tg= 0, ts=0.008;


    // vectors for Pos, Vel, Acc are using to store output trajectory pos, vel, acc for plotting purpose
    std::vector< std::vector<double>> T_vec, P_vec, V_vec, A_vec;
    T_vec.resize( n_jts);
    V_vec.resize( n_jts);
    A_vec.resize( n_jts);
    P_vec.resize( n_jts);

    //vector for Pos, Vel, Acc for eac joint at each sample time
    std::vector< std::vector<double> > TPVA_jt;
    TPVA_jt.resize(n_jts);
    for(int jt=0; jt<n_jts; jt++)
        TPVA_jt[jt].resize(4);

    // scurve coeffients, Tj: time with max jerk, Ta: time with max acc, Tv: time wit max vel,
    // Jrk: vector contains values of jerk for each joint to keep synchrounization
    double Ta=0, Tv=0, Tj=0;
    std::vector<double> Jrk;
    Jrk.resize(n_jts);


    // for each segment, it is represented by P_init and P_goal for each joint
     std::vector<double> P_jt_init, P_jt_goal;
     P_jt_init.resize(n_jts);
     P_jt_goal.resize(n_jts);

    // for each point/segement in trajectory between two successful points
    for (int pt=1; pt<n_pts; pt++) { //n_pts
//        std::cout<< " pt: "<< pt<<std::endl<<std::endl<<std::endl<<std::endl;
        // initialize P_init and P_goal for all joints
        for (int jt=0; jt<n_jts; jt++) {
             P_jt_init[jt] = P_jt_wpt[jt][pt-1];
             P_jt_goal[jt] = P_jt_wpt[jt][pt];
         }
//------------------------ print initial and goal waypoints-------------------------
         std::cout <<std::endl << "===============seg_" << pt<<" info ==================="<<std::endl;
//         std::cout<< " ref_T,  tg,  t_start: "<< ref_T<<"   "  <<  tg<<"   "<<  t_start <<std::endl;
         std::cout<<"P_jt_init: ";
         for(auto &p : P_jt_init)
             std::cout<<"  "<< p<< "  ";
         std::cout<< std::endl;
         std::cout<<"P_jt_goal: ";
         for(auto &p : P_jt_goal)
             std::cout<<"  "<< p<< "  ";
         std::cout<< std::endl;
//-------------------------------------------------
         // call function to compute jerk for all joits, taking into account 1.synchrounization and 2.ref_T (if exist)
         // when ref_T equals zero, then it calculate min optimal time
         compute_ndof_scurve_Coefs_per_waypt (P_jt_init,  P_jt_goal,  ref_T, Tj, Ta, Tv, Jrk,sm, vm, am, jm );

 //------------------------ print jerk for each joint to have synchrounization and ref_T -------------------------
         std::cout<<"jerk_for_all_joint: ";
         for(auto &j : Jrk)
             std::cout<<"  "<< j<< "  ";
         std::cout<< std::endl;
         std::cout<<"time_for_segment: " << ref_T<<  std::endl;
//-------------------------------------------------
        t_start=tg; // start of trajectory (first segement)
        // for each sample time insie the ref_T of the coressponding segment, calculate Pos, Vel, Acc using sample function
        while ( tg-t_start <= ref_T) {
            bool check = sample_scurve_ndof(P_jt_init,  Jrk, Ta, Tv,  Tj,  tg,  t_start,   TPVA_jt);
//            std::cout<< "check: " << check << std::endl;
            // store Pos, Vel, Acc for plotting later
             for (int jt =0; jt<n_jts; jt++) {
                T_vec[jt].push_back(TPVA_jt[jt][0]);
                P_vec[jt].push_back(TPVA_jt[jt][1]);
                V_vec[jt].push_back(TPVA_jt[jt][2]);
                A_vec[jt].push_back(TPVA_jt[jt][3]);
                }
            tg+=ts; // next sample
        }


        std::string jt_name;
        //for each joint plot Pos, Vel, Acc
        for (int jt =0; jt<n_jts; jt++){
            jt_name= "jt_" + std::to_string(jt);
             plt::named_plot(jt_name,T_vec[jt], P_vec[jt]);
        }
        plt::xlabel("time");
        plt::ylabel("pos");
        plt::title("position_scurve_0_vel_at_waypts");
        plt::legend();
        plt::show();


        for (int jt =0; jt<n_jts; jt++){
            jt_name= "jt_" + std::to_string(jt);
             plt::named_plot(jt_name, T_vec[jt], V_vec[jt]);
        }
        plt::xlabel("time");
        plt::ylabel("vel");
        plt::title("position_scurve_0_vel_at_waypts");
        plt::legend();
        plt::show();


        for (int jt =0; jt<n_jts; jt++){
            jt_name= "jt_" + std::to_string(jt);
             plt::named_plot(jt_name, T_vec[jt], A_vec[jt]);
        }
        plt::xlabel("time");
        plt::ylabel("acc");
        plt::title("position_scurve_0_vel_at_waypts");
        plt::legend();
        plt::show();



//        for (int jt =0; jt<n_jts; jt++)
//            plt::plot(T_vec[jt], V_vec[jt]);
//        plt::show();
//        for (int jt =0; jt<n_jts; jt++)
//            plt::plot(T_vec[jt], A_vec[jt]);
//        plt::show();


     }

}
