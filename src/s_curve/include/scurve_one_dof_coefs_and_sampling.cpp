#include"scurve_equations.cpp"
//============================== one_dof_scurve_coef ===========================
/* the main function that will be called to compute the scurve coeffients for specific trajectory with specific max limits of Pos, Vel, Acc, Jrk
 input argument:
    P_wpt: waypoints, ref_T: reference time for the trajectory, sm: max_Pos, vm: max_vel, am: max_acc, jm: max_jerk
 output_arg:
    jrk:  jrk per each segment in the path (whenever the direction changes we consider anew segment in the path)
    traj_T: vector of times of the inflection points for each segment in scurve the scurve
*/
void compute_1dof_scurve_coef(  std::vector<double> P_wpt, double ref_T,  std::vector< std::vector<double> > &traj_T, std::vector<double> & jrk,
                           std::vector<double> &seg_idx, double sm, double vm, double am, double jm){
//    std::cout<< " ====== one_dof_scurve_coef ============ "<< std::endl;
    int n_pts = P_wpt.size(); // number of waypoints
//calculate n_traj per joint, when the direction is reversed, a new segment is considered
    std::vector<double> trajs_idx;
    trajs_idx.push_back(0);
    int n_trajs=1;
    for(int pt=0; pt<n_pts-1; pt++){
        if( P_wpt[pt]*P_wpt[pt+1] < 0 ){
            n_trajs++;
            trajs_idx.push_back(pt);
        }
    }
    trajs_idx.push_back(n_pts-1);
    for(int seg=0; seg<trajs_idx.size() ; seg++)
        seg_idx.push_back( P_wpt[trajs_idx[seg]] );

//    std::cout<< "number of traj_segments (direction_change): "<< n_trajs << "  index wpts: " << std::endl;

    // variables
    traj_T[0].resize(n_trajs); // Tj
    traj_T[1].resize(n_trajs); // Ta
    traj_T[2].resize(n_trajs); // Tv
    traj_T[3].resize(n_trajs); // T
    jrk.resize(n_trajs); // Jrk

    std::vector<double>  Ds;
    std::vector<int> Ds_sgn;
    Ds.resize(n_trajs);
    Ds_sgn.resize(n_trajs);

// for each seg or change in dir
    for (int trj=0; trj< n_trajs; trj++){
        double p0 = P_wpt[ trajs_idx[trj] ];
        double pf = P_wpt[ trajs_idx[trj+1] ];
        Ds[trj]= pf - p0;
        Ds_sgn[trj] = (Ds[trj] > 0) ? 1 : ((Ds[trj] < 0) ? -1 : 0);

         if(fabs(Ds[trj])<=1e-5){
             Ds[trj] = 0;
             Ds_sgn[trj] = 0;
             traj_T[0][trj]= 0;
             traj_T[1][trj]= 0;
             traj_T[2][trj]= 0;
             traj_T[3][trj]= 0;
             jrk[trj]= 0;
//              std::cout<< "Ds["<< trj <<"] is zero ";
             continue;
         }
//          std::cout<< "Ds["<< trj <<"]= " << Ds[trj];
         double t = compute_time_for_jrk(fabs(Ds[trj]),traj_T[0][trj], traj_T[1][trj], traj_T[2][trj], sm, vm, am, jm );
         traj_T[3][trj] = t;
         jrk[trj] = Ds_sgn[trj]*jm;
     }

    if(ref_T > 0){
        double T_opt = 0;
        for (auto& t : traj_T[3]) //total optimal time for all phases (max_jek, max_vel, max_acc)
            T_opt += t;
        if( ref_T < T_opt)
            throw(std::invalid_argument("reference time is less than optimal time"));
        std::vector<double> DT;
        DT.resize(n_trajs);
        double  dt = ref_T - T_opt;
//        std::cout<< "ref_T, T_opt: "<< ref_T << "  "<< T_opt <<std::endl;

        for (int trj=0; trj< n_trajs; trj++){
          DT[trj] = traj_T[3][trj]*dt / T_opt;
//          std::cout<< "T, DT, Ds, Tj,Ta,Tv: "<< traj_T[3][trj]<<"  "<<  DT[trj]<<"  "<< Ds[trj]<<"  "<< traj_T[0][trj]<<"  "<< traj_T[1][trj]<<"  "<< traj_T[2][trj]<< std::endl;
          jrk[trj]= Ds_sgn[trj]*compute_jerk_for_time( traj_T[3][trj],  DT[trj], fabs(Ds[trj]), traj_T[0][trj], traj_T[1][trj], traj_T[2][trj], sm, vm, am, jm );
//          std::cout<< "new_jrk: "<< jrk[trj] << std::endl << std::endl ;
        }

    }
}



//========================= Sample function ==========================
/* to make this implementation similar to the joint_trajectory_controller
    it takes as input the coef computed by the one_dof_scurve_coef function
    and calculate vel, acc, jrk at each sample, for each sample we get vector TPVA has time, vel, acc, jrk at this sample
*/

bool sample_1dof_scurve( double tg,  std::vector<std::vector<double>> &traj_T, std::vector<double> &TPVA, std::vector<double> jrk, std::vector<double> seg_idx){
        // which segment
//    std::cout<< " ====== sample_scurve ============ "<< std::endl;
        int n_trajs = jrk.size();
        double T_tot = 0;
        for (auto& t : traj_T[3])
            T_tot += t;

        if( tg<0 || tg> T_tot )
            return false;

        std::vector<double> t_idx;
        t_idx.resize(n_trajs +1);
        t_idx[0]=0;
        double ts=0;
        for (int trj=0; trj<n_trajs; trj++) {
            ts += traj_T[3][trj];
            t_idx[trj+1]=ts;
        }

        int traj_idx =0;
        for (int trj=0; trj<n_trajs; trj++) {
            if( tg> t_idx[trj] && tg<= t_idx[trj+1] )
                traj_idx = trj;
        }

//        for(auto& p: t_idx)
//            std::cout<< p <<std::endl;

        double t = tg-t_idx[traj_idx];
//        std::cout<< "tg:  " << tg <<"  t: " << t <<"traj_idx   "<< traj_idx<<std::endl;
        double P0=0, V0=0, A0=0;
        // update ipts time
        std::vector<double> T_ipt;
        T_ipt.resize(8);
        update_ip_time(T_ipt, traj_T[0][traj_idx], traj_T[1][traj_idx], traj_T[2][traj_idx]);
//        std::cout<< "time was updated "<<std::endl;
//        for(auto& t: T_ipt)
//            std::cout<< t <<std::endl;
        // second update the ipt_time pos vel acc
        std::vector<double>  P_ipt, V_ipt, A_ipt;
        P_ipt.resize(8);
        V_ipt.resize(8);
        A_ipt.resize(8);
        update_ip_pos_vel_acc(P_ipt, V_ipt, A_ipt, T_ipt, jrk[traj_idx], P0, V0, A0);
//        for(auto& p: P_ipt)
//            std::cout<< p <<std::endl;
        // which phase this traj segment located
//        std::cout<< "PVA was updated "<<std::endl;
//        std::cout<< "t: "<<t<<  "   jrk: "<< jrk[traj_idx]<<std::endl;
        double p=0, v=0, a=0;
        if( t>=T_ipt[0] && t<T_ipt[1])
            phase_j1 (  traj_T[0][traj_idx], traj_T[1][traj_idx], traj_T[2][traj_idx],  P0,  V0,  A0,  jrk[traj_idx],  t,  a,  v,  p);
        else if(t<T_ipt[2])
            phase_a1 (  traj_T[0][traj_idx], traj_T[1][traj_idx], traj_T[2][traj_idx],  P0,  V0,  A0,  jrk[traj_idx],  t,  a,  v,  p);
        else if(t<T_ipt[3])
            phase_j2 (  traj_T[0][traj_idx], traj_T[1][traj_idx], traj_T[2][traj_idx],  P0,  V0,  A0,  jrk[traj_idx],  t,  a,  v,  p);
        else if(t<T_ipt[4])
            phase_v (  traj_T[0][traj_idx], traj_T[1][traj_idx], traj_T[2][traj_idx],  P0,  V0,  A0,  jrk[traj_idx],  t,  a,  v,  p);
        else if(t<T_ipt[5])
            phase_j3 (  traj_T[0][traj_idx], traj_T[1][traj_idx], traj_T[2][traj_idx],  P0,  V0,  A0,  jrk[traj_idx],  t,  a,  v,  p);
        else if(t<T_ipt[6])
            phase_a2 (  traj_T[0][traj_idx], traj_T[1][traj_idx], traj_T[2][traj_idx],  P0,  V0,  A0,  jrk[traj_idx],  t,  a,  v,  p);
        else if(t<=T_ipt[7])
            phase_j4 (  traj_T[0][traj_idx], traj_T[1][traj_idx], traj_T[2][traj_idx],  P0,  V0,  A0,  jrk[traj_idx],  t,  a,  v,  p);
        else
            return false;

        TPVA[0]= t;
        TPVA[1]= p + seg_idx[traj_idx];
        TPVA[2]= v;
        TPVA[3]= a;
        return true;
}

