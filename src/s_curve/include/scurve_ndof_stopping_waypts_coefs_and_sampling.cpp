#include "algorithm"
#include"scurve_equations.cpp"





//============================== compute_ndof_scurve_Coefs_per_waypt ===========================
/* the main function that will be called to compute the scurve coeffients for specific trajectory with specific max limits of Pos, Vel, Acc, Jrk
 input argument:
    P_jt_init: start_state, P_jt_goal: goal waypoints
    ref_T: reference time for the trajectory, sm: max_Pos, vm: max_vel, am: max_acc, jm: max_jerk
 output_arg:
    jrk:  jrk per each joint in the path
    Tj, Ta, and Tv: represents the time of the scurve with max jrk, acc, vel respectivily
*/
void compute_ndof_scurve_Coefs_per_waypt(std::vector<double> P_jt_init, std::vector<double> P_jt_goal, double &ref_T,
                            double &Tj, double &Ta, double &Tv,std::vector<double> & Jrk,
                           double sm, double vm, double am, double jm ){

    double n_jts = P_jt_goal.size();
//    Jrk.resize(n_jts);
    double P0=0, V0=0, A0=0, ds=0, syn_T=0;
    int ref_jt =0, ds_sgn=0;
    std::vector <double> Ds;
    std::vector <int> Ds_sgn;
    Ds.resize(n_jts);
    Ds_sgn.resize(n_jts);
//    std::cout<< "number of joints, and points: "<< n_jts <<std::endl;

    std::vector<double> T_jt_cpt;
    T_jt_cpt.resize(n_jts);

    for (int jt = 0; jt <n_jts ; jt++){ //n_jts
//        std::cout<< "\n ================= "<<jt << "================= " <<std::endl;
        // assign Ds for each joint for current point
        P0 = P_jt_init[jt];// [pt-1];
        ds = P_jt_goal[jt] - P0;
        ds_sgn = (ds > 0) ? 1 : ((ds < 0) ? -1 : 0);
        Ds[jt] =ds;
        Ds_sgn[jt] =ds_sgn;
//        std::cout<<"ds= "<< Ds[jt] <<std::endl;
         if(fabs(Ds[jt])<=1e-4){
             Ds[jt] = 0;
             Ds_sgn[jt] = 0;
//             std::cout<<"jt= "<< jt <<"  Ds is zero .. "<<std::endl;
             continue;
         }

         //compute min optimal time required by each joint to reach goal P_jt_goal
         double T = compute_time_for_jrk(fabs(Ds[jt]),Tj, Ta, Tv, sm, vm, am, jm );
         T_jt_cpt[jt] = T;
         Jrk[jt] = Ds_sgn[jt]*jm;
//         std::cout<< "initial total time= "<< T  <<"   initial jrk= "<< Jrk[jt] <<std::endl;
//         std::cout<< "Tj, Ta, Tv = "<< Tj << ", "<< Ta << ", "<< Tv <<std::endl;

    }

        // find syn_T
        syn_T = *max_element( T_jt_cpt.begin(), T_jt_cpt.end() );
        // get the ref_joit that has max time and assign values Tj, TA, Tv accordingly
        for(int i=0; i< n_jts; i++){
            if(T_jt_cpt [i] == syn_T){
                ref_jt = i;
                // compute new jrk time for each point to keep synchrounization bet them.
                syn_T = compute_time_for_jrk(fabs(Ds[ref_jt]),Tj, Ta, Tv, sm, vm, am, jm );
                break;
            }
        }

        //if ref_T < opt_T then cout using min opt time
        double DT = 0;
        if(ref_T < syn_T){
            std::cout<<" reference time is less than optimal min time, optimal time =  "<< syn_T << std::endl;
            ref_T = syn_T;
        }else {
             DT = ref_T - syn_T;
        }


        // assign jerk for all  joints, ref_jt takes max
        for (int jt=0; jt< n_jts; jt++) {
//            std::cout<< "jrk = "<< Jrk[jt] << "   Tj,Ta,Tv= "<< Tj<< ", " << Ta<< ", " << Tv<<std::endl;
            if(Ds[jt] == 0){
                Jrk[jt] = 0;
                continue;
            }
            // to ensure that the last assigned values for Tj, Ta, tv are corresponding to the syn_T
            double jk = compute_jerk_for_time(syn_T, DT, fabs(Ds[jt]),Tj, Ta, Tv, sm, vm, am, jm);
            Jrk[jt] = Ds_sgn[jt]*jk;
//            std::cout<< "jrk = "<< Jrk[jt] << "   Tj,Ta,Tv= "<< Tj<< ", " << Ta<< ", " << Tv<<std::endl;
        }


}

//===================================== sapming_function ========================================
bool sample_scurve_ndof(std::vector<double> P0_jt, std::vector<double> Jrk,
                   double Ta, double Tv, double Tj, double tg, double t_start,
                   std::vector< std::vector<double> > &TPVA_jt){

     unsigned long n_jts = P0_jt.size();
     std::vector<double> T_jt_ipt;  // vector with the time of the inflection points (when change from one phase to another j1, a2, j3, ... etc)
     for(int pt=0; pt<8; pt++)  //initialization, it has 8 inflection point inlcuding starting point
         T_jt_ipt.push_back( 0 );

     for (unsigned long jt=0; jt<n_jts; jt++) {
//          std::cout<< "=================== sample joint"<< jt <<" ==============" <<std::endl;
        double A0=0, V0=0, P0=P0_jt[jt];
        update_ip_time(T_jt_ipt, Tj, Ta, Tv); //update time for inflection point
//        std::cout<<"T_jt_ipt: ";
//        for(auto &j : T_jt_ipt)
//            std::cout<<"  "<< j<< "  ";
//        std::cout<< std::endl;
        double t = tg - t_start;
//        std::cout<<"tg: "<<tg<< std::endl;

        // which phase (of the seven phases) should be considered
        double p=0, v=0, a=0;
        if( t>=T_jt_ipt[0] && t<T_jt_ipt[1])
            phase_j1 ( Tj,  Ta,  Tv, P0,  V0, A0, Jrk[jt], t, a, v, p);
        else if(t<T_jt_ipt[2])
            phase_a1 ( Tj,  Ta,  Tv, P0,  V0, A0, Jrk[jt], t, a, v, p);
        else if(t<T_jt_ipt[3])
            phase_j2 ( Tj,  Ta,  Tv, P0,  V0, A0, Jrk[jt], t, a, v, p);
        else if(t<T_jt_ipt[4])
            phase_v ( Tj,  Ta,  Tv, P0,  V0, A0, Jrk[jt], t, a, v, p);
        else if(t<T_jt_ipt[5])
            phase_j3 ( Tj,  Ta,  Tv, P0,  V0, A0, Jrk[jt], t, a, v, p);
        else if(t<T_jt_ipt[6])
            phase_a2 ( Tj,  Ta,  Tv, P0,  V0, A0, Jrk[jt], t, a, v, p);
        else if(t<=T_jt_ipt[7])
            phase_j4 ( Tj,  Ta,  Tv, P0,  V0, A0, Jrk[jt], t, a, v, p);
        else
            return false;
//        std::cout<< "p,v,a= "<< p << ", "<< v << ", "<< a << std::endl;
        //return Pos, Vel, Acc at this time instant
        TPVA_jt[jt][0]= t +t_start;
        TPVA_jt[jt][1]= p ;
        TPVA_jt[jt][2]= v;
        TPVA_jt[jt][3]= a;

        }
        return true;
}

