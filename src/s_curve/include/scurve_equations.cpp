#include<math.h>
#include <iostream>
#include<vector>

double min_root(double r1, double r2);
double min_root(double r1, double r2,double r3);
int  cubic_eq_real_root (double a, double b, double c, double d, std::vector<double> &roots);
// equations that describe the motion in the seven phases of the s_curve
void phase_j1 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p);
void phase_a1 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p);
void phase_j2 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p);
void phase_v (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p);
void phase_j3 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p);
void phase_a2 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p);
void phase_j4 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p);


//======================  update_ip_time: ======================
/*this func updates times of the 7 inflection points at which the trajectory changes from one phase to another one
 * phases are j1,j2,j3,j4: coressponding to max_jrk, a1,a2: coressponding to motion with max_acc,  v: coressponding to motion with max_vel
*/
void update_ip_time ( std::vector<double> &T_jt_ipt, double Tj, double Ta, double Tv ){
    T_jt_ipt[0] = 0;
    T_jt_ipt[1] = Tj;
    T_jt_ipt[2] = Tj+Ta;
    T_jt_ipt[3] = 2*Tj+Ta;
    T_jt_ipt[4] = 2*Tj+Ta+Tv;
    T_jt_ipt[5] = 3*Tj+Ta+Tv;
    T_jt_ipt[6] = 3*Tj+2*Ta+Tv;
    T_jt_ipt[7] = 4*Tj+2*Ta+Tv; // or just T
}


//======================  update_ip_pos_vel_acc: ======================
void update_ip_pos_vel_acc(std::vector<double> &P_jt_ipt, std::vector<double> &V_jt_ipt,std::vector<double> &A_jt_ipt,
                           std::vector<double> &T_jt_ipt, double &J_jt_ipt,
                           double P0, double V0, double A0 ){
    double Tj = T_jt_ipt[1];
    double Ta = T_jt_ipt[2] - Tj;
    double Tv = T_jt_ipt[4] - 2* Tj - Ta;
    double jr = 0,a=0, v=0, p=0;

    //initials
    P_jt_ipt[0] = P0;
    V_jt_ipt[0] = V0;
    A_jt_ipt[0] = A0;

    jr = J_jt_ipt;
    phase_j1 ( Tj,  Ta,  Tv,  P0,  V0,  A0,  jr,  T_jt_ipt[1],  a,  v,  p);
    P_jt_ipt[1] = p;
    V_jt_ipt[1] = v;
    A_jt_ipt[1] = a;

    phase_a1 ( Tj,  Ta,  Tv,  P0,  V0,  A0,  jr,  T_jt_ipt[2],  a,  v,  p);
    P_jt_ipt[2] = p;
    V_jt_ipt[2] = v;
    A_jt_ipt[2] = a;

    phase_j2 ( Tj,  Ta,  Tv,  P0,  V0,  A0,  jr,  T_jt_ipt[3],  a,  v,  p);
    P_jt_ipt[3] = p;
    V_jt_ipt[3] = v;
    A_jt_ipt[3] = a;

    phase_v  ( Tj,  Ta,  Tv,  P0,  V0,  A0,  jr,  T_jt_ipt[4],  a,  v,  p);
    P_jt_ipt[4] = p;
    V_jt_ipt[4] = v;
    A_jt_ipt[4] = a;

    phase_j3 ( Tj,  Ta,  Tv,  P0,  V0,  A0,  jr,  T_jt_ipt[5],  a,  v,  p);
    P_jt_ipt[5] = p;
    V_jt_ipt[5] = v;
    A_jt_ipt[5] = a;

    phase_a2 ( Tj,  Ta,  Tv,  P0,  V0,  A0,  jr,  T_jt_ipt[6],  a,  v,  p);
    P_jt_ipt[6] = p;
    V_jt_ipt[6] = v;
    A_jt_ipt[6] = a;

    phase_j4 ( Tj,  Ta,  Tv,  P0,  V0,  A0,  jr,  T_jt_ipt[7],  a,  v,  p);
    P_jt_ipt[7] = p;
    V_jt_ipt[7] = v;
    A_jt_ipt[7] = a;

}

//======================  compute_time_for_jrk: ======================
//this func computes the trajectory time for specific jerk
// it will be called at the begining with the max jerk as an argument to minimum optimal time for th etrajectory
double compute_time_for_jrk(const double Ds, double &Tj, double &Ta, double &Tv, const double sm, const double vm, const double am, const double jm ){
//    std::cout<< " ====== compute_time_for_jrk ============ "<< std::endl;

    double Dthr1 = (am*vm)/jm + (vm*vm)/am ;
    double Dthr2 = 2*pow(am,3)/pow(jm,2);
    double Tjm= am/jm;
    double Tam= (vm/am) - (am/jm);// how to get theis one

    double D1 = (am/jm - vm/am)*(am*(am/jm - vm/am) - pow(am,2)/(2*jm)) + pow(am,3)/(2*pow(jm,2)) - (am*(am*(am/jm - vm/am) - pow(am,2)/jm))/jm - (am*(am*(am/jm - vm/am) - pow(am,2)/(2*jm)))/jm - (pow(am,2)*(am/jm - vm/am))/(2*jm);
    // D1 from matlab using symbolic tool, till now it gives same results of

//    std::cout<<" D1: "<< D1  <<std::endl;
//    std::cout<<" Dthr1: "<< Dthr1 <<"   Dthr2: "<< Dthr2 <<std::endl;
//    std::cout<<" Tjm: "<< Tjm <<"   Tam: "<< Tam <<std::endl;
    if (Ds>= Dthr1){
//        std::cout<<" compute_time_for_jrk case 1 ... "<<std::endl;
        Tj= Tjm;
        Ta= Tam;
        Tv= (Ds-Dthr1)/vm;
    }
    else if(Ds>=Dthr2){
//        std::cout<<" compute_time_for_jrk  case 2 ... "<<std::endl;
        Tj= Tjm;
//        Ta= sqrt( ((am*am)/(4*jm)) + (Ds/am) ) - ( (3*am)/(2*jm) );
        Tv= 0;
        double jr = jm;
        double Ta1 = -(3*pow(Tj,2)*jr - sqrt(Tj*jr*(jr*pow(Tj,3) + 2*jr*pow(Tj,2)*Tv + jr*Tj*pow(Tv,2) + 4*Ds)) + Tj*Tv*jr)/(2*Tj*jr);
        double Ta2 = -(3*pow(Tj,2)*jr + sqrt(Tj*jr*(jr*pow(Tj,3) + 2*jr*pow(Tj,2)*Tv + jr*Tj*pow(Tv,2) + 4*Ds)) + Tj*Tv*jr)/(2*Tj*jr);
//        std::cout<<" Ta1: "<< Ta1 <<"   Ta2: "<< Ta2 <<std::endl;
        Ta = min_root(Ta1, Ta2);
    }else {
//        std::cout<<" compute_time_for_jrk case 3 ... "<<std::endl;
        Tv=0;
        Ta=0;
        std::vector<double> rts;
        rts.resize(3);
        int sgn = (Ds>0)? 1:-1;
        int n_rts = cubic_eq_real_root (2*sgn*jm, 0.0, 0.0, -Ds, rts);
        Tj= rts[0];
    }
    return 4*Tj+2*Ta+Tv;
}


//======================  compute_jerk_for_time: ======================
//compute_jerk_for_time: computes the jerk of a trajectory for specific time
// it allows us to increase or decrease the total time of a trajectory (> min_optimal_time)
double compute_jerk_for_time( double &T, const double DT, const double Ds, double &Tj, double &Ta, double &Tv, const double sm, const double vm, const double am, const double jm ){
    // another way is to use same Tj, Ta, Tv as the max dof or syn_time in case of N_DoF
//    std::cout<< " ====== compute_jerk_for_time ============ "<< std::endl;

//    double T_opt =  compute_time_for_jrk( Ds, Tj, Ta, Tv,  sm, vm,  am,  jm );
////    std::cout<<" the optimal time is: " << T_opt <<std::endl;
//    if(T+DT < T_opt)
//        throw(std::invalid_argument("syn_time is less than min optimal time ... "));

    double tj= Tj, ta = Ta, tv=Tv;
    double kj= Tj/T, ka= Ta/T, kv= Tv/T;
//    std::cout<<" old_Tjav: "<< Tj << "  "<< Ta << "  "<< Tv << "  " << std::endl;
//    std::cout<<" T, DT, Ds: "<< T << "  "<< DT << "  "<< Ds << "  " << std::endl;
    //update TJ, Ta, Tv
    Tj = tj+kj*DT;
    Ta = ta+ka*DT;
    Tv = tv+kv*DT;
//    std::cout<<" new_Tjav: "<< Tj << "  "<< Ta << "  "<< Tv << "  " << std::endl;
    double jrk =0;
    jrk = Ds/(((ta + DT*ka)*(tj + DT*kj) + pow((tj + DT*kj),2) )*(tj + DT*kj) + ((ta + DT*ka)*(tj + DT*kj) +
              pow((tj + DT*kj),2) )*(tv + DT*kv) + 0.5*pow((tj + DT*kj),3) + 0.5*(ta + DT*ka)*pow((tj + DT*kj),2)+
              ((ta + DT*ka)*(tj + DT*kj) + 0.5*pow((tj + DT*kj),2) )*(ta + DT*ka) + ((ta + DT*ka)*(tj + DT*kj) +
              0.5*pow((tj + DT*kj),2) )*(tj + DT*kj));
    if(fabs(jrk) - jm > 1)
        throw(std::invalid_argument(" new_jrk > max_jrk "));
    if( fabs(Tj*jrk)-am > .1)
        throw(std::invalid_argument(" new_acc > max_acc "));
    if( fabs( jrk*pow((tj + DT*kj),2) + jrk*(ta + DT*ka)*(tj + DT*kj) )- vm > .1)
        throw(std::invalid_argument(" new_vel > max_vel "));
    T = 4*Tj +2*Ta +Tv;
    return jrk;
}


//======================  cubic_eq_real_root: ======================
//to find the real root/roots of a cubic equation
int  cubic_eq_real_root (double a, double b, double c, double d, std::vector<double> &roots)
{
        if (a == 0.000)
        {
            throw(std::invalid_argument("The coefficient of the cube of x is 0. Please use the utility for a SECOND degree quadratic. No further action taken."));
            return 0;
        } //End if a == 0

        if (d == 0.000)
        {
            throw(std::invalid_argument("One root is 0. Now divide through by x and use the utility for a SECOND degree quadratic to solve the resulting equation for the other two roots. No further action taken."));
            return 0;
        } //End if d == 0
//        std::cout<<"Eq info: "<<std::endl;
//        std::cout<< a <<", "<< b <<", "<< c <<", "<< d <<std::endl;
        b /= a;
        c /= a;
        d /= a;
        double disc, q, r, dum1, s, t, term1, r13;
        q = (3.0*c - (b*b))/9.0;
        r = -(27.0*d) + b*(9.0*c - 2.0*(b*b));
        r /= 54.0;
        disc = q*q*q + r*r;
        double root1 = 0; //The first root is always real.
        term1 = (b/3.0);
        if (disc > 1e-10) { // one root real, two are complex
//            std::cout<<"disc > 0, disc= "<< disc <<std::endl;
            s = r + sqrt(disc);
            s = ((s < 0) ? -pow(-s, (1.0/3.0)) : pow(s, (1.0/3.0)));
            t = r - sqrt(disc);
            t = ((t < 0) ? -pow(-t, (1.0/3.0)) : pow(t, (1.0/3.0)));
            double x1r= -term1 + s + t;
            term1 += (s + t)/2.0;
            double x3r = -term1,  x2r = -term1;
            term1 = sqrt(3.0)*(-t + s)/2;
            x2r = term1;
            double x3i = -term1;
            roots[0]= x1r;
            roots[1] = -100.0;
            roots[2]= -100.0;
            return 1;
        }
        // End if (disc > 0)
        // The remaining options are all real
        double x3i =0, x2r = 0;
        if (disc>=0.000 && disc< 1e-10){ // All roots real, at least two are equal.
//            std::cout<<"disc = 0 , disc= "<< disc <<std::endl;
            disc=0;
            r13 = ((r < 0) ? -pow(-r,(1.0/3.0)) : pow(r,(1.0/3.0)));
            double x1r= -term1 + 2.0*r13;
            double x3r = -(r13 + term1);
            x2r = -(r13 + term1);
            roots[0]= x1r;
            roots[1] = x2r;
            roots[2]= x3r;
            return 2;
        } // End if (disc == 0)
        // Only option left is that all roots are real and unequal (to get here, q < 0)
//        std::cout<<"disc < 0 , disc= "<< disc <<std::endl;
        q = -q;
        dum1 = q*q*q;
        dum1 = acos(r/sqrt(dum1));
        r13 = 2.0*sqrt(q);
        double x1r= -term1 + r13*cos(dum1/3.0);
        x2r = -term1 + r13*cos((dum1 + 2.0*M_PI)/3.0);
        double x3r = -term1 + r13*cos((dum1 + 4.0*M_PI)/3.0);
        roots[0]= x1r;
        roots[1] = x2r;
        roots[2]= x3r;
        return 3;
}  //End of cubicSolve


//======================  min_root: ======================
//to find the minimum root between two or three roots
double min_root(double r1, double r2){
    double rt=0;
if (r1>=0 && r2>=0 )
    rt= (r1 < r2) ? r1  :  r2;
else if (r1>=0 )
    rt= r1;
else if (r2>=0 )
    rt= r2;
return rt;
}

// -------------------------------------------
double min_root(double r1, double r2, double r3){
    double rt=0;
if (r1>=0 && r2>=0 && r3>=0)
    rt= (r1 < r2) ? ((r1 < r3) ? r1 : r3)   :   ((r2 < r3) ? r2 : r3);
else if (r1>=0 && r2>=0)
    rt= (r1 < r2) ? r1 : r2;
else if (r1>=0 && r3>=0)
    rt= (r1 < r3) ? r1 : r3;
else if (r2>=0 && r3>=0)
    rt= (r2 < r3) ? r2 : r3;
else if (r1>=0)
    rt= r1;
else if (r2>=0)
    rt= r2;
else if (r3>=0 )
    rt= r3;
return rt;
}

//-----------------------

//======================  seven phases of s_curve: ======================
//these functions describes the 7-sets of equations described s_curve
void phase_j1 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p){
    a  = jr*t + A0;  // jr should has correct sign
    v = jr*pow(t,2)/2 + A0*t + V0;
    p = jr*pow(t,3)/6 + A0*pow(t,2)/2+ V0*t + P0;  // this start after case 3 so both vi and si becomes v2 and s2
}

void phase_a1 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p){
    double A1=0, V1=0, P1=0;  // new inials for current phase
    double ts= t-Tj;
    phase_j1(Tj, Ta, Tv, P0, V0, A0, jr, t-ts, A1, V1, P1); // call previous phase to get inyials for current phase
    a  = A1;  // jr should has correct sign
    v =  A1*ts + V1;
    p =  A1*pow(ts,2)/2  + V1*ts + P1;  // this start after case 3 so both vi and si becomes v2 and s2
}

void phase_j2 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p){
    double A2=0, V2=0, P2=0;  // new inials for current phase
    double ts= t - (Tj+Ta);
    phase_a1(Tj, Ta, Tv, P0, V0, A0, jr, t-ts, A2, V2, P2); // call previous phase to get inyials for current phase
    a  = -jr*ts+ A2;  // jr should has correct sign
    v =  -jr*pow(ts,2)/2+ A2*ts + V2;
    p =  -jr*pow(ts,3)/6 + A2*pow(ts,2)/2  + V2*ts + P2;  // this start after case 3 so both vi and si becomes v2 and s2
}

void phase_v (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p){
    double A3=0, V3=0, P3=0;  // new inials for current phase
    double ts= t -(2*Tj+Ta);
    phase_j2(Tj, Ta, Tv, P0, V0, A0, jr, t-ts, A3, V3, P3); // call previous phase to get inyials for current phase
    a  = A3;  // jr should has correct sign
    v =  V3;
    p =  V3*ts + P3;  // this start after case 3 so both vi and si becomes v2 and s2
}

void phase_j3 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p){
    double A4=0, V4=0, P4=0;  // new inials for current phase
    double ts= t - (2*Tj+Ta+Tv);
    phase_v(Tj, Ta, Tv, P0, V0, A0, jr, t-ts, A4, V4, P4); // call previous phase to get inyials for current phase
    a  = -jr*ts+ A4;  // jr should has correct sign
    v =  -jr*pow(ts,2)/2+ A4*ts + V4;
    p =  -jr*pow(ts,3)/6 + A4*pow(ts,2)/2  + V4*ts + P4;  // this start after case 3 so both vi and si becomes v2 and s2
}

void phase_a2 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p){
    double A5=0, V5=0, P5=0;  // new inials for current phase
    double ts= t - (3*Tj+Ta+Tv);
    phase_j3(Tj, Ta, Tv, P0, V0, A0, jr, t-ts, A5, V5, P5); // call previous phase to get inyials for current phase
    a  = A5;  // jr should has correct sign
    v =  A5*ts + V5;
    p =  A5*pow(ts,2)/2  + V5*ts + P5;  // this start after case 3 so both vi and si becomes v2 and s2
}

void phase_j4 (double Tj, double Ta, double Tv, double P0, double V0, double A0, double jr, double t, double &a, double &v, double &p){
    double A6=0, V6=0, P6=0;  // new inials for current phase
    double ts= t - (3*Tj+2*Ta+Tv);
    phase_a2(Tj, Ta, Tv, P0, V0, A0, jr, t-ts, A6, V6, P6); // call previous phase to get inyials for current phase
    a  = jr*ts + A6;  // jr should has correct sign
    v = jr*pow(ts,2)/2 + A6*ts + V6;
    p = jr*pow(ts,3)/6 + A6*pow(ts,2)/2+ V6*ts + P6;  // this start after case 3 so both vi and si becomes v2 and s2
}





