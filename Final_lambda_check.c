//
//  Final_lambda_check.c
//  
//
//  Created by Viraj Sanghai on 10/02/2017.
//
//

//
//  Final_light_lambda.c
//
//
//  Created by Viraj Sanghai on 21/01/2017.
//
//



#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>


#define PI 3.14159265358979323846

//Lambda_only_code

int main (int argc, char *argv[])
{   clock_t begin = clock();
    
    FILE *f =fopen(argv[1],"w");
    int n_tot;
    for (n_tot=1250; n_tot>0;n_tot --){
        /*Initial size of scale factor-Typical current separation of galaxies of about 1 Mpc*/
        double X_0 =0.5;
        double X_initial= X_0;
        /*initial position of observer*/
        double x = -0.4;
        double y = 0.1;
        double z= 0.0;
        double x_mu[3] = {x,y,z};
        /* t is the lookback time. The inital time, t_0 =0.*/
        //double t_0 = 3000.0;
        double t = 0.0;
        /*speed of light*/
        double c = 1.0;
        /*Hubble constant*/
        double H_0 =  0.0002443023433231266;//0.7324/3000.0;
        /*Value of cosmological constant in 1/Mpc^2*/
        double lambda = (3.0)*H_0*H_0;
        double lam_div_3 = lambda/3.0;
        double lam_div_6 = lambda/6.0;
        double lam_div_12 = lambda/12.0;
        /*The edge of the cell is given by 2*X where X is the position from the
         %centre of the cell to the centre of the boundary and behaves like the scale factor. We choose the inital cell size to be about a 1 Mpc. The scale factor is then related to the looback time in the following way. Not cosmic time */
        /*Value of X*/
        double X=X_0*exp(-sqrt(lam_div_3)*t);
        /*printf("value of X is %.15f\n",X);*/
        
        /*%Choosing an end point for our cell size before the light ray diverges*/
        double X_end =0.2;
        
        /*Time derivative of scale factor set using the Friedmann equation at leading order*/
        double X_tder=X_0*sqrt(lam_div_3)*exp(-sqrt(lam_div_3)*t);
        double H = sqrt(lam_div_3) ;
        double H_sq = H*H;
        double r_sq = x*x + y*y + z*z;
        /*Initial total potential*/
        double Phi =(lam_div_6)*r_sq;
        double Psi = -(lam_div_12)*r_sq;
        double g_ab[16] = {-1+2.0*Phi,0,0,0,0,1+2.0*Psi,0,0,0,0,1+2.0*Psi,0,0,0,0,1+2.0*Psi};
        //Four velocity with index raised- Infact all vectors with _a are with their index raised
        double u_a[4] = {(1.0 + H_sq*r_sq/2.0 + Phi),H*x,H*y,H*z};
        
        /*Initial four velocity with index lowered - Infact all vectors with _undera are with their index lowered */
        double uinitial_undera[4] = {(-1 - H_sq*r_sq/2.0 + Phi),H*x,H*y,H*z};
        //double u_undera_prev[4] = {(-1 - H_sq*r_sq/2.0 + Phi),H*x,H*y,H*z};
        double u_undera[4] = {-1.0 - H_sq*r_sq/2.0 + Phi, H*x,H*y,H*z};
        double ubar_a[4] = {1,0,0,0};
        //initial frequency
        double omega_zero = 1.0;
        
        
        //Seed initial random directions of light ray - set using the product of time of clock in secs and time of clock in microsecs
        struct timeval tseed;
        //srand((unsigned)time(&tseed));
        gettimeofday(&tseed, NULL);
        srand(tseed.tv_usec*tseed.tv_sec);
        printf("tseed = %.15ld \n", tseed.tv_usec*tseed.tv_sec);
        printf("rand_max = %d \n", RAND_MAX);
        // Weighted by the value of angles on a sphere
        double theta= acos((rand())*2.0/(double)(RAND_MAX) - 1);
        printf("theta = %.15f \n", theta);
        double varphi = (rand())*2.0*PI/(double)(RAND_MAX);
        printf("varphi = %.15f \n ", varphi);
        
        /*Inital value of line of sight vector*/
        /*dbar is the background value.*/
        double dbar_t = 0.0;
        double dbar_x=-cos(varphi)*sin(theta);
        double dbar_y =-sin(varphi)*sin(theta);
        double dbar_z = -cos(theta);
        double dbar_mu[3] = {dbar_x, dbar_y, dbar_z};
        double dbar_a[4] = {0,dbar_x, dbar_y, dbar_z};
        double kbar_a[4];
        double kbar_mu[3];
        double kbar_t;
        double k_a[4];
        double k_mu[3];
        double k_t = 0;
        double d_a[4]={0,dbar_x, dbar_y, dbar_z};
        
        /*Inital value of wavevector of light*/
        /*kbar is the background value.*/
        kbar_a[0] = omega_zero*(ubar_a[0] + dbar_a[0]);
        kbar_a[1] = omega_zero*(ubar_a[1] + dbar_a[1]);
        kbar_a[2] = omega_zero*(ubar_a[2] + dbar_a[2]);
        kbar_a[3] = omega_zero*(ubar_a[3] + dbar_a[3]);
        
        kbar_t =kbar_a[0];
        kbar_mu[0] = kbar_a[1];
        kbar_mu[1] = kbar_a[2];
        kbar_mu[2] = kbar_a[3];
        
        double dbar_dot_x_mu;
        dbar_dot_x_mu = dbar_mu[0]*x_mu[0]+dbar_mu[1]*x_mu[1] + dbar_mu[2]*x_mu[2];
        
        //Set inital values of k
        k_t = omega_zero*(1.0 + Phi + H_sq*r_sq/2.0 + H*dbar_dot_x_mu);
        
        k_mu[0] = omega_zero*((1.0 - Psi)*dbar_mu[0] + x_mu[0]*(H_sq*dbar_dot_x_mu/2.0) + H*x_mu[0]);
        k_mu[1] = omega_zero*((1.0 - Psi)*dbar_mu[1] + x_mu[1]*(H_sq*dbar_dot_x_mu/2.0) + H*x_mu[1]);
        k_mu[2] = omega_zero*((1.0 - Psi)*dbar_mu[2] + x_mu[2]*(H_sq*dbar_dot_x_mu/2.0) + H*x_mu[2]);
        
        k_a[0] =k_t;
        k_a[1] = k_mu[0];
        k_a[2] = k_mu[1];
        k_a[3] = k_mu[2];
        
        double u_undera_dot_k_a = u_undera[0]*k_a[0] + u_undera[1]*k_a[1] + u_undera[2]*k_a[2]+ u_undera[3]*k_a[3];
        
        /*Inital value of Sachs basis which is taken orthogonal to the wavevector of light*/
        double sybar_t = 0;
        double sybar_x=sin(varphi);
        double sybar_y =-cos(varphi);
        double sybar_z = 0;
        double sybar_mu[3] = {sybar_x,sybar_y,sybar_z};
        double sybar_a[4] ={sybar_t,sybar_x,sybar_y,sybar_z};
        double sy_a[4];
        double sy_mu[3];
        double sy_t;
        
        //Initial value of frequency
        double omega_prev = -u_undera_dot_k_a;
        //printf("omega_prev %.15f \n", omega_prev);
        
        double sybar_dot_x_mu;
        sybar_dot_x_mu = sybar_mu[0]*x_mu[0] + sybar_mu[1]*x_mu[1] + sybar_mu[2]*x_mu[2];
        
        //Initial sy_a
        sy_t =  H*sybar_dot_x_mu;
        sy_mu[0] = (1.0 - Psi)*sybar_mu[0] + x_mu[0]*(H_sq*(sybar_dot_x_mu)/2.0);
        sy_mu[1] = (1.0 - Psi)*sybar_mu[1] + x_mu[1]*(H_sq*(sybar_dot_x_mu)/2.0);
        sy_mu[2] = (1.0 - Psi)*sybar_mu[2] + x_mu[2]*(H_sq*(sybar_dot_x_mu)/2.0);
        
        sy_a[0] = sy_t;
        sy_a[1] = sy_mu[0];
        sy_a[2] = sy_mu[1];
        sy_a[3] = sy_mu[2];
        
        double szbar_t = 0;
        double szbar_x=-cos(varphi)*cos(theta);
        double szbar_y =-sin(varphi)*cos(theta);
        double szbar_z = sin(theta);
        
        double szbar_mu[3] = {szbar_x,szbar_y,szbar_z};
        double szbar_a[4] ={szbar_t,szbar_x,szbar_y,szbar_z};
        double sz_a[4];
        double sz_mu[3];
        double sz_t;
        
        double szbar_dot_x_mu = 0;
        szbar_dot_x_mu = szbar_mu[0]*x_mu[0] + szbar_mu[1]*x_mu[1] + szbar_mu[2]*x_mu[2];
        
        //Initial sz_a
        sz_t =  H*szbar_dot_x_mu;
        sz_mu[0] = (1.0 - Psi)*szbar_mu[0] + x_mu[0]*(H_sq*(szbar_dot_x_mu)/2.0);
        sz_mu[1] = (1.0 - Psi)*szbar_mu[1] + x_mu[1]*(H_sq*(szbar_dot_x_mu)/2.0);
        sz_mu[2] = (1.0 - Psi)*szbar_mu[2] + x_mu[2]*(H_sq*(szbar_dot_x_mu)/2.0);
        
        sz_a[0] = sz_t;
        sz_a[1] = sz_mu[0];
        sz_a[2] = sz_mu[1];
        sz_a[3] = sz_mu[2];
        
        /*Initialize the Wronski matrix which is use to calculate the jacobi matrix*/
        double W[16] = {1.0,0,0,0,0,1.0,0,0,0,0,1.0,0,0,0,0,1.0};
        //The optical tidal matrix is null for a cosmological constant
        double R_yy = 0;
        double R_yz = 0;
        double R_zz = 0;
        
        /*time step*/
        double dt= 0.01;
        
        /*Counters to count the number of reflections and iterations*/
        int n =0;
        int n2 =0;
        int n3=0;
        int n4=0;
        int n5=0;
        //Initial affine parameter
        double lam = 0;
        
        //Initial derivatives of potentials
        double Phider_mu[3]={(lam_div_3)*x,(lam_div_3)*y,(lam_div_3)*z};
        double Psider_mu[3]={-(lam_div_6)*x,-(lam_div_6)*y,-(lam_div_6)*z};
        //Initial change in affine parameter
        double dlam = 0;
        
        //Dot products needed for evaluation of sachs basis evolution
        double Phider_mu_dot_k_mu;
        double k_mu_dot_k_mu;
        double Phider_mu_dot_sy_mu;
        double du_dot_sy_a;
        double Phider_mu_u_undera_mu;
        double Psider_mu_u_undera_mu;
        double sy_mu_dot_sy_mu;
        double sy_mu_dot_k_mu;
        double Psider_mu_dot_sy_mu;
        double sy_mu_u_undera_mu;
        double Phider_mu_dot_sz_mu;
        double du_dot_sz_a;
        double sz_mu_dot_sz_mu;
        double sz_mu_dot_k_mu;
        double Psider_mu_dot_sz_mu;
        double sz_mu_u_undera_mu;
        double Psider_mu_dot_k_mu;
        double k_mu_u_undera_mu;
        double sy_u_d_brackets;
        double sz_u_d_brackets;
        
        //Initial second derivatives of potentials
        double Phisecder_lambda[9]= {(lambda/3.0),0,0,0,(lambda/3.0),0,0,0,(lambda/3.0)};
        double Psisecder_lambda[9]= {-(lambda/6.0),0,0,0,-(lambda/6.0),0,0,0,-(lambda/6.0)};
        
        //Dot products needed for evaluation of Wronski matrix evolution
        double sy_Phisecder_sy;
        double kbar_mu_dot_kbar_mu;
        double sy_Psisecder_sy;
        double k_Phisecder_k;
        double k_Psisecder_k;
        double sybar_mu_dot_sybar_mu;
        double sy_Phisecder_sz;
        double sy_Psisecder_sz;
        double sz_Phisecder_sz;
        double sz_Psisecder_sz;
        double szbar_mu_dot_szbar_mu;
        //Deformation rate matrix
        double S[4];
        double inv_den;
        double n_a[4];
        double n_undera[4];
        double k_a_dot_n_undera;
        double sy_a_dot_n_undera;
        double sz_a_dot_n_undera;
        double k_a_dot_k_a;
        double sy_a_dot_sy_a;
        double sz_a_dot_sz_a;
        double u_undera_dot_sy_a;
        double u_undera_dot_sz_a;
        double u_a_dot_n_undera;
        //Angular diameter distance
        double D_A;
        double k_mu_dot_x_mu;
        double omega=1.0;
        double i,j,k;
        // Sum of time part in matter potentials for 11^3 point sources
        double den1= 71.418862580858601, den2, final_pot=0, X_two, x_num,y_num,z_num;
        double final_der[3];
        double M_den_cube;
        double final_pot_xx, final_pot_yy, final_pot_zz, final_pot_xy, final_pot_yz, final_pot_xz, den2_sq;
        /*double pot_array_xx, pot_array_yy, pot_array_zz, pot_array_xy, pot_array_yz, pot_array_xz;*/
        double den2_cube, den2_five,  M_den_five;
        double x_num_sq, y_num_sq, z_num_sq, r=sqrt(r_sq), final_m, totalz=1.0, det;
        double save_z=0.1;
      
        
        // Set radius of compact object at centre R
        double R=0.01;
        double R_sq, M_R_cube;
        double z_FRW;
        double dlam_corr, x_out, y_out, z_out;
        //Conformal time step
        double dtau=0.005;
        double u_undera_dot_u_a;
        double u_undera_dot_d_a;
        double d_a_dot_sy_a;
        double sy_a_dot_k_a;
        double sz_a_dot_k_a;
        double u_der[16]={0, -H_sq*x + Phider_mu[0], -H_sq*y + Phider_mu[1], -H_sq*z + Phider_mu[2],0,H,0,0,0,0,H,0,0, 0,0,H};
        double z_diff;
        double D_A_flrw;
        
        //For testing code
        /*double sy_a_prev[4],sz_a_prev[4];
         
         sy_a_prev[0] = sy_t;
         sy_a_prev[1] = sy_mu[0];
         sy_a_prev[2] = sy_mu[1];
         sy_a_prev[3] = sy_mu[2];
         
         sz_a_prev[0] = sz_t;
         sz_a_prev[1] = sz_mu[0];
         sz_a_prev[2] = sz_mu[1];
         sz_a_prev[3] = sz_mu[2];
         
         double LHS_y, LHS_z;*/
        
        
        /*d_a[0] = k_a[0]/(omega_zero) - u_a[0];
         d_a[1] = k_a[1]/(omega_zero) - u_a[1];
         d_a[2] = k_a[2]/(omega_zero) - u_a[2];
         d_a[3] = k_a[3]/(omega_zero) - u_a[3];
         
         u_undera_dot_sy_a = u_undera[0]*sy_a[0] + u_undera[1]*sy_a[1] + u_undera[2]*sy_a[2]+ u_undera[3]*sy_a[3];
         u_undera_dot_sz_a = u_undera[0]*sz_a[0] + u_undera[1]*sz_a[1] + u_undera[2]*sz_a[2]+ u_undera[3]*sz_a[3];
         u_undera_dot_u_a = u_undera[0]*u_a[0] + u_undera[1]*u_a[1] + u_undera[2]*u_a[2]+ u_undera[3]*u_a[3];
         u_undera_dot_d_a = u_undera[0]*d_a[0] + u_undera[1]*d_a[1] + u_undera[2]*d_a[2]+ u_undera[3]*d_a[3];
         d_a_dot_sy_a = (-1+2.0*Phi)*d_a[0]*sy_a[0]+(1+2.0*Psi)*d_a[1]*sy_a[1]+(1+2.0*Psi)*d_a[2]*sy_a[2]+(1+2.0*Psi)*d_a[3]*sy_a[3];
         sy_a_dot_k_a = (-1+2.0*Phi)*k_a[0]*sy_a[0]+(1+2.0*Psi)*k_a[1]*sy_a[1]+(1+2.0*Psi)*k_a[2]*sy_a[2]+(1+2.0*Psi)*k_a[3]*sy_a[3];
         sz_a_dot_k_a = (-1+2.0*Phi)*k_a[0]*sz_a[0]+(1+2.0*Psi)*k_a[1]*sz_a[1]+(1+2.0*Psi)*k_a[2]*sz_a[2]+(1+2.0*Psi)*k_a[3]*sz_a[3];
         
         //Inital constraint tests of code
         printf("check4 %.15f \n",u_undera_dot_sy_a);
         printf("check5 %.15f \n",u_undera_dot_sz_a);
         //printf("check6 %.15f \n",u_a_dot_n_undera);
         printf("check7 %.15f \n",u_undera_dot_u_a);
         printf("check8 %.15f \n",u_undera_dot_d_a);
         printf("check9 %.15f \n",d_a_dot_sy_a);
         printf("check10 %.15f \n",sy_a_dot_k_a);
         printf("check11 %.15f \n",sz_a_dot_k_a);*/
        
        double s_change_t;
        double s_change_mu[3];
        double sy_Phisecder_k, sz_Phisecder_k, sy_Psisecder_k, sz_Psisecder_k, sy_mu_dot_sz_mu;
        double sign=1.0;
        double dz_dlam, dlam_save, du_dot_k_a, W_save[16], t_save, lam_save, X_save;
        double totalzcheck = 1.0;
        double omega_prevcheck=1.0;
        double dt_x,dt_y,dt_z, dt_corr;
        double X_next, X_x, X_y, X_z;
        
        
        
        //While redshift is less than 1.5
        while((totalzcheck-1)<1.5){
            
            //very low redshift test
            //while((totalz-1)<0.00001){
            
            //Two reflection test
            //while(n2<2){
            
            /*Evaluate position of boundary*/
            X = X_0*exp(-sqrt(lam_div_3)*t);
            X_tder=X*sqrt(lam_div_3);
            H= X_tder/X;
            H_sq = H*H;
            r_sq = x*x + y*y + z*z;
            
            /*Potentials Phi and Psi for a universe with only a cosmological constant*/
            Phi = (lam_div_6)*r_sq;
            Psi = -(lam_div_12)*r_sq;
            
            /*Four velocity with index raised*/
            u_a[0] = (1.0 + H_sq*r_sq/2.0 + Phi);
            u_a[1]=H*x;
            u_a[2]=H*y;
            u_a[3]=H*z;
            
            /*Initial four velocity with index lowered*/
            u_undera[0] = -1.0 - H_sq*r_sq/2.0 + Phi;
            u_undera[1]=u_a[1];
            u_undera[2]=u_a[2];
            u_undera[3]=u_a[3];
            
            
            u_undera_dot_k_a = u_undera[0]*k_a[0] + u_undera[1]*k_a[1] + u_undera[2]*k_a[2]+ u_undera[3]*k_a[3];
            omega = -u_undera_dot_k_a;
            totalzcheck = omega;
            //omega_prevcheck = -u_undera_dot_k_a;
            
            /*d_a[0] = k_a[0]/(-u_undera_dot_k_a) - u_a[0];
             d_a[1] = k_a[1]/(-u_undera_dot_k_a) - u_a[1];
             d_a[2] = k_a[2]/(-u_undera_dot_k_a) - u_a[2];
             d_a[3] = k_a[3]/(-u_undera_dot_k_a) - u_a[3];*/
            
            u_undera_dot_sy_a = u_undera[0]*sy_a[0] + u_undera[1]*sy_a[1] + u_undera[2]*sy_a[2]+ u_undera[3]*sy_a[3];
            u_undera_dot_sz_a = u_undera[0]*sz_a[0] + u_undera[1]*sz_a[1] + u_undera[2]*sz_a[2]+ u_undera[3]*sz_a[3];
            Phider_mu_dot_sy_mu = Phider_mu[0]*sy_mu[0] + Phider_mu[1]*sy_mu[1] + Phider_mu[2]*sy_mu[2];
            
            //Short run test of code
            
            /*printf("check4 %.15f \n",u_undera_dot_sy_a);
             printf("check5 %.15f \n",u_undera_dot_sz_a);
             d_a_dot_sy_a = (-1+2.0*Phi)*d_a[0]*sy_a[0]+(1+2.0*Psi)*d_a[1]*sy_a[1]+(1+2.0*Psi)*d_a[2]*sy_a[2]+(1+2.0*Psi)*d_a[3]*sy_a[3];
             printf("check9 %.15f \n",d_a_dot_sy_a);*/
            /*printf("sz_t %.15f \n",sy_t);
             printf("sz_brack %.15f \n",(k_t/omega)*sy_u_d_brackets);
             printf("sz_du %.15f \n",(k_t/omega)*du_dot_sy_a);
             printf("sz_leftover %.15f \n",(Phider_mu_dot_sy_mu)*dlam*k_t + (Phider_mu_dot_k_mu)*dlam*sy_t);
             printf("kt_omega %.15f \n",(k_t/omega));
             printf("main %.15f \n", k_t*Phider_mu_dot_sz_mu);
             printf("dot %.15f \n", sz_mu_dot_sz_mu);*/
            /*printf("LHS_y %.15f \n", LHS_y);
             printf ("RHS_y %.15f \n", du_dot_sy_a);
             printf("y_sum %.15f \n", LHS_y + du_dot_sy_a);
             
             printf("LHS_z %.15f \n", LHS_z);
             printf ("RHS_z %.15f \n", du_dot_sz_a);
             printf("z_sum %.15f \n", LHS_z + du_dot_sz_a);
             
             printf("s_change_t %.15f \n", s_change_t);
             printf ("s_change_mu[0] %.15f \n", s_change_mu[0]);
             printf ("s_change_mu[1] %.15f \n", s_change_mu[1]);
             printf ("s_change_mu[2] %.15f \n", s_change_mu[2]);*/
            
            /*Derivative of Phi and Psi*/
            Phider_mu[0] = (lam_div_3)*x;
            Phider_mu[1] = (lam_div_3)*y;
            Phider_mu[2] = (lam_div_3)*z;
            Psider_mu[0] = -(lam_div_6)*x;
            Psider_mu[1] = -(lam_div_6)*y;
            Psider_mu[2] = -(lam_div_6)*z;
            
            //Derivatives of 4-velocity
            u_der[0]=0;
            u_der[1]= -H_sq*x + Phider_mu[0];
            u_der[2]= -H_sq*y+ Phider_mu[1];
            u_der[3]= -H_sq*z + Phider_mu[2];
            u_der[4]= 0;
            u_der[5]= H;
            u_der[6]=0;
            u_der[7]=0;
            u_der[8]=0;
            u_der[9]= 0;
            u_der[10]= H;
            u_der[11]= 0;
            u_der[12]= 0;
            u_der[13]= 0;
            u_der[14]= 0;
            u_der[15]= H;
                
            /*z_FRW = (X_initial/X) - 1.0;
            printf("z_FRW %.15f \n", z_FRW);
            printf("totalz %.15f \n", totalzcheck-1.0);
            printf("diff %.15f \n", (totalzcheck-1.0 - z_FRW));*/
            
            // Save values of redshift and angular diameter distance at each 0.1 redshift
            if  (totalzcheck-1 >= save_z){
                //printf("value of det is %.15f\n",det);
                
                z_FRW = (X_initial/X) - 1.0;
                /*printf("z_FRW %.15f \n", z_FRW);
                printf("totalz %.15f \n", totalzcheck-1.0);
                printf("diff %.15f \n", (totalzcheck-1.0 - z_FRW)/z_FRW);*/
                
                //Tests for code
                /*k_a_dot_k_a = (-1+2.0*Phi)*k_a[0]*k_a[0]+(1+2.0*Psi)*k_a[1]*k_a[1]+(1+2.0*Psi)*k_a[2]*k_a[2]+(1+2.0*Psi)*k_a[3]*k_a[3];
                 sy_a_dot_sy_a = (-1+2.0*Phi)*sy_a[0]*sy_a[0]+(1+2.0*Psi)*sy_a[1]*sy_a[1]+(1+2.0*Psi)*sy_a[2]*sy_a[2]+(1+2.0*Psi)*sy_a[3]*sy_a[3];
                 sz_a_dot_sz_a = (-1+2.0*Phi)*sz_a[0]*sz_a[0]+(1+2.0*Psi)*sz_a[1]*sz_a[1]+(1+2.0*Psi)*sz_a[2]*sz_a[2]+(1+2.0*Psi)*sz_a[3]*sz_a[3];
                 u_undera_dot_sy_a = u_undera[0]*sy_a[0] + u_undera[1]*sy_a[1] + u_undera[2]*sy_a[2]+ u_undera[3]*sy_a[3];
                 u_undera_dot_sz_a = u_undera[0]*sz_a[0] + u_undera[1]*sz_a[1] + u_undera[2]*sz_a[2]+ u_undera[3]*sz_a[3];
                 u_undera_dot_u_a = u_undera[0]*u_a[0] + u_undera[1]*u_a[1] + u_undera[2]*u_a[2]+ u_undera[3]*u_a[3];
                 u_undera_dot_d_a = u_undera[0]*d_a[0] + u_undera[1]*d_a[1] + u_undera[2]*d_a[2]+ u_undera[3]*d_a[3];
                 d_a_dot_sy_a = (-1+0*2.0*Phi)*d_a[0]*sy_a[0]+(1+0*2.0*Psi)*d_a[1]*sy_a[1]+(1+0*2.0*Psi)*d_a[2]*sy_a[2]+(1+0*2.0*Psi)*d_a[3]*sy_a[3];
                 //u_undera_dot_sz_a = (-1.0+2.0*Phi)*u_a[0]*sz_a[0] + (1.0+2.0*Psi)*u_a[1]*sz_a[1] + (1.0+2.0*Psi)*u_a[2]*sz_a[2]+ (1.0+2.0*Psi)*u_a[3]*sz_a[3];
                 sy_a_dot_k_a = (-1+2.0*Phi)*k_a[0]*sy_a[0]+(1+2.0*Psi)*k_a[1]*sy_a[1]+(1+2.0*Psi)*k_a[2]*sy_a[2]+(1+2.0*Psi)*k_a[3]*sy_a[3];
                 sz_a_dot_k_a = (-1+2.0*Phi)*k_a[0]*sz_a[0]+(1+2.0*Psi)*k_a[1]*sz_a[1]+(1+2.0*Psi)*k_a[2]*sz_a[2]+(1+2.0*Psi)*k_a[3]*sz_a[3];
                 
                 printf("check1 %.15f \n",k_a_dot_k_a);
                 printf("check2 %.15f \n",sy_a_dot_sy_a);
                 printf("check3 %.15f \n",sz_a_dot_sz_a);
                 printf("check4 %.15f \n",u_undera_dot_sy_a);
                 printf("check5 %.15f \n",u_undera_dot_sz_a);
                 printf("check6 %.15f \n",u_a_dot_n_undera);
                 printf("check7 %.15f \n",u_undera_dot_u_a);
                 printf("check8 %.15f \n",u_undera_dot_d_a);
                 printf("check9 %.15f \n",d_a_dot_sy_a);
                 printf("check10 %.15f \n",sy_a_dot_k_a);
                 printf("check11 %.15f \n",sz_a_dot_k_a);
                 printf("LHS_y %.15f \n", LHS_y);
                 printf ("RHS_y %.15f \n", du_dot_sy_a);
                 printf("y_sum %.15f \n", LHS_y + du_dot_sy_a);
                 
                 printf("LHS_z %.15f \n", LHS_z);
                 printf ("RHS_z %.15f \n", du_dot_sz_a);
                 printf("z_sum %.15f \n", LHS_z + du_dot_sz_a);
                 
                 printf("s_change_t %.15f \n", s_change_t);
                 printf ("s_change_mu[0] %.15f \n", s_change_mu[0]);
                 printf ("s_change_mu[1] %.15f \n", s_change_mu[1]);
                 printf ("s_change_mu[2] %.15f \n", s_change_mu[2]);*/
                
                
                
                Phider_mu_u_undera_mu = Phider_mu[0]*u_undera[1] + Phider_mu[1]*u_undera[2] + Phider_mu[2]*u_undera[3];
                Psider_mu_u_undera_mu = Psider_mu[0]*u_undera[1] + Psider_mu[1]*u_undera[2] + Psider_mu[2]*u_undera[3];
                k_mu_dot_k_mu = k_mu[0]*k_mu[0] + k_mu[1]*k_mu[1] + k_mu[2]*k_mu[2];
                Psider_mu_dot_k_mu = Psider_mu[0]*k_mu[0] + Psider_mu[1]*k_mu[1] + Psider_mu[2]*k_mu[2];
                k_mu_u_undera_mu = k_mu[0]*u_undera[1] + k_mu[1]*u_undera[2] + k_mu[2]*u_undera[3];
                Phider_mu_dot_k_mu = Phider_mu[0]*k_mu[0] + Phider_mu[1]*k_mu[1] + Phider_mu[2]*k_mu[2];
                
                
                du_dot_k_a = (k_a[0]*u_der[0] + k_a[1]*u_der[4] + k_a[2]*u_der[8] + k_a[3]*u_der[12])*k_t
                + (k_a[0]*u_der[1] + k_a[1]*u_der[5] + k_a[2]*u_der[9] + k_a[3]*u_der[13])*k_mu[0] + (k_a[0]*u_der[2] + k_a[1]*u_der[6] + k_a[2]*u_der[10] + k_a[3]*u_der[14])*k_mu[1] + (k_a[0]*u_der[3] + k_a[1]*u_der[7] + k_a[2]*u_der[11] + k_a[3]*u_der[15])*k_mu[2] ;
                
                dz_dlam = -( du_dot_k_a + (Phider_mu_u_undera_mu)*k_t*k_t + (Phider_mu_dot_k_mu)*k_t*u_undera[0] + (Phider_mu_dot_k_mu)*k_t*u_undera[0] + (Psider_mu_u_undera_mu)*(k_mu_dot_k_mu) - 2.0*(Psider_mu_dot_k_mu)*(k_mu_u_undera_mu));
                
                dlam_save = -(1.0/(dz_dlam))*(totalzcheck-1.0 -save_z);
                
                //printf("dlam %.15f \n", dlam_save);
                //printf("totalz %.15f \n", totalzcheck-1.0);
                
                sy_Phisecder_sy = (sy_mu[0]*Phisecder_lambda[0] + sy_mu[1]*Phisecder_lambda[3] + sy_mu[2]*Phisecder_lambda[6])*sy_mu[0] \
                + (sy_mu[0]*Phisecder_lambda[1] + sy_mu[1]*Phisecder_lambda[4] + sy_mu[2]*Phisecder_lambda[7])*sy_mu[1] \
                + (sy_mu[0]*Phisecder_lambda[2] + sy_mu[1]*Phisecder_lambda[5] + sy_mu[2]*Phisecder_lambda[8])*sy_mu[2];
                
                k_mu_dot_k_mu = k_mu[0]*k_mu[0]+ k_mu[1]*k_mu[1]+ k_mu[2]*k_mu[2] ;
                
                sy_Psisecder_sy = (sy_mu[0]*Psisecder_lambda[0] + sy_mu[1]*Psisecder_lambda[3] + sy_mu[2]*Psisecder_lambda[6])*sy_mu[0] \
                + (sy_mu[0]*Psisecder_lambda[1] + sy_mu[1]*Psisecder_lambda[4] + sy_mu[2]*Psisecder_lambda[7])*sy_mu[1] \
                + (sy_mu[0]*Psisecder_lambda[2] + sy_mu[1]*Psisecder_lambda[5] + sy_mu[2]*Psisecder_lambda[8])*sy_mu[2];
                
                k_Psisecder_k = (k_mu[0]*Psisecder_lambda[0] + k_mu[1]*Psisecder_lambda[3] + k_mu[2]*Psisecder_lambda[6])*k_mu[0] \
                + (k_mu[0]*Psisecder_lambda[1] + k_mu[1]*Psisecder_lambda[4] + k_mu[2]*Psisecder_lambda[7])*k_mu[1] \
                + (k_mu[0]*Psisecder_lambda[2] + k_mu[1]*Psisecder_lambda[5] + k_mu[2]*Psisecder_lambda[8])*k_mu[2];
                
                k_Phisecder_k = (k_mu[0]*Phisecder_lambda[0] + k_mu[1]*Phisecder_lambda[3] + k_mu[2]*Phisecder_lambda[6])*k_mu[0] \
                + (k_mu[0]*Phisecder_lambda[1] + k_mu[1]*Phisecder_lambda[4] + k_mu[2]*Phisecder_lambda[7])*k_mu[1] \
                + (k_mu[0]*Phisecder_lambda[2] + k_mu[1]*Phisecder_lambda[5] + k_mu[2]*Phisecder_lambda[8])*k_mu[2];
                
                sy_mu_dot_sy_mu = sy_mu[0]*sy_mu[0]+ sy_mu[1]*sy_mu[1]+ sy_mu[2]*sy_mu[2];
                
                sy_Phisecder_sz = (sy_mu[0]*Phisecder_lambda[0] + sy_mu[1]*Phisecder_lambda[3] + sy_mu[2]*Phisecder_lambda[6])*sz_mu[0] \
                + (sy_mu[0]*Phisecder_lambda[1] + sy_mu[1]*Phisecder_lambda[4] + sy_mu[2]*Phisecder_lambda[7])*sz_mu[1] \
                + (sy_mu[0]*Phisecder_lambda[2] + sy_mu[1]*Phisecder_lambda[5] + sy_mu[2]*Phisecder_lambda[8])*sz_mu[2];
                
                sy_Psisecder_sz = (sy_mu[0]*Psisecder_lambda[0] + sy_mu[1]*Psisecder_lambda[3] + sy_mu[2]*Psisecder_lambda[6])*sz_mu[0] \
                + (sy_mu[0]*Psisecder_lambda[1] + sy_mu[1]*Psisecder_lambda[4] + sy_mu[2]*Psisecder_lambda[7])*sz_mu[1] \
                + (sy_mu[0]*Psisecder_lambda[2] + sy_mu[1]*Psisecder_lambda[5] + sy_mu[2]*Psisecder_lambda[8])*sz_mu[2];
                
                sz_Phisecder_sz = (sz_mu[0]*Phisecder_lambda[0] + sz_mu[1]*Phisecder_lambda[3] + sz_mu[2]*Phisecder_lambda[6])*sz_mu[0] \
                + (sz_mu[0]*Phisecder_lambda[1] + sz_mu[1]*Phisecder_lambda[4] + sz_mu[2]*Phisecder_lambda[7])*sz_mu[1] \
                + (sz_mu[0]*Phisecder_lambda[2] + sz_mu[1]*Phisecder_lambda[5] + sz_mu[2]*Phisecder_lambda[8])*sz_mu[2];
                
                sz_Psisecder_sz = (sz_mu[0]*Psisecder_lambda[0] + sz_mu[1]*Psisecder_lambda[3] + sz_mu[2]*Psisecder_lambda[6])*sz_mu[0] \
                + (sz_mu[0]*Psisecder_lambda[1] + sz_mu[1]*Psisecder_lambda[4] + sz_mu[2]*Psisecder_lambda[7])*sz_mu[1] \
                + (sz_mu[0]*Psisecder_lambda[2] + sz_mu[1]*Psisecder_lambda[5] + sz_mu[2]*Psisecder_lambda[8])*sz_mu[2];
                
                sz_mu_dot_sz_mu = sz_mu[0]*sz_mu[0]+ sz_mu[1]*sz_mu[1]+ sz_mu[2]*sz_mu[2];
                
                sy_Phisecder_k = (sy_mu[0]*Phisecder_lambda[0] + sy_mu[1]*Phisecder_lambda[3] + sy_mu[2]*Phisecder_lambda[6])*k_mu[0] \
                + (sy_mu[0]*Phisecder_lambda[1] + sy_mu[1]*Phisecder_lambda[4] + sy_mu[2]*Phisecder_lambda[7])*k_mu[1] \
                + (sy_mu[0]*Phisecder_lambda[2] + sy_mu[1]*Phisecder_lambda[5] + sy_mu[2]*Phisecder_lambda[8])*k_mu[2];
                
                sz_Phisecder_k = (sz_mu[0]*Phisecder_lambda[0] + sz_mu[1]*Phisecder_lambda[3] + sz_mu[2]*Phisecder_lambda[6])*k_mu[0] \
                + (sz_mu[0]*Phisecder_lambda[1] + sz_mu[1]*Phisecder_lambda[4] + sz_mu[2]*Phisecder_lambda[7])*k_mu[1] \
                + (sz_mu[0]*Phisecder_lambda[2] + sz_mu[1]*Phisecder_lambda[5] + sz_mu[2]*Phisecder_lambda[8])*k_mu[2];
                
                sy_Psisecder_k = (sy_mu[0]*Psisecder_lambda[0] + sy_mu[1]*Psisecder_lambda[3] + sy_mu[2]*Psisecder_lambda[6])*k_mu[0] \
                + (sy_mu[0]*Psisecder_lambda[1] + sy_mu[1]*Psisecder_lambda[4] + sy_mu[2]*Psisecder_lambda[7])*k_mu[1] \
                + (sy_mu[0]*Psisecder_lambda[2] + sy_mu[1]*Psisecder_lambda[5] + sy_mu[2]*Psisecder_lambda[8])*k_mu[2];
                
                sz_Psisecder_k = (sz_mu[0]*Psisecder_lambda[0] + sz_mu[1]*Psisecder_lambda[3] + sz_mu[2]*Psisecder_lambda[6])*k_mu[0] \
                + (sz_mu[0]*Psisecder_lambda[1] + sz_mu[1]*Psisecder_lambda[4] + sz_mu[2]*Psisecder_lambda[7])*k_mu[1] \
                + (sz_mu[0]*Psisecder_lambda[2] + sz_mu[1]*Psisecder_lambda[5] + sz_mu[2]*Psisecder_lambda[8])*k_mu[2];
                
                sy_mu_dot_sz_mu = sy_mu[0]*sz_mu[0]+ sy_mu[1]*sz_mu[1]+ sy_mu[2]*sz_mu[2];
                
                /*optical tidal matrix*/
                R_yy = sy_Phisecder_sy*(k_t*k_t)+(sy_Psisecder_sy)*(k_mu_dot_k_mu) + k_Psisecder_k *(sy_mu_dot_sy_mu) + k_Phisecder_k*(sy_t*sy_t) - 2.0*sy_Phisecder_k*sy_t*k_t - 2.0*sy_Psisecder_k*sy_mu_dot_k_mu;
                
                R_yz = sy_Phisecder_sz*(k_t*k_t)+sy_Psisecder_sz*(k_mu_dot_k_mu)+ k_Phisecder_k*(sy_t*sz_t) + k_Psisecder_k*(sy_mu_dot_sz_mu) - sy_Phisecder_k*sz_t*k_t - sz_Phisecder_k*sy_t*k_t - sy_Psisecder_k*sz_mu_dot_k_mu - sz_Psisecder_k*sy_mu_dot_k_mu ;
                
                
                R_zz = sz_Phisecder_sz*(k_t*k_t)+sz_Psisecder_sz*(k_mu_dot_k_mu) + k_Psisecder_k*(sz_mu_dot_sz_mu) + k_Phisecder_k*(sz_t*sz_t)- 2.0*sz_Phisecder_k*sz_t*k_t- 2.0*sz_Psisecder_k*sz_mu_dot_k_mu;
                
                //printf("Ryy %.15f \n Ryz %.15f \n R_zz %.15f \n", R_yy,R_yz, R_zz);
                
                /*Calculate Wronski matrix in each cell*/
                /*W = W + [zeros(2) eye(2); R_AB zeros(2)]*W.*dlam;*/
                
                W_save[0]= W[0] + W[8]*dlam_save;
                W_save[1]= W[1] + W[9]*dlam_save;
                W_save[2]= W[2] + W[10]*dlam_save;
                W_save[3]= W[3] + W[11]*dlam_save;
                W_save[4]= W[4] + W[12]*dlam_save;
                W_save[5]= W[5] + W[13]*dlam_save;
                W_save[6]= W[6] + W[14]*dlam_save;
                W_save[7]= W[7] + W[15]*dlam_save;
                W_save[8]= W[8] + (R_yy*W[0] + R_yz*W[4])*dlam_save;
                W_save[9]= W[9] + (R_yy*W[1] + R_yz*W[5])*dlam_save;
                W_save[10]= W[10] + (R_yy*W[2] + R_yz*W[6])*dlam_save;
                W_save[11]= W[11] + (R_yy*W[3] + R_yz*W[7])*dlam_save;
                W_save[12]= W[12] + (R_yz*W[0] + R_zz*W[4])*dlam_save;
                W_save[13]= W[13] + (R_yz*W[1] + R_zz*W[5])*dlam_save;
                W_save[14]= W[14] + (R_yz*W[2] + R_zz*W[6])*dlam_save;
                W_save[15]= W[15] + (R_yz*W[3] + R_zz*W[7])*dlam_save;
                
                lam_save = lam + dlam_save;
                
                t_save = t - k_t*dlam_save;
                
                X_save = X_0*exp(-sqrt(lam_div_3)*t_save);
                
                z_FRW = (X_initial/X_save) - 1.0;
                
                //Determinant of jacobi matrix
                det= (W_save[2]*W_save[7] - W_save[3]*W_save[6]);
                inv_den =1.0/(det);
                S[0] = inv_den*(W_save[10]*W_save[7] - W_save[11]*W_save[6]);
                S[1] = inv_den*(W_save[11]*W_save[2] - W_save[10]*W_save[3]);
                S[2] = inv_den*(W_save[14]*W_save[7] - W_save[15]*W_save[6]);
                S[3] = inv_den*(W_save[15]*W_save[2] - W_save[14]*W_save[3]);
                
                
                //z_diff = (totalz-1 - z_FRW)/z_FRW;
                D_A = omega_zero*sqrt(det);
                //D_A_flrw = save_z*(1.0/sqrt(lambda/3.0))/(1.0+save_z);
                //printf("value of D_A is %.15f\n", D_A[n2]);
                fprintf(f," %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %d \n", theta, varphi, save_z, z_FRW, t_save, lam_save, D_A, W_save[2], W_save[3], W_save[6], W_save[7], S[0], S[1], S[2], S[3], 1);
                save_z += 0.1;
                
                
            }
            
            //change in conformal time
            dt = dtau*X/X_0;
            //dt=0.005;
            
            /*dlam is the affine parameter*/
            dlam = -dt/k_t;
            
            /*Change position of light ray*/
            x_mu[0] += dlam*k_mu[0];
            x_mu[1] += dlam*k_mu[1];
            x_mu[2] += dlam*k_mu[2];
            
            /*Store new position in x,y,z*/
            x = x_mu[0];
            y = x_mu[1];
            z= x_mu[2];
            
            X_next = X_0*exp(-sqrt(lam_div_3)*(t+dt));
          
            
            if (fabs(x)>=X_next || fabs(y)>=X_next || fabs(z)>=X_next)
            {   //Save outside position of light ray
                x_out = x;
                y_out = y;
                z_out = z;
                
                if (x>0)
                {//x=X;
                    dt_x = (x_out - X_next)/(-k_mu[0]/k_t + H*X); }
                
                else if (x<0)
                {//x=-X;
                    dt_x = (x_out + X_next)/(-k_mu[0]/k_t - H*X);}
                
                if (y>0)
                {//x=X;
                    dt_y = (y_out - X_next)/(-k_mu[1]/k_t + H*X); }
                
                else if (y<0)
                {//x=-X;
                    dt_y = (y_out + X_next)/(-k_mu[1]/k_t - H*X);}
                
                if (z>0)
                {//z=X;
                    dt_z = (z_out - X_next)/(-k_mu[2]/k_t + H*X);}
                
                else if (z<0)
                {//z=-X;
                    dt_z = (z_out + X_next)/(-k_mu[2]/k_t - H*X);}
                
                X_x = X_0*exp(-sqrt(lam_div_3)*(t+dt - dt_x));
                X_y = X_0*exp(-sqrt(lam_div_3)*(t+dt - dt_y));
                X_z = X_0*exp(-sqrt(lam_div_3)*(t+dt - dt_z));
                
                
                /*Due to numerical error when light ray hits boundary, reposition light ray on the boundary
                 In the x,y and z direction*/
                if (fabs(x_mu[0])>=X_next && fabs(x_mu[1])>=X_next && fabs(x_mu[2])>=X_next){
                    
                    if (fabs(dt_x)>=fabs(dt_y) && fabs(dt_x)>=fabs(dt_z)){
                        
                        X=X_x;
                        
                        if (x>0)
                        {//x=X;
                            x_mu[0]=X;
                        }
                        
                        else if (x<0)
                        {//x=-X;
                            x_mu[0]=-X;}
                        
                        dt_corr = dt_x;
                        dlam_corr = -dt_x/k_t;
                        //dlam_corr = -fabs(dt_x);
                        x_mu[1] = x_mu[1] - dlam_corr*k_mu[1];
                        x_mu[2] = x_mu[2] - dlam_corr*k_mu[2];
                    }
                    
                    else if (fabs(dt_y)>=fabs(dt_x) && fabs(dt_y)>=fabs(dt_z)){
                        
                        X=X_y;
                        
                        if (y>0)
                        {//y=X;
                            x_mu[1]=X;
                        }
                        
                        else if (y<0)
                        {//y=-X;
                            x_mu[1]=-X;
                        }
                        
                        dt_corr = dt_y;
                        dlam_corr = -dt_y/k_t;
                        //dlam_corr = -fabs(dt_y);
                        x_mu[0] = x_mu[0] - dlam_corr*k_mu[0];
                        x_mu[2] = x_mu[2] - dlam_corr*k_mu[2];
                    }
                    
                    else if (fabs(dt_z)>=fabs(dt_x) && fabs(dt_z)>=fabs(dt_y)){
                        
                        X=X_z;
                        
                        if (z>0)
                        {//z=X;
                            x_mu[2]=X;
                        }
                        
                        else if (z<0)
                        {//z=-X;
                            x_mu[2]=-X;
                        }
                        
                        dt_corr = dt_z;
                        dlam_corr = -dt_z/k_t;
                        //dlam_corr = -fabs(dt_z);
                        x_mu[0] = x_mu[0] - dlam_corr*k_mu[0];
                        x_mu[1] = x_mu[1] - dlam_corr*k_mu[1];
                    }
                    
                }
                
                else if (fabs(x_mu[0])>=X_next && fabs(x_mu[1])>=X_next){
                    
                    if (fabs(dt_x)>=fabs(dt_y)){
                        
                        X=X_x;
                        
                        if (x>0)
                        {//x=X;
                            x_mu[0]=X;
                        }
                        
                        else if (x<0)
                        {//x=-X;
                            x_mu[0]=-X;}
                        
                        dt_corr = dt_x;
                        dlam_corr = -dt_x/k_t;
                        //dlam_corr = -fabs(dt_x);
                        x_mu[1] = x_mu[1] - dlam_corr*k_mu[1];
                        x_mu[2] = x_mu[2] - dlam_corr*k_mu[2];
                    }
                    
                    else if (fabs(dt_y)>fabs(dt_x)){
                        
                        X=X_y;
                        
                        if (y>0)
                        {//y=X;
                            x_mu[1]=X;
                        }
                        
                        else if (y<0)
                        {//y=-X;
                            x_mu[1]=-X;
                        }
                        
                        dt_corr = dt_y;
                        dlam_corr = -dt_y/k_t;
                        //dlam_corr = -fabs(dt_y);
                        x_mu[0] = x_mu[0] - dlam_corr*k_mu[0];
                        x_mu[2] = x_mu[2] - dlam_corr*k_mu[2];
                    }
                    
                    
                }
                
                else if (fabs(x_mu[0])>=X_next && fabs(x_mu[2])>=X_next){
                    
                    if (fabs(dt_x)>fabs(dt_z)){
                        
                        X=X_x;
                        
                        if (x>0)
                        {//x=X;
                            x_mu[0]=X;
                        }
                        
                        else if (x<0)
                        {//x=-X;
                            x_mu[0]=-X;}
                        
                        dt_corr = dt_x;
                        dlam_corr = -dt_x/k_t;
                        //dlam_corr = -fabs(dt_x);
                        x_mu[1] = x_mu[1] - dlam_corr*k_mu[1];
                        x_mu[2] = x_mu[2] - dlam_corr*k_mu[2];
                    }
                    
                    
                    else if (fabs(dt_z)>=fabs(dt_x)){
                        
                        X=X_z;
                        
                        if (z>0)
                        {//z=X;
                            x_mu[2]=X;
                        }
                        
                        else if (z<0)
                        {//z=-X;
                            x_mu[2]=-X;
                        }
                        
                        dt_corr = dt_z;
                        dlam_corr = -dt_z/k_t;
                        //dlam_corr = -fabs(dt_z);
                        x_mu[0] = x_mu[0] - dlam_corr*k_mu[0];
                        x_mu[1] = x_mu[1] - dlam_corr*k_mu[1];
                    }
                    
                }
                
                else if (fabs(x_mu[1])>=X_next && fabs(x_mu[2])>=X_next){
                    
                    
                    if (fabs(dt_y)>=fabs(dt_z)){
                        
                        X=X_y;
                        
                        if (y>0)
                        {//y=X;
                            x_mu[1]=X;
                        }
                        
                        else if (y<0)
                        {//y=-X;
                            x_mu[1]=-X;
                        }
                        
                        dt_corr = dt_y;
                        dlam_corr = -dt_y/k_t;
                        //dlam_corr = -fabs(dt_y);
                        x_mu[0] = x_mu[0] - dlam_corr*k_mu[0];
                        x_mu[2] = x_mu[2] - dlam_corr*k_mu[2];
                    }
                    
                    else if (fabs(dt_z)>fabs(dt_y)){
                        
                        X=X_z;
                        
                        if (z>0)
                        {//z=X;
                            x_mu[2]=X;
                        }
                        
                        else if (z<0)
                        {//z=-X;
                            x_mu[2]=-X;
                        }
                        
                        dt_corr = dt_z;
                        dlam_corr = -dt_z/k_t;
                        //dlam_corr = -fabs(dt_z);
                        x_mu[0] = x_mu[0] - dlam_corr*k_mu[0];
                        x_mu[1] = x_mu[1] - dlam_corr*k_mu[1];
                    }
                    
                }
                
                
                
                else if (fabs(x_mu[0])>=X_next){
                    
                    X=X_x;
                    
                    if (x>0)
                    {//x=X;
                        x_mu[0]=X;
                        /*printf("kx %.15f  \n", k_mu[0]);
                         printf("dlam_corr %.15f \n", dlam_corr);
                         printf("x_out %.15f \n", x_out);
                         printf("X %.15f \n", X);*/
                    }
                    
                    else if (x<0)
                    {//x=-X;
                        x_mu[0]=-X;}
                    
                    dt_corr = dt_x;
                    dlam_corr = -dt_x/k_t;
                    //dlam_corr = -fabs(dt_x);
                    x_mu[1] = x_mu[1] - dlam_corr*k_mu[1];
                    x_mu[2] = x_mu[2] - dlam_corr*k_mu[2];
                }
                
                else if (fabs(x_mu[1])>=X_next){
                    
                    X=X_y;
                    
                    if (y>0)
                    {//y=X;
                        x_mu[1]=X;
                    }
                    
                    else if (y<0)
                    {//y=-X;
                        x_mu[1]=-X;
                    }
                    
                    dt_corr = dt_y;
                    dlam_corr = -dt_y/k_t;
                    //dlam_corr = -fabs(dt_y);
                    x_mu[0] = x_mu[0] - dlam_corr*k_mu[0];
                    x_mu[2] = x_mu[2] - dlam_corr*k_mu[2];
                }
                
                else if (fabs(x_mu[2])>=X_next){
                    
                    X=X_z;
                    
                    if (z>0)
                    {//z=X;
                        x_mu[2]=X;
                    }
                    
                    else if (z<0)
                    {//z=-X;
                        x_mu[2]=-X;
                    }
                    
                    dt_corr = dt_z;
                    dlam_corr = -dt_z/k_t;
                    //dlam_corr = -fabs(dt_z);
                    x_mu[0] = x_mu[0] - dlam_corr*k_mu[0];
                    x_mu[1] = x_mu[1] - dlam_corr*k_mu[1];
                }
                
                /*x_mu[0] = x_mu[0] - dlam_corr*k_mu[0];
                 x_mu[1] = x_mu[1] - dlam_corr*k_mu[1];
                 x_mu[2] = x_mu[2] - dlam_corr*k_mu[2];*/
                //printf("dlam_corr %.15f  \n", dlam_corr);
                x = x_mu[0];
                y = x_mu[1];
                z= x_mu[2];
                
                //printf("dlam_corr %.15f \n", dlam_corr);
                
                /*k_mu_dot_k_mu = k_mu[0]*k_mu[0]+ k_mu[1]*k_mu[1]+ k_mu[2]*k_mu[2] ;
                 
                 dlam_corr = sqrt(((x_out-x)*(x_out-x)+ (y_out-y)*(y_out-y)+ (z_out-z)*(z_out-z))/(k_mu_dot_k_mu));*/
                
                lam = lam + dlam - dlam_corr;
                
                t = t + dt - dt_corr;
                
                //printf("dlam %.15f  \n", dlam);
                dlam  = dlam - dlam_corr;
                //printf("dlam %.15f  \n", dlam);
                /*if (dlam>0)
                 {
                 dlam = -dlam;}*/
                //n2 =n2 +1;
            }
            
            
            
            
            /*Dot products we require*/
            
            Phider_mu_dot_k_mu = Phider_mu[0]*k_mu[0] + Phider_mu[1]*k_mu[1] + Phider_mu[2]*k_mu[2];
            
            /*Dot products we require*/
            Phider_mu_dot_sy_mu = Phider_mu[0]*sy_mu[0] + Phider_mu[1]*sy_mu[1] + Phider_mu[2]*sy_mu[2];
            
            //du_dot_sy_a = (sy_a[0]*u_der[0] + sy_a[1]*u_der[4] + sy_a[2]*u_der[8] + sy_a[3]*u_der[12])*k_t + (sy_a[0]*u_der[1]*k_a[1] + sy_a[0]*u_der[2]*k_a[2] + sy_a[0]*u_der[3]*k_a[3]) + H*(sy_a[1]*k_a[1] + sy_a[2]*k_a[2] + sy_a[3]*k_a[3]) ;
            
            //du_dot_sz_a = (sz_a[0]*u_der[0] + sz_a[1]*u_der[4] + sz_a[2]*u_der[8] + sz_a[3]*u_der[12])*k_t + (sz_a[0]*u_der[1]*k_a[1] + sz_a[0]*u_der[2]*k_a[2] + sz_a[0]*u_der[3]*k_a[3]) + H*(sz_a[1]*k_a[1] + sz_a[2]*k_a[2] + sz_a[3]*k_a[3]) ;
            
            du_dot_sy_a = (sy_a[0]*u_der[0] + sy_a[1]*u_der[4] + sy_a[2]*u_der[8] + sy_a[3]*u_der[12])*k_t
            + (sy_a[0]*u_der[1] + sy_a[1]*u_der[5] + sy_a[2]*u_der[9] + sy_a[3]*u_der[13])*k_mu[0] + (sy_a[0]*u_der[2] + sy_a[1]*u_der[6] + sy_a[2]*u_der[10] + sy_a[3]*u_der[14])*k_mu[1] + (sy_a[0]*u_der[3] + sy_a[1]*u_der[7] + sy_a[2]*u_der[11] + sy_a[3]*u_der[15])*k_mu[2] ;
            
            du_dot_sz_a = (sz_a[0]*u_der[0] + sz_a[1]*u_der[4] + sz_a[2]*u_der[8] + sz_a[3]*u_der[12])*k_t
            + (sz_a[0]*u_der[1] + sz_a[1]*u_der[5] + sz_a[2]*u_der[9] + sz_a[3]*u_der[13])*k_mu[0] + (sz_a[0]*u_der[2] + sz_a[1]*u_der[6] + sz_a[2]*u_der[10] + sz_a[3]*u_der[14])*k_mu[1] + (sz_a[0]*u_der[3] + sz_a[1]*u_der[7] + sz_a[2]*u_der[11] + sz_a[3]*u_der[15])*k_mu[2] ;
            
            
            Phider_mu_u_undera_mu = Phider_mu[0]*u_undera[1] + Phider_mu[1]*u_undera[2] + Phider_mu[2]*u_undera[3];
            Psider_mu_u_undera_mu = Psider_mu[0]*u_undera[1] + Psider_mu[1]*u_undera[2] + Psider_mu[2]*u_undera[3];
            sy_mu_dot_sy_mu = sy_mu[0]*sy_mu[0] + sy_mu[1]*sy_mu[1] + sy_mu[2]*sy_mu[2];
            sy_mu_dot_k_mu = sy_mu[0]*k_mu[0] + sy_mu[1]*k_mu[1] + sy_mu[2]*k_mu[2];
            Psider_mu_dot_sy_mu = Psider_mu[0]*sy_mu[0] + Psider_mu[1]*sy_mu[1] + Psider_mu[2]*sy_mu[2];
            sy_mu_u_undera_mu = sy_mu[0]*u_undera[1] + sy_mu[1]*u_undera[2] + sy_mu[2]*u_undera[3];
            Phider_mu_dot_sz_mu = Phider_mu[0]*sz_mu[0] + Phider_mu[1]*sz_mu[1]+ Phider_mu[2]*sz_mu[2];
            
            sz_mu_dot_sz_mu = sz_mu[0]*sz_mu[0] + sz_mu[1]*sz_mu[1] + sz_mu[2]*sz_mu[2];
            sz_mu_dot_k_mu = sz_mu[0]*k_mu[0] + sz_mu[1]*k_mu[1] + sz_mu[2]*k_mu[2];
            Psider_mu_dot_sz_mu = Psider_mu[0]*sz_mu[0]+ Psider_mu[1]*sz_mu[1] + Psider_mu[2]*sz_mu[2];
            sz_mu_u_undera_mu = sz_mu[0]*u_undera[1] + sz_mu[1]*u_undera[2] + sz_mu[2]*u_undera[3];
            Psider_mu_dot_k_mu = Psider_mu[0]*k_mu[0] + Psider_mu[1]*k_mu[1] + Psider_mu[2]*k_mu[2];
            k_mu_u_undera_mu = k_mu[0]*u_undera[1] + k_mu[1]*u_undera[2] + k_mu[2]*u_undera[3];
            
            //sy_u_d_brackets = k_t*(Phider_mu_dot_sy_mu + H*sy_t);
            sy_u_d_brackets = du_dot_sy_a + (Phider_mu_u_undera_mu)*sy_t*k_t + (Phider_mu_dot_k_mu)*sy_t*u_undera[0] + (Phider_mu_dot_sy_mu)*k_t*u_undera[0] + (Psider_mu_u_undera_mu)*(sy_mu_dot_k_mu) - (Psider_mu_dot_k_mu)*(sy_mu_u_undera_mu) - (Psider_mu_dot_sy_mu)*(k_mu_u_undera_mu);
            
            //sz_u_d_brackets = k_t*(Phider_mu_dot_sz_mu + H*sz_t);
            sz_u_d_brackets = du_dot_sz_a + (Phider_mu_u_undera_mu)*sz_t*k_t + (Phider_mu_dot_k_mu)*sz_t*u_undera[0] + (Phider_mu_dot_sz_mu)*k_t*u_undera[0] + (Psider_mu_u_undera_mu)*(sz_mu_dot_k_mu) - (Psider_mu_dot_k_mu)*(sz_mu_u_undera_mu) - (Psider_mu_dot_sz_mu)*(k_mu_u_undera_mu);
            
            
            sy_Phisecder_sy = (sy_mu[0]*Phisecder_lambda[0] + sy_mu[1]*Phisecder_lambda[3] + sy_mu[2]*Phisecder_lambda[6])*sy_mu[0] \
            + (sy_mu[0]*Phisecder_lambda[1] + sy_mu[1]*Phisecder_lambda[4] + sy_mu[2]*Phisecder_lambda[7])*sy_mu[1] \
            + (sy_mu[0]*Phisecder_lambda[2] + sy_mu[1]*Phisecder_lambda[5] + sy_mu[2]*Phisecder_lambda[8])*sy_mu[2];
            
            k_mu_dot_k_mu = k_mu[0]*k_mu[0]+ k_mu[1]*k_mu[1]+ k_mu[2]*k_mu[2] ;
            
            sy_Psisecder_sy = (sy_mu[0]*Psisecder_lambda[0] + sy_mu[1]*Psisecder_lambda[3] + sy_mu[2]*Psisecder_lambda[6])*sy_mu[0] \
            + (sy_mu[0]*Psisecder_lambda[1] + sy_mu[1]*Psisecder_lambda[4] + sy_mu[2]*Psisecder_lambda[7])*sy_mu[1] \
            + (sy_mu[0]*Psisecder_lambda[2] + sy_mu[1]*Psisecder_lambda[5] + sy_mu[2]*Psisecder_lambda[8])*sy_mu[2];
            
            k_Psisecder_k = (k_mu[0]*Psisecder_lambda[0] + k_mu[1]*Psisecder_lambda[3] + k_mu[2]*Psisecder_lambda[6])*k_mu[0] \
            + (k_mu[0]*Psisecder_lambda[1] + k_mu[1]*Psisecder_lambda[4] + k_mu[2]*Psisecder_lambda[7])*k_mu[1] \
            + (k_mu[0]*Psisecder_lambda[2] + k_mu[1]*Psisecder_lambda[5] + k_mu[2]*Psisecder_lambda[8])*k_mu[2];
            
            k_Phisecder_k = (k_mu[0]*Phisecder_lambda[0] + k_mu[1]*Phisecder_lambda[3] + k_mu[2]*Phisecder_lambda[6])*k_mu[0] \
            + (k_mu[0]*Phisecder_lambda[1] + k_mu[1]*Phisecder_lambda[4] + k_mu[2]*Phisecder_lambda[7])*k_mu[1] \
            + (k_mu[0]*Phisecder_lambda[2] + k_mu[1]*Phisecder_lambda[5] + k_mu[2]*Phisecder_lambda[8])*k_mu[2];
            
            sy_mu_dot_sy_mu = sy_mu[0]*sy_mu[0]+ sy_mu[1]*sy_mu[1]+ sy_mu[2]*sy_mu[2];
            
            sy_Phisecder_sz = (sy_mu[0]*Phisecder_lambda[0] + sy_mu[1]*Phisecder_lambda[3] + sy_mu[2]*Phisecder_lambda[6])*sz_mu[0] \
            + (sy_mu[0]*Phisecder_lambda[1] + sy_mu[1]*Phisecder_lambda[4] + sy_mu[2]*Phisecder_lambda[7])*sz_mu[1] \
            + (sy_mu[0]*Phisecder_lambda[2] + sy_mu[1]*Phisecder_lambda[5] + sy_mu[2]*Phisecder_lambda[8])*sz_mu[2];
            
            sy_Psisecder_sz = (sy_mu[0]*Psisecder_lambda[0] + sy_mu[1]*Psisecder_lambda[3] + sy_mu[2]*Psisecder_lambda[6])*sz_mu[0] \
            + (sy_mu[0]*Psisecder_lambda[1] + sy_mu[1]*Psisecder_lambda[4] + sy_mu[2]*Psisecder_lambda[7])*sz_mu[1] \
            + (sy_mu[0]*Psisecder_lambda[2] + sy_mu[1]*Psisecder_lambda[5] + sy_mu[2]*Psisecder_lambda[8])*sz_mu[2];
            
            sz_Phisecder_sz = (sz_mu[0]*Phisecder_lambda[0] + sz_mu[1]*Phisecder_lambda[3] + sz_mu[2]*Phisecder_lambda[6])*sz_mu[0] \
            + (sz_mu[0]*Phisecder_lambda[1] + sz_mu[1]*Phisecder_lambda[4] + sz_mu[2]*Phisecder_lambda[7])*sz_mu[1] \
            + (sz_mu[0]*Phisecder_lambda[2] + sz_mu[1]*Phisecder_lambda[5] + sz_mu[2]*Phisecder_lambda[8])*sz_mu[2];
            
            sz_Psisecder_sz = (sz_mu[0]*Psisecder_lambda[0] + sz_mu[1]*Psisecder_lambda[3] + sz_mu[2]*Psisecder_lambda[6])*sz_mu[0] \
            + (sz_mu[0]*Psisecder_lambda[1] + sz_mu[1]*Psisecder_lambda[4] + sz_mu[2]*Psisecder_lambda[7])*sz_mu[1] \
            + (sz_mu[0]*Psisecder_lambda[2] + sz_mu[1]*Psisecder_lambda[5] + sz_mu[2]*Psisecder_lambda[8])*sz_mu[2];
            
            sz_mu_dot_sz_mu = sz_mu[0]*sz_mu[0]+ sz_mu[1]*sz_mu[1]+ sz_mu[2]*sz_mu[2];
            
            sy_Phisecder_k = (sy_mu[0]*Phisecder_lambda[0] + sy_mu[1]*Phisecder_lambda[3] + sy_mu[2]*Phisecder_lambda[6])*k_mu[0] \
            + (sy_mu[0]*Phisecder_lambda[1] + sy_mu[1]*Phisecder_lambda[4] + sy_mu[2]*Phisecder_lambda[7])*k_mu[1] \
            + (sy_mu[0]*Phisecder_lambda[2] + sy_mu[1]*Phisecder_lambda[5] + sy_mu[2]*Phisecder_lambda[8])*k_mu[2];
            
            sz_Phisecder_k = (sz_mu[0]*Phisecder_lambda[0] + sz_mu[1]*Phisecder_lambda[3] + sz_mu[2]*Phisecder_lambda[6])*k_mu[0] \
            + (sz_mu[0]*Phisecder_lambda[1] + sz_mu[1]*Phisecder_lambda[4] + sz_mu[2]*Phisecder_lambda[7])*k_mu[1] \
            + (sz_mu[0]*Phisecder_lambda[2] + sz_mu[1]*Phisecder_lambda[5] + sz_mu[2]*Phisecder_lambda[8])*k_mu[2];
            
            sy_Psisecder_k = (sy_mu[0]*Psisecder_lambda[0] + sy_mu[1]*Psisecder_lambda[3] + sy_mu[2]*Psisecder_lambda[6])*k_mu[0] \
            + (sy_mu[0]*Psisecder_lambda[1] + sy_mu[1]*Psisecder_lambda[4] + sy_mu[2]*Psisecder_lambda[7])*k_mu[1] \
            + (sy_mu[0]*Psisecder_lambda[2] + sy_mu[1]*Psisecder_lambda[5] + sy_mu[2]*Psisecder_lambda[8])*k_mu[2];
            
            sz_Psisecder_k = (sz_mu[0]*Psisecder_lambda[0] + sz_mu[1]*Psisecder_lambda[3] + sz_mu[2]*Psisecder_lambda[6])*k_mu[0] \
            + (sz_mu[0]*Psisecder_lambda[1] + sz_mu[1]*Psisecder_lambda[4] + sz_mu[2]*Psisecder_lambda[7])*k_mu[1] \
            + (sz_mu[0]*Psisecder_lambda[2] + sz_mu[1]*Psisecder_lambda[5] + sz_mu[2]*Psisecder_lambda[8])*k_mu[2];
            
            sy_mu_dot_sz_mu = sy_mu[0]*sz_mu[0]+ sy_mu[1]*sz_mu[1]+ sy_mu[2]*sz_mu[2];
            
            /*optical tidal matrix*/
            R_yy = sy_Phisecder_sy*(k_t*k_t)+(sy_Psisecder_sy)*(k_mu_dot_k_mu) + k_Psisecder_k *(sy_mu_dot_sy_mu) + k_Phisecder_k*(sy_t*sy_t) - 2.0*sy_Phisecder_k*sy_t*k_t - 2.0*sy_Psisecder_k*sy_mu_dot_k_mu;
            
            R_yz = sy_Phisecder_sz*(k_t*k_t)+sy_Psisecder_sz*(k_mu_dot_k_mu)+ k_Phisecder_k*(sy_t*sz_t) + k_Psisecder_k*(sy_mu_dot_sz_mu) - sy_Phisecder_k*sz_t*k_t - sz_Phisecder_k*sy_t*k_t - sy_Psisecder_k*sz_mu_dot_k_mu - sz_Psisecder_k*sy_mu_dot_k_mu ;
            
            
            R_zz = sz_Phisecder_sz*(k_t*k_t)+sz_Psisecder_sz*(k_mu_dot_k_mu) + k_Psisecder_k*(sz_mu_dot_sz_mu) + k_Phisecder_k*(sz_t*sz_t)- 2.0*sz_Phisecder_k*sz_t*k_t- 2.0*sz_Psisecder_k*sz_mu_dot_k_mu;
            
            
            
            
            /*Calculate Wronski matrix in each cell*/
            /*W = W + [zeros(2) eye(2); R_AB zeros(2)]*W.*dlam;*/
            
            W[0]= W[0] + W[8]*dlam;
            W[1]= W[1] + W[9]*dlam;
            W[2]= W[2] + W[10]*dlam;
            W[3]= W[3] + W[11]*dlam;
            W[4]= W[4] + W[12]*dlam;
            W[5]= W[5] + W[13]*dlam;
            W[6]= W[6] + W[14]*dlam;
            W[7]= W[7] + W[15]*dlam;
            W[8]= W[8] + (R_yy*W[0] + R_yz*W[4])*dlam;
            W[9]= W[9] + (R_yy*W[1] + R_yz*W[5])*dlam;
            W[10]= W[10] + (R_yy*W[2] + R_yz*W[6])*dlam;
            W[11]= W[11] + (R_yy*W[3] + R_yz*W[7])*dlam;
            W[12]= W[12] + (R_yz*W[0] + R_zz*W[4])*dlam;
            W[13]= W[13] + (R_yz*W[1] + R_zz*W[5])*dlam;
            W[14]= W[14] + (R_yz*W[2] + R_zz*W[6])*dlam;
            W[15]= W[15] + (R_yz*W[3] + R_zz*W[7])*dlam;
            
            //Determinant of jacobi matrix
            det= (W[2]*W[7] - W[3]*W[6]);
            
            
            /*Evolve Sachs basis by parallely transporting it along the geodesic*/
            sy_t = sy_t + (Phider_mu_dot_sy_mu)*dlam*k_t + (Phider_mu_dot_k_mu)*dlam*sy_t + dlam*(k_t/omega)*(sy_u_d_brackets);
            
            sz_t = sz_t + (Phider_mu_dot_sz_mu)*dlam*k_t + (Phider_mu_dot_k_mu)*dlam*sz_t + dlam*(k_t/omega)*(sz_u_d_brackets);
            
            sy_mu[0] = sy_mu[0] + (-Psider_mu_dot_k_mu*sy_mu[0] - Psider_mu_dot_sy_mu*k_mu[0])*dlam + (Phider_mu[0]*sy_t)*(dlam*k_t) + Psider_mu[0]*sy_mu_dot_k_mu*(dlam)+ dlam*(k_mu[0]/omega)*(sy_u_d_brackets);
            
            sy_mu[1] = sy_mu[1] + (-Psider_mu_dot_k_mu*sy_mu[1] - Psider_mu_dot_sy_mu*k_mu[1])*dlam + (Phider_mu[1]*sy_t)*(dlam*k_t) + Psider_mu[1]*sy_mu_dot_k_mu*(dlam)+ dlam*(k_mu[1]/omega)*(sy_u_d_brackets);
            
            sy_mu[2] = sy_mu[2] + (-Psider_mu_dot_k_mu*sy_mu[2] - Psider_mu_dot_sy_mu*k_mu[2])*dlam + (Phider_mu[2]*sy_t)*(dlam*k_t) + Psider_mu[2]*sy_mu_dot_k_mu*(dlam) + dlam*(k_mu[2]/omega)*(sy_u_d_brackets);
            
            sz_mu[0] = sz_mu[0] + (-Psider_mu_dot_k_mu*sz_mu[0] - Psider_mu_dot_sz_mu*k_mu[0])*dlam + (Phider_mu[0]*sz_t)*(dlam*k_t) + Psider_mu[0]*sz_mu_dot_k_mu*(dlam)+ dlam*(k_mu[0]/omega)*(sz_u_d_brackets);
            
            sz_mu[1] = sz_mu[1] + (-Psider_mu_dot_k_mu*sz_mu[1] - Psider_mu_dot_sz_mu*k_mu[1])*dlam + (Phider_mu[1]*sz_t)*(dlam*k_t) + Psider_mu[1]*sz_mu_dot_k_mu*(dlam) + dlam*(k_mu[1]/omega)*(sz_u_d_brackets);
            
            sz_mu[2] = sz_mu[2] + (-Psider_mu_dot_k_mu*sz_mu[2] - Psider_mu_dot_sz_mu*k_mu[2])*dlam + (Phider_mu[2]*sz_t)*(dlam*k_t) + Psider_mu[2]*sz_mu_dot_k_mu*(dlam)+ dlam*(k_mu[2]/omega)*(sz_u_d_brackets);
            
            //Save values of sachs basis in 4-vector
            sy_a[0] = sy_t;
            sy_a[1] = sy_mu[0];
            sy_a[2] = sy_mu[1];
            sy_a[3] = sy_mu[2];
            
            sz_a[0] = sz_t;
            sz_a[1] = sz_mu[0];
            sz_a[2] = sz_mu[1];
            sz_a[3] = sz_mu[2];
            
            //For test of code
            /*s_change_t = (sy_a[0]-sy_a_prev[0])/dlam;
             s_change_mu[0] = (sy_mu[0]-sy_a_prev[1])/dlam;
             s_change_mu[1] = (sy_a[2]-sy_a_prev[2])/dlam;
             s_change_mu[2] = (sy_a[3]-sy_a_prev[3])/dlam;
             
             LHS_y = ((sy_a[0]-sy_a_prev[0])/dlam)*u_undera[0] + ((sy_a[1]-sy_a_prev[1])/dlam)*u_undera[1] + ((sy_a[2]-sy_a_prev[2])/dlam)*u_undera[2] + ((sy_a[3]-sy_a_prev[3])/dlam)*u_undera[3];
             
             LHS_z = ((sz_a[0]-sz_a_prev[0])/dlam)*u_undera[0] + ((sz_a[1]-sz_a_prev[1])/dlam)*u_undera[1] + ((sz_a[2]-sz_a_prev[2])/dlam)*u_undera[2] + ((sz_a[3]-sz_a_prev[3])/dlam)*u_undera[3];
             
             sy_a_prev[0] = sy_t;
             sy_a_prev[1] = sy_mu[0];
             sy_a_prev[2] = sy_mu[1];
             sy_a_prev[3] = sy_mu[2];
             
             sz_a_prev[0] = sz_t;
             sz_a_prev[1] = sz_mu[0];
             sz_a_prev[2] = sz_mu[1];
             sz_a_prev[3] = sz_mu[2];*/
            
            
            /*Evolve wavevector of light using the geodesic equation*/
            k_t = k_t + 2.0*(Phider_mu_dot_k_mu)*dlam*k_t;
            
            
            k_mu_dot_k_mu = k_mu[0]*k_mu[0]+ k_mu[1]*k_mu[1]+ k_mu[2]*k_mu[2] ;
            
            k_mu[0] = k_mu[0] + (Phider_mu[0]*k_t*k_t + Psider_mu[0]*(k_mu_dot_k_mu) - 2.0*(Psider_mu_dot_k_mu)*k_mu[0])*dlam;
            k_mu[1] = k_mu[1] + (Phider_mu[1]*k_t*k_t + Psider_mu[1]*(k_mu_dot_k_mu) - 2.0*(Psider_mu_dot_k_mu)*k_mu[1])*dlam;
            k_mu[2] = k_mu[2] + (Phider_mu[2]*k_t*k_t + Psider_mu[2]*(k_mu_dot_k_mu) - 2.0*(Psider_mu_dot_k_mu)*k_mu[2])*dlam;
            
            //Save values of wavevector in 4-vector
            k_a[0] = k_t;
            k_a[1] = k_mu[0];
            k_a[2] = k_mu[1];
            k_a[3] = k_mu[2];
            
            //Check if we've reached a caustic and stop code
            if (det <0) {
                z_FRW = (X_initial/X) - 1;
                inv_den =1.0/(W[2]*W[7] - W[3]*W[6]);
                S[0] = inv_den*(W[10]*W[7] - W[11]*W[6]);
                S[1] = inv_den*(W[11]*W[2] - W[10]*W[3]);
                S[2] = inv_den*(W[14]*W[7] - W[15]*W[6]);
                S[3] = inv_den*(W[15]*W[2] - W[14]*W[3]);
                omega = -u_undera_dot_k_a;
                totalz = totalz*omega/omega_prev;
                /*Evaluate angular diameter distance*/
                D_A = omega_zero*sqrt(det);
                printf(" %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %.15f \t %d \n", theta, varphi, totalz-1, z_FRW, t, lam, D_A, W[2], W[3], W[6], W[7], S[0], S[1], S[2], S[3],0);
                //printf(" caus %.15f \t %.15f \t %.15f \n", D_A[n2], totalz-1, det);
                break;}
            
            
            //Total number of steps
            n = n +1;
            
            
            /*Reflection code if light hits boundary*/
            /*if (x >=X || -x >=X || y>=X || -y >=X || z>=X || -z>=X)*/
            if (fabs(x)>=X || fabs(y)>=X || fabs(z)>=X)
            {
                
                /*Deformation rate matrix
                 S = ((D_new - D_prev)./dlam)/(D_new);
                 S = W(3:4,3:4)*inv((W(1:2,3:4)));*/
                
                X_tder=X*sqrt(lam_div_3);
                H= X_tder/X;
                H_sq = H*H;
                /*we need the potentials again*/
                r_sq = x*x + y*y + z*z;
                /*Potentials Phi and Psi for a universe with only a cosmological constant*/
                Phi = (lam_div_6)*r_sq;
                Psi = -(lam_div_12)*r_sq;
                
                /*Choosing the normal depending on whether we need to reflect in
                 %the x,y, or z direction*/
                if (fabs(x_mu[0])>=X)
                {/*Normal with index raised*/
                    if (x>0)
                    {sign = 1.0; }
                    
                    else if (x<0)
                    {sign =-1.0;}
                    
                    n_a[0] = sign*H*x;
                    n_a[1]= sign*(1.0-Psi+(H_sq)*(x*x)/2.0);
                    n_a[2]= 0;
                    n_a[3]= 0;
                    /*Normal with index lowered*/
                    n_undera[0] = -sign*H*x;
                    n_undera[1] = sign*(1.0+Psi+(H_sq)*(x*x)/2.0);
                    n_undera[2] = 0;
                    n_undera[3] = 0;
                    n3=n3+1;}
                else if (fabs(x_mu[1])>=X)
                {/*Normal with index raised*/
                    if (y>0)
                    {sign = 1.0; }
                    
                    else if (y<0)
                    {sign =-1.0;}
                    
                    n_a[0] = sign*H*y;
                    n_a[1]= 0;
                    n_a[2]= sign*(1.0-Psi+(H_sq)*(y*y)/2.0);
                    n_a[3]= 0;
                    /*Normal with index lowered*/
                    n_undera[0] = -sign*H*y;
                    n_undera[1] = 0;
                    n_undera[2] = sign*(1.0+Psi+(H_sq)*(y*y)/2.0);
                    n_undera[3] = 0;
                    n4=n4+1;}
                else if (fabs(x_mu[2])>=X)
                {/*Normal with index raised*/
                    if (z>0)
                    {sign = 1.0; }
                    
                    else if (z<0)
                    {sign =-1.0;}
                    
                    n_a[0] = sign*H*z;
                    n_a[1]= 0;
                    n_a[2]= 0;
                    n_a[3]= sign*(1.0-Psi+(H_sq)*(z*z)/2.0);
                    /*Normal with index lowered*/
                    n_undera[0] = -sign*H*z;
                    n_undera[1] = 0;
                    n_undera[2] = 0;
                    n_undera[3] = sign*(1.0+Psi+(H_sq)*(z*z)/2.0);
                    n5=n5+1;}
                
                
                
                /*Four velocity with index raised*/
                u_a[0] = (1.0 + H_sq*r_sq/2.0 + Phi);
                u_a[1]=H*x;
                u_a[2]=H*y;
                u_a[3]=H*z;
                
                
                /*Initial four velocity with index lowered*/
                u_undera[0] = -1.0 - H_sq*r_sq/2.0 + Phi;
                u_undera[1]=u_a[1];
                u_undera[2]=u_a[2];
                u_undera[3]=u_a[3];
                
                u_undera_dot_k_a = u_undera[0]*k_a[0] + u_undera[1]*k_a[1] + u_undera[2]*k_a[2]+ u_undera[3]*k_a[3];
                
                //For testing code
                /*d_a[0] = k_a[0]/(-u_undera_dot_k_a) - u_a[0];
                 d_a[1] = k_a[1]/(-u_undera_dot_k_a) - u_a[1];
                 d_a[2] = k_a[2]/(-u_undera_dot_k_a) - u_a[2];
                 d_a[3] = k_a[3]/(-u_undera_dot_k_a) - u_a[3];*/
                
                
                /*Store frequency of incoming light ray*/
                omega = -u_undera_dot_k_a;
                totalz = omega;
                //omega_prev = -u_undera_dot_k_a;
                
                k_a_dot_n_undera = n_undera[0]*k_a[0] + n_undera[1]*k_a[1] + n_undera[2]*k_a[2]+ n_undera[3]*k_a[3];
                sy_a_dot_n_undera = n_undera[0]*sy_a[0] + n_undera[1]*sy_a[1] + n_undera[2]*sy_a[2]+ n_undera[3]*sy_a[3];
                sz_a_dot_n_undera = n_undera[0]*sz_a[0] + n_undera[1]*sz_a[1] + n_undera[2]*sz_a[2]+ n_undera[3]*sz_a[3];
                
                
                /*Reflection vectors*/
                k_a[0] = k_a[0] - 2.0*(k_a_dot_n_undera)*n_a[0];
                k_a[1] = k_a[1] - 2.0*(k_a_dot_n_undera)*n_a[1];
                k_a[2] = k_a[2] - 2.0*(k_a_dot_n_undera)*n_a[2];
                k_a[3] = k_a[3] - 2.0*(k_a_dot_n_undera)*n_a[3];
                
                sy_a[0] = sy_a[0] - 2.0*(sy_a_dot_n_undera)*n_a[0];
                sy_a[1] = sy_a[1] - 2.0*(sy_a_dot_n_undera)*n_a[1];
                sy_a[2] = sy_a[2] - 2.0*(sy_a_dot_n_undera)*n_a[2];
                sy_a[3] = sy_a[3] - 2.0*(sy_a_dot_n_undera)*n_a[3];
                
                sz_a[0] = sz_a[0] - 2.0*(sz_a_dot_n_undera)*n_a[0];
                sz_a[1] = sz_a[1] - 2.0*(sz_a_dot_n_undera)*n_a[1];
                sz_a[2] = sz_a[2] - 2.0*(sz_a_dot_n_undera)*n_a[2];
                sz_a[3] = sz_a[3] - 2.0*(sz_a_dot_n_undera)*n_a[3];
                
                //Test code after reflection
                
                /*printf("check1 %.15f \n",k_a_dot_k_a);
                 printf("check2 %.15f \n",sy_a_dot_sy_a);
                 printf("check3 %.15f \n",sz_a_dot_sz_a);
                 printf("check4 %.15f \n",u_undera_dot_sy_a);
                 printf("check5 %.15f \n",u_undera_dot_sz_a);
                 printf("check6 %.15f \n",u_a_dot_n_undera);
                 printf("check7 %.15f \n",u_undera_dot_u_a);
                 printf("check8 %.15f \n",u_undera_dot_d_a);
                 printf("check9 %.15f \n",d_a_dot_sy_a);
                 printf("check10 %.15f \n",sy_a_dot_k_a);
                 printf("check11 %.15f \n",sz_a_dot_k_a);
                 printf("LHS_y %.15f \n", LHS_y);
                 printf ("RHS_y %.15f \n", du_dot_sy_a);
                 printf("y_sum %.15f \n", LHS_y + du_dot_sy_a);
                 
                 printf("LHS_z %.15f \n", LHS_z);
                 printf ("RHS_z %.15f \n", du_dot_sz_a);
                 printf("z_sum %.15f \n", LHS_z + du_dot_sz_a);
                 
                 printf("s_change_t %.15f \n", s_change_t);
                 printf ("s_change_mu[0] %.15f \n", s_change_mu[0]);
                 printf ("s_change_mu[1] %.15f \n", s_change_mu[1]);
                 printf ("s_change_mu[2] %.15f \n", s_change_mu[2]);*/
                
                
                
                /*Store frequency of outgoing light ray
                 Out(1,n2+1) = -u_undera*transpose(k_a);*/
                omega = -u_undera_dot_k_a;
                
                
                /*Evaluate angular diameter distance*/
                //D_A[n2] = omega_zero*sqrt(W[2]*W[7] - W[3]*W[6]);
                
                //fprintf(f," %.15f \t %.15f \n", D_A[n2], in[n2+1]);
                
                /*fwrite(D_A, sizeof(double), sizeof(D_A),f);*/
                
                /*Store expansion scalar;*/
                /*Expansion[n2] = -(S[0]*S[2])/2.0;*/
                
                /*Store shear rate*/
                /*Shear[n2] = sqrt((S[0]*S[0] + S[1]*S[2])*(S[1]*S[2] + S[3]*S[3])/2.0 - (S[0]*S[2]/2.0)*(S[0]*S[2]/2.0));*/
                
                /*Store new values of wavevector*/
                k_t = k_a[0];
                k_mu[0] = k_a[1];
                k_mu[1] = k_a[2];
                k_mu[2] = k_a[3];
                
                //k_mu_dot_x_mu = k_mu[0]*x_mu[0] + k_mu[1]*x_mu[1]+ k_mu[2]*x_mu[2];
                
                /*Store new values of sachs basis*/
                sy_t = sy_a[0];
                sy_mu[0] = sy_a[1];
                sy_mu[1] = sy_a[2];
                sy_mu[2] = sy_a[3];
                
                
                sz_t = sz_a[0];
                sz_mu[0] = sz_a[1];
                sz_mu[1] = sz_a[2];
                sz_mu[2] = sz_a[3];
                
                //Reflection counter
                n2 = n2+1;
                
              
                
                //printf("REFLECTION \n ");
                continue;
                
                
                
            }
            /*Increase time step and affine parameter*/
            t = t+ dt;
            lam += dlam;
            
            
        }
        
        printf("n2 %d  \n", n2);
        printf(" %.15f \n", totalzcheck-1.0);
        
        
    }
    
    
    
    fclose(f);
    
    //Time code
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time taken %.15f seconds \n", time_spent);
    
    
    return 0;
    
    
    
    
    
    
}







