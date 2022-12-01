
#ifdef RADIATIVE_2D_COOLING

#ifdef GPU
#define HOST   __host__
#define DEVICE __device__
#else
#define HOST  
#define DEVICE
#endif

// Bell & Lin type opacity from Zhaohuan Zhu (Zhu et al. 2012)
HOST DEVICE
static inline real opac_Zhu2012(real rho, real T)
{

  // Convert rho / T to cgs
  rho *= (MSTAR_CGS/MSTAR) / ((R0_CGS/R0)*(R0_CGS/R0)*(R0_CGS/R0)) ;
  T   *= (G_CGS/G) * (MSTAR_CGS/MSTAR) / ((R0_CGS/R0)*(R_MU_CGS/R_MU)) ; 

  if (T < 1) T = 1 ;

  real xlop, xlp, xlt, kappa,pre;
  pre=rho*T*R_MU_CGS;

  if (pre < 0. || T < 0.) {
#ifndef GPU
    fprintf(stderr, "error: pre or T negative\n");
    fprintf(stderr, "pre: %g T: %g\n", pre, T);
    prs_exit (1);
#endif
    return (1.);
  }
  if (pre == 0 || T == 0)
    return (1.);

  xlp = log10(pre);
  xlt = log10(T);

  if(xlt<2.23567+0.01899*(xlp-5.)){
    xlop=1.5*(xlt-1.16331)-0.736364;
  }
  else if (xlt<2.30713+0.01899*(xlp-5.)){
    xlop=-3.53154212*xlt+8.767726-(7.24786-8.767726)*(xlp-5.)/16.;
  }
  else if (xlt<2.79055){
    xlop=1.5*(xlt-2.30713)+0.62 ;
  }
  else if  (xlt<2.96931){
    xlop=-5.832*xlt+17.7;
  }
  else if (xlt<3.29105+(3.29105-3.07651)*(xlp-5.)/8.){
    xlop=2.129*xlt-5.9398;
  }
  else if (xlt<3.08+0.028084*(xlp+4)){
    xlop=129.88071-42.98075*xlt+(142.996475-129.88071)*0.1*(xlp+4);
  }
  else if (xlt<3.28+xlp/4.*0.12){
    xlop=-15.0125+4.0625*xlt;
  }
  else if (xlt<3.41+0.03328*xlp/4.){
    xlop=58.9294-18.4808*xlt+(61.6346-58.9294)*xlp/4.;
  }
  else if (xlt<3.76+(xlp-4)/2.*0.03){
    xlop=-12.002+2.90477*xlt+(xlp-4)/4.*(13.9953-12.002);
  }
  else if (xlt<4.07+(xlp-4)/2.*0.08){
    xlop=-39.4077+10.1935*xlt+(xlp-4)/2.*(40.1719-39.4077);
  }
  else if (xlt<5.3715+(xlp-6)/2.*0.5594){
    xlop=17.5935-3.3647*xlt+(xlp-6)/2.*(17.5935-15.7376);
  }
  else{
    xlop=-0.48;
  }
  
  if (xlop<3.586*xlt-16.85&&xlt<4.){
    xlop=3.586*xlt-16.85;
  }

  kappa=pow(10.,xlop);

  // convert kappa back to code units
  kappa *= (MSTAR_CGS/MSTAR) / ((R0_CGS/R0)*(R0_CGS/R0)) ;

  return(kappa);
}

#endif
