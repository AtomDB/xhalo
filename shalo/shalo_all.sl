variable Th_cm, dT_cm, dust_cm, HF_cm;
static variable Isca_calc_AA = Assoc_Type [];
static variable Isca_calc;
static variable E_cm=0;
static variable dust;

static variable IntDM, IntEn, IntHF;

variable EPS=6.0e-8;
variable EULER=0.57721566;
variable MAXIT=100;
variable FPMIN=1.0e-30;
variable TRUE=1;

define siiSingle( x ) {
   %
   %  Calculates SinIntegral(x), Si(x) = \int_0^x sin(t)/t dt
   % 
   variable i,k,odd;
   variable a,err,fact,sign1,sum,sumc,sums,t,term;
   variable h,b,c,d,del; % complex 

   t=abs(x);
   if (t == 0.0) {
      return 0.0;
   }
   if (t > 2.0) {
      b=1.0 + t*(1i); 
      c=1.0/FPMIN;
      d = 1.0/b;
      h = 1.0/b;
      i = 1;
      del = 2.0;
      while (abs(del-1.0) > EPS) {
	 i=i+1;
	 a = -(i-1)*(i-1);
	 b = b + 2;
	 d = 1./(a*d + b);
	 c =b + a/c; 
	 del= c * d;
	 h= h*del;
	 if (i > MAXIT) return NULL;
      }
      h = ( cos(t) -(1i)*sin(t))*h;
      if (x >= 0.0) {
	 return PI/2.0 + Imag(h);
      } else {
	 return -(PI/2.0 + Imag(h));
      }
   } else {
      if (t < sqrt(FPMIN)) {
	 if (x >= 0.0) {
	    return t;
	 } else {
	    return -t;
	 }
      } else {
	 sum=0.0;
	 sums=0.0;
	 sumc=0.0;
	 sign1=1.0;
	 fact=1.0;
	 odd=TRUE;
	 k = 0;
	 err = 1.0;
	 while (err > EPS) {
	    k = k+1;
	    fact = fact*t/k;
	    term=fact/k;
	    sum = sum + sign1*term;
	    err=term/abs(sum);
	    if (odd) {
	       sign1 = -sign1;
	       sums=sum;
	       sum=sumc;
	    } else {
	       sumc=sum;
	       sum=sums;
	    }
	    if (k > MAXIT) return NULL;
	    odd = not odd;
	 }
      }
      if (x < 0.0 ) {
	 return -sums;
      } else {
	 return sums;
      }
   }
   return NULL;
}

define sii( x ) {
   
   variable result, i;
   
   if (length(x) == 1) {
      result = siiSingle(x);
   } else {
      result = Double_Type[length(x)];
      for (i=0;i<length(x);i++) result[i] = siiSingle(x[i]);
   }

   return result;
}

define csevl(x, cs, n) {
   % april 1977 version.  w. fullerton, c3, los alamos scientific lab.
   %
   % evaluate the n-term chebyshev series cs at x.  adapted from
   % r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).  also see fox
   % and parker, chebyshev polys in numerical analysis, oxford press, p.56.
   %
   %             input arguments --
   % x      value at which the series is to be evaluated.
   % cs     array of n terms of a chebyshev series.  in eval-
   %        uating cs, only half the first coef is summed.
   % n      number of terms in array cs.
   %
   % Ported to S-lang 6/25/03 by Randall K Smith;

   if (n < 1) {
      message("csevl: number of terms le 0");
      return NULL;
   }
   if (n > 1000) {
      message("csevl: number of terms gt 1000");
      return NULL;
   }
   
   if ((x < -1.1) or  (x > 1.1)) {
      message("csevl: x outside (-1,+1)");
      return NULL;
   }
   variable i, ni, b2;
   variable b1 = 0.;
   variable b0 = 0.;
   variable twox = 2.*x;
   for (i=1;i<n+1;i++) {
      b2 = b1;
      b1 = b0;
      ni = n - i;
      b0 = twox*b1 - b2 + cs[ni];
   }

   return 0.5 * (b0-b2);
}

define besj1(x) {
   % sept 1983 edition.  w. fullerton, c3, los alamos scientific lab.
   %
   % series for bj1        on the interval  0.          to  1.60000d+01
   %                                        with weighted error   4.48e-17
   %                                         log weighted error  16.35
   %                               significant figures required  15.77
   %                                    decimal places required  16.89
   %
   
   variable bj1cs =[-.11726141513332787,-.25361521830790640,
		    .050127080984469569 ,    -.004631514809625081 ,
		    .000247996229415914 ,    -.000008678948686278 ,
		    .000000214293917143 ,    -.000000003936093079 ,
		    .000000000055911823 ,    -.000000000000632761 ,
		    .000000000000005840 ,    -.000000000000000044 ];
   %
   % series for bm1        on the interval  0.          to  6.25000d-02
   %                                        with weighted error   5.61e-17
   %                                         log weighted error  16.25
   %                               significant figures required  14.97
   %                                    decimal places required  16.91
   %
   variable bm1cs =[.1047362510931285,
		    .00442443893702345,      -.00005661639504035,
		    .00000231349417339,      -.00000017377182007,
		    .00000001893209930,      -.00000000265416023,
		    .00000000044740209,      -.00000000008691795,
		    .00000000001891492,      -.00000000000451884,
		    .00000000000116765,      -.00000000000032265,
		    .00000000000009450,      -.00000000000002913,
		    .00000000000000939,      -.00000000000000315,
		    .00000000000000109,      -.00000000000000039,
		    .00000000000000014,      -.00000000000000005];
   %
   % series for bth1       on the interval  0.          to  6.25000d-02
   %                                        with weighted error   4.10e-17
   %                                         log weighted error  16.39
   %                               significant figures required  15.96
   %                                    decimal places required  17.08
   %
   variable bth1cs =[ .74060141026313850,-.004571755659637690,
		     .000119818510964326,      -.000006964561891648,
		     .000000655495621447,      -.000000084066228945,
		     .000000013376886564,      -.000000002499565654,
		     .000000000529495100,      -.000000000124135944,
		     .000000000031656485,      -.000000000008668640,
		     .000000000002523758,      -.000000000000775085,
		     .000000000000249527,      -.000000000000083773,
		     .000000000000029205,      -.000000000000010534,
		     .000000000000003919,      -.000000000000001500,
		     .000000000000000589,      -.000000000000000237,
		     .000000000000000097,      -.000000000000000040];
     %
   variable pi4=0.78539816339744831;
   variable ntj1 = 12;
   variable ntm1 = 21;
   variable ntth1 = 24;
   variable xsml = sqrt (8.0*2.22E-16);
   variable xmin = 2.0*2.22e-308;
   variable xmax = 1.0/2.22E-16;
   variable result;

   variable y = abs(x);
   if (y < 4.0)  {
      result = 0.0;
      if (y==0.0) return result;
      if (y<=xmin) {
	 message("besj1: abs(x) so small j1 underflows");
	 return result;
      }
      if (y > xsml) return x * (.25 + csevl(.125*y*y-1., bj1cs, ntj1));
      if (y > xmin) return 0.5*x;
   }
%
   if (y>xmax) message("besj1: no precision because abs(x) is big");
   variable z = 32.0/y^2 - 1.0;
   variable ampl = (0.75 + csevl(z, bm1cs, ntm1)) / sqrt(y);
   variable theta = y - 3.0*pi4 + csevl(z, bth1cs, ntth1) / y;
   result = sign(x)*abs(ampl)*cos(theta);

   return result;
}

define besj0(X) {
   %***BEGIN PROLOGUE  BESJ0
   %***DATE WRITTEN   770401   (YYMMDD)
   %***REVISION DATE  820801   (YYMMDD)
   %***CATEGORY NO.  C10A1
   %***KEYWORDS  BESSEL FUNCTION,FIRST KIND,ORDER ZERO,SPECIAL FUNCTION
   %***AUTHOR  FULLERTON, W., (LANL)
   %***PURPOSE  Computes the Bessel function of the first kind of order
   %            zero
   %***DESCRIPTION
   %
   % BESJ0(X) calculates the Bessel function of the first kind of
   % order zero for real argument X.
   %
   % Series for BJ0        on the interval  0.          to  1.60000D+01
   %                                        with weighted error   7.47E-18
   %                                         log weighted error  17.13
   %                               significant figures required  16.98
   %                                    decimal places required  17.68
   %
   % Series for BM0        on the interval  0.          to  6.25000D-02
   %                                        with weighted error   4.98E-17
   %                                         log weighted error  16.30
   %                               significant figures required  14.97
   %                                    decimal places required  16.96
   %
   % Series for BTH0       on the interval  0.          to  6.25000D-02
   %                                        with weighted error   3.67E-17
   %                                         log weighted error  16.44
   %                               significant figures required  15.53
   %                                    decimal places required  17.13
   %***REFERENCES  (NONE)
   %***ROUTINES CALLED  CSEVL,INITS,R1MACH,XERROR
   %***END PROLOGUE  BESJ0

   variable BJ0CS = [.100254161968939137, -.665223007764405132,
		     .248983703498281314, -.0332527231700357697,
		     .0023114179304694015,-.0000991127741995080,
		     .0000028916708643998,-.0000000612108586630,
		     .0000000009838650793,-.0000000000124235515,
		     .0000000000001265433 ,-.0000000000000010619,
		     .0000000000000000074 ];
   variable BM0CS = [.09284961637381644,-.00142987707403484,
		     .00002830579271257,-.00000143300611424,
		     .00000012028628046,-.00000001397113013, 
		     .00000000204076188,-.00000000035399669,
		     .00000000007024759,-.00000000001554107, 
		     .00000000000376226,-.00000000000098282, 
		     .00000000000027408,-.00000000000008091, 
		     .00000000000002511,-.00000000000000814, 
		     .00000000000000275,-.00000000000000096, 
		     .00000000000000034,-.00000000000000012,
		     .00000000000000004];
   variable BTH0CS = [-.24639163774300119,.001737098307508963 ,
		   -.000062183633402968, .000004368050165742 ,
		   -.000000456093019869, .000000062197400101 ,
		   -.000000010300442889, .000000001979526776 ,
		   -.000000000428198396, .000000000102035840 ,
		   -.000000000026363898, .000000000007297935 ,
		   -.000000000002144188, .000000000000663693 ,
		   -.000000000000215126, .000000000000072659 ,
		   -.000000000000025465, .000000000000009229 ,
		   -.000000000000003448, .000000000000001325 ,
		   -.000000000000000522, .000000000000000210 ,
		   -.000000000000000087, .000000000000000036 ];

   variable result;
   variable PI4=0.78539816339744831;
   variable NTJ0 = 13;
   variable NTM0 = 21;
   variable NTTH0 = 24;
   variable XSML = sqrt(4.0*2.22E-16);; % sqrt(4.0*R1MACH(3));
   variable XMAX = 1.0/2.22E-16; % 1.0/R1MACH(4);
      
   variable Y = abs(X);
   if (Y <= 4.0) {
      result = 1.0;
      if (Y > XSML) result = csevl(.125*Y*Y-1., BJ0CS, NTJ0);
      return result;
   }

   if (Y > XMAX) message("BESJ0: NO PRECISION BECAUSE ABS(X) IS TOO BIG");
   variable Z = 32.0/Y^2 - 1.0;
   variable AMPL = (0.75 + csevl(Z, BM0CS, NTM0)) / sqrt(Y);
   variable THETA = Y - PI4 + csevl(Z, BTH0CS, NTTH0) / Y;
   result = AMPL * cos(THETA);

   return result;
}

define erfc() { %
 
   variable x, z, t, ans, neg, i;
   
   if (_NARGS != 1) {
      message("Call as:");
      message("erfc(x) where");
      message("x is equal to or larger than zero.");
      message("---");
      message("Code returns (2/sqrt{pi})\int_x^infy exp(-t^2) dt");
      return NULL;
   }
   
   x = ();
   z=abs(x);
   t=1.0/(1.0+0.5*z);
   ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
                t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
                t*(-0.82215223+t*0.17087277)))))))));

   if (typeof(x) == Array_Type) {
      neg = where(x < 0);
      if (length(neg) != 0) {
	 for (i=0;i<length(neg);i++) ans[neg[i]] = 2.0-ans[neg[i]];
      }
   } else {
      if (x < 0) ans = 2.0-ans;
   }
	
   return ans;
}

define interpol() { % x, y, z_in

   variable x, y, z_in; 
   variable grad,d1,df,f, result;   % double 
   variable inc, jl, jm, ju; % int
   variable n, iZ, z;
   
   if (_NARGS != 3) {
      message("Call as:");
      message("result = interpol(x, y, x_int)");
      message("where x and y are the input independent and dependant variables, respectively");
      message("and x_int is the array (or scalar) to interpolate y values for.");
      if (_NARGS != 0) _pop_n(_NARGS);
      return NULL;
   } 

   (x, y, z_in) = ();
   
   n = length(x);

   if (typeof(z_in) != Array_Type) {
      z = [z_in];
   } else {
      z = @z_in;
   }
   result = 0.0*@z;
   
   inc = (x[n-1] > x[0]);

   for (iZ=0;iZ<length(z);iZ++) {
      if (( inc == 1 and ((z[iZ] < x[0]) or (z[iZ] > x[n-1]))) or
	  ( inc == 0 and ((z[iZ] > x[0]) or (z[iZ] < x[n-1])))) {
	 vmessage("interpol: Asking for %e, min is %e, max is %e",z[iZ],
		  x[0],x[n-1]);
	 message("interpol: Returning 0 for this value.");
	 result[iZ] = 0;
      } else {

	 jl = 0;
	 ju = n;
	 while (ju - jl > 1) {
	    jm = (ju + jl) / 2;
	    if ((z[iZ] > x[jm]) == inc) {
	       jl = jm;
	    } else {
	       ju = jm;
	    }
	 }

	 % ------ Z is sandwiched between JL and JU ------

	 if ((x[jl] > 0. and x[ju] > 0.) and (y[jl] > 0. and y[ju] > 0.)) {
	    grad = (log10(y[ju]) - log10(y[jl]))/(log10(x[ju]) - log10(x[jl]));
	    df = grad * (log10(z[iZ]) - log10(x[jl]));
	    result[iZ] = 10.0^(log10(y[jl]) + df);
	 } else {
	    result[iZ] = y[jl]+(y[ju]-y[jl])*((z[iZ]-x[jl])/(x[ju]-x[jl]));
	 }
      }
   }

   if (typeof(z_in) != Array_Type) result = result[0];

   return result;
}

define halo_func() { % a

   variable a;
   variable SrPerArcmin2 = 8.4616e-8;
   variable K1, K2, y1, y2;
   variable c, beta1, beta2, n_a, fdst_henke, Isca;
   c = 1.1*SrPerArcmin2; % cm^2 / arcmin^2
   
   if (_NARGS == 1) { 
      a = ();
   } else {
      message("Not a callable routine.");
      return NULL;
   }

   if (dust_cm.ZDAtype == -1) {
      n_a = dust_cm.norm*(@dust_cm.sizefunc)(a); % um^-1
   } else {
      n_a = dust_cm.norm*(@dust_cm.sizefunc)(a,dust_cm.ZDAtype); % um^-1
   }
   
   fdst_henke = interpol(dust_cm.henke_E, dust_cm.henke_F, E_cm);

   %
   %  Now, which function we use depends on HF_cm
   %
   
   if (HF_cm == 1) {  % Gaussian RG approximation, with accurate Henke factors
      beta1 = E_cm*a*(Th_cm/60.)*sqrt(0.4575);

      Isca = c * n_a * a^6 * sqrt(PI) * fdst_henke * (dust_cm.rho/3.0)^2 *
	(erfc(beta1)/(2*beta1));
   
   }
   if (HF_cm == 2) { % RG with infinite cylinder dust grains
      variable v = 1.474*a*(Th_cm/60.)*E_cm;
      variable J0v = besj0(v);
      variable J1v = besj1(v);
      y1 = 4*(1/(3*PI*v))*
	(4*v-2*PI*v^2*J0v^2 +2*PI*v*J0v*J1v + (PI-2*PI*v^2)*J1v^2);
      Isca = c * n_a * a^6 *fdst_henke * (dust_cm.rho/3.0)^2 * y1;
   }
   if (HF_cm == 3) { % Exact RG approximation, with accurate Henke factors
      K1 = 1.474*a*(Th_cm/60.)*E_cm;
      y1 = 9*(3 + 5*K1^2  + 2*K1^5*PI + (-3 + K1^2  - 2*K1^4)*cos(2*K1) -
	      6*K1*sin(2*K1) - sin(2*K1)*K1^3 - 4*K1^5*sii(2*K1)) / (30*K1^6 );
      Isca = c * n_a * a^6 *fdst_henke * (dust_cm.rho/3.0)^2 * y1;
   }
   if (HF_cm == 4) { % Mie table solution, or Exact RG.
      if (dust_cm.tableload != 1) load_dsdO_table(dust_cm);
   
      if (check_dsdO_table(a, E_cm, Th_cm, dust_cm) == 1) {
	 Isca = n_a * dsdO_val(a, E_cm, Th_cm, dust_cm);
      } else {
	 K1 = 1.474*a*(Th_cm/60.)*E_cm;
	 y1 = 9*(3 + 5*K1^2  + 2*K1^5*PI + (-3 + K1^2  - 2*K1^4)*cos(2*K1) -
		 6*K1*sin(2*K1) - sin(2*K1)*K1^3 - 4*K1^5*sii(2*K1)) / (30*K1^6 );
	 Isca = c * n_a * a^6 *fdst_henke * (dust_cm.rho/3.0)^2 * y1;
      }
   }

   if (dT_cm != 0) Isca = Isca*(dT_cm/60.); % dT_cm is in arcsec.
   
   return Isca;
}

%
%  A general catch-all routine for calculating X-ray halos.  Assumes
%  (for the moment) that Xmin = 0, Xmax = 1, and f(x) = 1.
%
define shalo() {

   variable E, a, NH, Estr, DustModel, HaloModel, SrcRad;
   variable c, iT, iD, dVec;
   variable Theta, Theta_Hi, Isca, Th, Th_hi;
   variable dqag_result, dqag_abserr, dqag_neval, dqag_ier;
   
   if ((_NARGS !=6)and(_NARGS !=7)) {
      message("Call as:");
      message("halo(NH, E, dustmodel, halomodel, SrcRad, Theta, [Theta_Hi]) where");
      message("NH, the column density in cm^-2");
      message("E, the energy of the X-ray, is in keV,");
      message("dustmodel, an integer from 1-xxx");
      message("halomodel, an integer from 1-xxx");
      message("SrcRad, the source radius in arcminutes");
      message("theta, the (observed) scattering angle is in arcseconds");
      message("and is passed in directly from sherpa.");
      _pop_n(_NARGS);
      return NULL;
   }

   if (_NARGS == 6) { 
      (NH, E, DustModel, HaloModel, SrcRad, Theta) = ();
      Theta_Hi = Theta; 
   }
   if (_NARGS == 7) (NH, E, DustModel, HaloModel, SrcRad, Theta, Theta_Hi)=();
   
   if (__is_initialized(&dust) == 0) dust = setdustparms();

   DustModel = int(DustModel);
   HaloModel = int(HaloModel);
   
   if ((HaloModel < 1) or (HaloModel > 4)) {
      vmessage("HaloModel must be 1,2,3 or 4, %d is not legal.",HaloModel);
      return NULL;
   }
   if ((DustModel < 1) or (DustModel > 52)) {
      vmessage("DustModel must be between 1-52, %d is not legal.",DustModel);
      return NULL;
   }

   %
   %  If SrcRad > 0, we need larger theta coverage, at least in to
   %  SrcRad distance.  However, we can space the theta's logarithmically,
   %  with perhaps 30 per dex.
   %
   if (SrcRad > 0.0) {
      Th = @Theta; 
      Th_hi = @Theta_Hi;
      variable loTheta = log10(SrcRad*60.);
      variable hiTheta = log10(max(Theta_Hi*1.05)); % Some breathing room
      variable NTheta  = (hiTheta - loTheta)*30;
      Theta = 10^[loTheta:hiTheta:1./NTheta];
      Theta_Hi = @Theta;
   }
      
   HF_cm = HaloModel;
   Estr = string(E)+","+string(DustModel)+","+string(HaloModel);

   if ((1==0) and (assoc_key_exists (Isca_calc_AA, Estr) != 0)) {
      Isca_calc = NH*Isca_calc_AA[Estr];
   } else {
      E_cm = E;
      switch (DustModel)
	{ case 1 : dVec = [0,1]; }      % MRN 
	{ case 2 : dVec = [2,3]; }      % WD01
	{ case 3 : dVec = [4,5]; }      % Witt, Smith, Dwek
	{ case 4 : dVec = [6,7,8]; }    % ZDA BARE-GR-S
	{ case 5 : dVec = [11,12,13]; } % ZDA BARE-GR-FB
	{ case 6 : dVec = [16,17,18]; } % ZDA BARE-GR-B
	{ case 7 : dVec = [36,37,38,39]; } % ZDA BARE-AC-S
	{ case 8 : dVec = [41,42,43,44]; } % ZDA BARE-AC-FG
	{ case 9 : dVec = [46,47,48,49]; } % ZDA BARE-AC-B
	{ case 10: dVec = [21,22,23,25]; } % ZDA COMP-GR-S
	{ case 11: dVec = [26,27,28,30]; } % ZDA COMP-GR-FG
	{ case 12: dVec = [31,32,33,35]; } % ZDA COMP-GR-B
	{ case 13: dVec = [51,52,53,54,55]; } % ZDA COMP-AC-S
	{ case 14: dVec = [56,57,58,59,60]; } % ZDA COMP-AC-FG
	{ case 15: dVec = [61,62,63,64,65]; } % ZDA COMP-AC-B
	{ case 16: dVec = [66,68,69,70]; } % ZDA COMP-NC-S
	{ case 17: dVec = [71,73,74,75]; } % ZDA COMP-NC-FG
	{ case 18: dVec = [76,78,79,80]; } % ZDA COMP-NC-B
	{ case 19: dVec = [81,82]; }       % WD01 SMC
	{ case 20: dVec = [0,1]; }  % MRN, repeated
	{ case 21: dVec = [83,84]; }  % WD01 type 1
	{ case 22: dVec = [85,86]; }  % WD01 type 2
	{ case 23: dVec = [87,88]; }  % WD01 type 3
	{ case 24: dVec = [89,90]; }  % WD01 type 4
	{ case 25: dVec = [91,92]; }  % WD01 type 5
	{ case 26: dVec = [93,94]; }  % WD01 type 6
	{ case 27: dVec = [95,96]; }  % WD01 type 7
	{ case 28: dVec = [97,98]; }  % WD01 type 8
	{ case 29: dVec = [99,100]; }  % WD01 type 9
	{ case 30: dVec = [101,102]; }  % WD01 type 10
	{ case 31: dVec = [103,104]; }  % WD01 type 11
	{ case 32: dVec = [105,106]; }  % WD01 type 12
	{ case 33: dVec = [107,108]; }  % WD01 type 13
	{ case 34: dVec = [109,110]; }  % WD01 type 14
	{ case 35: dVec = [111,112]; }  % WD01 type 15
	{ case 36: dVec = [113,114]; }  % WD01 type 16
	{ case 37: dVec = [115,116]; }  % WD01 type 17
	{ case 38: dVec = [117,118]; }  % WD01 type 18
	{ case 39: dVec = [119,120]; }  % WD01 type 19
	{ case 40: dVec = [121,122]; }  % WD01 type 20
	{ case 41: dVec = [123,124]; }  % WD01 type 21
	{ case 42: dVec = [125,126]; }  % WD01 type 22
	{ case 43: dVec = [127,128]; }  % WD01 type 23
	{ case 44: dVec = [129,130]; }  % WD01 type 24
	{ case 45: dVec = [131,132]; }  % WD01 type 25
	{ case 46: dVec = [133,134]; }  % WD01 type 26 LMC avg
	{ case 47: dVec = [135,136]; }  % WD01 type 27 LMC avg
	{ case 48: dVec = [137,138]; }  % WD01 type 28 LMC avg
	{ case 49: dVec = [139,140]; }  % WD01 type 29 LMC 2
	{ case 50: dVec = [141,142]; }  % WD01 type 30 LMC 2
	{ case 51: dVec = [143,144]; }  % WD01 type 31 LMC 2
	{ case 52: dVec = [145,146]; }  % WD01 type 32  (SMC bar)

      Isca_calc = Double_Type [length(Theta)];
      for (iT=0;iT<length(Theta);iT++) {

	 Th_cm = (Theta[iT] + Theta_Hi[iT])/2.0;
	 dT_cm = Theta_Hi[iT] - Theta[iT];
	 Isca_calc[iT] = 0.0;

	 for (iD=0;iD<length(dVec);iD++) {
	    dust_cm = dust[dVec[iD]];
	    (dqag_result, dqag_abserr, dqag_neval, dqag_ier) =
	      dqag(&halo_func, dust[dVec[iD]].Amin, dust[dVec[iD]].Amax);
	    if (dqag_ier != 0) vmessage("dqag: %f, %d, %d",
					dqag_abserr, dqag_neval, dqag_ier);
	    %vmessage("%d, I[%d], %d evals",dVec[iD], iT, dqag_neval);
	    Isca_calc[iT] += dqag_result;
	 }
      }
      Isca_calc_AA[Estr] = Isca_calc;
   }
   %
   %  Now if needed do approximate convolution to handle extended sources.
   %  See notes in "Misc Calcs and Notes" folder, dated Jan 7, 2005.
   %  Note we'll need to re-interpolate results here as well.
   %
   if (SrcRad > 0.0) {
      Isca_calc = extendSrc(SrcRad, Theta, Isca_calc);
      Isca_calc = interpol(Theta, Isca_calc, Th);
   }
   Isca = NH*Isca_calc;
   
   if (Float_Type == _typeof(Theta)) Isca = typecast(Isca, Float_Type);
   
   return Isca;
}

define inthalo_func() { % theta
   
   % haloex_tabfunc returns result in arcmin^-2.  Need to rescale
   
   variable NH = 1.e22;  % baseline value
   variable theta,result, srcrad;

   theta = (); % is in arcmin
   srcrad = 0.0;
   
   result = PI*(theta)*shalo(NH, IntEn, IntDM, IntHF, srcrad, [theta*60.]);
   return result[0];
}

define intshalo() { % E, theta0, theta1, dustmodel
   
   %
   %  Integrate shalo over some limits; ie, do
   %    \int I(theta) pi theta dtheta
   %   
  
   if (_NARGS !=5) {
      message("Call as:");
      message("intshalo(E, theta0, theta1, dustmodel, HaloModel) where");
      message("E, the energy of the X-ray, is in keV,");
      message("theta{0,1}, the (observed) scattering angles are in arcseconds");
      _pop_n(_NARGS);
      return NULL;
   }
   
   variable En, DustModel, HaloModel, theta0, theta1;
   (En, theta0, theta1, DustModel, HaloModel) = (); 
   IntDM = DustModel;
   IntEn = En;
   IntHF = HaloModel;
   
   variable result = quadint2(theta0/60., theta1/60., &inthalo_func);
   return result/1.e22;
   
}

() = register_model("shalo",["NH","E","DustModel","HaloModel","SrcRad"], 1,
		    [1.e22, 1.0  ,1 ,3, 0.0],  % default values
		    [0    , 1.e-3,1 ,1, 0.0],  % Minimum values
		    [1.e25, 10   ,52,4, 30.],  % Maximum values
		    [1    , 0    ,0 ,0, 0]); % Fix everything but NH
