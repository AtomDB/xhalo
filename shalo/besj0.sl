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

