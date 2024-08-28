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
