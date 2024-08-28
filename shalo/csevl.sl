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

