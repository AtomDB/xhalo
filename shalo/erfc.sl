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
