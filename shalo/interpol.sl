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
