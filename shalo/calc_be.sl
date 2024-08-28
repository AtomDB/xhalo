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
   } else {
      z_in = ();
      y = ();
      x = ();
   }
   
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

define num_models_mdl(modelstr, mdl) {
   
   variable iM, result;

   result = 0;
   for (iM=0;iM<length(mdl.model);iM++) {
      if (strncmp(modelstr, mdl.model[iM], strlen(modelstr)) == 0) result++;
   }
   return result;
}

define get_model_parmname(modelstr, Nmodel, mdl) {

   variable iM, result, nummatch, model, parmname;
   result = 0;
   variable parsestr = "%[a-zA-Z0-9][%[a-zA-Z0-9]]";

   for (iM=0;iM<length(mdl.model);iM++) {
      if (strncmp(modelstr, mdl.model[iM], strlen(modelstr)) == 0) result++;
      if (result == Nmodel) {
	 nummatch = sscanf(mdl.model[iM], parsestr, &model, &parmname);
      }
   }
   return parmname;
}

define get_parm_value(parmname, item, mdl) {

   variable result, iM;
   
   variable parsestr = Sprintf("%s.%s",parmname,item,2);
   for (iM=0;iM<length(mdl.model);iM++) {
      if (strncmp (parsestr, mdl.parname[iM], strlen(parsestr)) == 0) {
	 result = mdl.parvalue[iM];
      }
   }
   return result;
   
}

define calc_be() { % theta, R, Elo, Ehi (optional)

   variable term1, term2, term3, B, Ee;
   
   variable g_vec = [1,1.5,2,2.5,3,4,5];
   variable a_vec = [0.283,0.147,0.103,0.0852,0.0742,0.0725,0.0922];
   variable y1_vec= [0.8,1.3,1.8,2.2,2.7,3.4,4.0];
   variable y2_vec= [0.00045,0.011,0.032,0.1,0.18,0.38,0.65];

   variable t = 1.;     % Assume theta = 1'
   variable R = 1.0;    % Assume R = 1kpc
   variable Elo = 0.1;  % keV
   variable Ehi = 10.0; % keV

   if (_NARGS == 0) {   
      message("Assuming object size = 1', distance = 1 kpc, E = 0.1-10 keV");
      message("Include input parameters if you want different values, eg");
      message("calc_be(1.7, 4.0, 0.1, 2.4) for size = 1.7', D=4 kpc, E=0.1-2.4");
   } else {
      if (_NARGS == 1) t = ();
      if (_NARGS == 2) (t, R) = ();
      if (_NARGS == 3) (t, R, Elo) = ();
      if (_NARGS == 4) (t, R, Elo, Ehi) = ();
   }

   %
   %  Scan the mdl, seaching for power laws (matching to powlaw1d[XXX]
   %  where XXX is the model name
   %
   variable mdl = sherpa_get_mdl();

   variable Npowlaw = num_models_mdl("powlaw1d",mdl);
   
   if (Npowlaw/4 > 1) {
      vmessage("Found %d powlaw1d models; using first.",Npowlaw/4);
   }
   
   variable parmname = get_model_parmname("powlaw1d", 1, mdl);
   vmessage("Using parameter %s",parmname);
   
   variable ampl = get_parm_value(parmname, "ampl", mdl);
   variable alpha = get_parm_value(parmname, "gamma", mdl);
   variable gamma = 2*alpha-1.0;
   
   variable a = interpol(g_vec, a_vec, gamma);
   variable y1= interpol(g_vec, y1_vec, gamma);
   variable y2= interpol(g_vec, y2_vec, gamma);
   
   term1 = 1.15e-4*(t^(-6/7.))*(R^(-2/7.))*(ampl/a)^(2./7)*(25.87)^((1-gamma)/7.);
   term2 = (0.197)^(2*(2-gamma)/7.)*(gamma-2)^(-2/7.);
   term3 = (- (Ehi/y2)^(-gamma/2. + 1) + (Elo/y1)^(-gamma/2. + 1))^(2/7.);
   B = term1*term2*term3;

   Ee = B*B*t*t*t*R*R*R*2.878e52;
   
   vmessage("B (uG) = %f",B*1.e6);
   vmessage("E (erg) = %e",Ee);
}
