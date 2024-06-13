variable runonce;
variable dust_cm, di;
%variable ZDA = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/ZDA_BGF_0.37_0.37.fits");
variable ZDA = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/ZDA.fits");


if (__is_initialized(&runonce) == 0) {
   typedef struct {
      name, abppM, Amin, Amax, sizefunc, ap, Z, rho, norm, atomweight, 
	henke_N, henke_E, henke_F, WDtype, ZDAtype,
	file, tableload, Tvec, Avec, Evec, dsdO,
	Tmin, Tmax, Amin_tab, Amax_tab, Emin, Emax,
      	fileAct, tableloadAct, TvecAct, AvecAct, EvecAct, dsdA,
	TminAct, TmaxAct, Amin_tabAct, Amax_tabAct, EminAct, EmaxAct
   } dustparmS;
   variable runonce=1;
}

define ZDA_graindist() { % dtype, a

%  From Zubko, Dwek & Arendt 2004, ApJS, 152, 211
%  Output is dN/dA, in [um^-1 H^-1].  Input:   DTYPE = Grain type, 1-75, A = grain radius (um)

   variable dtype, a, gg;
   variable Norm, c0, b0, b1, a1, m1, b2, a2, m2, b3, a3, m3, b4, a4, m4;

   if ((_NARGS != 1) and (_NARGS !=2)) {
      message("Call as:");
      message("result = ZDA_graindist(a,dtype)");
      message("where");
      message(" dtypeindex = 0-74, the dust grain model to use");
      message(" a = dust grain radius, in um.");
      return NULL;
   }
	      
   if (_NARGS == 1) {
      a = ();
      dtype = dust_cm.ZDAtype;
   } else {
      (a, dtype) = ();
   }

   if ((dtype < 0) or (dtype >= length(ZDA.A))) {
      vmessage("Illegal Dustype %d",dtype);
      return NULL;
   }

   %
   %  Get values for this particular index
   %
   Norm = ZDA.A[dtype];  c0= ZDA.C0[dtype];   b0= ZDA.B0[dtype];
   b1= ZDA.B1[dtype];    a1= ZDA.A1[dtype];   m1= ZDA.M1[dtype];
   b2= ZDA.B2[dtype];    a2= ZDA.A2[dtype];   m2= ZDA.M2[dtype];
   b3= ZDA.B3[dtype];    a3= ZDA.A3[dtype];   m3= ZDA.M3[dtype];
   b4= ZDA.B4[dtype];    a4= ZDA.A4[dtype];   m4= ZDA.M4[dtype];
   %
   %  Calculate sizes
   %
   if (Norm == 0.0) {
      vmessage("Warning: Called with norm = 0.");
      return 0.0;
   }

   % Fixed ZDA.fits to have uniform values.
   gg = c0 + b0*log10(a) - b1*(abs(log10(a/a1)))^m1 - b2*(abs(log10(a/a2)))^m2 -
     b3*(abs(a-a3))^m3 - b4*(abs(a-a4))^m4;
   
   %if (strcmp(ZDA.NAME[dtype],"COMP") != 1) {
   %} else {
   %   gg = c0 + b0*log10(a) - b1*(abs(log(a/a1)))^m1 - b2*(abs(log(a/a2)))^m2 -
   %	b3*(abs(a-a3))^m3 - b4*(abs(a-a4))^m4;
   %}

   %vmessage("f(%d, %f) = %e (%e)", dtype, a, Norm*10^gg,gg);
   return Norm*10^gg;
}

%
%  Goal is to return the unnormalized dust size distribution, 
%  n(a)da with units of grains/H atom with sizes between a, a+da.
%  Input is dust grain size, in um; output is grains/um.
%
define dustsize0(Agr) { return Agr^(-3.5); } 
define dustsize1(Agr) { return Agr^(-3.5); }
define dustsize2(Agr) { return 1.e-4*WD_graindist(7,1,Agr/1.e4); } % 
define dustsize3(Agr) { return 1.e-4*WD_graindist(7,2,Agr/1.e4); } % 
define dustsize4(Agr) { if (Agr < 0.5) { return Agr^(-3.5);} else { return Agr^(-4);} } 
define dustsize5(Agr) { if (Agr < 0.5) { return Agr^(-3.5);} else { return Agr^(-4);} }
define dustsizeSMCs(Agr) { return 1.e-4*WD_graindist(32,1,Agr/1.e4); } % 
define dustsizeSMCg(Agr) { return 1.e-4*WD_graindist(32,2,Agr/1.e4); } % 

define dustsizeWDslc(Agr) { return 1.e-4*WD_graindist(dust_cm.WDtype,1,Agr/1.e4); } %
define dustsizeWDgra(Agr) { return 1.e-4*WD_graindist(dust_cm.WDtype,2,Agr/1.e4); } %

define dust_mass_dist(Agr) {
   %
   % Agr is in um.  
   %
   variable dnda = (@dust_cm.sizefunc)(Agr);  % grains/H atom/um
   variable dustmass = (4.*PI*dust_cm.rho/3.)*(Agr/1e4)^3; % g/grain

   return dust_cm.norm*dnda*dustmass; % g/H atom/um
}

define setdustparms() {
   
   variable dustparm = dustparmS [147];
   variable i, j, data;
   variable ppM=1.E-6;
   
   dustparm[0].name = "Silicate";
   dustparm[0].abppM = 33. * ppM;
   dustparm[0].Amin = 0.0050; % um
   dustparm[0].Amax = 0.250;  % um
   dustparm[0].sizefunc = &dustsize0;
   dustparm[0].atomweight = 172.0;
   dustparm[0].Z = dustparm[0].abppM * dustparm[0].atomweight/1.4;
   dustparm[0].rho = 3.3; % g/cm^3
   dustparm[0].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_3.30.fits";
   dustparm[0].tableload=0;
   dustparm[0].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Silicate.fits";
   dustparm[0].tableloadAct=0;
   data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Silicate_f.dat");
   dustparm[0].henke_N = data._nrows;
   dustparm[0].henke_E = Float_Type [data._nrows];
   dustparm[0].henke_F = Float_Type [data._nrows];
   
   for (i=0;i<dustparm[0].henke_N;i++) dustparm[0].henke_E[i] = data.col1[i];
   for (i=0;i<dustparm[0].henke_N;i++) dustparm[0].henke_F[i] = data.col2[i];
   dustparm[0].ZDAtype = -1;
   dustparm[0].WDtype = -1;
   %
   %  Now to calculate normalization factor.  Should be equal to 
   %  Z*1.4/(TotMass*N_A), where TotMass is the integral over the dust mass 
   %  distribution, N_A is Avagadro's number.  
   %
   dustparm[0].norm = 1.0;
   dust_cm = dustparm[0];
   variable TotMass = quadint(dust_cm.Amin,dust_cm.Amax,&dust_mass_dist);

   dustparm[0].norm = dustparm[0].Z*1.4/(TotMass*6.023e23); % unitless
%   TotMass = quadint(dust.Amin,dust.Amax,&dust_mass_dist);
%   vmessage("Silicate DustMass/GasMass = %e, Correct = 0.0044",TotMass/2.12e-24);


   %
   %  Next dust element.
   %
   dustparm[1].name = "Graphite";
   dustparm[1].abppM = 270. * ppM;
   dustparm[1].Amin = 0.0050;
   dustparm[1].Amax = 0.250;
   dustparm[1].sizefunc = &dustsize1;
   dustparm[1].ap = [3.5,  0.0,  0.0,  0.0];
   dustparm[1].atomweight = 12.0;
   dustparm[1].Z = dustparm[1].abppM * dustparm[1].atomweight/1.4;
   dustparm[1].rho = 2.2;
   dustparm[1].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_2.20.fits";
   dustparm[1].tableload=0;
   dustparm[1].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Graphite.fits";
   dustparm[1].tableloadAct=0;
   data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Graphite_f.dat");
   dustparm[1].henke_N = data._nrows;
   dustparm[1].henke_E = Float_Type [data._nrows];
   dustparm[1].henke_F = Float_Type [data._nrows];
   
   for (i=0;i<dustparm[1].henke_N;i++) dustparm[1].henke_E[i] = data.col1[i];
   for (i=0;i<dustparm[1].henke_N;i++) dustparm[1].henke_F[i] = data.col2[i];
   dustparm[1].ZDAtype = -1;
   dustparm[1].WDtype = -1;

   dustparm[1].norm = 1.0;
   dust_cm = dustparm[1];
   TotMass = quadint(dustparm[1].Amin,dustparm[1].Amax,&dust_mass_dist);

   dustparm[1].norm = dustparm[1].Z*1.4/(TotMass*6.023e23);

%   TotMass = quadint(dustparm[1].Amin,dustparm[1].Amax,&dust_mass_dist);
%   vmessage("Graphite DustMass/GasMass = %e, Correct = 0.0022",TotMass/2.12e-24);


   %  New dust model: Silicate
   %  Weingartner \& Draine, ApJ, 2001, 548, 296
   %  Model 3, with R_V = 3.1 and b_C = 6e-5 (see pg 301 for discussion)
   dustparm[2].name = "WD_Silicate";
   dustparm[2].abppM = 36.3 * ppM;
   dustparm[2].Amin = 0.0003; % From Fig 2
   dustparm[2].Amax = 0.5;
   dustparm[2].sizefunc = &dustsize2;
   dustparm[2].atomweight = 172.0;
   dustparm[2].Z = dustparm[2].abppM * dustparm[2].atomweight/1.4;
   dustparm[2].rho = 3.5;
   dustparm[2].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_3.30.fits";
   dustparm[2].tableload=0;
   dustparm[2].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Silicate.fits";
   dustparm[2].tableloadAct=0;
   data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Silicate_f.dat");
   dustparm[2].henke_N = data._nrows;
   dustparm[2].henke_E = Float_Type [data._nrows];
   dustparm[2].henke_F = Float_Type [data._nrows];
   
   for (i=0;i<dustparm[2].henke_N;i++) dustparm[2].henke_E[i] = data.col1[i];
   for (i=0;i<dustparm[2].henke_N;i++) dustparm[2].henke_F[i] = data.col2[i];
   dustparm[2].ZDAtype = -1;
   dustparm[2].WDtype = -1;  % really is, but old-style WD

   % This size function pre-normalized.
   dustparm[2].norm = 1.0;
   dust_cm = dustparm[2];
   TotMass = quadint(dustparm[2].Amin,dustparm[2].Amax,&dust_mass_dist);
 
   dustparm[2].norm = dustparm[2].Z*1.4/(TotMass*6.023e23); % unitless
%   TotMass = quadint(dustparm[2].Amin,dustparm[2].Amax,&dust_mass_dist);
%   vmessage("WD Silicate DustMass/GasMass = %e, Correct = %e",TotMass/2.12e-24, 2.98e-27*3.5/2.12e-24);
 
   
   %  New dust model: Graphite
   %  Weingartner \& Draine, ApJ, 2001, 548, 296
   %  Model 3, with R_V = 3.1 and b_C = 6e-5 (see pg 301 for discussion)
   dustparm[3].name = "WD_Graphite";
   dustparm[3].abppM = 330 * 0.7 * ppM; % 30% in gas phase
   dustparm[3].Amin = 0.0003; % um, From Fig 2
   dustparm[3].Amax = 1.250;  % um, From Fig 2
   dustparm[3].sizefunc = &dustsize3;
   dustparm[3].ap = [0.0,  0.0,  0.0,  0.0];
   dustparm[3].atomweight = 12.0;
   dustparm[3].Z = dustparm[3].abppM * dustparm[3].atomweight/1.4;
   dustparm[3].rho = 2.24; % g/cm^3
   dustparm[3].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_2.20.fits";
   dustparm[3].tableload=0;
   dustparm[3].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Graphite.fits";
   dustparm[3].tableloadAct=0;
   data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Graphite_f.dat");
   dustparm[3].henke_N = data._nrows;
   dustparm[3].henke_E = Float_Type [data._nrows];
   dustparm[3].henke_F = Float_Type [data._nrows];
   
   for (i=0;i<dustparm[3].henke_N;i++) dustparm[3].henke_E[i] = data.col1[i];
   for (i=0;i<dustparm[3].henke_N;i++) dustparm[3].henke_F[i] = data.col2[i];
   dustparm[3].ZDAtype = -1;
   dustparm[3].WDtype = -1; % really is, but old-style WD

   % This size function pre-normalized.
   dustparm[3].norm = 1.0;
   dust_cm = dustparm[3];
   TotMass = quadint(dustparm[3].Amin,dustparm[3].Amax,&dust_mass_dist);
 
   dustparm[3].norm = dustparm[3].Z*1.4/(TotMass*6.023e23); % cm^-1
 
 % TotMass = quadint(dustparm[3].Amin,dustparm[3].Amax,&dust_mass_dist);
%  vmessage("WD Graphite DustMass/GasMass = %e, Correct = %e",TotMass/2.12e-24, 2.07e-27*2.24/2.12e-24);

   %  New dust model: SMC Silicate
   %  Weingartner \& Draine, ApJ, 2001, 548, 296
   %  Model 32 for the SMC 
   dustparm[81].name = "WD_SMCSilicate";
   dustparm[81].abppM = 36.3 * ppM /4.0; % abundance reduced by 4x for SMC
   dustparm[81].Amin = 0.0003; % From Fig 2
   dustparm[81].Amax = 0.5;
   dustparm[81].sizefunc = &dustsizeSMCs;
   dustparm[81].atomweight = 172.0;
   dustparm[81].Z = dustparm[81].abppM * dustparm[81].atomweight/1.4;
   dustparm[81].rho = 3.5;
   dustparm[81].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_3.30.fits";
   dustparm[81].tableload=0;
   dustparm[81].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Silicate.fits";
   dustparm[81].tableloadAct=0;
   data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Silicate_f.dat");
   dustparm[81].henke_N = data._nrows;
   dustparm[81].henke_E = Float_Type [data._nrows];
   dustparm[81].henke_F = Float_Type [data._nrows];
   
   for (i=0;i<dustparm[81].henke_N;i++) dustparm[81].henke_E[i] = data.col1[i];
   for (i=0;i<dustparm[81].henke_N;i++) dustparm[81].henke_F[i] = data.col2[i];
   dustparm[81].ZDAtype = -1;
   dustparm[81].WDtype = -1; % really is, but old-style WD

   % This size function pre-normalized.
   dustparm[81].norm = 1.0;
   dust_cm = dustparm[81];
   TotMass = quadint(dustparm[81].Amin,dustparm[81].Amax,&dust_mass_dist);
 
   dustparm[81].norm = dustparm[81].Z*1.4/(TotMass*6.023e23); % unitless
%   TotMass = quadint(dustparm[81].Amin,dustparm[81].Amax,&dust_mass_dist);
%   vmessage("WD Silicate DustMass/GasMass = %e, Correct = %e",TotMass/2.12e-24, 2.98e-27*3.5/2.12e-24);
 
   
   %  New dust model: SMC Graphite
   %  Weingartner \& Draine, ApJ, 2001, 548, 296
   %  Model 32 for the SMC 
   dustparm[82].name = "WD_SMCGraphite";
   dustparm[82].abppM = 330 * 0.7 * ppM/4.0; % 30% in gas phase, abundance reduced by 4x for SMC
   dustparm[82].Amin = 0.0003; % um, From Fig 2
   dustparm[82].Amax = 1.250;  % um, From Fig 2
   dustparm[82].sizefunc = &dustsizeSMCg;
   dustparm[82].ap = [0.0,  0.0,  0.0,  0.0];
   dustparm[82].atomweight = 12.0;
   dustparm[82].Z = dustparm[82].abppM * dustparm[82].atomweight/1.4;
   dustparm[82].rho = 2.24; % g/cm^3
   dustparm[82].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_2.20.fits";
   dustparm[82].tableload=0;
   dustparm[82].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Graphite.fits";
   dustparm[82].tableloadAct=0;
   data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Graphite_f.dat");
   dustparm[82].henke_N = data._nrows;
   dustparm[82].henke_E = Float_Type [data._nrows];
   dustparm[82].henke_F = Float_Type [data._nrows];
   
   for (i=0;i<dustparm[82].henke_N;i++) dustparm[82].henke_E[i] = data.col1[i];
   for (i=0;i<dustparm[82].henke_N;i++) dustparm[82].henke_F[i] = data.col2[i];
   dustparm[82].ZDAtype = -1;
   dustparm[82].WDtype = -1; % really is, but old-style WD

   % This size function pre-normalized.
   dustparm[82].norm = 1.0;
   dust_cm = dustparm[82];
   TotMass = quadint(dustparm[82].Amin,dustparm[82].Amax,&dust_mass_dist);
 
   dustparm[82].norm = dustparm[82].Z*1.4/(TotMass*6.023e23); % cm^-1
 
 %  TotMass = quadint(dustparm[3].Amin,dustparm[3].Amax,&dust_mass_dist);
 %  vmessage("WD Graphite DustMass/GasMass = %e, Correct = %e",TotMass/2.12e-24, 2.07e-27*2.24/2.12e-24);
   % Load all the WD models 
   for (i=0;i<32;i++) {
      variable indx = 83+2*i;
      dustparm[indx].name = "WD_Silicate";
      dustparm[indx].abppM = 36.3 * ppM;
      % For LMC models, reduce abundance by 1.6 (WD01 page 304)
      if (indx > 132 and indx < 145 ) dustparm[indx].abppM /= 1.6;
      % For SMC models, reduce abundance by 4 (WD01 page 304)
      if (indx > 144 and indx < 147 ) dustparm[indx].abppM /= 4.0;
      dustparm[indx].Amin = 0.0003; % From Fig 2
      dustparm[indx].Amax = 0.5;
      dustparm[indx].sizefunc = &dustsizeWDslc;
      dustparm[indx].atomweight = 172.0;
      dustparm[indx].Z = dustparm[indx].abppM * dustparm[indx].atomweight/1.4;
      dustparm[indx].rho = 3.5;
      dustparm[indx].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_3.30.fits";
      dustparm[indx].tableload=0;
      dustparm[indx].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Silicate.fits";
      dustparm[indx].tableloadAct=0;
      data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Silicate_f.dat");
      dustparm[indx].henke_N = data._nrows;
      dustparm[indx].henke_E = Float_Type [data._nrows];
      dustparm[indx].henke_F = Float_Type [data._nrows];

      for (j=0;j<dustparm[indx].henke_N;j++) dustparm[indx].henke_E[j] = data.col1[j];
      for (j=0;j<dustparm[indx].henke_N;j++) dustparm[indx].henke_F[j] = data.col2[j];

      dustparm[indx].ZDAtype = -1;
      dustparm[indx].WDtype = i+1;

      % This size function pre-normalized.
      dustparm[indx].norm = 1.0;
      dust_cm = dustparm[indx];
      TotMass = quadint(dustparm[indx].Amin,dustparm[indx].Amax,&dust_mass_dist);
 
      dustparm[indx].norm = dustparm[indx].Z*1.4/(TotMass*6.023e23); % unitless
      %   TotMass = quadint(dustparm[81].Amin,dustparm[81].Amax,&dust_mass_dist);
      %   vmessage("WD Silicate DustMass/GasMass = %e, Correct = %e",TotMass/2.12e-24, 2.98e-27*3.5/2.12e-24);

      dustparm[indx+1].name = "WD_Graphite";
      dustparm[indx+1].abppM = 330 * 0.7 * ppM;
      % For LMC models, reduce abundance by 1.6 (WD01 page 304)
      if (indx > 132 and indx < 145 ) dustparm[indx+1].abppM /= 1.6;
      % For SMC models, reduce abundance by 4 (WD01 page 304)
      if (indx > 144 and indx < 147 ) dustparm[indx+1].abppM /= 4.0;
      
      dustparm[indx+1].Amin = 0.0003; % um, From Fig 2
      dustparm[indx+1].Amax = 1.250;  % um, From Fig 2
      if (indx == 134) dustparm[indx+1].Amax = 2.5;  % um, see WD01 Fig 18
      dustparm[indx+1].sizefunc = &dustsizeWDgra;
      dustparm[indx+1].ap = [0.0,  0.0,  0.0,  0.0];
      dustparm[indx+1].atomweight = 12.0;
      dustparm[indx+1].Z = dustparm[indx+1].abppM * dustparm[indx+1].atomweight/1.4;
      dustparm[indx+1].rho = 2.24; % g/cm^3
      dustparm[indx+1].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_2.20.fits";
      dustparm[indx+1].tableload=0;
      dustparm[indx+1].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Graphite.fits";
      dustparm[indx+1].tableloadAct=0;
      data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Graphite_f.dat");
      dustparm[indx+1].henke_N = data._nrows;
      dustparm[indx+1].henke_E = Float_Type [data._nrows];
      dustparm[indx+1].henke_F = Float_Type [data._nrows];
      
      for (j=0;j<dustparm[indx+1].henke_N;j++) dustparm[indx+1].henke_E[j] = data.col1[j];
      for (j=0;j<dustparm[indx+1].henke_N;j++) dustparm[indx+1].henke_F[j] = data.col2[j];
      dustparm[indx+1].ZDAtype = -1;
      dustparm[indx+1].WDtype = i+1;
      
      % This size function pre-normalized.
      dustparm[indx+1].norm = 1.0;
      dust_cm = dustparm[indx+1];
      TotMass = quadint(dustparm[indx+1].Amin,dustparm[indx+1].Amax,&dust_mass_dist);
 
      dustparm[indx+1].norm = dustparm[indx+1].Z*1.4/(TotMass*6.023e23); % cm^-1
   }

   %
   %  New dust model: Silicate with large grain sizes
   %    Specifically from Witt, Smith, Dwek (2001)
   %
   dustparm[4].name = "Lrg_Silicate";
   dustparm[4].abppM = 33. * ppM;
   dustparm[4].Amin = 0.0050; 
   dustparm[4].Amax = 2.0;
   dustparm[4].sizefunc = &dustsize4;
   dustparm[4].atomweight = 172.0;
   dustparm[4].Z = dustparm[4].abppM * dustparm[4].atomweight/1.4;
   dustparm[4].rho = 3.3;
   dustparm[4].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_3.30.fits";
   dustparm[4].tableload=0;
   dustparm[4].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Silicate.fits";
   dustparm[4].tableloadAct=0;
   data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Silicate_f.dat");
   dustparm[4].henke_N = data._nrows;
   dustparm[4].henke_E = Float_Type [data._nrows];
   dustparm[4].henke_F = Float_Type [data._nrows];
   
   for (i=0;i<dustparm[4].henke_N;i++) dustparm[4].henke_E[i] = data.col1[i];
   for (i=0;i<dustparm[4].henke_N;i++) dustparm[4].henke_F[i] = data.col2[i];
   dustparm[4].ZDAtype = -1;
   dustparm[4].WDtype = -1; 
   
   % Normalize size function
   dustparm[4].norm = 1.0;
   dust_cm = dustparm[4];
   TotMass = quadint(dustparm[4].Amin,0.5,&dust_mass_dist) +
     quadint(0.5,dustparm[4].Amax,&dust_mass_dist);
   dustparm[4].norm = dustparm[4].Z*1.4/(TotMass*6.023e23); % cm^-1

%   TotMass = quadint(dustparm[4].Amin,0.5,&dust_mass_dist) +
%     quadint(0.5,dustparm[4].Amax,&dust_mass_dist);
%   vmessage("Lrg Silicate DustMass/GasMass = %e, Correct = 0.0044",
%	    (TotMass/2.12e-24));

   
   %  New dust model: Large Graphite
   %  Specifically from Witt, Smith, Dwek (2001)
   %  
   dustparm[5].name = "Lrg_Graphite";
   dustparm[5].abppM = 270 * ppM; 
   dustparm[5].Amin = 0.0050; % um
   dustparm[5].Amax = 2.0;  % um
   dustparm[5].sizefunc = &dustsize5;
   dustparm[5].atomweight = 12.0;
   dustparm[5].Z = dustparm[5].abppM * dustparm[5].atomweight/1.4;
   dustparm[5].rho = 2.2; % g/cm^3
   dustparm[5].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_2.20.fits";
   dustparm[5].tableload=0;
   dustparm[5].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Graphite.fits";
   dustparm[5].tableloadAct=0;
   data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Graphite_f.dat");
   dustparm[5].henke_N = data._nrows;
   dustparm[5].henke_E = Float_Type [data._nrows];
   dustparm[5].henke_F = Float_Type [data._nrows];
   
   for (i=0;i<dustparm[5].henke_N;i++) dustparm[5].henke_E[i] = data.col1[i];
   for (i=0;i<dustparm[5].henke_N;i++) dustparm[5].henke_F[i] = data.col2[i];
   dustparm[5].ZDAtype = -1;
   dustparm[5].WDtype = -1; 

   % Normalize size function 
   dustparm[5].norm = 1.0;
   dust_cm = dustparm[5];
   TotMass = quadint(dustparm[5].Amin,0.5,&dust_mass_dist) +
     quadint(0.5,dustparm[5].Amax,&dust_mass_dist);

   dustparm[5].norm = dustparm[5].Z*1.4/(TotMass*6.023e23); % cm^-1

%   TotMass = quadint(dustparm[5].Amin,0.5,&dust_mass_dist) +
%     quadint(0.5,dustparm[5].Amax,&dust_mass_dist);
%   vmessage("Lrg Graphite DustMass/GasMass = %e, Correct = 0.0022",
%	    (TotMass/2.12e-24));

   %
   %  Load all the Zubko, Dwek & Arendt dust models.
   %
   for (i=0;i<75;i++) {
      dustparm[i+6].name = ZDA.NAME[i]+","+ZDA.TYPE[i];
      dustparm[i+6].abppM = ppM*(ZDA.ABPPMC[i] + ZDA.ABPPMO[i] + ZDA.ABPPMSI[i] +
				 ZDA.ABPPMMG[i] +ZDA.ABPPMFE[i]+ ZDA.ABPPMN[i]);
      dustparm[i+6].Amin = ZDA.AMIN[i];
      dustparm[i+6].Amax = ZDA.AMAX[i];
      dustparm[i+6].ZDAtype = i;
      dustparm[i+6].WDtype = -1;
      dustparm[i+6].sizefunc = &ZDA_graindist;
      if (ZDA.A[i] != 0) {
	 dustparm[i+6].atomweight = ppM*(ZDA.ABPPMC[i]*12.0 + 
					 ZDA.ABPPMO[i]*16.0 + 
					 ZDA.ABPPMSI[i]*28.09 +
					 ZDA.ABPPMMG[i]*24.31 +
					 ZDA.ABPPMFE[i]*55.85 + 
					 ZDA.ABPPMN[i]*14.01)/dustparm[i+6].abppM;
	 dustparm[i+6].Z = dustparm[i+6].abppM * dustparm[i+6].atomweight/1.4;
	 
	 if (ZDA.ITYPE[i] == 1) dustparm[i+6].rho = 2.24; % PAHs; g/cm^3
	 if (ZDA.ITYPE[i] == 2) {
	    if (ZDA.INAME[i] <=6) dustparm[i+6].rho = 2.24; % Graphite; g/cm^3
	    if (ZDA.INAME[i] > 6) dustparm[i+6].rho = 1.83; % Amorphous Carbon; g/cm^3
	 }
	 if (ZDA.ITYPE[i] == 3) dustparm[i+6].rho = 3.5; % Silicates; g/cm^3
	 if (ZDA.ITYPE[i] == 4) dustparm[i+6].rho = 3.5; % Silicates; g/cm^3
	 if (ZDA.ITYPE[i] == 5) {
	    if (ZDA.INAME[i] == 4) dustparm[i+6].rho = 1.05; % COMP-GR-S; g/cm^3
	    if (ZDA.INAME[i] == 5) dustparm[i+6].rho = 1.05; % COMP-GR-FG; g/cm^3
	    if (ZDA.INAME[i] == 6) dustparm[i+6].rho = 1.89; % COMP-GR-B; g/cm^3
	    if (ZDA.INAME[i] == 10)dustparm[i+6].rho = 0.84; % COMP-AC-S; g/cm^3
	    if (ZDA.INAME[i] == 11)dustparm[i+6].rho = 1.05; % COMP-AC-FG; g/cm^3
	    if (ZDA.INAME[i] == 12)dustparm[i+6].rho = 1.05; % COMP-AC-B; g/cm^3
	    if (ZDA.INAME[i] == 13)dustparm[i+6].rho = 1.05; % COMP-NC-S; g/cm^3
	    if (ZDA.INAME[i] == 14)dustparm[i+6].rho = 1.05; % COMP-NC-FG; g/cm^3
	    if (ZDA.INAME[i] == 15)dustparm[i+6].rho = 1.05; % COMP-NC-B; g/cm^3
	 }
	 %
	 %  This is a placeholder until I get the real optical constants file.
	 %  At most energies, this value is within 5% of 1, so not too important.
	 %
	 % Default value
	 data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Graphite_f.dat");
	 % Graphite
	 if ((dustparm[i+6].rho > 2.0) and (dustparm[i+6].rho < 3.0)) {
	    dustparm[i+6].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_2.20.fits";
	    dustparm[i+6].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Graphite.fits";
	    data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Graphite_f.dat");
	    
	 }
	 % Silicate
	 if ((dustparm[i+6].rho > 3.0) and (dustparm[i+6].rho < 4.0)) {
	    dustparm[i+6].file = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/DsdO_3.30.fits";
	    dustparm[i+6].fileAct = "/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Scatter_Silicate.fits";
	    data = readfile("/local/data/orcrist2/valencic/home_halo_fitting/ciao34/dat/Silicate_f.dat");
	 }
	 
	 dustparm[i+6].tableload=0;
	 dustparm[i+6].tableloadAct=0;
	 
	 dustparm[i+6].henke_N = data._nrows;
	 dustparm[i+6].henke_E = Float_Type [data._nrows];
	 dustparm[i+6].henke_F = Float_Type [data._nrows];
	 for (j=0;j<dustparm[i+6].henke_N;j++) dustparm[i+6].henke_E[j] = data.col1[j];
	 for (j=0;j<dustparm[i+6].henke_N;j++) dustparm[i+6].henke_F[j] = data.col2[j];
	 
	 % This size function pre-normalized.
	 dustparm[i+6].norm = 1.0;
	 
	 %
	 %  Eliminate this, except for checking.
	 %  Note: Would need to do something special for Silicate #2 cases.
	 %
	%dust_cm = dustparm[i+6];
	%TotMass = quadint(dustparm[i+6].Amin,dustparm[i+6].Amax,&dust_mass_dist);
	%(TotMass,,,) = dqag(&dust_mass_dist,
	%		     dustparm[i+6].Amin,dustparm[i+6].Amax);

	% if (TotMass > 0)  {
	%    dustparm[i+6].norm = dustparm[i+6].Z*1.4/(TotMass*6.023e23); % cm^-1
	%    vmessage("Norm[%d] = %f",i, dustparm[i+6].norm);
	% } else {
	%    dustparm[i+6].norm = 1.0;
	%    vmessage("WARNING: %s(%s) was everywhere zero.",ZDA.NAME[i],ZDA.TYPE[i]);
	% }
      } else {
	 dustparm[i+6].norm = 0.0;
      }
   }
   return dustparm;
}
