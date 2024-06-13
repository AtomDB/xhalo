variable Th_cm, dT_cm, dust_cm, x0_cm, F_cm;
static variable Isca_calc_AA = Assoc_Type [];
static variable Isca_calc;
static variable E_cm=0;
static variable dust;

static variable IntDM, IntEn, IntX0, IntHF;

private define halox_func() { % a

   variable a;
   variable SrPerArcmin2 = 8.4616e-8;
   variable K1, K2, y1, y2, u, y;
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
   %  Now, which function we use depends on F_cm
   %
   
   if (HF_cm == 1) {  % Gaussian RG approximation, with accurate Henke factors
      beta1 = E_cm*a*(Th_cm/60.)*sqrt(0.4575);
      Isca = c * n_a * (dust_cm.rho/3.0)^2 * a^6 * fdst_henke *
	exp(-((beta1/(1-x0_cm))^2))/((1-x0_cm)^2);
   }
   if (HF_cm == 2) { % RG with infinite cylinder dust grains; NOT YET AVAILABLE
      return NULL;
   }
   if (HF_cm == 3) { % Exact RG approximation, with accurate Henke factors
      u = 1.474*a*(Th_cm/60.)*E_cm/(1-x0_cm);
      y = (9/(1-x0_cm)^2) * (sin(u) - u*cos(u))^2/u^6;
      Isca = c * n_a * a^6 *fdst_henke * (dust_cm.rho/3.0)^2 * y;
   }
   if (HF_cm == 4) { % Table model newly available
      if (dust_cm.tableloadAct != 1) load_dsdA_table(dust_cm);
      if (check_dsdA_table(a, E_cm, Th_cm/(1-x0_cm), dust_cm) == 1) {
	 Isca = n_a * dsdA_val(a, E_cm, Th_cm/(1-x0_cm), dust_cm);
      } else {
	 u = 1.474*a*(Th_cm/60.)*E_cm/(1-x0_cm);
	 y = (9/(1-x0_cm)^2) * (sin(u) - u*cos(u))^2/u^6;
	 Isca = c * n_a * a^6 *fdst_henke * (dust_cm.rho/3.0)^2 * y;
      }

   }

   if (dT_cm != 0) Isca = Isca*(dT_cm/60.); % dT_cm is in arcsec.
   
   return Isca;
}

%
%  A general catch-all routine for calculating X-ray halos at any point.
%
define shalox() {

   variable E, x0, a, NH, Estr, DustModel, HaloModel, SrcRad;
   variable c, iT, iD, dVec;
   variable Theta, Theta_Hi, Isca, Th, Th_hi;
   variable dqag_result, dqag_abserr, dqag_neval, dqag_ier;
   
   if ((_NARGS !=7)and(_NARGS !=8)) {
      message("Call as:");
      message("shalox(NH, E, x0, dustmodel, halomodel, SrcRad, Theta, [Theta_Hi]) where");
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

   if (_NARGS == 7) { 
      (NH, E, x0, DustModel, HaloModel, SrcRad, Theta) = ();
      Theta_Hi = Theta; 
   }
   if (_NARGS == 8) (NH, E, x0, DustModel, HaloModel, SrcRad, Theta, Theta_Hi)=();
   
   if (__is_initialized(&dust) == 0) dust = setdustparms();

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

   if (x0 > 1.0) {  % Return zero halo if source is between dust and observer
      Isca = 0.*Theta;
      if (Float_Type == _typeof(Theta)) Isca = typecast(Isca, Float_Type);
      return Isca;
   }
   
   x0_cm = x0;
   HF_cm = HaloModel;
   Estr = string(E)+","+string(x0)+","+string(DustModel)+","+string(HaloModel);

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
	      dqag(&halox_func, dust[dVec[iD]].Amin, dust[dVec[iD]].Amax);
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

define inthalox_func() { % theta
   
   % haloex_tabfunc returns result in arcmin^-2.  Need to rescale
   
   variable NH = 1.e22;  % baseline value
   variable theta,result, srcrad;

   theta = (); % is in arcmin
   srcrad = 0.0;
   
   result = PI*(theta)*shalox(NH, IntEn, IntX0, IntDM, IntHF, srcrad, [theta*60.]);
   return result[0];
}

define intshalox() { % E, theta0, theta1, dustmodel
   
   %
   %  Integrate shalo over some limits; ie, do
   %    \int I(theta) pi theta dtheta
   %   
  
   if (_NARGS !=6) {
      message("Call as:");
      message("intshalo(E, x0, theta0, theta1, dustmodel, HaloModel) where");
      message("E, the energy of the X-ray, is in keV,");
      message("x0, the position of the dust cloud relative to source");
      message("theta{0,1}, the (observed) scattering angles are in arcseconds");
      _pop_n(_NARGS);
      return NULL;
   }
   
   variable En, x0, DustModel, HaloModel, theta0, theta1;
   (En, x0, theta0, theta1, DustModel, HaloModel) = (); 
   IntDM = DustModel;
   IntEn = En;
   IntHF = HaloModel;
   IntX0 = x0;
   
   variable result = quadint2(theta0/60., theta1/60., &inthalox_func);
   return result/1.e22;
   
}

() = register_model("shalox",["NH","E","x0","DustModel","HaloModel","SrcRad"], 1,
		    [1.e22, 1.0  ,0.5, 1 ,3, 0.0],  % default values
		    [0    , 1.e-3,0.0, 1 ,1, 0.0],  % Minimum values
		    [1.e25, 10   ,1-1.e-9,52,4, 30.],  % Maximum values
		    [1    , 0    ,1  ,0 ,0, 0]); % Fix everything but NH and x0
