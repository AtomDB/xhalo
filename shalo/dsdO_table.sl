private define GetHdrKey() { % header ,keyword, integer
   
   variable i, Nel, result=NULL;
   variable Key, Value, matched=0, len;
   variable header, keyword, intval;
   
   if (_NARGS != 3) {
      message("Bad call");
      return NULL;
   }
   (header, keyword, intval) = ();
        
   for (i=0;i<length(header);i++) {
      Nel = sscanf(header[i],"%[a-zA-Z0-9_] = %[^/]",&Key, &Value);
      if (Nel == 2) {
         Key = strcompress(Key," ");
         Value = strcompress(Value," ");
         len = max([strlen(keyword),strlen(Key)]);
	 if (strncmp(Key,keyword,len) == 0) {
	    if (matched == 1) message("Odd.");
	    matched = 1;
	    result = Value;
         }
      }
   }
   
   %
   %  Check for a result
   %
   if (result == NULL) return NULL;
   %
   %  Return appropriate type
   %
   result = (atof(result));
   if (intval == 1) {  % integer
      %result = atof(result);
      result = int(result);
   }

   return result;
}


define load_dsdO_table(dust)  {
   %
   %  Data are in dust.file.  Order in file is 
   %  data.Energy, data.Size, data.ThetaObs (fastest)
   %
   variable data = readfile(dust.file);
   %
   %  Need to extract Tvec, Avec, Evec, and Intensity(E, A, T)
   %  Can derive vectors from header if necessary.
   variable NTheta = GetHdrKey(data._header, "INUMT", 1);
   variable NSize  = GetHdrKey(data._header, "INUMA", 1);
   variable NEnergy= GetHdrKey(data._header, "INUME", 1);
   
   variable Tvec = data.ThetaObs[[0:NTheta:1]];
   variable Avec = data.Size[[0:NSize*NTheta:NTheta+1]];
   variable Evec = data.Energy[[0:NEnergy*NSize*NTheta:NSize*(NTheta+1)]];
   
   dust.Tvec = @Tvec;
   dust.Tmin = min(Tvec);
   dust.Tmax = max(Tvec);
   dust.Avec = @Avec;
   dust.Amin_tab = min(Avec);
   dust.Amax_tab = max(Avec);
   dust.Evec = @Evec;
   dust.Emin = min(Evec);
   dust.Emax = max(Evec);
   dust.dsdO = _reshape(data.Intensity,[NEnergy,NSize,NTheta+1]);

   vmessage("Loaded %s",dust.file);
   dust.tableload=1;
   
}

define check_dsdO_table(a, En, Theta, dust) {

   if ((Theta >= dust.Tmin) and (Theta < dust.Tmax) and
       (a >= dust.Amin_tab) and (a < dust.Amax_tab) and
       (En >= dust.Emin) and (En < dust.Emax)) { 
      return 1;
   } else {
      return 0;
   }
   return NULL; 
}

define table_norm(a,dust) {
   
%
% a is in um, rho is in g/cc
%
%    norm = 1.20d-15/((rho/3.0)*a^3)
 
   variable N_A = 6.023e23;
   return (dust.atomweight*dust.abppM)/((4*PI*a^3/3)*dust.rho*1.e-12*N_A);

}

define tabinv(x, v) {
   %
   %  Given a value x and a vector v in increasing order, find the 
   %  index i such that v(i) <= x < v(i+1).  
   
   variable gp = max(where(x >= v));

   if (length(gp) == 0) {
      printf("tabinv bad at %e",x);
      return NULL;
   }
   
   if (gp+1 == length(v)) return gp-1;

   return gp;
     
   
}

define dsdO_val(a, En, Theta, dust) {

   %
   %  Return interpolated best value from table.
   %  (ugh, writing this code again).  
   %
   variable i,j,k, result;
   variable wt_i = Float_Type[2];
   variable wt_j = Float_Type[2];
   variable wt_k = Float_Type[2];
   
   variable i_Tlo = tabinv(Theta, dust.Tvec);
   variable j_Alo = tabinv(a, dust.Avec);
   variable k_Elo = tabinv(En, dust.Evec);

   if ((i_Tlo == NULL) or (j_Alo == NULL) or (k_Elo == NULL)) return 0;
   %
   %  There are 8 points of interest.  Each has a different weight.
   %
   wt_i[1] = (Theta - dust.Tvec[i_Tlo])/(dust.Tvec[i_Tlo+1]-dust.Tvec[i_Tlo]);
   wt_i[0] = 1 - wt_i[1];

   wt_j[1] = (a - dust.Avec[j_Alo])/(dust.Avec[j_Alo+1]-dust.Avec[j_Alo]);
   wt_j[0] = 1 - wt_j[1];
   
   wt_k[1] = (En - dust.Evec[k_Elo])/(dust.Evec[k_Elo+1]-dust.Evec[k_Elo]);
   wt_k[0] = 1 - wt_k[1];

   result = 0.0;
   for (i=0;i<2;i++) {
      for (j=0;j<2;j++) {
	 for (k=0;k<2;k++) {
	    result = result + wt_i[i]*wt_j[j]*wt_k[k]*
	      dust.dsdO[k_Elo+k,j_Alo+j,i_Tlo+i];
	 }
      }
   }
   return result/table_norm(a,dust);
}


define load_dsdA_table(dust)  {
   %
   %  Data are in dust.file.  Order in file is 
   %  data.Energy, data.Size, data.ThetaObs (fastest)
   %
   variable data = readfile(dust.fileAct);
   %
   %  Need to extract Tvec, Avec, Evec, and Intensity(E, A, T)
   %  Can derive vectors from header if necessary.
   variable NTheta = GetHdrKey(data._header, "NUMT", 1);
   variable NSize  = GetHdrKey(data._header, "NUMA", 1);
   variable NEnergy= GetHdrKey(data._header, "NUME", 1);

   variable Tvec = data.THETAACT[0,*];
   variable Avec = data.CORESIZE[[0:NSize-1:1]];
   variable Evec = data.ENERGY[[0:(NEnergy-1)*NSize:NSize]];

   dust.TvecAct = @Tvec;
   dust.TminAct = min(Tvec);
   dust.TmaxAct = max(Tvec);
   dust.AvecAct = @Avec;
   dust.Amin_tabAct = min(Avec);
   dust.Amax_tabAct = max(Avec);
   dust.EvecAct = @Evec;
   dust.EminAct = min(Evec);
   dust.EmaxAct = max(Evec);
   dust.dsdA = _reshape(data.INTENSIT,[NEnergy,NSize,NTheta+1]);

   vmessage("Loaded %s",dust.fileAct);
   dust.tableloadAct=1;
   
}

define check_dsdA_table(a, En, Theta, dust) {

   if ((Theta >= dust.TminAct) and (Theta < dust.TmaxAct) and
       (a >= dust.Amin_tabAct) and (a < dust.Amax_tabAct) and
       (En >= dust.EminAct) and (En < dust.EmaxAct)) { 
      return 1;
   } else {
      return 0;
   }
   return NULL; 
}

define dsdA_val(a, En, Theta, dust) {

   %
   %  Return interpolated best value from table.
   %  (ugh, writing this code again).  
   %
   variable i,j,k, result;
   variable wt_i = Float_Type[2];
   variable wt_j = Float_Type[2];
   variable wt_k = Float_Type[2];
   
   variable i_Tlo = tabinv(Theta, dust.TvecAct);
   variable j_Alo = tabinv(a, dust.AvecAct);
   variable k_Elo = tabinv(En, dust.EvecAct);

   if ((i_Tlo == NULL) or (j_Alo == NULL) or (k_Elo == NULL)) return 0;
   %
   %  There are 8 points of interest.  Each has a different weight.
   %
   wt_i[1] = (Theta - dust.TvecAct[i_Tlo])/
     (dust.TvecAct[i_Tlo+1]-dust.TvecAct[i_Tlo]);
   wt_i[0] = 1 - wt_i[1];

   wt_j[1] = (a - dust.AvecAct[j_Alo])/(dust.AvecAct[j_Alo+1]-dust.AvecAct[j_Alo]);
   wt_j[0] = 1 - wt_j[1];
   
   wt_k[1] = (En - dust.EvecAct[k_Elo])/(dust.EvecAct[k_Elo+1]-dust.EvecAct[k_Elo]);
   wt_k[0] = 1 - wt_k[1];

   result = 0.0;
   for (i=0;i<2;i++) {
      for (j=0;j<2;j++) {
	 for (k=0;k<2;k++) {
	    result = result + wt_i[i]*wt_j[j]*wt_k[k]*
	      dust.dsdA[k_Elo+k,j_Alo+j,i_Tlo+i];
	 }
      }
   }
   return result/table_norm(a,dust);
}
