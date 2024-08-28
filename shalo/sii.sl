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
