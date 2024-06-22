#include<math.h>
/**  
 *  @brief Contains function to find roots of cubic equation 
 *  
 */

int cubic_solve(double a,double b,double c,
        double d, double *roots){

	 ///  @brief Solve the roots of equation ax^3 +  bx^2 + cx + d = 0
  	 ///
	 ///  @param a  coefficient of x^3
	 ///  @param b  coefficient of x^2
	 ///  @param c  coefficient of x^1
	 ///  @param d  coefficient of x^0
	 ///  @param roots  Roots of the equation;   
     
	 ///  @return   number of real roots. 
	 ///
	 ///  @details
	 /// 
 	 ///  @note     The array roots is 1 dimensional size 6. The roots are orgainzed
     /// as roots[0], roots[1] will be real and imaginary component of first root.
     ///         

    int nroot;
    double  pi = 3.141592654;
    double DD, p, q, phi, temp1, temp2, y1,y2,y3, u, v, y2r, y2i;
    DD=0e0;  p=0e0;  q=0e0;  phi=0e0;  temp1=0e0;  temp2=0e0;  
    y1=0e0; y2=0e0; y3=0e0;  u=0e0; v=0e0;  y2r=0e0;  y2i=0e0;

    p  = c/a - b*b/a/a/3.;
    q  = (2.*b*b*b/a/a/a - 9.*b*c/a/a + 27.*d/a) / 27.;
    DD = p*p*p/27. + q*q/4.;
    nroot = 0;
    if(DD < 0.){
        nroot = 3;
        phi = acos(-q/2./sqrt(fabs(p*p*p)/27.));
        temp1 = 2.*sqrt(fabs(p)/3.);
        y1 =  temp1*cos(phi/3.);
        y2 = -temp1*cos((phi+pi)/3.);
        y3 = -temp1*cos((phi-pi)/3.);
    } else{
        nroot = 1;
        temp1 = -q/2. + sqrt(DD);
        temp2 = -q/2. - sqrt(DD);
        u = pow(fabs(temp1),(1./3.));
        v = pow(fabs(temp2),(1./3.));
        if(temp1 < 0.) u=-u;
        if(temp2 < 0.) v=-v;
        y1  = u + v;
        y2r = -(u+v)/2.;
        y2i =  (u-v)*sqrt(3.)/2.;
    }
    temp1 = b/a/3.;
    y1 = y1-temp1;
    y2 = y2-temp1;
    y3 = y3-temp1;
    y2r=y2r-temp1;

    if(DD < 0.){
        roots[0] = y1; roots[1] = 0e0;

        roots[2] = y2; roots[3] = 0e0;

        roots[4] = y3; roots[5] = 0e0;

    }else if(DD == 0.0e0) {
        roots[0] = y1; roots[1] = 0e0;
        roots[2] = y2r; roots[3] = 0e0;
        roots[4] = y2r; roots[5] = 0e0;
    }else{

        roots[0] = y1; roots[1] = 0e0;
        roots[2] = y2r; roots[3] = y2i;
        roots[4] = y2r; roots[5] = -y2i;
    }

return nroot;
}
