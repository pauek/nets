// 
//  Copyright (c) 2007, Pau Fernández
//

#include <cmath>

double sinc(double x)
{
  static const double half_pi=0.5*3.141592653589;
  return (x==0.0?1.0:sin(half_pi*x)/(half_pi*x));
}
double sinc_range() { return 15.0; }

double bspline(double x)
{
  if (x < -2.0) return 0.0;
  if (x < -1.0) return 0.5*(4.0+x*(8.0+x*(5.0+x)));
  if (x <  0.0) return 0.5*(2.0+x*x*(-5.0+x*-3.0));
  if (x <  1.0) return 0.5*(2.0+x*x*(-5.0+x*3.0));
  if (x <  2.0) return 0.5*(4.0+x*(-8.0+x*(5.0-x)));
  return 0.0;
}
double bspline_range() { return 2.0; }

double cosinus(double x)
{
  if (x < -1.0) return 0.0;
  if (x > 1.0) return 0.0;
  return .5*(1+cos(3.141592*x));
}
double cosinus_range() { return 2.0; }


void interp_enlarge(uint factor,
		    const vector< vector<double> >& M,
		    vector< vector<double> >& res,
		    double (*filter)(double),
		    double range)
{
  const uint XSZ=M[0].size(),YSZ=M.size();
  const uint NXSZ=XSZ*factor,NYSZ=YSZ*factor;
  
  res=vector< vector<double> >(NYSZ,vector<double>(NXSZ,0.0));

  for (uint j=0;j<YSZ;j++) {
    for (uint i=0;i<XSZ;i++) {
      for (uint _j=0;_j<NYSZ;_j++) {
	for (uint _i=0;_i<NXSZ;_i++) {
	  double x=double(_i)/double(factor)-double(i);
	  double y=double(_j)/double(factor)-double(j);
	  double f_x=(*filter)(x);
	  double f_y=(*filter)(y);
	  res[_j][_i]+=M[j][i]*f_x*f_y;
	}
      }
    }
  }
}

void interp_enlarge_centered(uint factor,
			     const vector< vector<double> >& M,
			     vector< vector<double> >& res,
			     double (*filter)(double),double range)
{
  const uint XSZ=M[0].size(),YSZ=M.size();
  const uint NXSZ=XSZ*factor,NYSZ=YSZ*factor;
  
  res=vector< vector<double> >(NYSZ,vector<double>(NXSZ,0.0));

  for (int j=-1;j<int(YSZ)+1;j++) {
    for (int i=-1;i<int(XSZ)+1;i++) {
      int _x=i,_y=j;
      if (_x < 0) _x=0; if (_x >= int(XSZ)) _x=XSZ-1;
      if (_y < 0) _y=0; if (_y >= int(YSZ)) _y=YSZ-1;
      double value=M[_y][_x];
      int irange=int(double(factor)*range);
      int xfrom=-irange+i*factor;
      int xto=i*factor+irange;
      int yfrom=-irange+j*factor;
      int yto=j*factor+irange;
      if (xfrom < 0) xfrom=0; if (xfrom >= int(NXSZ)) xfrom=NXSZ-1;
      if (xto < 0)   xto=0;   if (xto   >= int(NXSZ)) xto=NXSZ-1;
      if (yfrom < 0) yfrom=0; if (yfrom >= int(NYSZ)) yfrom=NYSZ-1;
      if (yto < 0)   yto=0;   if (yto   >= int(NYSZ)) yto=NYSZ-1;
      for (uint _j=yfrom;_j<=uint(yto);_j++) {
	for (uint _i=xfrom;_i<=uint(xto);_i++) {
	  double x=double(_i)/double(factor)-double(i)-.5;
	  double y=double(_j)/double(factor)-double(j)-.5;
	  double f_x=(*filter)(x);
	  double f_y=(*filter)(y);
	  res[_j][_i]+=value*f_x*f_y;
	}
      }
    }
  }
}

