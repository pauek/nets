// 
//  Copyright (c) 2007, Pau Fernández
//

#ifndef _GRAPHICSGENERAL_H_
#define _GRAPHICSGENERAL_H_

#include <iostream>
#include <string>
#include <assert.h>
#include <math.h>
// #include "../utils/coord.H"
using namespace std;

namespace graphics {

//////////////////////////////////////////////////////////////////
// long RGB packing

inline long RGB(uint r,uint g,uint b) { return ((r&0xFF) << 16) | ((g&0xFF) << 8) | (b&0xFF); }
inline uint getR(long rgb) { return ((rgb >> 16) &  0xFF); }
inline uint getG(long rgb) { return ((rgb >> 8) & 0xFF); }
inline uint getB(long rgb) { return (rgb & 0xFF); }

////////////////////////////////////////////////////////////////////
// point

  struct point { 
    double x,y;
    point(double _x=0.0,double _y=0.0):x(_x),y(_y) {}
    
    point operator+(const point& p) const { return point(x+p.x,y+p.y); }
    point& operator+=(const point& p) { x+=p.x; y+=p.y; return *this; }
    point operator-(const point& p) const { return point(x-p.x,y-p.y); }
    point& operator-=(const point& p) { x-=p.x; y-=p.y; return *this; }
  };
  
inline ostream& operator<<(ostream& o,point p)
{ return o << p.x << ' ' << p.y; }

// typedef coord2D point;

///////////////////////////////////////////////////////////
// Color con indices entre 0.0 y 1.0

struct color;
color HSVtoRGB(color c);

struct color
{ 
  double x,y,z; 
  bool bHSV;
  color(double _x=0.0,double _y=0.0,double _z=0.0,bool HSV=false)
    :x(_x),y(_y),z(_z),bHSV(HSV) { }

  color operator+(const color& c) const 
  {
    color tmp=c;
    tmp.x+=x,tmp.y+=y,tmp.z+=z;
    return tmp;
  }

  color toRGB() const { assert(bHSV); return HSVtoRGB(*this); }
  long toLong() const { return RGB(uint(x*255.0),uint(y*255.0),uint(z*255.0)); }
};

inline ostream& operator<<(ostream& o,color c)
{ return o << c.x << ' ' << c.y << ' ' << c.z; }

/////////////////////////////////////////////////////////
// Especificacion de fuente

struct font
{
  string typeface;
  double size;
  font(string tf="Times-Roman",double s=10)
    :typeface(tf),size(s) {}
};

//////////////////////////////////////////////////////////////
// Rutina HSV -> RGB

inline color HSVtoRGB(color c)
{
  assert(c.bHSV);
  double h=c.x*360.0,s=c.y,v=c.z;
  double f,p,q,t,i;
  int a;
    
  if (s==0.0) {
    // achromatic (grey)
    return color(v,v,v);
  }
    
  h /= 60;                        // sector 0 to 5
  i=floor( h );
  a=(int)i;
  f=h-i;                      // factorial part of h
  p=v*(1-s);
  q=v*(1-s*f);
  t=v*(1-s*(1-f));

  double r,g,b;
  switch( a ) {
  case 0: r = v, g = t, b = p; break;
  case 1: r = q, g = v, b = p; break;
  case 2: r = p, g = v, b = t; break;
  case 3: r = p, g = q, b = v; break;
  case 4: r = t, g = p, b = v; break;
  case 5: r = v, g = p, b = q; break;
  default: r = g = b = 0; break;
  }
  return color(r,g,b);
}

}

#endif
