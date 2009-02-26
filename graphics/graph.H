// 
//  Copyright (c) 2007, Pau Fernández
//

#ifndef _GRAPHICSGRAPH_H_
#define _GRAPHICSGRAPH_H_

#include "general.H"
#include "postscript.H"

using namespace graphics;

struct axis_prms
{
  bool reverse;
  double offset;

  font labelfont,ticksfont;
  string label;

  double start,stop;
  double sep;
  uint precision,minor_ticks;
  bool fixedpoint;

  double linewidth;
  color  rgb;
  double ticklength,major_minor_ratio;

  axis_prms() {
    // default values
    reverse=false;
    offset=0;
    labelfont=font("Times-Roman",8);
    ticksfont=font("Helvetica",5);
    label="";
    fixedpoint=true;
    start=0.0,stop=1.0,sep=0.1;
    precision=0;
    minor_ticks=0;
    linewidth=0.2;
    rgb=color(0,0,0);
    ticklength=1;
    major_minor_ratio=0.5;
  }
};


// Generica axes
void draw_axis(EPSarea& a,double width,double height,
	       point orig,const axis_prms& prms,bool blog=false,
	       bool mirror=false,bool swap=false);
// void draw_vertical_axis(EPSarea& a,double width,double height,
// 			point orig,const axis_prms& prms,bool blog=false,bool north=false);

// The four different axes
inline void draw_south_axis(EPSarea& a,double width,double height,point orig,const axis_prms& prms,bool blog=false)
{ draw_axis(a,width,height,orig,prms,blog,false,false); }
inline void draw_north_axis(EPSarea& a,double width,double height,point orig,const axis_prms& prms,bool blog=false)
{ draw_axis(a,width,height,orig,prms,blog,true,false); }
inline void draw_west_axis(EPSarea& a,double width,double height,point orig,const axis_prms& prms,bool blog=false)
{ draw_axis(a,width,height,orig,prms,blog,false,true); }
inline void draw_east_axis(EPSarea& a,double width,double height,point orig,const axis_prms& prms,bool blog=false)
{ draw_axis(a,width,height,orig,prms,blog,true,true); }

// Colorbar
void draw_colorbar(EPSarea& a,double width,double height,
		   point orig,const vector<color>& colormap,bool reverse=false,bool swap=false);

// Two different colorbars
inline void draw_vertical_colorbar(EPSarea& a,double width,double height,point orig,
			    const vector<color>& colormap,bool reverse=false)
{ draw_colorbar(a,width,height,orig,colormap,reverse,false); }
inline void draw_horizontal_colorbar(EPSarea& a,double width,double height,point orig,
			      const vector<color>& colormap,bool reverse=false)
{ draw_colorbar(a,width,height,orig,colormap,reverse,true); }

#endif

