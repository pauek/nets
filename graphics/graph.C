// 
//  Copyright (c) 2007, Pau Fernández
//

#include <iostream>
#include <iomanip>
#include "graph.H"

using namespace std;

void ticks(double start,double stop,double sep,uint minor_ticks,uint precision,bool bfixed,
	   vector<double>& major,vector<double>& minor,vector<string>& labels)
{
  double v=floor(start/sep)*sep;
  while (v <= stop) {
    if (v >= start) {
      // major
      major.push_back(v);
      ostringstream ostr;
      if (bfixed) ostr.setf(ios::fixed,ios::floatfield);
      ostr.precision(precision);
      ostr << v;
      labels.push_back(ostr.str());
    }
    // minor
    double minsep=sep/double(minor_ticks+1);
    for (uint k=1;k<=minor_ticks;k++) {
      double vv=v+minsep*double(k);
      if (vv > start && vv < stop) minor.push_back(vv);
    }
    v+=sep;
  }
}

void log_ticks(double start,double stop,double sep,uint minor_ticks,uint precision,bool bfixed,
	       vector<double>& major,vector<double>& minor,vector<string>& labels)
{
  double _start=log10(start),_stop=log10(stop);
  double _sep=log10(sep);
  double v=floor(_start/_sep)*_sep;
  while (v <= _stop+1e-5) {
    if (v >= _start-1e-5) {
      // major
      major.push_back(v);
      ostringstream ostr;
      if (bfixed) ostr.setf(ios::fixed,ios::floatfield);
      ostr.precision(precision);
      ostr << pow(10.0,v);
      labels.push_back(ostr.str());
    }
    // minor
    double _v=pow(10.0,v);
    double minsep=_v*(sep-1)/double(minor_ticks+1);
    for (uint k=1;k<=minor_ticks;k++) {
      double vv=log10(_v+minsep*double(k));
      if (vv > _start && vv < _stop) minor.push_back(vv);
    }
    v+=_sep;
  }
}

inline double translate(double val,double x,double w,
			double start,double stop,bool blog)
{ 
  if (!blog) 
    return (val-start)/(stop-start)*w+x; 
  else {
    double _start=log10(start),_stop=log10(stop);
    return (val-_start)/(_stop-_start)*w+x; 
  }
}

class mirror_t
{
  double _x1,_x2;
  double _y1,_y2;
  bool _xact,_yact,_swap;
public:
  mirror_t() {}
  mirror_t(point p1,point p2,bool xact,bool yact,bool swap)
    :_x1(p1.x),_x2(p2.x),_y1(p1.y),_y2(p2.y),
     _xact(xact),_yact(yact),_swap(swap) {}
  
  inline point operator()(point p)
  { 
    point tmp=(_swap?point(p.y,p.x):p);
    return point(_xact ? _x2-(tmp.x-_x1) : tmp.x,
		 _yact ? _y2-(tmp.y-_y1) : tmp.y); 
  }
};

///////////////////////////////////////////////////////////////////////////////
// General axis

void draw_axis(EPSarea& a,double width,double height,
	       point orig,const axis_prms& prms,bool blog,
	       bool mirr,bool swap)
{
  bool _mirrx=(swap?!mirr:prms.reverse),_mirry=(swap?prms.reverse:mirr);
  point  _orig=(swap?point(orig.y,orig.x):orig);
  double _height=(swap?width:height),_width=(swap?height:width);
  mirror_t mirror(orig,point(orig.x+width,orig.y+height),_mirrx,_mirry,swap);
  double _curr=prms.offset+_orig.y;
  
  // axis bar
  a.setcolor(prms.rgb);
  a.setlinewidth(prms.linewidth);
  a.line(mirror(point(_orig.x,_curr)),
	 mirror(point(_orig.x+_width,_curr)));

  vector<double> M,m;
  vector<string> L;
  if (!blog) {
    ticks(prms.start,prms.stop,
	  prms.sep,prms.minor_ticks,
	  prms.precision,prms.fixedpoint,
	  M,m,L);
  }
  else {
    log_ticks(prms.start,prms.stop,
	      prms.sep,prms.minor_ticks,
	      prms.precision,prms.fixedpoint,
	      M,m,L);
  }

  // Ticks
  vector<double>::const_iterator Mi=M.begin(),Mend=M.end();
  for (;Mi!=Mend;Mi++) {
    double _c=translate(*Mi,_orig.x,_width,prms.start,prms.stop,blog);
    a.line(mirror(point(_c,_curr)),
	   mirror(point(_c,_curr+prms.ticklength)));
  }
  vector<double>::const_iterator mi=m.begin(),mend=m.end();
  for (;mi!=mend;mi++) {
    double _c=translate(*mi,_orig.x,_width,prms.start,prms.stop,blog);
    a.line(mirror(point(_c,_curr)),
	   mirror(point(_c,_curr+prms.ticklength*prms.major_minor_ratio)));
  }
  _curr+=prms.ticklength*1.50;

  // Tick labels
  //  a.setcolor(color(0,0,0));
  a.setfont(prms.ticksfont);
  double fontthick=prms.ticksfont.size*(swap?1.5:1.0);
  vector<string>::const_iterator Li=L.begin(),Lend=L.end();
  for (Mi=M.begin();Li!=Lend;Li++,Mi++) {
    double _c=translate(*Mi,_orig.x,_width,prms.start,prms.stop,blog);
    double corr_m=0.0,corr_s=0.0;
    if (!mirr && !swap) corr_m=fontthick*.75;
    if (swap) corr_s=-0.19*fontthick;
    point tmp=mirror(point(_c+corr_s,_curr+corr_m));
    if (swap)  
      if (mirr) a.text(tmp,*Li);
      else      a.righttext(tmp,*Li);
    else
      a.centertext(tmp,*Li);
  }
  _curr+=prms.ticksfont.size*1.1;
  
  // Title
  a.setfont(prms.labelfont);
  double remaining=(_orig.y+_height)-_curr;
  if (remaining < prms.labelfont.size) {
    cerr << "draw_horizontal_axis: space for axis label too thin" << endl;

    point tmp=mirror(point(_orig.x+_width/2.0,_curr));
    a.centertext(tmp,prms.label);
  }
  else {
    remaining-=prms.labelfont.size;
    double off=remaining/2.0;
    double corr_m=0.0;
    if (!mirr && !swap) corr_m=fontthick;
    if (!swap) {
      point tmp=mirror(point(_orig.x+_width/2.0,_curr+off+corr_m));
      a.centertext(tmp,prms.label);
    }
    else {
      // TODO!!!
    }
  }
}

void draw_colorbar(EPSarea& a,double width,double height,point orig,
		   const vector<color>& colormap,bool reverse,bool swap)
{
  bool _mirrx=(swap?reverse:false),_mirry=(swap?false:!reverse);
  point  _orig=(swap?point(orig.y,orig.x):orig);
  double _height=(swap?width:height),_width=(swap?height:width);
  mirror_t mirror(orig,point(orig.x+width,orig.y+height),_mirrx,_mirry,swap);

  const uint SZ=colormap.size();
  double curr=_orig.y,sep=_height/double(SZ);
  vector<color>::const_iterator cmi=colormap.begin(),cmend=colormap.end();
  for (uint k=0;cmi!=cmend;++cmi,++k) {
    a.setcolor(*cmi);
    a.fillbox(mirror(point(_orig.x,curr)),
	      mirror(point(_orig.x+_width,(k!=(SZ-1)?curr+sep*1.5:curr+sep))));
    curr+=sep;
  }
}

