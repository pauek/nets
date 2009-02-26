// 
//  Copyright (c) 2007, Pau Fernández
//

#ifndef _POSTSCRIPT_H_
#define _POSTSCRIPT_H_

#include <fstream>
#include <sstream>
#include <vector>
#include "general.H"

using namespace std;

namespace graphics {

class EPSarea
{
  double _xcm,_ycm;
  double _xsz,_ysz;
//   ostringstream _tmp;
  ofstream _out;
  string _filename,_tmpfile;

  typedef enum { LINE,CIRCLE,BOX,IMAGE,MAXPARTS } parts;

  vector<bool> _active;
  int  _currpage;
  bool _multipage;
  static string _preamble_part[];
  
public:
  EPSarea(double xcm,double ycm,double xsz,double ysz,string filename,bool multipage=false)
    :_xcm(xcm),_ycm(ycm),_xsz(xsz),_ysz(ysz),_filename(filename),
     _active(MAXPARTS,false),_currpage(1),_multipage(multipage)
  { start(); }

  ~EPSarea() 
  { end(); }

  double xsize() const { return _xsz; }
  double ysize() const { return _ysz; }

  void setgray(double g) { _out << g << " sg "; }
  void setcolor(color c) { _out << c << (c.bHSV?" sh ":" sr "); }
  void setlinewidth(double w) { _out << w << " sw "; }

  void setclipbox(point p1,point p2);
  void removeclipbox();

  void line(point p1,point p2);
  void line(const vector<point>& p);

  void curve(point p1,point p2,point p3,point p4);
  void curve(const vector<point>& p);

  void box(point p1,point p2);
  void fillbox(point p1,point p2);

  void circle(point p,double r);
  void fillcircle(point p,double r);

  void bnimage(point p1,point p2,
	       const vector< vector<double> >& M);
  void colorimage(point p1,point p2,
		  const vector< vector<double> >& M,
		  const vector<color>& colormap);

  void setfont(font f);
  void text(point p,string s);
  void centertext(point p,string s);
  void righttext(point p,string s);

  void start();
  void end();
  
  void next_page();
};

}

#endif
