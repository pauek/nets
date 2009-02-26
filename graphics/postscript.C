// 
//  Copyright (c) 2007, Pau Fernández
//

#include <assert.h>
#include <iomanip>
#include "postscript.H"

// S = stroke  F = fill
// d = def     c = circle  C = "fillcircle"
// m = moveto  l = lineto  gr = gsave   gr = grestore
// sr = setrgbcolor sh = sethsbcolor b = "box"
// sc = scale   tr = translate   rl = rlineto
// sb = "setclipbox"
// rs = rectstroke   rf = rectfill    rc = rectclip

using namespace graphics;

string EPSarea::_preamble_part[]={
  "/m{moveto}d/l{lineto}d/rl{rlineto}d/k{curveto}d", // LINE
  "/cc{0 360 arc}d/c{cc S}bind d/C{cc F}bind d", // CIRCLE
  "/sq{np m exch dup neg 0 0 5 3 roll 0 rl rl rl cp}d"
  "/rs{rectstroke}d/rf{rectfill}d/rc{rectclip}d", // BOX
  "/cf{currentfile}d/rhs{readhexstring}d/ci{colorimage}d/im{image}d", //IMAGE
};

void EPSarea::line(point p1,point p2)
{
  _active[LINE]=true;
  _out << p1.x << ' ' << p1.y << ' ' << "m "
       << p2.x << ' ' << p2.y << " l S ";
}

void EPSarea::line(const vector<point>& p)
{
  _active[LINE]=true;
  vector<point>::const_iterator pi=p.begin(),piend=p.end();
  _out << pi->x << ' ' << pi->y << ' ' << "m ";
  for (++pi;pi!=piend;pi++) {
    _out << pi->x << ' ' << pi->y << ' ' << "l ";
  }
  _out << "S ";
}

void EPSarea::curve(point p1,point p2,point p3,point p4)
{
  _active[LINE]=true;
  _out << p1.x << ' ' << p1.y << ' ' << "m "
       << p2.x << ' ' << p2.y << ' '
       << p3.x << ' ' << p3.y << ' '
       << p4.x << ' ' << p4.y << " k S ";
}

void EPSarea::curve(const vector<point>& p)
{
  _active[LINE]=true;
  assert(p.size() == 4);
  vector<point>::const_iterator pi=p.begin(),piend=p.end();
  _out << pi->x << ' ' << pi->y << ' ' << "m ";
  for (++pi;pi!=piend;pi++) {
    _out << pi->x << ' ' << pi->y << ' ';
  }
  _out << " k S ";
}


void EPSarea::circle(point p, double r)
{
  _active[CIRCLE]=true;
  _out << p.x << ' ' << p.y << ' ' << r << " c " ;
}

void EPSarea::fillcircle(point p, double r)
{
  _active[CIRCLE]=true;
  _out << p.x << ' ' << p.y << ' ' << r << " C ";
}

void EPSarea::box(point p1,point p2)
{
  _active[LINE]=true,_active[BOX]=true;
  _out << p1.x << ' ' << p1.y << ' ' 
       << p2.x-p1.x << ' ' << p2.y-p1.y << " rs ";
}

void EPSarea::fillbox(point p1,point p2)
{
  _active[LINE]=true,_active[BOX]=true;
  _out << p1.x << ' ' << p1.y << ' '
       << p2.x-p1.x << ' ' << p2.y-p1.y << " rf ";
}

void EPSarea::setclipbox(point p1,point p2)
{
  _active[LINE]=true,_active[BOX]=true;
  _out << "clipsave "
       << p1.x << ' ' << p1.y  << ' ' 
       << p2.x-p1.x << ' ' << p2.y-p1.y << " rc ";
}

void EPSarea::removeclipbox()
{ _out << "cliprestore "; }

char int2hex[16]={'0','1','2','3','4','5','6','7',
		  '8','9','a','b','c','d','e','f'};

void print_hex(ofstream& o,const vector<uint>& v)
{ 
  vector<uint>::const_iterator vi=v.begin(),viend=v.end();
  string s;
  for (;vi!=viend;++vi) {
    assert(*vi >= 0 && *vi < 256);
    s+=int2hex[(*vi)/16];
    s+=int2hex[*vi % 16];
  }
  o << s;
}

void EPSarea::bnimage(point p1,point p2,
		      const vector< vector<double> >& M)
{
  _active[IMAGE]=true;
  const uint YSZ=M.size(),XSZ=M[0].size();
  _out << endl << "% B/N Image" << endl;
  _out << "gs " << endl
       << p1.x << ' ' << p1.y << " tr "
       << p2.x-p1.x << ' ' << p2.y-p1.y << " sc "
       << XSZ << ' ' << YSZ << " 8 "
       << '[' << XSZ << " 0 0 " << YSZ << " 0 0] "
       << "{ cf " << XSZ << " string rhs pop } im " << endl;
  vector< vector<double> >::const_iterator mi=M.begin(),miend=M.end();
  for (;mi!=miend;++mi) {
    vector<double>::const_iterator vi=mi->begin(),viend=mi->end();
    for (;vi!=viend;++vi) {
      vector<uint> v(1,uint((*vi)*255.0));
      print_hex(_out,v);
    }
    _out << endl;
  }
  _out << "gr ";
}

void EPSarea::colorimage(point p1,point p2,
			 const vector< vector<double> >& M,
			 const vector<color>& colormap)
{
  _active[IMAGE]=true;
  const uint YSZ=M.size(),XSZ=M[0].size();
  _out << endl << "% Color Image" << endl;
  _out << "gs " << endl
       << p1.x << ' ' << p1.y << " tr "
       << p2.x-p1.x << ' ' << p2.y-p1.y << " sc "
       << XSZ << ' ' << YSZ << " 8 "
       << '[' << XSZ << " 0 0 " << YSZ << " 0 0] "
       << "{ cf " << 3*XSZ << " string rhs pop } false 3 ci" << endl;
  vector< vector<double> >::const_iterator mi=M.begin(),miend=M.end();
  for (;mi!=miend;++mi) {
    vector<double>::const_iterator vi=mi->begin(),viend=mi->end();
    for (;vi!=viend;++vi) {
      if (!(*vi <= 1.0 && *vi >= 0.0)) {
	cerr << *vi << endl;
	assert(*vi <= 1.0 && *vi >= 0.0);
      }
      uint idx=uint((*vi)*255.0);
      vector<uint> v;
      v.push_back(uint(colormap[idx].x*255.0));
      v.push_back(uint(colormap[idx].y*255.0));
      v.push_back(uint(colormap[idx].z*255.0));
      print_hex(_out,v);
    }
    _out << endl;
  }
  _out << "gr ";
}

void EPSarea::setfont(font f)
{
  _out << "/" << f.typeface << " findfont " 
       << f.size << " scalefont setfont ";
}

void EPSarea::text(point p,string s)
{
  _out << p << " m gs 1 -1 sc (" << s << ") show gr ";
}

void EPSarea::centertext(point p,string s)
{
  _out << p << " m gs 1 -1 sc (" << s << ") dup stringwidth pop 2 div neg 0 rmoveto show gr ";
}

void EPSarea::righttext(point p,string s)
{
  _out << p << " m gs 1 -1 sc (" << s << ") dup stringwidth pop neg 0 rmoveto show gr ";
}

void EPSarea::start()
{
  // 1cm = 28.3464 pt
  double xptsize=_xcm*28.3464,yptsize=_ycm*28.3464;

  _out.open(_filename.c_str());
  // Header
  _out << "%!PS-Adobe-3.0" << endl
       << "%%BoundingBox: " << 0 << ' ' << 0 << ' ' 
       << int(xptsize) << ' ' << int(yptsize) << endl
       << "%%EndComments" << endl;
  // Procedure Definitions
  _out << "%%BeginProlog" << endl
       << "/d{def}def/S{stroke}d/sr{setrgbcolor}d/sh{sethsbcolor}d/sg{setgray}d" << endl
       << "/gs{gsave}d/gr{grestore}d/sc{scale}d/tr{translate}d/F{fill}d" << endl
       << "/np{newpath}d/cp{closepath}d/sw{setlinewidth}d" << endl;

  for (uint k=0;k<MAXPARTS;k++) {
    //     if (_active[k]) _out << _preamble_part[k] << endl;
    _out << _preamble_part[k] << endl;
  }

  // Change coordinate mapping
  _out << "0 " << yptsize << " tr " 
       << xptsize/_xsz << ' ' << -yptsize/_ysz << " sc" << endl;

  _out << "%%EndProlog" << endl;
  
  if (_multipage) {
    _out << "%%Page: " << _currpage << ' ' << _currpage << endl;
  }
}

void EPSarea::next_page()
{
  assert(_multipage);
  ++_currpage;
  _out << "%%PageTrailer" << endl;
  _out << "%%Page: " << _currpage << ' ' << _currpage << endl;
}

void EPSarea::end()
{
  _out << endl;
  if (_multipage) _out << "%%PageTrailer" << endl;
  _out << "%%Trailer" << endl
       << "%%EOF" << endl;
  _out.close();
}
