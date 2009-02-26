// 
//  Copyright (c) 2007, Pau Fernández
//

#include <fstream>
#include <math.h>
#include <utils/param.H>
#include <utils/stl.H>
#include <graphics/postscript.H>
using namespace graphics;

#include "adj_list.H"
#include "stats.H"
#include "hierclust.H"
#include "io.H"

using namespace net;

bool read_similarity_measure(istream& i,
			     vector< map<uint,double> >& meas,
			     map<uint,string>& names)
{
  string n1,n2;
  double v;
  map<string,uint> _names;
  int curr_idx=0;

  while (i) {
    i >> n1 >> n2 >> v;
    if (!i) {
      if (i.eof()) break;
      else {
	cerr << "read_similarity_measure: Format must be '<name1> <name2> <value>'" << endl;
	return false;
      }
    }
    else {
      // first string
      uint idx1;
      map<string,uint>::iterator it=_names.find(n1);
      if (it!=_names.end()) idx1=it->second;
      else {
	idx1=curr_idx++;
	_names[n1]=idx1;
      }

      // second string
      uint idx2;
      it=_names.find(n2);
      if (it!=_names.end()) idx2=it->second;
      else {
	idx2=curr_idx++;
	_names[n2]=idx2;
      }

      if (meas.size() < 1+max(idx1,idx2)) meas.resize(1+max(idx1,idx2));
      map<uint,double>::iterator mi=meas[idx1].find(idx2);
      if (mi!=meas[idx1].end()) {
	cerr << "read_similarity_measure: " 
	     << "warning: repeated pair in similarity measure... " 
	     << "using most recent." << endl;
      }
      meas[idx1][idx2]=v;
    }
  }
  
  // reverse mapping
  names.clear(); 
  map<string,uint>::const_iterator mi=_names.begin(),miend=_names.end();
  for (;mi!=miend;++mi) {
    names[mi->second]=mi->first;
  }

  return true;
}

void prepare_image(const vector< map<uint,double> >& meas,
		   const vector<int>& ord,
		   vector< vector<double> >& image)
{
  const uint SZ=ord.size();

  image.clear();
  image.resize(ord.size());
  for (uint i=0;i<SZ;i++) image[i].resize(SZ);

  for (uint i=0;i<SZ;i++) {
    for (uint j=0;j<SZ;j++) {
      double val=0.0;
      map<uint,double>::const_iterator mi=meas[ord[i]].find(ord[j]);
      if (mi != meas[ord[i]].end()) val=mi->second;
      image[i][j]=1-val;
    }
  }
}

void print_pgm(ostream& o,const vector< map<uint,double> >& meas,
	       const vector<int>& ord,bool binverse)
{
  const uint SZ=ord.size();
  
  o << "P2" << endl << SZ << ' ' << SZ << endl << "256" << endl;
  for (uint i=0;i<SZ;i++) {
    for (uint j=0;j<SZ;j++) {
      double val=0.0;
      map<uint,double>::const_iterator mi=meas[ord[i]].find(ord[j]);
      if (mi != meas[ord[i]].end()) val=mi->second;
      o << int(255*(binverse?val:1.0-val)) << ' ';
    }
    o << endl;
  }
}

void print_pnm(ostream& o,const vector< map<uint,double> >& meas,
	       const vector<int>& ord,const vector<color>& colormap,
	       bool binverse)
{
  const uint SZ=ord.size();
  assert(colormap.size() >= 256);
  o << "P3" << endl << SZ << ' ' << SZ << endl << "256" << endl;
  for (uint i=0;i<SZ;i++) {
    for (uint j=0;j<SZ;j++) {
      color c(1,1,1);
      if (binverse) c=color(0,0,0);
      map<uint,double>::const_iterator mi=meas[ord[i]].find(ord[j]);
      if (mi != meas[ord[i]].end()) c=colormap[255-int(mi->second*255)];
      uint r=uint(c.x*255),g=uint(c.y*255),b=uint(c.z*255);
      o << r << ' ' << g << ' ' << b << ' ';
    }
    o << endl;
  }
}

void draw_dendrogram(branch<int>* b,const map<uint,string>& names,EPSarea& eps,
		     point topleft,point bottomright)
{
  const uint N=b->size();
  assert(N > 0);
  if (N == 1) {
    map<uint,string>::const_iterator it=names.find(b->data());
    if (it!=names.end()) {
      string s=it->second;
      eps.text(point(.7*bottomright.x+.3*topleft.x,
		     bottomright.y+0.15*(topleft.y-bottomright.y)),s);
    }
  }
  else {
    const double height=topleft.y-bottomright.y;
    const double width=bottomright.x-topleft.x;

    //    const point mid(bottomright.x,(topleft.y+bottomright.y)/2.0);
    const uint rsz=b->right()->size(),lsz=b->left()->size();
    const double rfrac=double(rsz)/double(rsz+lsz),lfrac=double(lsz)/double(rsz+lsz);
    const double rh=rfrac*height,lh=lfrac*height;
    const double rw=rfrac*width,lw=lfrac*width;
    
    const point rmid(bottomright.x-lw,topleft.y-rh/2.0);
    const point lmid(bottomright.x-rw,bottomright.y+lh/2.0);

    vector<point> pp;
    pp.push_back(rmid);
    pp.push_back(point(bottomright.x,rmid.y));
    pp.push_back(point(bottomright.x,lmid.y));
    pp.push_back(lmid);
    eps.line(pp);

    draw_dendrogram(b->right(),names,eps,topleft,point(rmid.x,topleft.y-rh));
    draw_dendrogram(b->left(),names,eps,point(topleft.x,bottomright.y+lh),
		    point(lmid.x,bottomright.y));
  }
}

void DoubleToOcean(double& r,double& g,double& b, double val)
{
  r=1-val*3;   if (r < 0.0) r=0.0;
  g=1-val*3/2; if (g < 0.0) g=0.0;
  b=1-val;
}

void ocean_colormap(vector<color>& colormap,bool binverse=false)
{
  for (double v=0.0;v<1.0;v+=1.0/256.0) {
    color c;
    DoubleToOcean(c.x,c.y,c.z,1-v);
    if (binverse) c=color(1-c.x,1-c.y,1-c.z);
    colormap.push_back(c);
  }
}

int main(int argc,char** argv)
{
  double WF=0.3,linewidth=0.01,textwidth=10.0;
  double eps_height=20.0;
  bool binverse=false;
  bool boceancolormap=false;

  vector<param*> prms;
  prms.push_back(make_param('I',"inverse_colors",binverse));
  prms.push_back(make_param('d',"dendro_width_factor",WF));
  prms.push_back(make_param('l',"linewidth",linewidth));
  prms.push_back(make_param('t',"textwidth",textwidth));
  prms.push_back(make_param('e',"eps_height_in_cm",eps_height));
  prms.push_back(make_param('o',"oceancolormap",boceancolormap));

  vector<string> args;
  string usage =
    "Usage: Ghierclust <command> <infile(.mtx)> <outfile>\n"
    "Copyright (c) 2007, Pau Fernandez\n\n"
    "  Compute the hierarchical clustering from a matrix\n"
    "  command:\n"
    "    order      Print the order\n"
    "    image      Generate an image (.pnm)\n"
    "    dendro     Print the dendrogram\n"
    "    diagram    Print the whole diagram";

  parse_params_ex(prms,argc,argv,usage,"Ghierclust",args,3);
  string infile=args[1],outfile=args[2],command=args[0];

  // Read similarity measure
  map<uint,string> names;
  vector< map<uint,double> > meas;
  ifstream fin(infile.c_str());
  if (!fin.is_open()) {
    cerr << "Ghierclust: File " << infile << " doesn't seem to exist..." << endl;
    return -1;
  }
  if (!read_similarity_measure(fin,meas,names)) {
    cerr << "Problems reading similarity measure." << endl;
    return -1;
  }

  // Execute hierarchical clustering
  vector<branch<int>*> roots;
  hierarchical_clustering(roots,meas);
  
  vector<int> order;
  uint max_size=0;
  for (uint k=0;k<roots.size();k++) {
    vector<int> ord;
    extract_order(roots[k],ord);
    copy(ord.begin(),ord.end(),back_inserter(order));
    max_size=std::max(max_size,uint(roots[k]->size()));
  }

  ofstream fout(outfile.c_str());
  if (!fout.is_open()) {
    cerr << "Coudln't open '" << outfile << "' for writing." << endl;
    return -1;
  }

  // Carry out commands
  if (command=="image") {
    if (boceancolormap) {
      vector<color> cm;
      ocean_colormap(cm,binverse);
      print_pnm(fout,meas,order,cm,binverse);
    }
    else {
      print_pgm(fout,meas,order,binverse);
    }
    return 0;
  }
  
  if (command=="order") {
    vector<int>::const_iterator oi=order.begin(),oiend=order.end();
    for (;oi!=oiend;++oi) {
      fout << names[*oi] << endl;
    }
    return 0;
  }

  if (command=="dendro") {
    uint totalsz=names.size();
    double fontsize=0.85/double(totalsz);
    double wtext=fontsize*textwidth;
    double curr=0.005;
    double width_factor=WF*double(max_size)/double(totalsz);
    double height_in_cm=eps_height;
    
    EPSarea eps((width_factor+wtext)*height_in_cm,height_in_cm,
		width_factor+wtext,1.0,outfile);
    eps.setlinewidth(linewidth/height_in_cm);
    eps.setfont(font("Times-Roman",fontsize));
    for (uint k=0;k<roots.size();k++) {
      uint sz=roots[k]->size();
      double frac=double(sz)/double(totalsz)*0.99;
      draw_dendrogram(roots[k],names,eps,
		      point(0.99*width_factor,curr),
		      point(0.01*width_factor+(1-frac)*0.98*width_factor,curr+frac));
      curr+=frac;
    }
  }

  if (command=="diagram") {
    uint totalsz=names.size();
    double curr=0.0;
    double dendro_width_factor=WF*double(max_size)/double(totalsz);
    double height_in_cm=eps_height;
    double width_in_cm=(1+dendro_width_factor)*height_in_cm;
    
    EPSarea eps(width_in_cm,height_in_cm,
		1+dendro_width_factor+0.01,1.0,outfile);
    eps.setlinewidth(linewidth/height_in_cm);
    if (binverse) {
      eps.setcolor(color(0,0,0));
      eps.fillbox(point(0,0),point(1+dendro_width_factor+0.01,1.0));
      eps.setcolor(color(1,1,1)); // white lines & text
    }
    for (uint k=0;k<roots.size();k++) {
      uint sz=roots[k]->size();
      double frac=double(sz)/double(totalsz);
      map<uint,string> nonames;
      draw_dendrogram(roots[k],nonames,eps,
		      point(1.0,curr),
		      point(dendro_width_factor*frac+1.0,curr+frac));
      curr+=frac;
    }

    vector< vector<double> > _image;
    prepare_image(meas,order,_image);
    
    if (boceancolormap) {
      vector<color> cm;
      ocean_colormap(cm,binverse);
      eps.colorimage(point(1.0,0.0),point(0.0,1.0),_image,cm);
    }
    else {
      eps.bnimage(point(1.0,0.0),point(0.0,1.0),_image);
    }
  }

  if (command == "sizes") {
    for (uint k=0;k<roots.size();k++) 
      {
	list<int> sizes;
	extract_sizes(roots[k],sizes);
	list<int>::const_iterator li=sizes.begin(),liend=sizes.end();
	for (;li!=liend;++li) {
	  fout << *li << endl;
	}
	fout << endl;
      }
  }
}
