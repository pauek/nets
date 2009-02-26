// 
//  Copyright (c) 2007, Pau Fernández
//

#include <fstream>
#include <graphics/postscript.H>
#include <graphics/graph.H>
#include <utils/param.H>
#include <utils/interp.H>

#include "adj_list.H"
#include "stats.H"
#include "io.H"

using namespace net;

int main(int argc,char** argv)
{
  uint bpd=4;
  uint IT=1000;
  uint factor=5;
  bool bzscore=false,bsymm=false;
  bool bnolog=false,bcolor=false,binverse=false;
  string epsfile="<none>";

  vector<param*> prms;
  prms.push_back(make_param('l',"no_logbin",bnolog));
  prms.push_back(make_param('s',"symmetric",bsymm));
  prms.push_back(make_param('b',"bins_per_decade",bpd));
  prms.push_back(make_param('I',"iterations",IT));
  prms.push_back(make_param('z',"zscore",bzscore));
  prms.push_back(make_param('f',"enlarge_factor",factor));
  prms.push_back(make_param('e',"eps_file",epsfile));
  prms.push_back(make_param('c',"eps_color",bcolor));
  prms.push_back(make_param('i',"inverse",binverse));

  string usage = 
    "Usage: Gmaslov <input_graph_file> <output_matrix>\n"
    "Copyright (c) 2007, Pau Fernandez";
  vector<string> args;
  parse_params_ex(prms,argc,argv,usage,"Gmaslov",args,2);

  string grfile=args[0],mtxfile=args[1];

  typedef adj_list<std::string,std::string> graph;
  graph g;
  if (!read_graph_from_file(g,grfile)) {
    cerr << "Couldn't read graph " << grfile << endl;
    return -1;
  }

  vector< vector<double> > P,Z;
  progress_bar pb(0,IT);
  corr_matrix_vs_random(g,bpd,IT,P,Z,!bnolog,bsymm,&pb);

  vector< vector<double> > M=(bzscore?Z:P);
  if (epsfile!="<none>") {
    // enlarge + min/max
    vector< vector<double> > M2(M.size()*factor,vector<double>(M[0].size()*factor,0.0));
    interp_enlarge_centered(factor,M,M2,cosinus,cosinus_range());
    double max=M2[0][0],min=M2[0][0];
    for (uint j=0;j<M2.size();j++) {
      for (uint i=0;i<M2[0].size();i++) {
	if (max < M2[j][i]) max=M2[j][i];
	if (min > M2[j][i]) min=M2[j][i];
      }
    }

    // normalize
    vector< vector<double> > M3=M2;
    for (uint j=0;j<M2.size();j++) {
      for (uint i=0;i<M2[0].size();i++) {
	if (max-min > 0.0) M3[M2.size()-1-j][i]=1.0-(M2[j][i]-min)/(max-min);
      }
    }
    M2=M3;
    
    double waxes=15.0,img_size=85.0,wbar=10.0,sep=5.0;
    double mrg=10.0;
    double img_xsz=img_size,img_ysz=img_size;
    if (M[0].size() > M.size()) 
      img_ysz*=double(M.size())/double(M[0].size());
    else if (M[0].size() < M.size())
      img_xsz*=double(M[0].size())/double(M.size());
    
    double _xsz=waxes+img_xsz+sep+wbar+waxes;
    double _ysz=mrg+img_ysz+waxes;
    EPSarea a(_xsz/10.0,_ysz/10.0,_xsz,_ysz,epsfile);
    if (binverse) {
      a.setcolor(color(0,0,0));
      a.fillbox(point(0,0),point(_xsz,_ysz));
    }
    vector<color> cm;
    for (uint k=0;k<256;k++) {
      if (bcolor)
	cm.push_back(HSVtoRGB(color(double(k)/255.0*.666,1.0,1.0,true)));
      else {
	double v=double(k)/255.0;
	cm.push_back(color(v,v,v));
      }
    }
    if (bcolor) 
      a.colorimage(point(waxes,mrg),point(waxes+img_xsz,mrg+img_ysz),M2,cm);
    else 
      a.bnimage(point(waxes,mrg),point(waxes+img_xsz,mrg+img_ysz),M2);

    axis_prms ap;
    ap.offset=0.5;
    ap.ticklength=2.0;
    //    ap.label="K";
    ap.labelfont=font("Times-Roman",5.5);
    ap.ticksfont=font("Helvetica",4);
    ap.fixedpoint=true;
    ap.precision=1;
    double _bpd=bpd;
    while (_bpd > 3) _bpd/=2.0;
    ap.sep=pow(10.0,1/double(_bpd));
    ap.minor_ticks=0;
    ap.start = 1;
    ap.stop  = pow(10.0,(double(M2[0].size())/double(factor))/double(bpd));
    if (binverse) ap.rgb=color(1.0,1.0,1.0);

    draw_south_axis(a,img_xsz,waxes,point(waxes,mrg+img_ysz),ap,true);

    ap.reverse=true;
    ap.start = 1;
    ap.stop  = pow(10.0,(double(M2.size())/double(factor))/double(bpd));
    draw_west_axis(a,waxes,img_ysz,point(0,mrg),ap,true);

    draw_vertical_colorbar(a,wbar,img_ysz,point(waxes+img_xsz+sep,mrg),cm,true);

    ap.sep=2.0;
    while (ap.sep/(max-min) > .2) ap.sep/=2.0;
    ap.start=min,ap.stop=max;
    ap.minor_ticks=1;
    ap.precision=(ap.sep < 1.0 ? 2:1);
    draw_east_axis(a,waxes,img_ysz,point(waxes+img_xsz+sep+wbar,mrg),ap,false);
  }

  {
    ofstream mtx(mtxfile.c_str());
    for (uint j=0;j<M.size();j++) {
      for (uint i=0;i<M[0].size();i++) {
	mtx << M[j][i] << ' ';
      }
      mtx << endl;
    }
  }
}
