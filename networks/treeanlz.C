// 
//  Copyright (c) 2007, Pau Fernández
//

#include <fstream>
#include <cmath>
using namespace std;

#include <utils/param.H>
#include <utils/stl.H>
#include <graphics/postscript.H>
using namespace graphics;

#include "adj_list.H"
#include "stats.H"
#include "io.H"
#include "community.H"
using namespace net;

typedef adj_list<string,string> graph;
typedef adj_list<pair<string,double>,double> bgraph;
typedef adj_list<tree_tag,null>   resgraph;

bool read_measure(istream& i,
		  vector< map<int,double> >& meas,
		  map<string,int>& _names)
{
  string n1,n2;
  double v;
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
      int idx1;
      map<string,int>::iterator it=_names.find(n1);
      if (it!=_names.end()) idx1=it->second;
      else {
	idx1=curr_idx++;
	_names[n1]=idx1;
      }

      // second string
      int idx2;
      it=_names.find(n2);
      if (it!=_names.end()) idx2=it->second;
      else {
	idx2=curr_idx++;
	_names[n2]=idx2;
      }

      if (int(meas.size()) < 1+max(idx1,idx2)) meas.resize(1+max(idx1,idx2));
      map<int,double>::iterator mi=meas[idx1].find(idx2);
      if (mi!=meas[idx1].end()) {
	cerr << "read_similarity_measure: " 
	     << "warning: repeated pair in similarity measure... " 
	     << "using most recent." << endl;
      }
      meas[idx1][idx2]=v;
    }
  }
  
  return true;
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

void print_pnm(ostream& o,const vector< map<int,double> >& meas,
	       const map<string,int>& names,
	       const vector<string>& ord,const vector<color>& colormap,
	       bool binverse)
{
  const uint SZ=ord.size();
  assert(colormap.size() >= 256);
  o << "P3" << endl << SZ << ' ' << SZ << endl << "256" << endl;
  vector<string>::const_iterator oi=ord.begin(),oiend=ord.end();
  for (;oi!=oiend;oi++) {
    map<string,int>::const_iterator ni=names.find(*oi);
    assert( ni != names.end() );
    int i=ni->second;
    vector<string>::const_iterator oj=ord.begin(),ojend=ord.end();
    for (;oj!=ojend;++oj) {
      map<string,int>::const_iterator nj=names.find(*oj);
      assert( nj != names.end() );
      int j=nj->second;
      color c(1,1,1);
      if (binverse) c=color(0,0,0);
      map<int,double>::const_iterator mi=meas[i].find(j);
      if (mi != meas[i].end()) c=colormap[255-int(mi->second*255)];
      uint r=uint(c.x*255),g=uint(c.y*255),b=uint(c.z*255);
      o << r << ' ' << g << ' ' << b << ' ';
    }
    o << endl;
  }
}

void print_pgm(ostream& o,const vector< map<int,double> >& meas,
	       const map<string,int>& names,
	       const vector<string>& ord,bool binverse)
{
  const uint SZ=ord.size();
  
  o << "P2" << endl << SZ << ' ' << SZ << endl << "256" << endl;
  vector<string>::const_iterator oi=ord.begin(),oiend=ord.end();
  for (;oi!=oiend;oi++) {
    map<string,int>::const_iterator ni=names.find(*oi);
    assert( ni != names.end() );
    int i=ni->second;
    vector<string>::const_iterator oj=ord.begin(),ojend=ord.end();
    for (;oj!=ojend;++oj) {
      map<string,int>::const_iterator nj=names.find(*oj);
      assert( nj != names.end() );
      int j=nj->second;

      double val=0.0;
      map<int,double>::const_iterator mi=meas[i].find(j);
      if (mi != meas[i].end()) val=mi->second;
      o << int(255*(binverse?val:1.0-val)) << ' ';
    }
    o << endl;
  }
}

struct _dendr {
  int idx,size;
  point tl,br;
  _dendr(int i,int sz,point t,point b)
    :idx(i),size(sz),tl(t),br(b) {}
};

void draw_dendrogram(const resgraph& R,int root,EPSarea& eps,
		     point tl,point br)
{
  stack<_dendr> S;
  S.push(_dendr(root,R.tag(root)._size,tl,br));

  while (!S.empty()) 
    {
      _dendr curr=S.top(); S.pop();
      point topleft=curr.tl,bottomright=curr.br;

      const uint N=curr.size;
      assert(N > 0);
      if (N == 1) {
	string s=R.tag(curr.idx)._name;
	eps.text(point(.7*bottomright.x+.3*topleft.x,
		       bottomright.y+0.15*(topleft.y-bottomright.y)),s);
      }
      else {
	const double height=topleft.y-bottomright.y;
	const double width=bottomright.x-topleft.x;
	
	resgraph::edge_const_iterator ei=R.nbrs_const_iterate(curr.idx);
	const int _r=ei.index(); ++ei;
	const int _l=ei.index(); ++ei; assert( ei.end() );

	const uint rsz=R.tag(_r)._size,lsz=R.tag(_l)._size;
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

	S.push(_dendr(_r,R.tag(_r)._size,
		      topleft,point(rmid.x,topleft.y-rh)));
	S.push(_dendr(_l,R.tag(_l)._size,
		      point(topleft.x,bottomright.y+lh),point(lmid.x,bottomright.y)));
      }
    } 
}


void prepare_image(const vector< map<int,double> >& meas,
		   const map<string,int>& names,
		   const vector<string>& ord,
		   vector< vector<double> >& image)
{
  const uint SZ=ord.size();

  image.resize(SZ);
  for (uint i=0;i<SZ;i++) image[i].resize(SZ);

  vector<string>::const_iterator oi=ord.begin(),oiend=ord.end();
  for (uint _i=0;oi!=oiend;oi++,++_i) {
    map<string,int>::const_iterator ni=names.find(*oi);
    assert( ni != names.end() );
    int i=ni->second;
    vector<string>::const_iterator oj=ord.begin(),ojend=ord.end();
    for (uint _j=0;oj!=ojend;++oj,++_j) {
      map<string,int>::const_iterator nj=names.find(*oj);
      assert( nj != names.end() );
      int j=nj->second;

      double val=0.0;
      map<int,double>::const_iterator mi=meas[i].find(j);
      if (mi != meas[i].end()) val=mi->second;
      image[SZ-1-_i][SZ-1-_j]=1-val;
    }
  }
}

void extract_order(resgraph& R,int root,vector<string>& ord)
{
  stack<int> S;
  S.push(root);
  while (!S.empty()) 
    {
      int curr=S.top(); S.pop();      
      assert(R.indegree(curr) == 2);
      resgraph::edge_const_iterator ei=R.nbrs_const_iterate(curr);
      for (;!ei.end();++ei) {
	if (R.tag(ei.index())._size==1) 
	  ord.push_back(R.tag(ei.index())._name);
	else {
	  S.push(ei.index());
	}
      }
    }
}

bool is_tree_name(string str,char c)
{
  if (str.empty()) return false;
  if (str[0]!='<') return false;
  if (str[str.size()-1]!='>') return false;
  if (str[1]!=c)   return false;

  istringstream iss(str.substr(2,str.size()-1));
  int idx;
  iss >> idx; // the rest must be an integer
  if (!iss) return false;

  return true;
}

void clean_tree(resgraph& R,resgraph& clean,
		const vector<int>& roots)
{
  // label the "level"
  vector<int> level(R.num_nodes(),1);
  bool onechange=true;
  while (onechange) {
    onechange=false;
    for (uint k=0;k<R.num_nodes();k++) {
      if (R.indegree(k) > 0) {
	resgraph::edge_const_iterator ei=R.nbrs_const_iterate(k);
	uint right=ei.index(); ++ei;
	uint left=ei.index();  ++ei; assert(ei.end());
	  
	int old=level[k];
	if (level[right]==level[left]) level[k]=level[right]+1;
	else {
	  level[k]=max(level[right],level[left]);
	}
	if (old != level[k]) onechange=true;
      }
    }
  }

  // collapse whenever the same level
  stack<int> Q;
  for (uint k=0;k<roots.size();k++) Q.push(roots[k]);
  while (!Q.empty()) {
    int curr=Q.top(); Q.pop();
    if (R.indegree(curr) > 0) {
      resgraph::edge_const_iterator ei=R.nbrs_const_iterate(curr);	
      int to_collapse=-1;
      for (;!ei.end();++ei) {
	if (level[curr]==level[ei.index()]) to_collapse=ei.index();
      }
      if (to_collapse < 0) { // none found
	ei=R.nbrs_const_iterate(curr);
	while (!ei.end()) Q.push(ei.index()),++ei;
	continue;
      }
      
      ei=R.nbrs_const_iterate(to_collapse);
      vector<int> nbrscoll(2);
      for (uint k=0;!ei.end();++ei,++k) {
	assert(k < 2);
	nbrscoll[k]=ei.index();
      }
      R.remove_edge(to_collapse,nbrscoll[0]);
      R.remove_edge(to_collapse,nbrscoll[1]);
      R.add_edge(curr,nbrscoll[0]);
      R.add_edge(curr,nbrscoll[1]);
      R.remove_edge(curr,to_collapse);
      Q.push(curr);
    }
  }  

  // mask of degree > 0
  uint total=0;
  vector<bool> mask(R.num_nodes(),false);
  for (uint k=0;k<R.num_nodes();k++) {
    if (R.degree(k) > 0) {
      mask[k]=true;
      total++;
    }
  }

  clean.resize_and_clear(total);
  R.subgraph(mask,clean);
}

bool is_root_name(string str) { return is_tree_name(str,'r'); }
bool is_branch_name(string str) { return is_tree_name(str,'b'); }

int main(int argc,char** argv)
{
  double WF=0.3,linewidth=0.01,textwidth=10.0;
  double eps_height=20.0;
  bool binverse=false;
  bool boceancolormap=false;
  string measfile="<none>";

  vector<param*> prms;
  prms.push_back(make_param('m',"measure_file",measfile));
  prms.push_back(make_param('I',"inverse_colors",binverse));
  prms.push_back(make_param('d',"dendro_width_factor",WF));
  prms.push_back(make_param('l',"linewidth",linewidth));
  prms.push_back(make_param('t',"textwidth",textwidth));
  prms.push_back(make_param('e',"eps_height_in_cm",eps_height));
  prms.push_back(make_param('o',"oceancolormap",boceancolormap));

  vector<string> args;
  string usage =
    "Usage: Gtreeanlz <command> <infile> <outfile>\n"
    "Copyright (c) 2007, Pau Fernandez\n\n"
    "command:\n"
    "    order      Print the order\n"
    "    tree       Print the clean tree\n"
    "    image      Generate an image (.pnm)\n"
    "    dendro     Print the dendrogram\n"
    "    diagram    Print the whole diagram";

  parse_params_ex(prms,argc,argv,usage,"Gtreeanlz",args,3);
  string infile=args[1],outfile=args[2],command=args[0];

  // Read the tree
  resgraph R;
  if (!read_graph_from_file(R,infile)) {
    cerr << "Couldn't read graph " << infile << endl;
    return -1;
  }
  
  assert(!R.is_undirected()); // simple check
  vector<int> roots;
  int totalsz=0;
  for (uint i=0;i<R.num_nodes();i++) {
    if (is_root_name(R.tag(i)._name)) {
      roots.push_back(i);
      totalsz+=R.tag(i)._size;
    }
  }
  
  // read measure
  vector< map<int,double> > meas;
  map<string,int> names;
  if (measfile!="<none>") {
    ifstream fin(measfile.c_str());
    read_measure(fin,meas,names);
  }
  
  // order
  vector<string> order;
  int max_size=0;
  for (uint k=0;k<roots.size();k++) {
    vector<string> ord;
    extract_order(R,roots[k],ord);
    copy(ord.begin(),ord.end(),back_inserter(order));
    max_size=max(max_size,R.tag(roots[k])._size);
  }

  ofstream fout(outfile.c_str());
  if (!fout.is_open()) {
    cerr << "Coudln't open '" << outfile << "' for writing." << endl;
    return -1;
  }

  // Carry out commands
  if (command=="image") {
    if (measfile=="<none>") {
      cerr << "To print an image a measure file is needed." << endl;
      return -1;
    }
    if (boceancolormap) {
      vector<color> cm;
      ocean_colormap(cm,binverse);
      print_pnm(fout,meas,names,order,cm,binverse);
    }
    else {
      print_pgm(fout,meas,names,order,binverse);
    }
    return 0;
  }

  if (command=="tree") {
    resgraph C;
    clean_tree(R,C,roots);
    write_graph_to_file(C,outfile);
    return 0;
  }
  
  if (command=="order") {
    vector<string>::const_iterator oi=order.begin(),oiend=order.end();
    for (;oi!=oiend;++oi) {
      fout << *oi << endl;
    }
    return 0;
  }

  if (command=="dendro") {
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
      int sz=R.tag(roots[k])._size;
      double frac=double(sz)/double(totalsz)*0.99;
      draw_dendrogram(R,roots[k],eps,
		      point(0.99*width_factor,curr),
		      point(0.01*width_factor+(1-frac)*0.98*width_factor,curr+frac));
      curr+=frac;
    }
  }

  if (command=="diagram") {
    if (measfile=="<none>") {
      cerr << "To print a diagram a measure file is needed." << endl;
      return -1;
    }
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
      uint sz=R.tag(roots[k])._size;
      double frac=double(sz)/double(totalsz);
      map<uint,string> nonames;
      draw_dendrogram(R,roots[k],eps,
		      point(1.0,curr),
		      point(dendro_width_factor*frac+1.0,curr+frac));
      curr+=frac;
    }

    vector< vector<double> > _image;
    prepare_image(meas,names,order,_image);
    
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
    for (uint k=0;k<R.num_nodes();k++) 
      fout << R.tag(k)._size << endl;
  }
}
