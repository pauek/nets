// 
//  Copyright (c) 2007, Pau Fernández
//

#include <cmath>
#include <fstream>
#include <list>

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

pair<int,int> max_edge_btwns(bgraph& B)
{
  const uint N=B.num_nodes();
  pair<int,int> idxs(-1,-1);
  double mx=-1.0;
  for (uint i=0;i<N;i++) {
    bgraph::edge_const_iterator ei=B.nbrs_const_iterate(i);
    for (;!ei.end();++ei) {
      if (ei.tag() > mx) {
	idxs=make_pair(i,ei.index());
	mx=ei.tag();
      }
    }
  }
  return idxs;
}

template<class A,class B>
ostream& operator<<(ostream& o,pair<A,B> p)
{ return o << p.first << ' ' << p.second; }

inline string new_name(string pre,int idx,string suff)
{
  ostringstream oss;
  oss << pre << idx << suff;
  return oss.str();
}

inline string new_root_name(int idx) { return new_name("<r",idx,">"); }
inline string new_branch_name(int idx) { return new_name("<b",idx,">"); }

void communities(graph& G,resgraph& R,bool boutputtimes=false)
{
  assert(G.is_undirected());
  const uint N=G.num_nodes(),E=G.num_edges();

  // Prepare initial structure
  vector<int> comp(N);
  connected_components(G,comp);
  int mx=0;
  for (uint i=0;i<N;i++) mx=max(mx,comp[i]);
  int numcomp=mx+1;

  int curr_root_idx=0,curr_branch_idx=0;
  
  vector<int> roots(numcomp);
  for (int k=0;k<numcomp;k++) roots[k]=R.add_node();
  vector< set<int> > members(numcomp);
  for (uint i=0;i<N;i++) members[comp[i]].insert(i);

  std::list< pair<set<int>,int> > curr_roots;
  for (int i=0;i<numcomp;i++) {
    R.tag(roots[i])._size=members[i].size();
    R.tag(roots[i])._name=new_root_name(curr_root_idx++);
    curr_roots.push_back(make_pair(members[i],roots[i]));
  }

  // Main loop
  bgraph B;
  betweenness_init_result_network(G,B);
  betweenness_centrality(B);

  progress_bar P(0,B.num_edges(),cerr);
  while (B.num_edges() > 0) 
    {
      // max btwns edge?
      pair<int,int> meb=max_edge_btwns(B);
      B.remove_edge(meb.first,meb.second);      

      // determine the root affected
      list< pair<set<int>,int> >::iterator ri=curr_roots.begin(),riend=curr_roots.end();
      list< pair<set<int>,int> >::iterator iroot=curr_roots.end();
      uint max_size=0;
      for (;ri!=riend;++ri) {
	if (ri->first.find(meb.first) != ri->first.end()) {
	  assert(ri->first.find(meb.second) != ri->first.end());
	  assert(iroot == curr_roots.end());
	  iroot=ri;
	}
	max_size = max((size_t)max_size,ri->first.size());
      }

      // Accumulate number of edges necessary to break module
      tree_tag& t=R.tag(iroot->second);
      t._num_edges++;
      
      // components
      const set<int> rootset=iroot->first;
      uint num_comp=connected_components(B,comp,rootset);

      if (num_comp > 1) {
	set<int> sright,sleft;
	vector<int>::const_iterator ci=comp.begin(),ciend=comp.end();
	for (int k=0;ci!=ciend;++ci,++k) {
	  assert(*ci < 2 && *ci > -2);
	  switch (*ci) {
	  case 0: sright.insert(k); break;
	  case 1: sleft.insert(k);  break;
	  };
	}
	assert(!sleft.empty());

	// subdivide
	if (boutputtimes) 
	  cout << E-B.num_edges() << ' ' 
	       << 1-double(B.num_edges())/double(E) << ' ' 
	       << rootset.size() << ' ' << double(max_size)/double(N) << ' '
	       << sright.size() << ' ' << sleft.size() << endl;

	// right branch
	int right=R.add_node();
	tree_tag& rtag=R.tag(right);
	if (sright.size()==1) {
	  rtag._name=B.tag(*sright.begin()).first;
	  rtag._size=1;
	}
	else {
	  rtag._name=new_branch_name(curr_branch_idx++);
	  rtag._size=sright.size();
	}
	
	// left branch
	int left=R.add_node();
	tree_tag& ltag=R.tag(left);
	if (sleft.size()==1) {
	  ltag._name=B.tag(*sleft.begin()).first;
	  ltag._size=1;
	}
	else {
	  ltag._name=new_branch_name(curr_branch_idx++);
	  ltag._size=sleft.size();
	}
	
	R.add_edge(iroot->second,right);
	R.add_edge(iroot->second,left);
	
	curr_roots.erase(iroot);
	curr_roots.push_back(make_pair(sright,right));
	curr_roots.push_back(make_pair(sleft,left));
      }

      betweenness_reset_result_network(B,rootset);
      betweenness_centrality(B,rootset); // refresh only locally

      ++P;
    }
}

int main(int argc,char** argv)
{
  bool boutputtimes=false;

  vector<param*> prms;
  prms.push_back(make_param('t',"outputtimes",boutputtimes));

  vector<string> args;
  string usage =
    "Usage: community <infile> <outfile>\n"
    "Copyright (c) 2007, Pau Fernandez\n\n"
    "   Find the community structure of a network (Girvan & Newman)";
  parse_params_ex(prms,argc,argv,usage,"community",args,2);
  string infile=args[0],outfile=args[1];

  graph G;
  resgraph R;

  if (!read_graph_from_file(G,infile)) {
    cerr << "Couldn't read graph " << infile << endl;
    return -1;
  }
  if (!G.is_undirected()) {
    cerr << "warning: making graph undirected";
    G.to_undirected();
  }
  communities(G,R,boutputtimes);   // modularize
  write_graph_to_file(R,outfile);
}
