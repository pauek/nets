// 
//  Copyright (c) 2007, Pau Fernández
//

#include <fstream>

#include <utils/param.H>

#include "adj_list.H"
#include "stats.H"
#include "io.H"

using namespace net;

int main(int argc,char** argv)
{
  bool bnormalize=false,bnode=false,bsort=false,breverse=false;

  vector<param*> prms;
  vector<string> args;
  string usage = 
    "Usage: Gbtwns <infile> <outfile>\n"
    "Copyright (c) 2007, Pau Fernandez";

  prms.push_back(make_param('m',"normalize",bnormalize));
  prms.push_back(make_param('n',"node_betweenness",bnode));
  prms.push_back(make_param('s',"sort",bsort));
  prms.push_back(make_param('r',"reverse",breverse));
  parse_params_ex(prms,argc,argv,usage,"Gbtwns",args,2);
  string infile=args[0];
  string outfile=args[1];

  typedef adj_list<std::string,std::string> graph;
  graph G;
  if (!read_graph_from_file(G,infile)) {
    cerr << "Couldn't read file " << infile;
  }

  typedef adj_list<pair<string,double>,double> resgraph;
  resgraph resG;
  betweenness_init_result_network(G,resG);
  vector<bool> mask(resG.num_nodes(),true);
  betweenness_centrality(resG,mask,bnormalize);

  ofstream out(outfile.c_str());
  const uint N=resG.num_nodes();
  map<double,string> sorted;
  if (bnode) {
    for (uint k=0;k<N;k++) {
      if (!bsort) out << resG.tag(k).first << ' ' << resG.tag(k).second << endl;
      else {
	sorted.insert(make_pair(resG.tag(k).second,resG.tag(k).first));
      }
    }
  }
  else {
    for (uint k=0;k<N;k++) {
      resgraph::edge_const_iterator ei=resG.nbrs_const_iterate(k);
      for (;!ei.end();++ei) {
	if (!bsort) {
	  out << resG.tag(ei.index()).first << ' ' 
	      << resG.tag(k).first <<  ' ' 
	      << ei.tag() << endl;
	}
	else {
	  string s1=resG.tag(ei.index()).first;
	  string s2=resG.tag(k).first;
	  sorted.insert(make_pair(ei.tag(),s1+" "+s2));
	}
      }
    }
  }
  if (bsort) {
    if (breverse) {
      map<double,string>::const_reverse_iterator si=sorted.rbegin(),siend=sorted.rend();
      for (;si!=siend;++si) {
	out << si->second << ' ' << si->first << endl;
      }
    }
    else {
      map<double,string>::const_iterator si=sorted.begin(),siend=sorted.end();
      for (;si!=siend;++si) {
	out << si->second << ' ' << si->first << endl;
      }
    }
  }
}
  
