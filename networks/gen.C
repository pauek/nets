// 
//  Copyright (c) 2007, Pau Fernández
//

#include <fstream>
#include <unistd.h>
#include <utils/param.H>
#include <utils/stl.H>

#include "adj_list.H"
#include "stats.H"
#include "io.H"
#include "gen.H"

using namespace net;


int main(int argc,char** argv)
{
  uint N=1000,E=2000;
  double p1=2.0,p2=5.0;
  long seed=-1;

  vector<param*> prms;
  string usage =
    "Usage: Ggen [options] <type> <outfile>\n"
    "Copyright (c) 2007, Pau Fernandez\n\n"
    "   Generate a random graph"
    "\n   <type>: \"er\"    - Erdos-Renyi"
    "\n           \"ba\"    - Barabasi-Albert (Preferential Attachment)"
    "\n           \"sf\"    - Scale-Free with cutoff (exp - p, cutoff - q)"
    "\n           \"sfin\"  - Directed Scale-Free with cutoff (only indegrees) (exp - p, cutoff - q)";

  vector<string> args;
  prms.push_back(make_param('N',"num_nodes",N));
  prms.push_back(make_param('E',"num_edges",E));
  prms.push_back(make_param('p',"param1",p1));
  prms.push_back(make_param('q',"param2",p2));
  prms.push_back(make_param('s',"seed",seed));
  parse_params_ex(prms,argc,argv,usage,"Ggen",args,2);

  const string command=args[0],outfile=args[1];

  Uniform<double> rng;
  rng.seed(seed < 0 ? long(time(0)+getpid()) : seed);

  typedef adj_list<std::string,std::string> graph;
  graph G(N);

  map<string,uint> menu;
  menu["er"]=0;
  menu["ba"]=1;
  menu["sf"]=2;
  menu["sfin"]=3;

  map<string,uint>::const_iterator comm=menu.find(command);
  if (comm==menu.end()) {
    cerr << "Didn't understand command " << command << endl;
    return -1;
  }

  switch (comm->second) {
  case 0:{ // Erdos-Renyi
    erdos_renyi(G,E);
    break;
  }
  case 1:{ // pref. attach.
    G.to_undirected();
    barabasi_albert(G);
    break;
  }
  case 2:{
    G.to_undirected();
    scale_free_with_cutoff(G,p1,p2);
    break;
  }
  case 3:{
    scale_free_indegree(G,p1,p2);
    break;
  }
  };

  if (!write_graph_to_file(G,outfile)) {
    cerr << "Couldn't write graph " << outfile << endl;
    return -1;
  }
}
