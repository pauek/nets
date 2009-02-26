// 
//  Copyright (c) 2007, Pau Fernández
//

#include <fstream>

#include <utils/param.H>
#include <utils/stl.H>

#include "adj_list.H"
#include "stats.H"
#include "io.H"

using namespace net;

int main(int argc,char** argv)
{
  bool bdegree_dist=false;
  bool bindegree_dist=false;
  bool boutdegree_dist=false;
  bool bcum_degree_dist=false;
  bool bclustering_dist=false;
  bool bcorrelations=false;

  string filename("stdout");
  vector<param*> prms;
  string usage = 
    "Usage: Gdstats [option(just one)] <graph_file(.ladj)>\n"
    "Copyright (c) 2007, Pau Fernandez\n\n"
    "  Degree statistics of a network";

  prms.push_back(make_param('d',"degree_dist",bdegree_dist));
  prms.push_back(make_param('i',"indegree_dist",bindegree_dist));
  prms.push_back(make_param('o',"outdegree_dist",boutdegree_dist));
  prms.push_back(make_param('D',"cum_degree_dist",bcum_degree_dist));
  prms.push_back(make_param('c',"clustering_dist",bclustering_dist));
  prms.push_back(make_param('r',"correlations",bcorrelations));

  vector<string> args;
  parse_params_ex(prms,argc,argv,usage,"ndstats",args,1);

  uint count=0;
  count+=(bdegree_dist?1:0);
  count+=(bindegree_dist?1:0);
  count+=(boutdegree_dist?1:0);
  count+=(bcum_degree_dist?1:0);
  count+=(bclustering_dist?1:0);
  count+=(bcorrelations?1:0);

  if (count!=1 || args.size() > 1) {
    print_usage(prms,usage);
    exit(0);
  }

  string graph_file(args[0]);

  typedef adj_list<std::string,std::string> graph;
  graph g;
  read_graph_from_file(g,graph_file);

  ostream* poutput=&cout;
  if (filename != "stdout") 
    poutput=new ofstream(filename.c_str());
  ostream& o=*poutput;

  if (bdegree_dist || bindegree_dist || boutdegree_dist) {
    vector<double> dd;
    if (bdegree_dist) 
      degree_distribution(g,dd);
    else if (bindegree_dist) 
      indegree_distribution(g,dd);
    else if (boutdegree_dist) 
      outdegree_distribution(g,dd);
    
    for (uint i=1;i<dd.size();i++) 
      if (dd[i] != 0) o << i << ' ' << dd[i] << endl;
  }

  if (bcum_degree_dist) {
    map<uint,double> cdd;
    cum_degree_distribution(g,cdd);
    map<uint,double>::const_iterator ci=cdd.begin(),ciend=cdd.end();
    for (;ci!=ciend;ci++) 
      if (ci->first!=0 && ci->second!=0.0) o << ci->first << ' ' << ci->second << endl;
  }

  if (bcorrelations) {
    vector<double> nd;
    graph gcopy=g;
    gcopy.to_undirected();
    avg_nbrs_degree(gcopy,nd);
    for (uint i=1;i<nd.size();i++) 
      if (nd[i]!=0)
	o << i << ' ' << nd[i] << endl;
  }

  if (bclustering_dist) {
    vector<double> cd;
    graph gcopy=g;
    gcopy.to_undirected();
    clustering_dist(gcopy,cd);
    for (uint i=1;i<cd.size();i++) 
      if (cd[i]!=0) o << i << ' ' << cd[i] << endl;
  }
}
