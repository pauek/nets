// 
//  Copyright (c) 2007, Pau Fernández
//

#include <fstream>

#include <utils/param.H>
#include <utils/stl.H>

#include "adj_list.H"
#include "io.H"
#include "stats.H"

using namespace net;

int main(int argc,char** argv)
{
  bool bclustering=false,bavgmindist=false;
  bool bglobal_eff=false,bassortmix=false;
  bool blocal_eff=false,bundirected=false;
  bool bnumvertices=false,bnumedges=false,bavgconn=false;
  bool bperc_thres=false;
  bool bcomponents=false;
  bool bhorizontal=false;

  string filename("stdout");
  vector<param*> prms;
  string usage = 
    "Usage: Gstats [option]... <graph_file(.ladj)>...\n"
    "Copyright (c) 2007, Pau Fernandez\n\n"
    "  Compute network statistics";

  prms.push_back(make_param('U',"undirected",bundirected));
  prms.push_back(make_param('N',"numvertices",bnumvertices));
  prms.push_back(make_param('E',"numedges",bnumedges));
  prms.push_back(make_param('K',"avgconnectivity",bavgconn));
  prms.push_back(make_param('L',"avgmindist",bavgmindist));
  prms.push_back(make_param('C',"clustering",bclustering)); 
  prms.push_back(make_param('m',"components",bcomponents)); 
  prms.push_back(make_param('x',"assortative_mixing",bassortmix));
  prms.push_back(make_param('g',"global_effic",bglobal_eff));
  prms.push_back(make_param('l',"local_effic",blocal_eff));
  prms.push_back(make_param('p',"perc_thres",bperc_thres));
  prms.push_back(make_param('o',"outputfile",filename));
  prms.push_back(make_param('H',"horizontal",bhorizontal));

  vector<string> args;
  parse_params_ex(prms,argc,argv,usage,"stats",args,1);

  uint count=0;
  count+=(bundirected  ?1:0);
  count+=(bnumvertices ?1:0);
  count+=(bnumedges    ?1:0);
  count+=(bavgconn     ?1:0);
  count+=(bavgmindist  ?1:0);
  count+=(bclustering  ?1:0);
  count+=(bcomponents  ?1:0);
  count+=(bassortmix   ?1:0);
  count+=(bglobal_eff  ?1:0);
  count+=(blocal_eff   ?1:0);
  count+=(bperc_thres  ?1:0);

  bool ball=false;
  if (count==0) {
    count=11; // none selected -> all
    ball=true;
  }

  for (uint i=0;i<args.size();i++) {
    string graph_file(args[i]);

    ostream* poutput=&cout;
    if (filename != "stdout") 
      poutput=new ofstream(filename.c_str());
    ostream& o=*poutput;

    string prefix="";
    if (args.size() > 1) {
      o << args[i] << ": ";
      if (!bhorizontal) o << endl;
      prefix="   ";
    }

    typedef adj_list<std::string,std::string> graph;
    graph g;
    if (!read_graph_from_file(g,graph_file)) {
      cerr << "Couldn't read graph " << graph_file << endl;
      return -1;
    }
    
    bool isU=g.is_undirected();
    if (bundirected || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Undirected:       " << '\t';
      o << (isU?"true":"false");
      if (bhorizontal) o << ' '; else o << endl;
    }

    if (bnumvertices || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Number of vertices:" << '\t';
      o << g.num_nodes();
      if (bhorizontal) o << ' '; else o << endl;
    }

    uint num_edges=g.num_edges();
    if (bnumedges || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Number of edges:   " << '\t';
      o << num_edges ;
      if (bhorizontal) o << ' '; else o << endl;
    }

    if (bavgconn || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Avg. Connectivity: " << '\t';
      o << 2.0*double(num_edges)/double(g.num_nodes());
      if (bhorizontal) o << ' '; else o << endl;
    }
    
    if (bavgmindist || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Avg. Min. Distance:" << '\t';
      double avgmindist=avg_minimum_distance(g);
      if (avgmindist < 0.0) o << "inf";
      else o << avgmindist;
      if (bhorizontal) o << ' '; else o << endl;
    }
    
    if (bclustering || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Avg. Clust. Coeff.:" << '\t';
      graph gcopy=g;
      gcopy.to_undirected();
      o << avg_clustering_coeff(gcopy);
      if (bhorizontal) o << ' '; else o << endl;
    }

    if (bcomponents || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Components:        " << '\t';
      vector< pair<uint,uint> > comp;
      graph gcopy=g;
      gcopy.to_undirected();
      num_component_sizes(gcopy,comp);
      if (!comp.empty()) {
	o << comp[0].second << "(x" << comp[0].first << ')';
	for (uint k=1;k<comp.size();k++) {
	  o << " " << comp[k].second << "(x" << comp[k].first << ')';
	}
	if (bhorizontal) o << ' '; else o << endl;
      }
    }

    if (bassortmix || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Assortative Mixing:" << '\t';
      o << assortative_mixing_coeff(g); 
      if (bhorizontal) o << ' '; else o << endl;
    }
    
    if (bglobal_eff || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Global Efficiency:" << '\t';
      o << global_efficiency(g);
      if (bhorizontal) o << ' '; else o << endl;
    }
    
    if (blocal_eff || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Local Efficiency:" << '\t';
      o << avg_local_efficiency(g);
      if (bhorizontal) o << ' '; else o << endl;
    }
    
    if (bperc_thres || ball) {
      if (count > 1 && !bhorizontal) o << prefix << "Perc. threshold: " << '\t';
      o << percolation_threshold(g);
      if (bhorizontal) o << ' '; else o << endl;
    }

    if (bhorizontal) o << endl;
  }
}
