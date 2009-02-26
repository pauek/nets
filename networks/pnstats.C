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
  bool bclustering=false,bdegree=false,boutdegree=false;
  bool blocal_eff=false,bindegree=false,bpercolation=false;
  bool ball=false,bprint_index=true;

  string filename("stdout");
  vector<param*> prms;
  string usage = 
    "Usage: Gstats [options] <graph_file(.ladj)>...\n"
    "Copyright (c) 2007, Pau Fernandez\n\n"
    "  Per node statistics";

  prms.push_back(make_param('K',"degree",bdegree));
  prms.push_back(make_param('I',"indegree",bindegree));
  prms.push_back(make_param('O',"outdegree",boutdegree));
  prms.push_back(make_param('C',"clustering",bclustering)); 
  prms.push_back(make_param('l',"local_effic",blocal_eff));
  prms.push_back(make_param('p',"percolation",bpercolation));
  
  prms.push_back(make_param('a',"all",ball));
  prms.push_back(make_param('i',"print_index",bprint_index));
  prms.push_back(make_param('o',"outputfile",filename));

  vector<string> args;
  parse_params_ex(prms,argc,argv,usage,"stats",args,1);

  uint count=0;
  count+=(bdegree     ?1:0);
  count+=(bindegree   ?1:0);
  count+=(boutdegree  ?1:0);
  count+=(bclustering ?1:0);
  count+=(blocal_eff  ?1:0);
  count+=(bpercolation?1:0);

  if (ball) count=6;
  if (count==0 && !ball) {
    print_usage(prms,usage);
    exit(0);
  }

  string graph_file(args[0]);

  typedef adj_list<std::string,std::string> graph;
  graph g;
  if (!read_graph_from_file(g,graph_file)) {
    cerr << "Couldn't read graph " << graph_file << endl;
    return -1;
  }
  graph gu=g;
  gu.to_undirected();
    
  ostream* poutput=&cout;
  if (filename != "stdout") 
    poutput=new ofstream(filename.c_str());
  ostream& o=*poutput;
    
  bool isU=g.is_undirected();
  vector<int> lp;
  if (bpercolation) local_percolation(g,lp);
  for (uint i=0;i<g.num_nodes();i++) {
    if (bprint_index) o << i;

    ostringstream os;
    os << g.tag(i);
    std::string itag=os.str();
    if (itag.size() > 0) o << ' ' << itag;

    uint out=g.outdegree(i),in=g.indegree(i);
    if (bdegree || ball) o << ' ' << (out + in)/(isU?2:1);
    if (bindegree || ball) o << ' ' << in;
    if (boutdegree || ball) o << ' ' << out;
    if (blocal_eff || ball) o << ' ' << local_efficiency(g,i);
    if (bclustering || ball) o << ' ' << clustering_coeff(gu,i);
    if (bpercolation || ball) o << ' ' << lp[i];
    o << endl;
  }
}
