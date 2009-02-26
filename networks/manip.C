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
  string parameter="";

  vector<param*> prms;
  string usage = 
    "Usage: Gmanip [options] <command> <infile> <outfile>\n"
    "Copyright (c) 2007, Pau Fernandez\n\n"
    "   Manipulate a network\n"
    "\n   <command>: \"toundir\"  - To undirected graph"
    "\n              \"comp\"     - Extract component number <parameter>"
    "\n              \"subgraph\" - Extract subgraph with nodes <parameter>"
    "\n              \"remove\"   - Remove nodes <parameter>"
    "\n              \"rand\"     - Randomize keeping distribution (<parameter> = iterations)" ;

  vector<string> args;
  prms.push_back(make_param('p',"parameter",parameter));
  parse_params_ex(prms,argc,argv,usage,"Gcomp",args,3);

  const string command=args[0],infile=args[1],outfile=args[2];

  typedef adj_list<std::string,std::string> graph;
  graph G;
  if (!read_graph_from_file(G,infile)) {
    cerr << "Couldn't read graph " << infile << endl;
    return -1;
  }

  map<string,uint> menu;
  menu["t"]=0; menu["toundir"]=0;
  menu["c"]=1; menu["comp"]=1;
  menu["s"]=2; menu["subgraph"]=2;
  menu["r"]=3; menu["remove"]=3;
  menu["n"]=3; menu["rand"]=4;


  map<string,uint>::const_iterator comm=menu.find(command);
  if (comm==menu.end()) {
    cerr << "Didn't understand command " << command << endl;
    return -1;
  }
  switch (comm->second) {
  case 0:{ // to undirected
    G.to_undirected();
    break;
  }
  case 1:{ // components
    uint n;
    if (parameter.size() > 0) {
      istringstream sin(parameter);
      sin >> n;
      if (!sin) { 
	n=1; 
	cerr << "warning: didn't understand <parameter>, assuming first component" 
	     << endl; 
      }
    }
    else {
      n=1;
      cerr << "warning: assuming component number 1" << endl;
    }
    graph subG;
    extract_component(G,subG,n-1);
    G=subG;
    break;
  }
  case 2:{ // subgraph
    vector<bool> mask(G.num_nodes(),false);
    uint count=0;
    istringstream sin(parameter);
    string node;
    sin >> node;
    while (sin) {
      set<uint> idxs;
      G.index_from_tag(node,idxs);
      set<uint>::const_iterator ii=idxs.begin(),iiend=idxs.end();
      for (;ii!=iiend;++ii) mask[*ii]=true;
      ++count;
      sin >> node;
    }
    graph subG;
    subG.resize_and_clear(count);
    G.subgraph(mask,subG);
    G=subG;
    break;
  }
  case 3:{ // remove nodes
    vector<bool> mask(G.num_nodes(),true);
    istringstream sin(parameter);
    uint count=0;
    string node;
    sin >> node;
    while (sin) {
      set<uint> idxs;
      G.index_from_tag(node,idxs);
      set<uint>::const_iterator ii=idxs.begin(),iiend=idxs.end();
      for (;ii!=iiend;++ii) mask[*ii]=false;
      ++count;
      sin >> node;
    }
    graph subG;
    subG.resize_and_clear(G.num_nodes()-count);
    G.subgraph(mask,subG);
    G=subG;
    break;
  }
  case 4:{ // randomize connectivity
    istringstream sin(parameter);
    uint times=G.num_edges()/2;
    sin >> times;
    randomize_connectivity(G,times);
    break;
  }
  };

  if (!write_graph_to_file(G,outfile)) {
    cerr << "Couldn't write graph " << outfile << endl;
    return -1;
  }
}
