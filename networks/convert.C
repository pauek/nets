// 
//  Copyright (c) 2007, Pau Fernández
//

#include <fstream>
#include <utils/param.H>

#include "adj_list.H"
#include "io.H"

using namespace net;

int main(int argc,char** argv)
{
  vector<param*> prms;
  vector<string> args;
  string usage=
    "Usage: Gconvert <infile> <outfile>\n"
    "Copyright (c) 2007, Pau Fernandez\n\n"
    "  Convert a network to another format\n"
    "  (based on the filename extension)";

  parse_params_ex(prms,argc,argv,usage,"Gconvert",args,2);
  string infile=args[0];
  string outfile=args[1];

  typedef adj_list<std::string,std::string> graph;
  graph G;
  if (!read_graph_from_file(G,infile)) {
    cerr << "Couldn't read file " << infile;
  }
  if (!write_graph_to_file(G,outfile)) {
    cerr << "Couldn't write file " << outfile;
  }
}
  
