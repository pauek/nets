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
  vector<param*> prms;
  vector<string> args;
  string usage = 
    "Usage: Gtopovrlp <infile> <outfile(.mtx)>\n"
    "Copyright (c) 2007, Pau Fernandez";
    
  parse_params_ex(prms,argc,argv,usage,"Gtopovrlp",args,2);
  string infile=args[0];
  string outfile=args[1];

  typedef adj_list<std::string,std::string> graph;
  graph G;
  if (!read_graph_from_file(G,infile)) {
    cerr << "Couldn't read file " << infile;
  }

  G.to_undirected();
  vector< map<int,double> > topovrlp(G.num_nodes());
  topological_overlap(G,topovrlp);

  ofstream out(outfile.c_str());
  vector< map<int,double> >::const_iterator ti=topovrlp.begin(),tiend=topovrlp.end();
  for (uint k=0;ti!=tiend;++ti,++k) {
    map<int,double>::const_iterator iti=ti->begin(),itiend=ti->end();
    for (;iti!=itiend;++iti) {
      {
	ostringstream sout;
	sout << G.tag(k);
	string etag=sout.str();
	if (!etag.empty()) out << etag;
	else out << k;
      }
      out << ' ';
      {
	ostringstream sout;
	sout << G.tag(iti->first);
	string etag=sout.str();
	if (!etag.empty()) out << etag;
	else out << iti->first;
      }
      out << ' ' << iti->second << endl;
    }
  }
}
  
