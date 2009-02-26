// 
//  Copyright (c) 2007, Pau Fernández
//

#include <string>
#include <map>
#include <sstream>
#include <fstream>
// #include "adjlist.H"

using namespace std;

namespace net 
{
  class wrong_format {
    string _msg;
  public:
    wrong_format(string msg):_msg(msg) {}
    string msg() const { return _msg; }
  };

  inline void eat_white(istream& i)
  {
    char c=i.get();
    while (c==' ' || c=='\t') c=i.get();
    i.unget();
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Edge list format:
  // - edges are specified by their tags.
  // - as an optional third column, the tags for the edges can be specified
  // 
  // node_u0_tag node_v0_tag edge0_tag
  // node_u1_tag node_v1_tag edge
  // ...
  // node_uE_tag node_vE_tag

  template<class Network>
  void write_edgl(Network& g,ostream& o)
  {
    const uint SZ=g.num_nodes();
    for (uint k=0;k<SZ;k++) {
      typename Network::edge_const_iterator ei=g.nbrs_const_iterate(k);
      for (;!ei.end();++ei) {
	o << g.tag(k) << ' ' << g.tag(ei.index());
	ostringstream sout;
	sout << ei.tag();
	string etag=sout.str();
	if (etag.size()>0) o << ' ' << etag;
	o << endl;
      }
    }
  }

  template<class Network>
  void read_edgl(Network& g,istream& i)
  {
    // We assume the node tag to be able to read itself from an istream
    map<string,uint> tag2idx;
    g.clear();
    while (i) {
      char line[100000]; // to specify two nodes and the tag for the edge...
      i.getline(line,100000);
      if (strlen(line)==0) continue;

      istringstream si(line);
      string str1,str2,str3;

      uint u,v;
      si >> str1;
      if (!si) {
	cerr << str1 << endl;
	throw wrong_format("Tag for the first vertex expected");
      }
      map<std::string,uint>::iterator ui=tag2idx.find(str1);
      if (ui==tag2idx.end()) {
	istringstream is(str1);
	typename Network::node_tag nt;
	is >> nt;
	u=g.add_node(nt);
	tag2idx[str1]=u;
      }
      else u=ui->second;

      si >> str2;
      if (!si) throw wrong_format("Tag for the second vertex expected");
      map<std::string,uint>::iterator vi=tag2idx.find(str2);
      if (vi==tag2idx.end()) {
	istringstream is(str2);
	typename Network::node_tag nt;
	is >> nt;
	v=g.add_node(nt);
	tag2idx[str2]=v;
      }
      else v=vi->second;

      if (si) {
	typename Network::edge_tag et;
	si >> et;
	g.add_edge(u,v,et);
      }
      else {
	g.add_edge(u,v);
      }
    }

    if (g.check_undirected()) g.to_undirected();
  }

  ///////////////////////////////////////////////////////////////////////////////////
  // Adjacency list format:
  // - Optional tags. If '[]' parenthesis not found, the default tag is constructed.
  // - Mandatory consecutive nodes and all appearing
  // -------------------------------------------------------
  // numnodes
  // node0[tag]  nbr0[tag] nbr1[tag] nbr2[tag] ... nbrK[tag]
  // node1[tag]  ...
  // ...
  // nodeN[tag]  ...
  // -------------------------------------------------------

  template<class Network>
  void write_ladj(Network& g,ostream& o)
  {
    const uint SZ=g.num_nodes();
    o << SZ << endl;
    for (uint k=0;k<SZ;k++) {
      o << k;
      ostringstream nout;
      nout << g.tag(k);
      string ntag=nout.str();
      if (ntag.size() > 0) o << '[' << ntag << ']';
      typename Network::edge_const_iterator ei=g.nbrs_const_iterate(k);
      for (;!ei.end();++ei) {
	o << ' ' << ei.index();
	ostringstream eout;
	eout << ei.tag();
	string etag=eout.str();
	if (etag.size() > 0) o << '[' << etag << ']';
      }
      o << endl;
    }
  }

  template<class Network>
  void read_ladj(Network& g,istream& i)
  {
    // We assume the node tag to be able to read itself from a istream
    uint N;
    i >> N;
    if (!i) throw wrong_format("Expecting graph size at the beginning");
    g.resize_and_clear(N);

    int last=-1;
    while (i) {
      // vertex index
      uint u;
      i >> u;
      if (!i) break;
      if (last == -1) {
	if (u != 0) throw wrong_format("Adjacency list doesn't start with 0");
      }
      else {
	if (int(u) != last+1) throw wrong_format("Vertices must be consecutive");
	if (u >= N) throw wrong_format("Indices cannot exceed graph size");
      }
      last=int(u);

      // vertex tag (if any)
      typename Network::node_tag nt;
      char c=i.get();
      i.unget();
      if (c=='[') {
	c=i.get(); // the '['
	string in;
	c=i.get();
	while (c!=']') { 
	  in+=c; 
	  c=i.get(); 
	  if (!i) throw wrong_format("Parenthesis '[' not closed");
	}
      
	istringstream sin(in);
	sin >> nt;
      }
      g.set_node_tag(u,nt);

      eat_white(i);

      // edges
      eat_white(i);
      bool eol=(i.get()=='\n' || i.eof());
      i.unget();
      while (!eol) 
	{
	  uint v;
	  i >> v;
	  if (!i) {
	    string str;
	    i >> str;
	    throw wrong_format(string("Expected an integer: got \"")+str+string("\""));
	  }
	  else {
	    if (v >= N) throw wrong_format("Indices cannot exceed graph size");
	  }

	  // edge tag (if any)
	  typename Network::edge_tag et;
	  char c=i.get();
	  i.unget();
	  if (c=='[') {
	    c=i.get(); // the '['
	    string in;
	    c=i.get();
	    while (c!=']') { 
	      in+=c; 
	      c=i.get();
	      if (!i) throw wrong_format("Parenthesis '[' not closed");
	    }
	  
	    istringstream sin(in);
	    sin >> et;
	  }
	  g.add_edge(u,v,et);

	  // eol?
	  eat_white(i);
	  eol=(i.get()=='\n' || i.eof());
	  i.unget();
	}
    }

    if (g.check_undirected()) {
      g.to_undirected();
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////
  // Dot (neato) format

  template<class Network>
  void write_dot(Network& g,ostream& o)
  {
    const uint SZ=g.num_nodes();
    string name="digraph",link="->";
    if (g.is_undirected()) name="graph",link="--";

    o << name << " G {" << endl;
    for (uint k=0;k<SZ;k++) {
      bool any_to_write=false;
      typename Network::edge_const_iterator ei=g.nbrs_const_iterate(k);
      for (;!ei.end();++ei) { 
	if ((g.is_undirected() && (ei.index()>= k)) || 
	    !g.is_undirected() ) any_to_write=true; 
      }
      ei=g.nbrs_const_iterate(k);
      if (any_to_write) 
	{
	  // the links
	  o << "  { ";
	  for (;!ei.end();++ei) 
	    {
	      if ((g.is_undirected() && (ei.index() >= k)) || 
		  !g.is_undirected()) {
		ostringstream nout;
		nout << g.tag(ei.index());
		string ntag=nout.str();
		if (ntag.size() > 0) {
		  o << '\"' << ntag << "\" ";
		}
		else {
		  o << ei.index() << ' ';
		}
	      }
	    }
	  o << "} " << link << ' ';

	  // the node
	  ostringstream nout;
	  nout << g.tag(k);
	  string ntag=nout.str();
	  if (ntag.size() > 0) o << '\"' << ntag << '\"';
	  else o << k;
	  o << ';' << endl;
	}
    }
    o << "}" << endl;
  }

  ///////////////////////////////////////////////////////////////////////////////////
  // Pajek format

  template<class Network>
  void write_net(const Network& g,ostream& o)
  {
    const uint N=g.num_nodes();
    // nodes
    o << "*Vertices\t" << g.num_nodes() << endl;
    for (uint k=0;k<N;k++) {
      o << "  " << k+1 << " \"" << g.tag(k) << "\"" << endl;
    }
    
    // edges
    o << "*Arcs" << endl;
    for (uint k=0;k<N;k++) {
      typename Network::edge_const_iterator ei=g.nbrs_const_iterate(k);
      for (;!ei.end();++ei) {
	o << "  " << k+1 << " " << ei.index()+1 << " 1" << endl;
      }
    }    
  }

  /*
  template<class Network>
  read_net(const Network& g,istream& i)
  {
    g.clear();

    // nodes
    string str;
    i >> str;
    if (str!="*Vertices") {
      throw wrong_format("Section '*Vertices' not found");
    }

    uint n;
    i >> n;
    if (!i) {
      throw wrong_format("Number of Vertices not present");
    }

    g.resize(n);
    for (uint k=0;k<N;k++) 
      {
	uint u;
	i >> u; 
	if (!i) throw wrong_format("Integer expected");
	
	eat_white(i);
	
	char c;
	i >> c;
	if (c!='\"') throw wrong_format("Quote '\"' expected");
	
	string in;
	while (c!='\"') { 
	  in+=c; 
	  c=i.get();
	  if (!i) throw wrong_format("Quote '\"' not closed");
	}

	NodeTag nt
	istringstream sin(in);
	sin >> nt;

	g.set_node_tag(u,nt);
      }
    
    // edges
    // ...
  }
  */

  ///////////////////////////////////////////////////////////////////////////////////
  // Adjacency matrix format

  template<class Network>
  void write_madj(const Network& g,ostream& o)
  {
    const uint SZ=g.num_nodes();
    o << SZ << endl;
    for (uint i=0;i<SZ;i++) {
      for (uint j=0;j<SZ;j++) {
	if (g.connected(i,j)) o << "1 ";
	else cout << "0 ";
      }
      o << endl;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////
  // Chaco format

  ///////////////////////////////////////////////////////////////////////////////////
  // Tulip format

  ///////////////////////////////////////////////////////////////////////////////////
  // General read

  template<class Network>
  bool try_all_readers(Network& g,istream& i)
  {
    typedef void (*reader_func)(Network&,istream&);
    reader_func reader[]={read_edgl,read_ladj}; //,read_net};
    int k = 0;
    bool failed=true;
    while (failed) {
      failed=false;
      try {
	(*reader[k])(g,i);
      }
      catch (wrong_format e) {
	failed=true;
      }
      if (!failed) return true;
      ++k;
    }
    return false;
  }

  inline string file_extension(string filename)
  {
    uint idx=filename.size();
    while (filename[--idx]!='.');
    string ext;
    while (++idx < filename.size()) ext+=filename[idx];
    return ext;
  }

  template<class Network>
  bool read_graph_from_file(Network& g,string filename) 
  {
    // extension?
    string ext=file_extension(filename);

    ifstream fin(filename.c_str());
    if (!fin.is_open()) {
      cerr << "read_graph_from_file: File \"" << filename << "\" doesn't exist" << endl;
      return false;
    }
    if (ext=="edgl") {
      try { read_edgl(g,fin); }
      catch (wrong_format& e) {
	cerr << "read_graph_from_file: [edge_list] " << e.msg() << endl;
	return false;
      }
      return true;
    }
    if (ext=="ladj") {
      try { read_ladj(g,fin); }
      catch (wrong_format& e) {
	cerr << "read_graph_from_file: [adj_list] " << e.msg() << endl;
	return false;
      }
      return true;
    }

    return false;
  }

  template<class Network>
  bool write_graph_to_file(Network& g,string filename)
  {
    // extension?
    string ext=file_extension(filename);

    ofstream fout(filename.c_str());
    if (!fout.is_open()) {
      cerr << "write_graph_to_file: Cannot open file \"" << filename << "\"" << endl;
      return false;
    }

    if (ext=="edgl") {
      write_edgl(g,fout);
      return true;
    }
    if (ext=="ladj") {
      write_ladj(g,fout);
      return true;
    }
    if (ext=="madj") {
      write_madj(g,fout);
      return true;
    }
    if (ext=="dot") {
      write_dot(g,fout);
      return true;
    }
    if (ext=="net") {
      write_net(g,fout);
      return true;
    }

    cerr << "write_graph_to_file: Unknown format with extension '." << ext << "'" << endl;
    return false;
  }
}
