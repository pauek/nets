// 
//  Copyright (c) 2007, Pau Fern√°ndez
//

#include <fstream>
#include <math.h>
using namespace std;

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

void clean_tree(resgraph& R,resgraph& clean,
		const vector<int>& roots)
{
  // label the "level"
  vector<int> level(R.num_nodes(),1);
  bool onechange=true;
  while (onechange) {
    onechange=false;
    for (uint k=0;k<R.num_nodes();k++) {
      if (R.indegree(k) > 0) {
	resgraph::edge_const_iterator ei=R.nbrs_const_iterate(k);
	uint right=ei.index(); ++ei;
	uint left=ei.index();  ++ei; assert(ei.end());
	  
	int old=level[k];
	if (level[right]==level[left]) level[k]=level[right]+1;
	else {
	  level[k]=max(level[right],level[left]);
	}
	if (old != level[k]) onechange=true;
      }
    }
  }

  // collapse whenever the same level
  stack<int> Q;
  for (uint k=0;k<roots.size();k++) Q.push(roots[k]);
  while (!Q.empty()) {
    int curr=Q.top(); Q.pop();
    if (R.indegree(curr) > 0) {
      resgraph::edge_const_iterator ei=R.nbrs_const_iterate(curr);	
      int to_collapse=-1;
      for (;!ei.end();++ei) {
	if (level[curr]==level[ei.index()]) to_collapse=ei.index();
      }
      if (to_collapse < 0) { // none found
	ei=R.nbrs_const_iterate(curr);
	while (!ei.end()) Q.push(ei.index()),++ei;
	continue;
      }
      
      ei=R.nbrs_const_iterate(to_collapse);
      vector<int> nbrscoll(2);
      for (uint k=0;!ei.end();++ei,++k) {
	assert(k < 2);
	nbrscoll[k]=ei.index();
      }
      R.remove_edge(to_collapse,nbrscoll[0]);
      R.remove_edge(to_collapse,nbrscoll[1]);
      R.add_edge(curr,nbrscoll[0]);
      R.add_edge(curr,nbrscoll[1]);
      R.remove_edge(curr,to_collapse);
      Q.push(curr);
    }
  }  

  // mask of degree > 0
  uint total=0;
  vector<bool> mask(R.num_nodes(),false);
  for (uint k=0;k<R.num_nodes();k++) {
    if (R.degree(k) > 0) {
      mask[k]=true;
      total++;
    }
  }

  clean.resize_and_clear(total);
  R.subgraph(mask,clean);
}

bool is_tree_name(string str,char c)
{
  if (str.empty()) return false;
  if (str[0]!='<') return false;
  if (str[str.size()-1]!='>') return false;
  if (str[1]!=c)   return false;

  istringstream iss(str.substr(2,str.size()-1));
  int idx;
  iss >> idx; // the rest must be an integer
  if (!iss) return false;

  return true;
}

bool is_root_name(string str) { return is_tree_name(str,'r'); }
bool is_branch_name(string str) { return is_tree_name(str,'b'); }

void add_leafs(const resgraph& R,map<string,int>& name_map,int root,set<int>& s)
{
  queue<int> Q;
  Q.push(root);
  assert(s.empty());
//   cerr << "leaf " << R.tag(root)._name << endl;
  while (!Q.empty()) {
    uint curr=Q.front();
    resgraph::edge_const_iterator ei=R.nbrs_const_iterate(curr);
    for (;!ei.end();++ei) {
//       cerr << R.tag(ei.index())._name << ' ';
      if (R.indegree(ei.index()) == 0) {
	map<string,int>::iterator it=name_map.find(R.tag(ei.index())._name);
	assert(it!=name_map.end());
	int idx=it->second;
	s.insert(idx);
      }
      else Q.push(ei.index());
    }
    Q.pop();
  }
//   cerr << endl << endl;
}

inline int which_module(int k,vector< set<int> >& modules)
{
  vector< set<int> >::const_iterator mi=modules.begin(),miend=modules.end();
  for (int m=0;mi!=miend;++mi,++m) {
    if (mi->find(k) != mi->end()) return m;
  }
  return -1;
}

double modularity_measure(const graph& G,vector< set<int> >& modules)
{
  const uint NM=modules.size();
  vector<double> e(NM,0.0),a(NM,0.0);
  long int total=0;
  vector< set<int> >::const_iterator mi=modules.begin(),miend=modules.end();
  for (int Mi=0;mi!=miend;++mi,++Mi) {
//     cerr << '{';print_all(*mi,cerr,","); cerr << '}' << endl;
    set<int>::const_iterator si=mi->begin(),siend=mi->end();
    for (;si!=siend;++si) {
      int i=*si;
      graph::edge_const_iterator ei=G.nbrs_const_iterate(i);
      for (;!ei.end();++ei) {
	int j=ei.index();
	int Mj=which_module(j,modules);
// 	cerr << i << ' ' << j << "  " << mi << ' ' << mj << endl;
	if (Mj != -1) { 
	  // indegree only
	  if (Mi == Mj) e[Mi]+=1.0;
	  a[Mi]+=1.0;
	  total++;
	}
      }
    }
  }
//   cerr << "total " << total << endl;
  for (uint k=0;k<NM;k++) {
//     cerr << e[k] << ' ' << a[k] << endl;
    e[k]/=double(total);
    a[k]/=double(total);
  }

  double acum=0.0;
  for (uint k=0;k<NM;k++) acum+=e[k]-a[k]*a[k];

  // Formulaka del Guimer‡... (statistically significant??)
  // arXiv:cond-mat/0403660
  uint SZ=0;
  for (mi=modules.begin();mi!=miend;++mi) SZ+=mi->size();
  double random_case=(1-2.0/sqrt(double(SZ)))*pow(double(SZ-1)/double(total),2.0/3.0);
  
//   cout << SZ << ' ' << random_case << ' ' << acum << ' ';

  return acum/random_case;
}

int main(int argc,char** argv)
{
  vector<param*> prms;
  //  prms.push_back(make_param('o',"oceancolormap",boceancolormap));

  vector<string> args;
  string usage = 
    "Usage: Gmodmeas <graphfile> <treefile>\n"
    "Copyright (c) 2007, Pau Fernandez";

  parse_params_ex(prms,argc,argv,usage,"Gmodmeas",args,2);
  string grfile=args[0],treefile=args[1];

  // Read the tree
  resgraph R;
  if (!read_graph_from_file(R,treefile)) {
    cerr << "Couldn't read graph " << treefile << endl;
    return -1;
  }

  // Read the graph
  graph G;
  if (!read_graph_from_file(G,grfile)) {
    cerr << "Couldn't read graph " << grfile << endl;
    return -1;
  }

  map<string,int> name_map;
  for (uint k=0;k<G.num_nodes();k++) {
    name_map[G.tag(k)]=k;
  }

  // Calculate
  for (uint k=0;k<R.num_nodes();k++) {
    if (R.indegree(k) > 0) {
      uint root=k;

//       cerr << "Calculating for " << R.tag(k)._name << endl;

      // collect submodules
      vector< set<int> > submodules;
      resgraph::edge_const_iterator ei=R.nbrs_const_iterate(root);
      bool all_leafs=true;
      uint total_sz=0;
      for (;!ei.end();++ei) {
	if (R.indegree(ei.index()) > 0) {
	  all_leafs=false;
	  vector< set<int> >::iterator si=submodules.insert(submodules.end(),set<int>());
	  add_leafs(R,name_map,ei.index(),*si);
	  total_sz+=si->size();
	}
	else {
	  set<int> s;
	  map<string,int>::iterator it=name_map.find(R.tag(ei.index())._name);
	  assert(it!=name_map.end());
	  int idx=it->second;
	  s.insert(idx);
	  submodules.push_back(s);
	  total_sz++;
	}
      }

      // modularity measure
      if (!all_leafs && total_sz > 4) {
	double sz=double(R.tag(root)._size)/double(R.tag(0)._size);
	cout << R.tag(root)._name << ' ' << sz << ' ';
	double mm=modularity_measure(G,submodules);
	cout << mm << endl;
      }
    }
  }
}
