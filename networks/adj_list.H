// 
//  Copyright (c) 2007, Pau FernÃ¡ndez
//

#ifndef _ADJ_H_
#define _ADJ_H_

/////////////////////////////////////////////////////////////////////////
// Adjacency list

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include <assert.h>

// Blitz++
#include <random/discrete-uniform.h> 

namespace net 
{
  /////////////////////////////////////////
  // NodeTag: DefaultConstructible
  // EdgeTag: DefaultConstructible 

  struct null {};
  std::ostream& operator<<(std::ostream& o,const null& n) { return o; }
  std::istream& operator>>(std::istream& i,null& n) { return i; }

  template<class NodeTag=null,class EdgeTag=null>
  class adj_list
  {
    // Node structure
    struct node {
      uint _id,_od; // indegree + outdegree
      NodeTag _tag;
      node(uint i=0,uint o=0):_id(i),_od(o) {}
      node(uint i,uint o,NodeTag t):_id(i),_od(o),_tag(t) {}
    };

    // Edge structure
//     struct edge {
//       uint _idx;
//       EdgeTag _tag;
//       edge(uint i):_idx(i) {}
//       edge(uint i,EdgeTag t):_idx(i),_tag(t) {}
//       bool operator<(const edge& e) const { return _idx < e._idx; }
//       bool operator>(const edge& e) const { return _idx > e._idx; }
//     };

    struct _edge_iterator {
      typename std::map<uint,EdgeTag>::iterator _pos,_end;
      _edge_iterator(std::map<uint,EdgeTag>& m):_pos(m.begin()),_end(m.end()) {}
      bool end() const { return _pos==_end; }
      uint index() const { return _pos->first; }
      const EdgeTag& tag() const { return _pos->second; }
      void set_tag(const EdgeTag& et) { _pos->second=et; }
      void operator++(int) { _pos++; }
      void operator++() { ++_pos; }
    };

    struct _edge_const_iterator {
      typename std::map<uint,EdgeTag>::const_iterator _pos,_end;
      _edge_const_iterator(const std::map<uint,EdgeTag>& m):_pos(m.begin()),_end(m.end()) {}
      bool end() const { return _pos==_end; }
      uint index() const { return _pos->first; }
      EdgeTag tag() const { return _pos->second; }
      void operator++(int) { _pos++; }
      void operator++() { ++_pos; }
    };
  
    // data members
    uint         _N,_E;
    bool         _bundir;
    std::vector< std::map<uint,EdgeTag> > _adj;
    std::vector<node> _nod; // keep track of the degree to know what nodes are isolated
    std::set<uint>    _isolated;  // list of isolated nodes

    bool add_dir_edge(uint u,uint v,EdgeTag t);    // returns true if success
    bool remove_dir_edge(uint u,uint v); // returns true if success
  public:
    typedef _edge_iterator edge_iterator;
    typedef _edge_const_iterator edge_const_iterator;
    typedef NodeTag node_tag;
    typedef EdgeTag edge_tag;

    // Constructors
    adj_list(bool bundir=false):_N(0),_E(0),_bundir(bundir) {}

    adj_list(uint SZ,bool bundir=false)
      :_N(SZ),_E(0),_bundir(bundir),
       _adj(SZ),_nod(SZ)
    { for (uint k=0;k<SZ;k++) _isolated.insert(k),_nod[k]._tag=node_tag(); }

    // Getting info
    uint num_nodes() const { return _N; }
    bool is_undirected() const { return _bundir; }
    uint num_non_isolated() const { return _N-_isolated.size(); }
    uint num_edges() const { return (_bundir?_E/2:_E); }
    double average_connectivity() const 
    { return double(num_edges()*2)/double(num_nodes()); }

    bool connected(uint u,uint v) const
    { assert(bounds_ok(u) && bounds_ok(v)); return _adj[u].find(v)!=_adj[u].end(); }

    uint indegree(uint k) const { assert(bounds_ok(k)); return _nod[k]._id; }
    uint outdegree(uint k) const { assert(bounds_ok(k)); return _nod[k]._od; }
    uint degree(uint k) const 
    { assert(bounds_ok(k)); return (_bundir?_nod[k]._id:_nod[k]._id+_nod[k]._od); }

    pair<uint,uint> max_degree() const; // max (in,out)-degree

    // Tags
    void set_node_tag(uint k,const NodeTag& t) { assert(bounds_ok(k)); _nod[k]._tag=t; }
    NodeTag& tag(uint k) { assert(bounds_ok(k)); return _nod[k]._tag; }
    NodeTag tag(uint k) const { assert(bounds_ok(k)); return _nod[k]._tag; }
    bool empty_tag(uint k) const {
      std::ostringstream o;
      o << _nod[k]._tag;
      string s=o.str();
      return (s.empty());
    }

    void index_from_tag(const NodeTag& nt,std::set<uint>& idxs) const
    {
      typename std::vector<node>::const_iterator ni=_nod.begin(),niend=_nod.end();
      for (uint k=0;ni!=niend;++ni,++k) 
	if (nt==ni->_tag) idxs.insert(k);
    }

    // PROBLEMÓN: Y los grafos undirected qué?? Cómo se refrescan los tags de los 2 links???
    EdgeTag& get_edge_tag(uint u,uint v)
    { 
      assert(connected(u,v));
      if (is_undirected()) {
	typename std::map<uint,EdgeTag>::iterator ai1=_adj[u].find(v),ai2=_adj[v].find(u);
	assert(ai1!=_adj[u].end() && ai2!=_adj[v].end());
	assert(ai1->second == ai2->second);
	return ai1->second;
      }
      else {
	typename std::map<uint,EdgeTag>::iterator ai=_adj[u].find(v);
	assert(ai!=_adj[u].end());
	return ai->second;
      }
    }

    void set_edge_tag(uint u,uint v,const EdgeTag& t)
    { 
      if (is_undirected()) {
	typename std::map<uint,EdgeTag>::iterator ai1=_adj[u].find(v),ai2=_adj[v].find(u);
	assert(ai1!=_adj[u].end() && ai2!=_adj[v].end());
	assert(ai1->second == ai2->second);
	ai1->second=ai2->second=t;
      }
      else {
	typename std::map<uint,EdgeTag>::iterator ai=_adj[u].find(v);
	assert(ai!=_adj[u].end());
	ai->second=t;
      }
    }

    bool isolated(uint k) const 
    { assert(bounds_ok(k)); return (_nod[k]._id==0 && _nod[k]._od==0); }

    uint first_free() const {
      if (_isolated.empty()) return -1;
      else return *(_isolated.begin());
    }

    uint common_neighbours(uint i,uint j) const {
      assert(bounds_ok(i) && bounds_ok(j));
      uint count=0;
      const std::map<uint,EdgeTag>& nbrsi=_adj[i],nbrsj=_adj[j];
      typename std::map<uint,EdgeTag>::const_iterator ui=nbrsi.begin(),uiend=nbrsi.end();
      for (;ui!=uiend;++ui) {
	if (nbrsj.find(ui->first) != nbrsj.end()) ++count;
      }
      return count;
    }

    // Iterate neighbours
    edge_iterator nbrs_iterate(uint i)
    { assert(bounds_ok(i)); return edge_iterator(_adj[i]); }
    edge_const_iterator nbrs_const_iterate(uint i) const
    { assert(bounds_ok(i)); return edge_const_iterator(_adj[i]); }

    uint nbr_index(edge_iterator& ei) { return ei._pos->_idx; }
    EdgeTag nbr_tag(edge_iterator& ei) { return ei._pos->_tag; }
    const uint& nbr_index(edge_const_iterator& ei) const { return ei->_idx; }
    const EdgeTag& nbr_tag(edge_const_iterator& ei) const { return ei->_tag; }

    // check
    bool bounds_ok(uint i) const { return (i>=0 && i<_N); }

    // Modify operations
    void clear() {  // remove all connections
      typename std::vector< std::map<uint,EdgeTag> >::iterator a=_adj.begin(),aend=_adj.end();
      for (;a!=aend;a++) a->clear();
      // fill(_nod.begin(),_nod.end(),node(0,0)); // clear tags?
      for (uint k=0;k<_N;k++) _isolated.insert(k);
      _E=0;
    }

    void resize_and_clear(uint N) {
      _N=N;
      _E=0;
      _nod.resize(N);
      _adj.resize(N);
      _isolated.clear();
      typename std::vector< std::map<uint,EdgeTag> >::iterator a=_adj.begin(),aend=_adj.end();
      for (;a!=aend;++a) a->clear();
      fill(_nod.begin(),_nod.end(),node(0,0));
      for (uint k=0;k<N;k++) _isolated.insert(k);
    }

    uint add_node(NodeTag tag) 
    { 
      _adj.push_back(std::map<uint,EdgeTag>());
      _nod.push_back(node(0,0,tag)),++_N;
      const uint newidx=_N-1;
      _isolated.insert(newidx);
      return newidx;
    }

    uint add_node() { return add_node(NodeTag()); }

    bool add_edge(uint u,uint v,EdgeTag e); // returns true if success
    bool add_edge(uint u,uint v) { return add_edge(u,v,EdgeTag()); } // returns true if success
    bool remove_edge(uint u,uint v); // returns true if success

    // Random edge
    pair<uint,uint> random_edge() const;
    pair<uint,uint> random_free_edge(bool ballowselflinks=false) const;

    // Subgraph
    void subgraph(const std::vector<bool>& mask,adj_list& adj) const;

    // To undirected
    bool check_undirected() const;
    void to_undirected();
    void transpose();
  };

  /////////////////////////////////////
  // new edge:     u <---- v
  template<class NodeTag,class EdgeTag>
  bool adj_list<NodeTag,EdgeTag>::add_dir_edge(uint u,uint v,EdgeTag t)
  {
    assert(bounds_ok(u) && bounds_ok(v));
    pair<typename std::map<uint,EdgeTag>::iterator,bool> res=_adj[u].insert(make_pair(v,t));
    _isolated.erase(u),_isolated.erase(v);
    if (res.second) {
      ++_nod[u]._id,++_nod[v]._od;
      ++_E;
    }
    return res.second;
  }

  template<class NodeTag,class EdgeTag>
  bool adj_list<NodeTag,EdgeTag>::add_edge(uint u,uint v,EdgeTag t)
  {
    if (_bundir) { 
      add_dir_edge(u,v,t); 
      return add_dir_edge(v,u,t); 
    }
    else 
      return add_dir_edge(u,v,t);
  }

  template<class NodeTag,class EdgeTag>
  bool adj_list<NodeTag,EdgeTag>::remove_dir_edge(uint u,uint v)
  {
    assert(bounds_ok(u) && bounds_ok(v));
    uint res=_adj[u].erase(v);
    if (res > 0) {
      --_nod[u]._id,--_nod[v]._od;
      _E -= res;
      if (isolated(u)) _isolated.insert(u);
      if (isolated(v)) _isolated.insert(v);
    }
    return  res > 0;
  }

  template<class NodeTag,class EdgeTag>
  bool adj_list<NodeTag,EdgeTag>::remove_edge(uint u,uint v)
  {
    assert(connected(u,v));
    if (_bundir) { 
      remove_dir_edge(u,v); 
      return remove_dir_edge(v,u); 
    }
    else 
      return remove_dir_edge(u,v);
  }

  template<class NodeTag,class EdgeTag>
  pair<uint,uint> adj_list<NodeTag,EdgeTag>::max_degree() const
  {
    pair<uint,uint> mx;
    typename std::vector<node>::const_iterator ni=_nod.begin(),niend=_nod.end();
    for (;ni!=niend;++ni) {
      if (mx.first < ni->_id) mx.first=ni->_id;
      if (mx.second < ni->_od) mx.second=ni->_od;
    }
    return mx;
  }

  template<class NodeTag,class EdgeTag>
  pair<uint,uint> adj_list<NodeTag,EdgeTag>::random_free_edge(bool ballowselflinks) const
  {
    // find a random edge uniformly among the free ones
    uint free=(ballowselflinks?_N*_N-_E:(_N*(_N-1))-_E);
    assert(free > 0); // something to choose from
    int r;
    if (free>=2) {
      ranlib::DiscreteUniform<uint> drng(free);
      r=drng.random();
    }
    else {
      // no need to decide randomly
      r=0;
    }

    int to=-1;
    while (r >= 0) r-=_N-indegree(++to)-(ballowselflinks?0:1);
    int r2=-r-1;
    int from=-1;
    while (r2 >= 0) {
      ++from;
      while (_adj[to].find(from)!=_adj[to].end() || 
	     (!ballowselflinks && to==from)) ++from;
      --r2;
    }
    return make_pair(to,from);
  }

  template<class NodeTag,class EdgeTag>
  pair<uint,uint> adj_list<NodeTag,EdgeTag>::random_edge() const
  {
    const uint E=_E;
    assert(E > 0);
    int r;
    if (E == 1) r=0;
    else {
      ranlib::DiscreteUniform<uint> drng(E);
      r=drng.random();
    }
    int to=0;
    r-=indegree(to);
    while (r >= 0) r-=indegree(++to);
    int r2=-r;

    _edge_const_iterator ei=nbrs_const_iterate(to);
    while (--r2>0) { ++ei; assert(!ei.end()); }
    return make_pair(to,ei.index());
  }

  template<class NodeTag,class EdgeTag>
  void adj_list<NodeTag,EdgeTag>::subgraph(const std::vector<bool>& mask,adj_list& subg) const
  {
    assert(int(subg.num_nodes()) == count(mask.begin(),mask.end(),true));
    subg.clear();

    std::vector<uint> _map(mask.size(),0);
    uint idx=0;
    for (uint i=0;i<mask.size();i++) 
      if (mask[i]) {
	_map[i]=idx++;
	subg.set_node_tag(_map[i],tag(i));
      }

    for (uint k=0;k<_N;k++) {
      _edge_const_iterator ei=nbrs_const_iterate(k);
      for (;!ei.end();++ei) {
	uint u=k,v=ei.index();
	if (mask[u] && mask[v]) subg.add_edge(_map[u],_map[v],ei.tag());
      }
    }
  }

  template<class NodeTag,class EdgeTag>
  void adj_list<NodeTag,EdgeTag>::to_undirected()
  {
    // repeat necessary edges
    for (uint k=0;k<_N;k++) {
      edge_const_iterator ei=nbrs_const_iterate(k);
      for (;!ei.end();++ei) {
	add_edge(ei.index(),k,ei.tag());
      }
    }
    // switch to undir...
    _bundir=true;
  }

  // change the links' direction
  template<class NodeTag,class EdgeTag>
  void adj_list<NodeTag,EdgeTag>::transpose() 
  {
    // reverse adjlist
    std::vector< std::map<uint,EdgeTag> > _tadj(_adj.size()); 
    typename std::vector< std::map<uint,EdgeTag> >::const_iterator ai=_adj.begin(),aiend=_adj.end();
    for (uint i=0;ai!=aiend;++ai,++i) {
      typename std::map<uint,EdgeTag>::const_iterator ni=ai->begin(),niend=ai->end();
      for (;ni!=niend;++ni) {
	uint j=ni->first;
	_tadj[j].insert(make_pair(i,ni->second));
      }
    }
    _adj=_tadj;

    // reverse indegree + outdegree
    const uint N=_adj.size();
    for (uint k=0;k<N;k++) {
      swap(_nod[k]._id,_nod[k]._od);
    }
  }

  template<class NodeTag,class EdgeTag>
  bool adj_list<NodeTag,EdgeTag>::check_undirected() const
  {
    for (uint k=0;k<_N;k++) {
      edge_const_iterator ei=nbrs_const_iterate(k);
      for (;!ei.end();++ei) {
	if (!connected(ei.index(),k)) return false;
      }
    }
    return true;
  }

  //////////////////////////////////////////////////////////////
  // Undirected Adjacency List

  template<class NodeTag=null,class EdgeTag=null>
  class adj_list_un : public adj_list<NodeTag,EdgeTag>
  {
  public:
    adj_list_un():adj_list<NodeTag,EdgeTag>(true) {}
    adj_list_un(uint N):adj_list<NodeTag,EdgeTag>(N,true) {}
  };
}

#endif
