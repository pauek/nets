// 
//  Copyright (c) 2007, Pau Fern√°ndez
//

#ifndef _STATS_H_
#define _STATS_H_

#include <vector>
#include <map>
#include <set>
#include <queue>
#include <stack>

#include <utils/progress.H>
#include <utils/stl.H>

#include "adj_list.H"

namespace net 
{
  ///////////////////////////////////////////////////////////////////////////////
  // Degree distributions

  template<class Network>
  pair<double,double> avg_degree(const Network& g)
  {
    const uint N=g.num_nodes();

    uint D = 0;
    for (uint i=0;i<N;i++) 
      D += g.degree(i);

    return double(D)/double(N);
  }

  template<class Network>
  pair<double,double> avg_indegree(const Network& g)
  {
    const uint N=g.num_nodes();

    uint D = 0;
    for (uint i=0;i<N;i++) 
      D += g.indegree(i);

    return double(D)/double(N);
  }

  template<class Network>
  pair<double,double> avg_outdegree(const Network& g)
  {
    const uint N=g.num_nodes();

    uint D = 0;
    for (uint i=0;i<N;i++) 
      D += g.outdegree(i);

    return double(D)/double(N);
  }

  template<class Network>
  void _degree_distribution(const Network& g,vector<double>& dd,uint kind)
  {
    const uint N=g.num_nodes();
    pair<uint,uint> mx=g.max_degree();

    switch (kind) {
    case 0: dd.resize(mx.first+mx.second+1); break;
    case 1: dd.resize(mx.first+1); break;
    case 2: dd.resize(mx.second+1); break;
    }
    fill(dd.begin(),dd.end(),0.0);
    for (uint i=0;i<N;i++) {
      switch (kind) {
      case 0: dd[g.degree(i)]+=1.0; break;
      case 1: dd[g.indegree(i)]+=1.0; break;
      case 2: dd[g.outdegree(i)]+=1.0; break;
      default: assert(false);
      }
    }
  }

  template<class Network>
  void degree_distribution(const Network& g,vector<double>& dd)
  { _degree_distribution(g,dd,0); }

  template<class Network>
  void indegree_distribution(const Network& g,vector<double>& dd)
  { _degree_distribution(g,dd,1); }

  template<class Network>
  void outdegree_distribution(const Network& g,vector<double>& dd)
  { _degree_distribution(g,dd,2); }

  // Cumulative...

  void make_cumulative(const vector<double>& dd,map<uint,double>& cdd)
  {
    double acum=0;
    cdd.clear();
    vector<double>::const_reverse_iterator it=dd.rbegin(),itend=dd.rend();
    for (uint k=dd.size()-1;it!=itend;++it,--k) {
      if ((*it)!=0) {
	acum+=(*it);
	cdd[k]=acum;
      }
    }
  }

  template<class Network>
  void cum_degree_distribution(const Network& g,map<uint,double>& cdd)
  {
    vector<double> dd;
    degree_distribution(g,dd);
    make_cumulative(dd,cdd);
  }

  template<class Network>
  void cum_indegree_distribution(const Network& g,map<uint,double>& cdd)
  {
    vector<double> dd;
    indegree_distribution(g,dd);
    make_cumulative(dd,cdd);
  }

  template<class Network>
  void cum_outdegree_distribution(const Network& g,map<uint,double>& cdd)
  {
    vector<double> dd;
    outdegree_distribution(g,dd);
    make_cumulative(dd,cdd);
  }

  ///////////////////////////////////////////////////////////////////////////////
  // Clustering

  template<class Network>
  double clustering_coeff(const Network& g,uint i)
  {
    assert(g.is_undirected());
    if (g.degree(i) < 2) return 0.0;
    const uint degr=g.degree(i);
    typename Network::edge_const_iterator ni=g.nbrs_const_iterate(i);
    uint acum=0;
    for(;!ni.end();++ni) {
      if (ni.index() != i) {
	typename Network::edge_const_iterator ei=g.nbrs_const_iterate(ni.index());
	for (;!ei.end();++ei) 
	  if (g.connected(i,ei.index())) ++acum;
      }
    }
    return double(acum)/double(degr*(degr-1));
  }

  template<class Network>
  double avg_clustering_coeff(const Network& g)
  {
    assert(g.is_undirected());
    vector<double> cc(g.num_nodes());
    for (uint i=0;i<g.num_nodes();i++) cc[i]=clustering_coeff(g,i);
    return average(cc);
  }

  template<class Network>
  void clustering_dist(const Network& g,vector<double>& cd)
  {
    assert(g.is_undirected());
    vector<double> dd;
    degree_distribution(g,dd);
    cd.resize(dd.size());
    fill(cd.begin(),cd.end(),0.0);
    const uint N=g.num_nodes();
    for (uint i=0;i<N;i++) {
      const uint d=g.degree(i);
      if (d > 0) cd[d]+=clustering_coeff(g,i)/dd[d];
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Minimum distance

  template<class Network>
  void minimum_distance(const Network& g,uint from_node,
			const vector<bool>& to_mask,vector<int>& min_dist)
  {
    assert(g.num_nodes() == min_dist.size());
    assert(g.num_nodes() == to_mask.size());

    vector<bool> processed(g.num_nodes(),false);

    // Breadth First Search
    std::queue<uint>   tbp; // _to _be _processed
    min_dist[from_node]=0;
    processed[from_node]=true;
    tbp.push(from_node);
    while (!tbp.empty()) {
      const uint curr=tbp.front();
      typename Network::edge_const_iterator ei=g.nbrs_const_iterate(curr);
      for (;!ei.end();++ei) {
	const uint nbr=ei.index();
	if (to_mask[nbr] && !processed[nbr]) {
	  min_dist[nbr]=min_dist[curr]+1; // not processed yet
	  processed[nbr]=true;
	  tbp.push(nbr);
	}
      }
      tbp.pop();
    }
  }

  template<class Network>
  double avg_minimum_distance(const Network& g)
  {
    const uint N=g.num_nodes();
    vector<bool> mask_all(N,true);
    vector<int> mindist(N,-1);
    double acum=0.0;
    for (uint i=0;i<N;i++) {
      minimum_distance(g,i,mask_all,mindist);
      for (uint j=0;j<N;j++) 
	if (j!=i) {
	  acum+=mindist[j];
	  if (mindist[j]<0) return -1.0; // some nodes are not accessible from others
	}
    }
    return acum/double(N*(N-1));
  }

  template<class Network>
  void minimum_distance(const Network& g,vector< vector<int> >& mindist)
  {
    assert(g.num_nodes() == mindist.size());

    const uint N=g.num_nodes();
    vector<bool> all(N,true);
    for (uint i=0;i<N;i++) {
      mindist[i]=vector<int>(N,-1);
      minimum_distance(g,i,all,mindist[i]);
    }
  }


  ////////////////////////////////////////////////////////////////////
  // Efficiency

  template<class Network>
  double efficiency(const Network& g,const set<uint>& subset)
  {
    // isolated nodes with one or no neighbours shoud have local eff. = 0
    if (subset.size() <= 1) return 0.0; 
  
    // prepare mask
    vector<bool> is_in(g.num_nodes(),false);
    set<uint>::const_iterator i=subset.begin(),iend=subset.end();
    for (;i!=iend;i++) is_in[*i]=true;

    const uint SZ=subset.size();
    double acum=0.0;
    i=subset.begin();
    for (;i!=iend;i++)  {
      const uint u=*i;
      vector<int> min_dist(g.num_nodes(),0);
      minimum_distance(g,u,is_in,min_dist);
    
      set<uint>::const_iterator k=subset.begin(),kend=subset.end();
      for (;k!=kend;k++)
	if (u != *k && min_dist[*k] != 0) 
	  acum+=1.0/double(min_dist[*k]);
    }

    return acum/double(SZ*(SZ-1));
  }

  template<class Network>
  double local_efficiency(const Network& g,uint node)
  {
    assert(g.is_undirected());
    assert(node >= 0 && node < g.num_nodes());
    set<uint> subset;
    typename Network::edge_const_iterator ei=g.nbrs_const_iterate(node);
    for (;!ei.end();++ei) subset.insert(ei.index());
    return efficiency(g,subset);
  }

  template<class Network>
  double global_efficiency(const Network& g)
  {
    assert(g.is_undirected());
    set<uint> all_nodes;
    for (uint i=0;i<g.num_nodes();i++) all_nodes.insert(i);
    return efficiency(g,all_nodes);
  }

  template<class Network>
  double avg_local_efficiency(const Network& g)
  {
    assert(g.is_undirected());
    double acum=0.0;
    const uint N=g.num_nodes();
    for (uint i=0;i<N;i++)
      acum+=local_efficiency(g,i);
    
    return acum/double(N);
  }

  /////////////////////////////////////////////////////////////////////////
  // Components

  inline int find_first_untouched(const vector<int>& comp,
				  const vector<bool>& mask)
  {
    vector<int>::const_iterator ti=comp.begin(),tiend=comp.end();
    vector<bool>::const_iterator mi=mask.begin();
    for (uint idx=0;ti!=tiend;++ti,++mi,++idx) 
      if ((*ti)==-1 && (*mi)) return idx;
    return -1;
  }

  template<class Network>
  uint connected_components(const Network& g,vector<int>& comp)
  {
    vector<bool> mask(g.num_nodes(),true);
    return connected_components(g,comp,mask);
  }

  template<class Network>
  uint connected_components(const Network& g,vector<int>& comp,const set<int>& subset)
  {
    // prepare mask
    vector<bool> mask(g.num_nodes(),false);
    set<int>::const_iterator ssi=subset.begin(),ssiend=subset.end();
    for (;ssi!=ssiend;++ssi) mask[*ssi]=true;

    return connected_components(g,comp,mask);
  }

  template<class Network>
  uint connected_components(const Network& g,vector<int>& comp,
			    const vector<bool>& mask)
  {
    assert(g.num_nodes() == comp.size());
    assert(g.num_nodes() == mask.size());
    fill(comp.begin(),comp.end(),-1);

    uint currcomp=0;
    queue<uint> Q;
    int idx=find_first_untouched(comp,mask);
    while (idx!=-1) {
      Q.push(idx);
      comp[idx]=currcomp;
      while (!Q.empty()) {
	uint curr=Q.front();
	typename Network::edge_const_iterator ei=g.nbrs_const_iterate(curr);
	for (;!ei.end();++ei) {
	  const uint u=ei.index();
	  if (comp[u]<0) {
	    comp[u]=currcomp;
	    Q.push(u);
	  }
	}
	Q.pop();
      }
      idx=find_first_untouched(comp,mask);
      currcomp++;
    }

    return currcomp;
  }

  template<class Network>
  void component_sizes(const Network& g,vector<uint>& comp_sz)
  {
    vector<int> comp(g.num_nodes());
    uint ncc=connected_components(g,comp);

    comp_sz.resize(ncc);
    fill(comp_sz.begin(),comp_sz.end(),0);
    for (uint i=0;i<comp.size();i++) {
      comp_sz[comp[i]]++;
    }
    sort(comp_sz.begin(),comp_sz.end(),greater<int>());
  }

  template<class Network>
  void num_component_sizes(const Network& g,vector< pair<uint,uint> >& n_comp_sz)
  {
    vector<uint> comp_sz;
    component_sizes(g,comp_sz);

    n_comp_sz.clear();
    uint count=1,SZ=comp_sz.size();
    for (uint i=1;i<SZ;i++) {
      if (comp_sz[i-1] != comp_sz[i]) {
	n_comp_sz.push_back(make_pair(count,comp_sz[i-1]));
	count=1;
      }
      else count++;
    }
    n_comp_sz.push_back(make_pair(count,comp_sz[SZ-1]));
  }

  template<class Network>
  void extract_component(const Network& g,Network& subg,uint num)
  {
    assert(g.is_undirected());
    vector<int> comp(g.num_nodes());
    uint ncc=connected_components(g,comp);
    if (num > ncc) {
      cerr << "extract_component: The graph has only " << ncc << " components" << endl;
      exit(-1);
    }

    // Organize component sizes in order to choose one
    vector<int> comp_sz(ncc,0);
    for (uint k=0;k<comp.size();k++) comp_sz[comp[k]]++;
    typedef multimap<uint,uint, greater<uint> > my_map;
    my_map sz_and_idx;
    for (uint k=0;k<comp_sz.size();k++) sz_and_idx.insert(make_pair(comp_sz[k],k));
    my_map::const_iterator mi=sz_and_idx.begin();
    for (uint k=0;k<num;k++) mi++;
    uint which=mi->second,sz=mi->first;
  
    // Prepare mask
    vector<bool> mask(g.num_nodes(),false);
    for (uint i=0;i<g.num_nodes();i++) {
      if (comp[i]==int(which)) mask[i]=true;
    }

    subg=Network(sz);
    g.subgraph(mask,subg);
  }


  /////////////////////////////////////////////////////////////////
  // Topological Overlap
  
  template<class Network>
  void topological_overlap(const Network& g,vector< map<int,double> >& t_ov)
  {
    assert(g.is_undirected());
    assert(g.num_nodes() == t_ov.size());
    const uint N=g.num_nodes();
    for (uint i=0;i<N;i++) {
      for (uint j=0;j<N;j++) {
	if (i!=j) {
	  int cn=g.common_neighbours(i,j);
	  int md=min(g.degree(i),g.degree(j));
	  if (g.connected(i,j)) cn++;
	  if (cn!=0) t_ov[i][j]=double(cn)/double(md);
	}
      }
    }
  }

  /////////////////////////////////////////////////////////////////
  // Correlations

  template<class Network>
  void avg_nbrs_degree(const Network& g,vector<double>& nd)
  {
    vector<double> dd;
    degree_distribution(g,dd);
    nd.resize(dd.size());
    fill(nd.begin(),nd.end(),0.0);
    const uint N=g.num_nodes();
    for (uint i=0;i<N;i++) {
      const uint d=g.degree(i);
      if (d > 0) {
	double acum=0.0;
	typename Network::edge_const_iterator ei=g.nbrs_const_iterate(i);
	for (;!ei.end();++ei) 
	  acum+=double(g.degree(ei.index()));
	acum/=double(g.degree(i));
	nd[d]+=acum/dd[d];
      }
    }
  }

  template<class Network>
  double assortative_mixing_coeff(const Network& g)
  {
    const uint N=g.num_nodes();
    double sumjk=0.0,sumj=0.0,sumk=0.0,sumj2=0.0,sumk2=0.0;

    // excess-degree
    vector<uint> eindegr(N,0),eoutdegr(N,0);
    for (uint v=0;v<N;v++) {
      eindegr[v]=g.indegree(v);
      eoutdegr[v]=g.outdegree(v);
    }
    // mixing
    ulong l=0;
    for (uint v=0;v<N;v++) {
      typename Network::edge_const_iterator nv=g.nbrs_const_iterate(v);
      for (;!nv.end();++nv) {
	const uint j=eindegr[v],k=eoutdegr[nv.index()];
	sumjk+=double(j*k);
	sumj+=double(j),sumj2+=double(j*j);
	sumk+=double(k),sumk2+=double(k*k);
	l++;
      }
    }
    double L=double(l); // num edges
    return (sumjk-sumj*sumk/L)/
      sqrt((sumj2-sumj*sumj/L)*(sumk2-sumk*sumk/L));
  }

  ////////////////////////////////////////////////////////////////////////
  // Correlation Matrix + logBin

  template<class Network>
  void correlation_matrix(const Network& g,double bins_per_decade,
			  vector< vector<double> >& M,bool blog,bool bsymm)
  {
    // asumimos que M est· dimensionada y que est· a 0

    const uint N=g.num_nodes();
    const double log10=log(10.0);
    for (uint k=0;k<N;k++) {
      typename Network::edge_const_iterator ei=g.nbrs_const_iterate(k);
      const uint indegrk=g.indegree(k);
      for (;!ei.end();++ei) {
	uint i,j;
	if (blog) {
	  i=uint(floor(bins_per_decade*log(double(indegrk))/log10));
	  j=uint(floor(bins_per_decade*log(double(g.outdegree(ei.index())))/log10));
	}
	else {
	  i=indegrk;
	  j=g.outdegree(ei.index());
	}
      
	assert(i < M[0].size() && j < M.size());
	if (bsymm) {
	  M[j][i]+=0.5;
	  M[i][j]+=.5;
	}
	else {
	  M[j][i]+=1.0;
	}
      }
    }
  }

  template<class Network>
  void randomize_connectivity(Network& g,ulong times)
  {
    uint k=0;
    while (k<times) {
      pair<uint,uint> e1=g.random_edge();
      pair<uint,uint> e2=g.random_edge();
      // No pueden ser autolinks. El link que se anyadira no puede existir.
      // Debe haber vertices diferentes en e1 y e2.
      if (e1.first  != e2.first  && e1.first  != e2.second &&
	  e1.second != e2.first  && e1.second != e2.second &&
	  !g.connected(e1.first,e2.second) &&
	  !g.connected(e2.first,e1.second)) {
	g.remove_edge(e1.first,e1.second);
	g.remove_edge(e2.first,e2.second);
	g.add_edge(e1.first,e2.second);
	g.add_edge(e2.first,e1.second);
	++k;
      } 
    }
  }

//   inline void _print(vector< vector<double> >& M)
//   {
//     cerr << endl;
//     cerr.setf(ios::fixed);
//     cerr.precision(6);
//     vector< vector<double> >::const_iterator i=M.begin(),iend=M.end();
//     for (;i!=iend;++i) {
//       print_all(*i,cerr," ");
//       cerr << endl;
//     }
//     cerr << endl;
//   }

  template<class Network>
  void corr_matrix_vs_random(const Network& g,double bins_per_decade,uint IT,
			     vector< vector<double> >& P,vector< vector<double> >& Z,
			     bool blog,bool bsymm,progress_bar* p)
  {
    pair<uint,uint> _mx=g.max_degree();
    uint MXsz=(blog?(1+uint(floor(bins_per_decade*log(double(_mx.first))/log(10.0)))):_mx.first+1);
    uint MYsz=(blog?(1+uint(floor(bins_per_decade*log(double(_mx.second))/log(10.0)))):_mx.second+1);
    if (bsymm) {
      if (MXsz > MYsz) MYsz=MXsz;
      if (MYsz > MXsz) MXsz=MYsz;
    }
    vector< vector<double> > _M(MYsz,vector<double>(MXsz,0.0));
    vector< vector<double> > _Mr=_M,_Mr_avg=_M,_Mr_avg2=_M;
  
    correlation_matrix(g,bins_per_decade,_M,blog,bsymm);
    {
      for (uint k=0;k<IT;k++) {
	Network gcopy=g;
	randomize_connectivity(gcopy,g.num_edges()*3);
	correlation_matrix(gcopy,bins_per_decade,_Mr,blog,bsymm);
	for (uint j=0;j<MYsz;j++) {
	  for (uint i=0;i<MXsz;i++) {
	    _Mr_avg[j][i]+=_Mr[j][i];
	    _Mr_avg2[j][i]+=_Mr[j][i]*_Mr[j][i];
	    _Mr[j][i]=0; // for next time
	  }
	}
	if (p!=0) p->incr();
      }
    }

    for (uint j=0;j<MYsz;j++) {
      for (uint i=0;i<MXsz;i++) {
	_Mr_avg[j][i]/=double(IT);
	_Mr_avg2[j][i]/=double(IT);
      }
    }

    vector< vector<double> > _Mr_stddev=_Mr_avg;
    for (uint j=0;j<MYsz;j++) {
      for (uint i=0;i<MXsz;i++) 
	_Mr_stddev[j][i]=sqrt(_Mr_avg2[j][i]-_Mr_avg[j][i]*_Mr_avg[j][i]);
    }

    Z=P=_Mr_avg;
    for (uint j=0;j<MYsz;j++) 
      for (uint i=0;i<MXsz;i++) { 
	if (_Mr_avg[j][i]) P[j][i]=_M[j][i]/_Mr_avg[j][i];
	double diff=_M[j][i]-_Mr_avg[j][i];
	if (_Mr_stddev[j][i]==0.0 || fabs(diff) < 1e-20) 
	  Z[j][i]=0.0;
	else
	  Z[j][i]=diff/_Mr_stddev[j][i];
      }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Percolation threshold

  template<class Network>
  void local_percolation(const Network& g,vector<int>& v)
  {
    const uint N=g.num_nodes();
    v.resize(N);
    if (g.is_undirected()) {
      vector<int>::iterator vi=v.begin(),viend=v.end();
      for (uint i=0;vi!=viend;++vi,++i) {
	const uint k=g.degree(i);
	*vi=k*(k-2);
      }
    }
    else {
      vector<int>::iterator vi=v.begin(),viend=v.end();
      for (uint i=0;vi!=viend;++vi,++i) {
	const int j=g.indegree(i);
	const int k=g.outdegree(i);
	*vi=2*j*k-j-k;
      }
    }
  }  

  template<class Network>
  double percolation_threshold(const Network& g)
  {
    vector<int> locperc;
    local_percolation(g,locperc);
    vector<int>::const_iterator lpi=locperc.begin(),lpiend=locperc.end();
    double acum=0.0;
    for (;lpi!=lpiend;++lpi) acum+=double(*lpi);
    return acum/double(locperc.size());
  }

  //////////////////////////////////////////////////////////////////////////////
  // Betweenness centrality
  //
  // Based on: "A Faster Algorithm for Betweenness Centrality", Ulrik Brandes, J.Math.Sociol. 25(2):163-177 (2001)

  template<class Network>
  void betweenness_init_result_network(const Network& g,
				       adj_list<pair<typename Network::node_tag,double>,double>& res)
  {
    const uint N=g.num_nodes();
    res.resize_and_clear(N);
    for (uint i=0;i<N;i++) {
      res.tag(i)=make_pair(g.tag(i),0.0);
      typename Network::edge_const_iterator ei=g.nbrs_const_iterate(i);
      for (;!ei.end();++ei) res.add_edge(ei.index(),i,0.0); // transpose also...
    }
    if (g.is_undirected()) res.to_undirected();
  }
  
  template<class Network>
  void betweenness_reset_result_network(Network& res,const set<int>& subset)
  {
    const uint N=res.num_nodes();

    // precumpute mask
    vector<bool> mask(N,false);
    set<int>::const_iterator ssi=subset.begin(),ssiend=subset.end();
    for (;ssi!=ssiend;++ssi) mask[*ssi]=true;

    // erase links' weights
    for (uint i=0;i<N;i++) 
      if (mask[i]) {
	res.tag(i).second=0.0;
	typename Network::edge_iterator ei=res.nbrs_iterate(i);
	for (;!ei.end();++ei) ei.set_tag(0.0);
      }
  }  

  template<class Network>
  void betweenness_centrality(Network& res,const set<int>& subset,bool bnormalize=false)
  {
    // precumpute mask
    vector<bool> mask(res.num_nodes(),false);
    set<int>::const_iterator ssi=subset.begin(),ssiend=subset.end();
    for (;ssi!=ssiend;++ssi) mask[*ssi]=true;

    betweenness_centrality(res,mask,bnormalize);
  }

  template<class Network>
  void betweenness_centrality(Network& res,bool bnormalize=false)
  {
    vector<bool> mask(res.num_nodes(),true);
    betweenness_centrality(res,mask,bnormalize);
  }
    
  template<class Network>
  void betweenness_centrality(Network& res,vector<bool>& mask,bool bnormalize=false)
  {
    const uint N=res.num_nodes();
    for (uint i=0;i<N;i++) 
      if (mask[i]) {
	stack<uint> S;
	vector< set<uint> > pred(N);
	vector<uint> s(N,0);   s[i]=1;
	vector<int>  d(N,-1);  d[i]=0;

	// breadth first search 
	queue<uint> Q;
	Q.push(i);
	while (!Q.empty()) {
	  uint v=Q.front(); Q.pop();
	  S.push(v);
	  typename Network::edge_const_iterator ei=res.nbrs_const_iterate(v);
	  for (;!ei.end();++ei) {
	    const uint w=ei.index();
	    if (d[w] < 0) {
	      Q.push(w);
	      d[w]=d[v]+1;
	    }
	    if (d[w]==(d[v]+1)) {
	      s[w]+=s[v];
	      pred[w].insert(v);
	    }
	  }
	}

	// back-propagation
	const double startval=(bnormalize?1.0/double(N*(N-1)):1.0);
	vector<double> dd(N,0.0);
	while (!S.empty()) {
	  uint w=S.top(); S.pop();
	  set<uint>::const_iterator v=pred[w].begin(),vend=pred[w].end();
	  for (;v!=vend;++v) {
	    double flux=double(s[*v])/double(s[w])*(startval+dd[w]);
	    dd[*v]+=flux; 
	    res.get_edge_tag(*v,w)+=flux;
	    if (res.is_undirected()) res.get_edge_tag(w,*v)+=flux;
	    if (w != i) res.tag(w).second+=dd[w];
	  }
	}
      }
  } 
}

#endif
