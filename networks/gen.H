// 
//  Copyright (c) 2007, Pau Fernández
//

#ifndef _GEN_H_
#define _GEN_H_

#include <utils/random.H>

namespace net 
{

  template<class Network>
  void erdos_renyi(Network& g,uint L)
  {
    g.clear(); // remove edges
    
    for (uint l=0;l<L;l++) {
      pair<uint,uint> e=g.random_free_edge(true);
      g.add_edge(e.first,e.second);
    }
  }


  template<class Network>
  void barabasi_albert(Network& g)
  {
    g.to_undirected();
    const uint N=g.num_nodes();
    g.clear();
    g.add_edge(0,1);
    
    for (uint i=2;i<N;i++) 
      {
	// choose preferentially
	uint idx=0;
	DiscreteUniform<uint> drng(2*g.num_edges());
	int which=drng.random()-g.degree(idx);
	while (which >= 0) which-=g.degree(++idx);
	
	// new edge
	g.add_edge(i,idx);
      }
  }

  inline int sample_dist(const vector<double>& p)
  {
#ifndef NDEBUG
    double total=0.0;
    for (uint i=0;i<p.size();i++) total+=p[i];
    assert(fabs(total-1.0) < 1e-7);
#endif
    static Uniform<double> rng;
    int idx=0;
    double r=rng.random()-p[idx];
    while (r > 0.0) r-=p[++idx];
    return idx;
  }

  inline int sample_int_dist(vector<int>& p,int total)
  {
#ifndef NDEBUG
    int tot=0;
    for (uint i=0;i<p.size();i++) tot+=p[i];
    assert(total == tot);
#endif
    uint idx=0;
    DiscreteUniform<uint> drng(total);
    int which=drng.random()-p[idx];
    while (which >= 0) which-=p[++idx];
    return idx;
  }

  void scale_free_dist(vector<double>& probs,uint N,
		       double alpha,double kappa)
  {
    probs.resize(N);
    double total=0.0;
    for (uint i=1;i<N;i++) {
      probs[i]=pow(double(i),-alpha)*exp(-double(i)/double(kappa));
      total+=probs[i];
    }
    for (uint i=0;i<N;i++) probs[i]/=total;
  }

  template<class Network>
  void scale_free_with_cutoff(Network& g,double alpha,double kappa)
  {
    g.to_undirected();
    const uint N=g.num_nodes();

    // compute distribution
    vector<double> probs;
    scale_free_dist(probs,N,alpha,kappa);

    // stubs
    vector<int> stubs(N,0);
    int tot=0;
    for (uint i=0;i<N;i++) tot+=stubs[i]=sample_dist(probs);
    while (tot % 2 == 1) {
      tot-=stubs[0];
      tot+=stubs[0]=sample_dist(probs);
    }

    // edges
    while (tot > 0) {
      uint v1=sample_int_dist(stubs,tot);
      stubs[v1]--; tot--;
      uint v2=sample_int_dist(stubs,tot);
      stubs[v2]--; tot--;

      g.add_edge(v1,v2);
    }
  }

  template<class Network>
  void scale_free_indegree(Network& g,double alpha,double kappa)
  {
    const uint N=g.num_nodes();

    // compute distribution
    vector<double> probs;
    scale_free_dist(probs,N,alpha,kappa);

    // stubs
    vector<int> instubs(N,0);
    int tot=0;
    for (uint i=0;i<N;i++) {
      tot+=instubs[i]=sample_dist(probs);
    }

    // edges
    while (tot > 0) {
      uint v1=sample_int_dist(instubs,tot);
      instubs[v1]--; tot--;
      DiscreteUniform<uint> drng(N);
      uint v2=drng.random();
      g.add_edge(v1,v2);
    }
  }
}

#endif
