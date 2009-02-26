// 
//  Copyright (c) 2007, Pau Fernández
//

#ifndef _HIERCLUST_H_
#define _HIERCLUST_H_

////////////////////////////////////////////////////////////////////////
// Hierarchical clustering

#include <list>

template<class T>
class branch 
{
  T _data;
  int _size;
  branch *_left,*_right;

public:

  typedef T data_type;

  branch(T t):_data(t),_size(1),_left(0),_right(0) {}
  branch(branch* l,branch* r):_left(l),_right(r) {
    _size=_left->size()+_right->size();
  }
  branch():_size(0),_left(0),_right(0) {}
  bool leaf() const { return _left==0 && _right==0; }
  int size() const { return _size; }
  T& data() { return _data; }
  branch* left() const { return _left; }
  branch* right() const { return _right; }
  int depth() { 
    if (_size == 1) return 1;
    else return 1+max(_left->depth(),_right->depth());
  }
  double graphical_depth(double factor) {
    if (_size == 1) return 1;
    else {
      double k=double(min(_left->_size,_right->_size))*factor;
      return k+max(_left->graphical_depth(factor),
		   _right->graphical_depth(factor));
    }
  }

  void set_data(T t)  { _data=t; }
  void set_size(int sz)  { _size=sz; }
  void set_leaf(T t=0) { _data=t;_left=0;_right=0; }
  void set_left(branch* l)  { _left=l; }
  void set_right(branch* r) { _right=r; }
};

pair<int,int> highest_weight(vector< map<uint,double> >& weights)
{
  double max=0.0;
  pair<int,int> res(-1,-1);
  
  vector< map<uint,double> >::iterator i=weights.begin(),i_end=weights.end();
  for (uint id=0;i!=i_end;++i,++id) {
    map<uint,double>::iterator i2=i->begin(),i2_end=i->end();
    while (i2!=i2_end && i2->first < id) {
      if (i2->second > max) {
	res=make_pair(id,i2->first);
	max=i2->second;
      }
      ++i2;
    }
  }
  return res;
}

void collapse_rows_columns(pair<int,int> h,pair<double,double> sizes,
			   vector< map<uint,double> >& weights)
{
  const int SZ=weights.size();
  const double& sz1=sizes.first;
  const double& sz2=sizes.second;

  // Erase references to h.second from others 
  // (we will use the symmetric ones to restore their values)
  vector< map<uint,double> >::iterator w,w_end=weights.end();
  for (w=weights.begin();w!=w_end;++w) w->erase(h.second);

  // Calculate average weights
  for (int k=0;k<SZ;k++) {
    if (k!=h.first && k!=h.second) 
      {
	double w1=0.0,w2=0.0;
	
	map<uint,double>::iterator f1=weights[h.first].find(k);
	if (f1!=weights[h.first].end()) w1=f1->second;
	map<uint,double>::iterator f2=weights[h.second].find(k);
	if (f2!=weights[h.second].end()) w2=f2->second;

	if (w1 > 0.0 || w2 > 0.0) {
	  // const double w=(w1+w2)/2.0;
	  const double w=(w1*sz1+w2*sz2)/(sz1+sz2); // ???  <----
	  weights[h.first][k]=w;
	  weights[k][h.first]=w; // "symmetrize"
	}
      }
    }

  // Erase now that they are integrated with h.first
  weights[h.second].clear();
}

void hierarchical_clustering(vector<branch<int>*>& _roots,
			     vector< map<uint,double> >& _weights)
{
  vector< map<uint,double> > weights=_weights;

  const int SZ=weights.size();
  vector<branch<int>*> tree;
  set<int> roots;

  for (int i=0;i<SZ;++i) {
    tree.push_back(new branch<int>(i));
    roots.insert(i);
  }

  pair<int,int> h=highest_weight(weights);
  while (h.first != -1 && h.second != -1)
    {
      double rsz=double(tree[h.first]->size());
      double lsz=double(tree[h.second]->size());      
      
      // Update weights
      collapse_rows_columns(h,make_pair(rsz,lsz),weights);

      // Update tree
      branch<int>* bnew=new branch<int>(tree[h.first],tree[h.second]);
      tree[h.first]=bnew;
      tree[h.second]=0;
      roots.erase(h.second);

      h=highest_weight(weights);
    }

  _roots.clear();
  set<int>::const_iterator ri=roots.begin(),riend=roots.end();
  for (;ri!=riend;++ri) {
    _roots.push_back(tree[*ri]);
  }
}

template<class Branch>
void extract_leafs(Branch* b, std::list<typename Branch::data_type>& l)
{
  if (b->leaf()) {
    l.push_back(b->data());
  }
  else {
    extract_leafs(b->right(),l);
    extract_leafs(b->left(),l);
  }
}

template<class Branch>
void extract_order(Branch* root,vector<int>& order)
{
  std::list<int> l;
  extract_leafs(root,l);
  order.resize(l.size());
  copy(l.begin(),l.end(),order.begin());
}

template<class Branch>
void extract_sizes(Branch* root,list<int>& sizes,bool jump=false)
{
  queue<Branch*> Q;
  Q.push(root);

  while (!Q.empty()) 
    {
      branch<int>* curr=Q.front();      
      if (!curr->leaf()) {
	branch<int> *r=curr->right(),*l=curr->left();
	if (!jump) sizes.push_back(curr->size());
	if (l->size() > r->size()) swap(r,l);
	Q.push(r);
	Q.push(l);	
      }
      Q.pop();
    }
}

#endif
