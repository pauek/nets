// 
//  Copyright (c) 2007, Pau Fernández
//

#ifndef _COMMUNITY_H_
#define _COMMUNITY_H_

struct tree_tag
{
  string _name;
  int    _size;
  int    _num_edges;

  tree_tag(string n="",int sz=0)
    :_name(n),_size(sz),_num_edges(0) {}
};

ostream& operator<<(ostream& o,const tree_tag& t)
{ return o << t._size << ' ' << t._num_edges << ' ' << t._name; }

istream& operator>>(istream& i,tree_tag& t)
{ return i >> t._size >> t._num_edges >> t._name; }

#endif
