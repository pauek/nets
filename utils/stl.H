// 
//  Copyright (c) 2007, Pau Fernández
//

#ifndef _STL_UTILS_H_
#define _STL_UTILS_H_

#include <assert.h>
#include <math.h>

#include <iterator>
#include <string>
#include <iosfwd>
#include <utility>

using namespace std;

template<class Container>
inline void print_all(const Container& c,std::ostream& o,char* sep)
{
  typename Container::const_iterator i=c.begin(),iend=c.end();
  if (i!=iend) {
    o << *i;
    for (++i;i!=iend;++i) o << sep << *i;
  }
}

template<class Container>
inline void print_all_endl(const Container& c,std::ostream& o,char* sep)
{
  print_all(c,o,sep);
  o << endl;
}

template<class Container>
inline double average(const Container& c)
{
  typename Container::value_type acum=0;
  typename Container::const_iterator i=c.begin(),iend=c.end();
  for (;i!=iend;i++)
    acum+=*i;
  return double(acum)/double(c.size());
}

template<class Container>
inline double entropy(const Container& c)
{
  assert(c.size() > 1);
  const uint SZ=c.size();
  
  // average
  typename Container::value_type acum=0;
  typename Container::const_iterator i=c.begin(),iend=c.end();
  for (;i!=iend;i++) acum+=*i;

  // entropy (probs = divide by sum)
  double E=0.0;
  for (i=c.begin();i!=iend;i++) 
    if (*i > 0) {
      double pr=double(*i)/double(acum);
      E-=pr*log(pr)/log(double(SZ)); // normalized
    }

  return E;
}

/* -> an STL equivalent exists: max_element(c) <- 
template<class Container>
inline typename Container::value_type maximum(const Container& c)
{
  typename Container::value_type max;
  typename Container::const_iterator i=c.begin(),iend=c.end();

  max=*i++;
  for (;i!=iend;i++) if (*i > max) max=*i;

  return max;
}
*/

template<class Container>
inline pair<double,double> variance(const Container& c)
{
  double acum=0.0;
  typename Container::const_iterator i=c.begin(),iend=c.end();
  double avg=average(c);
  for (;i!=iend;i++) acum+=pow(*i-avg,2.0);
  
  if (c.size()>1)
    return make_pair(avg,double(acum)/double(c.size()-1));
  else 
    return make_pair(avg,0.0);
}

template<class Container>
inline pair<double,double> standard_deviation(const Container& c)
{
  pair<double,double> res=variance(c);
  return make_pair(res.first,sqrt(res.second));
}

inline std::string pf_format(const std::string& str)
{
  std::string tmp = str.substr(0,str.find("("));
  return tmp.substr(tmp.find(" ")+1);
}

#define die(str) \
 { cerr << pf_format(__PRETTY_FUNCTION__) << ": " << str << endl , exit(0); }

////////////////////////////////////////////////////////////////////
// Print pairs

template<class First,class Second>
ostream& operator<<(ostream& o,const pair<First,Second>& p) 
{ return o << '[' << p.first << ',' << p.second << ']'; }

//////////////////////////////////////////////////////////////////////
// Comparison of subsequences

template<class Seq>
inline bool starts_with(const Seq& _seq,const Seq& _subseq)
{
  assert(_seq.size() >= _subseq.size());
  typename Seq::const_iterator si=_subseq.begin(),siend=_subseq.end();
  typename Seq::const_iterator mi=_seq.begin();
  for (;si!=siend;++si,++mi) if (*si != *mi) return false;
  return true;
}

template<class Seq>
inline bool ends_with(const Seq& _seq,const Seq& _subseq)
{
  assert(_seq.size() >= _subseq.size());
  typename Seq::const_reverse_iterator si=_subseq.rbegin(),siend=_subseq.rend();
  typename Seq::const_reverse_iterator mi=_seq.rbegin();
  for (;si!=siend;++si,++mi) if (*si != *mi) return false;
  return true;
}

#endif
