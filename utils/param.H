// 
//  Copyright (c) 2007, Pau Fernández
//

#ifndef _PARAM_H_
#define _PARAM_H_

#include <getopt.h>
#include <assert.h>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;
typedef unsigned int uint;

////////////////////////////////////////////////////////////////
// base class

class param
{
  char   _short;
  string _long;
public:
  typedef enum { OK,HELP,UNKNOWN_OPTION } condition;

  param(char s,string l):_short(s),_long(l) {}
  virtual ~param() {}
  char shortname() const { return _short; }
  string longname() const { return _long; }

  virtual bool has_args() const=0;
  virtual bool parse(string arg)=0;
  virtual string print() const=0;
};

////////////////////////////////////////////////////////////////
// template derived _param class
//
// T type requirements: operator>> + operator<<

template<class T>
class _param : public param
{
  T& _var;
public:
  _param(char s,string l,T& v)
    :param(s,l),_var(v) {}

  bool has_args() const { return true; }

  bool parse(string arg) {
    istringstream i(arg.c_str());
    i >> _var;
    return bool(i); // state ok ?
  }

  string print() const {
    ostringstream o;
    o << _var;
    string tmp(o.str());
    return tmp;
  }
};

/////////////////////////////////////////////////////////////////////
// Specialization for bool params (no need for option arguments)

template<>
class _param<bool> : public param {
  bool& _var;
public:
  _param(char s,string l,bool& v):param(s,l),_var(v) {}
  bool has_args() const { return false; }
  bool parse(string arg) { _var=true; return true;  }
  string print() const { return ""; } // no need to print in print_params(...)
};

/////////////////////////////////////////////////////////////////////
// Specialization for string params (assign the string directly)

template<>
class _param<string> : public param {
  string& _var;
public:
  _param(char s,string l,string& v):param(s,l),_var(v) {}
  bool has_args() const { return true; }
  bool parse(string arg) { _var=arg; return true; }
  string print() const { return _var; } 
};

////////////////////////////////////////////////////////////////////
// Parameter "factory"

template<class T>
inline _param<T>* make_param(char s,string l,T& v)
{ return new _param<T>(s,l,v); }

template<class T>
inline _param<T>* make_param(string l,T& v)
{ return new _param<T>(0,l,v); }


/////////////////////////////////////////////////////////////////////
// Helper functions

inline pair<param::condition,int> parse_params(vector<param*>& v,int argc,char** argv)
{
  typedef struct option toption;
    
  string short_options;
  map<int,int> _short_to_idx;
  toption* plong_options=new toption[v.size()+2];

  // parameters
  for (uint i=0;i<v.size();i++) {
    int c=v[i]->shortname();
    if (c!=0) {
      short_options+=c;
      if (v[i]->has_args()) short_options+=':';
    }
    else {
      c=-i-2;
    }
    _short_to_idx[c]=i;
    char* ptmp=new char[v[i]->longname().size()+1];
    strcpy(ptmp,v[i]->longname().c_str());
    toption optnew={ptmp,v[i]->has_args()?1:0,0,c};
    plong_options[i]=optnew;
  }

  // help option
  short_options+='h';
  toption help_option={"help",0,0,'h'};
  plong_options[v.size()]=help_option;
  _short_to_idx['h']=-1;

  // terminal
  toption lastopt={0,0,0,0};
  plong_options[v.size()+1]=lastopt;
  
  int idx;
  while (true) 
    {
      int c=getopt_long(argc,argv,
			short_options.c_str(),
			plong_options,
			&idx);
      if (c==-1) break;
      else {
	if (c!=0) {
	  map<int,int>::iterator it=_short_to_idx.find(c);
	  if (it == _short_to_idx.end()) {
	    return make_pair(param::UNKNOWN_OPTION,0);
	  }
	  idx=it->second;
	}
	else assert(false);

	if (idx < 0) {
	  if (idx == -1) return make_pair(param::HELP,0);
	  idx=-idx-2;
	}
	assert (uint(idx) < v.size());

	string stroptarg="";
	if (optarg) stroptarg=string(optarg);
	v[idx]->parse(stroptarg);
      }
    }
  return make_pair(param::OK,optind);
}

inline void print_params(ostream& o,const vector<param*>& v)
{
  o << endl << "options:" << endl;
  for (uint i=0;i<v.size();i++) {
    if (v[i]->shortname() != 0) 
      o << "  -" << v[i]->shortname() << ",";
    else
      o << "     ";
    o << " --" << v[i]->longname();
    if (v[i]->has_args()) 
      o << " [" << v[i]->print() << "]";
    o << endl;
  }
  o << "  -h, --help" << endl << endl;
}

inline void print_usage(vector<param*>& prms,string usage_first_line)
{
  cout << "usage: " << usage_first_line << endl;
  print_params(cout,prms);
}

inline void parse_params_ex(vector<param*>& prms,int argc,char** argv,
			    string usage_first_line,string progname,
			    vector<string>& args,uint minargs)
{
  pair<param::condition,int> res=parse_params(prms,argc,argv);

  switch (res.first) {
  case param::HELP: 
    cerr << usage_first_line << endl;
    print_params(cerr,prms);
    exit(0);
  case param::UNKNOWN_OPTION:
    cerr << "Try `" << progname << " --help' for more information" << endl;
    exit(0);
  case param::OK: {
    const uint numargs=argc-res.second;
    if (numargs < minargs) {
      cerr << "Too few arguments." << endl;
      cerr << "Try `" << progname << " --help' for more information" << endl;
      exit(0);
    }
    else {
      for (uint i=res.second;i<uint(argc);i++) {
	args.push_back(argv[i]);
      }
    }
    break;
  }
  default:
    assert(false);
  }

}

#endif
