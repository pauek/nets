// 
//  Copyright (c) 2007, Pau Fernández
//

#ifndef _PROGRESS_H_
#define _PROGRESS_H_

#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>
using namespace std;

class progress_bar
{
  string _task_name;
  double _start,_end,_curr;
  ostream& _o;
  int _length;
  
  int _perc;

  void update() {
    double perc=(_curr-_start)/(_end-_start)*100.0;
    if (int(floor(perc))!=_perc) {
      _o << "\r[";
      for (int i=0;i<_length;i++) 
	_o << (double(i)/double(_length)*100.0 < floor(perc) ? '*' : ' ');
      _o << "] " << setw(3) << uint(floor(perc)) << "%" << flush;
    }
    _perc=int(floor(perc));
  }

public:
  progress_bar(double start,double end,
	       ostream& o=cout,int length=50)
    :_start(start),_end(end),_curr(start),
     _o(o),_length(length),_perc(-1) 
  { update(); }

  progress_bar(string task_name,double start,double end,
	       ostream& o=cout,int length=50)
    :_task_name(task_name),_start(start),_end(end),
     _curr(start),_o(o),_length(length),_perc(-1) 
  { _o << task_name << ":" << endl; update(); }

  ~progress_bar() { _o << endl; }

  void end() const { _o << endl; }

  void operator++(int) { _curr+=1.0; update(); }
  void operator++() { _curr+=1.0; update(); }
  void operator+=(double amount) { _curr+=amount; update(); }

  void incr(double amount=1.0) { _curr+=amount; update(); }
};

#endif
