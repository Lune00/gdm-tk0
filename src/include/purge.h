// purge.h, from Thinking in C++, Bruce Eckel

#ifndef purge_h
#define purge_h
#include <algorithm>

template<class Seq> 
void purge(Seq& c) 
{
  typename Seq::iterator i;
  for(i = c.begin(); i != c.end(); i++) 
    {
      delete *i;
      *i = NULL;
    }
}

// Iterator version:
template<class InpIt>
void purge(InpIt begin, InpIt end) 
{
  while(begin != end) 
    {
      delete *begin;
      *begin = NULL;
      begin++;
    }
}
#endif // purge_h
