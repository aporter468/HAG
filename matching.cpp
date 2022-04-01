//for use by HAG algorithms
//input: list of aggregation nodes + their receivers
//output: score, i.e. size of best matching on the request graph



#include <string>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#include <unordered_map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <assert.h>
#include <algorithm>    // std::make_heap, std::pop_heap, std::push_heap, std::sort_heap
#include <iterator>     // std::front_inserter
#include <vector>
#include <chrono>
#include <iostream>
#include <list>
using namespace boost;
typedef adjacency_list<vecS,vecS,undirectedS> req_graph;
 inline size_t key(int i,int j) {return (size_t) i << 32 | (unsigned int) j;}
namespace std{
  template<>
  struct hash<std::pair<int, int> >
  {
    size_t operator()(const std::pair<int, int>& p) const
    {
      size_t res = 17;
      res = res * 31 + hash<int>()(p.first);
      res = res * 31 + hash<int>()(p.second);
      return res;
    }
  };
}

int getMaxMatch(std::set <int>& verts, std::set<std::pair<int,int> >& edges )
{
	int n =verts.size();
	std::cout<<"n="<<n<<std::endl;
	req_graph g(n);

	//node set: list out all requester,data pairs that are relevant and map them to a vertex number
	
	std::set<std::pair<int,int> >::iterator edgeIter;
	for (edgeIter = edges.begin();edgeIter!=edges.end(); edgeIter++)
	{
		std::cout<<"add edge: "<<edgeIter->first<<" "<<edgeIter->second<<std::endl;
		add_edge(edgeIter->first,edgeIter->second,g);
	}
 	//add edge: each pair of requests on same node for a (u,v) in aggnode set, verts in endstovert with same 2nd end and first ends are an aggregation node
  std::vector<graph_traits<req_graph>::vertex_descriptor> mat(n);
  bool success = checked_edmonds_maximum_cardinality_matching(g, &mat[0]);
//  std::cout << std::endl << "Found a matching of size " << matching_size(g, &mat[0]) << std::endl;

	
	return matching_size(g,&mat[0]);
}



