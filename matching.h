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


int getMaxMatch(std::set<int>& verts,std::set<std::pair<int,int> >& edges);
