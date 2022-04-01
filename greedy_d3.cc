#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <unordered_map>
#include <set>
#include <assert.h>
#include <chrono>
#include <iostream>
#include <tuple>
#include <utility>



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


namespace N {
	typedef std::tuple<int, char, char> key_t;
	 
	struct key_hash : public std::unary_function<key_t, std::size_t>
	{
		std::size_t operator()(const key_t& k) const
		{
			return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k);
		}
	};
	 
	struct key_equal : public std::binary_function<key_t, key_t, bool>
	{
		bool operator()(const key_t& v0, const key_t& v1) const
		{
			return (
			std::get<0>(v0) == std::get<0>(v1) &&
			std::get<1>(v0) == std::get<1>(v1) &&
			std::get<2>(v0) == std::get<2>(v1)
			);
		}
	};


	typedef std::unordered_map<const key_t,int,key_hash,key_equal> map_t;
	
}
struct PairCount {
  PairCount(int u, int v, int _cnt) : fst(u), snd(v), cnt(_cnt) {}
  int fst, snd, cnt;
};

struct TripCount {
  TripCount(int u, int v,int w, int _cnt) : fst(u), snd(v),thr(w), cnt(_cnt) {}
  int fst, snd, thr, cnt;
};


struct pair_count_compare {
  bool operator()(const PairCount& lhs, const PairCount& rhs) const {
    if (lhs.cnt != rhs.cnt) return (lhs.cnt > rhs.cnt);
    if (lhs.fst != rhs.fst) return (lhs.fst < rhs.fst);
    return (lhs.snd < rhs.snd);
  }
};

struct trip_count_compare {
  bool operator()(const TripCount& lhs, const TripCount& rhs) const {
    if (lhs.cnt != rhs.cnt) return (lhs.cnt > rhs.cnt);
    if (lhs.fst != rhs.fst) return (lhs.fst < rhs.fst);
	if (lhs.snd != rhs.snd) return (lhs.snd < rhs.snd);
    return (lhs.thr < rhs.thr);
  }
};


void add_pair_count(int u, int v,
                    std::unordered_map<std::pair<int, int>, int>& counter,
                    std::set<PairCount, pair_count_compare>& heap)
{
  if (u > v) {int w = u; u = v; v = w;}
  if (counter.find(std::make_pair(u, v)) == counter.end()) {
    counter[std::make_pair(u, v)] = 1;
    PairCount pc(u, v, 1);
    heap.insert(pc);
  } else {
    int oldVal = counter[std::make_pair(u, v)];
    PairCount pc(u, v, oldVal);
    heap.erase(pc);
    counter[std::make_pair(u, v)] = oldVal + 1;
    pc.cnt = oldVal + 1;
    heap.insert(pc);
  }
}

void sub_pair_count(int u, int v,
                    std::unordered_map<std::pair<int, int>, int>& counter,
                    std::set<PairCount, pair_count_compare>& heap)
{
  if (u > v) {int w = u; u = v; v = w;}
  int oldVal = counter[std::make_pair(u, v)];
  PairCount pc(u, v, oldVal);
  printf("new pc: %d %d old val of %d\n",u,v,oldVal);
  heap.erase(pc);
  counter[std::make_pair(u, v)] = oldVal - 1;
  pc.cnt = oldVal - 1;
  heap.insert(pc);
}


void sub_trip_count(int u, int v,int w,
                    N::map_t& counter,
                    std::set<TripCount, trip_count_compare>& heap)
{
  int u2=u;
  int v2 = v;
  int w2 = w;
  printf("order into sub: %d %d %d\n",u,v,w);
 if(u<v && u<w && w<v)
 {
	 v2 = w;
	 w2 = v;
 }
 if(w<u && u<v && w<v)
 {
	 u2 = w;
	 v2 = u;
	 w2 = v;
 }
 if(w<v && w<u && v<u)
 {
	 u2 = w; v2 = v; w2 = u;
 }
 if(v<u && u< w&& v<w )
 {
	 u2 = v;
	 v2 = u;
	 w2 = w;
 }
 if(v< w&& w<u && v<u)
 {
	 u2 = v;
	 v2=w;
	 w2 =u;
 }
 u = u2;
 v = v2;
 w = w2;


  int oldVal = counter[std::make_tuple(u, v,w)];
  TripCount tc(u, v, w,oldVal);
  printf("new tc: %d %d %d old val of %d\n",u,v,w,oldVal);
  heap.erase(tc);
  counter[std::make_tuple(u, v,w)] = oldVal - 1;
  tc.cnt = oldVal - 1;
  heap.insert(tc);
}



int main(int argc, char** argv)
{
	using namespace N;
  //FILE* file = fopen("IMDB-MULTI/IMDB-MULTI_A.txt", "r");
 //FILE* file = fopen("REDDITBINARY/REDDITBINARY_A.txt", "r");
//FILE* file = fopen("PROTEINS/PROTEINS_A.txt", "r");
	//FILE* file = fopen("COLLAB/COLLAB_A.txt", "r");
  //FILE* file = fopen("BZR_MD/BZR_MD_A.txt", "r");
//FILE* file = fopen("amazon0302.txt","r");
// FILE* file = fopen("smalltest.txt","r");
// //FILE* file = fopen("../syntheticGraphs/testset1/graph_1.txt","r");
	FILE* file = fopen(argv[1],"r");	
	int nsteps = atoi(argv[2]);

	
	int u, v,w;
  int nv = 0;
   int ne =0;
  std::map<int, std::set<int>* > inEdges;
  bool readAsUndirected = true; 
  
  map_t counter_t;
  
   std::set<TripCount, trip_count_compare> heap_t;
  printf("start file read\n");

  while (fscanf(file, "%d %d", &u, &v) != EOF) {
	  ne++;
    if (std::max(u, v) >= nv)
      nv = std::max(u, v) + 1;
//	printf("edge: %d %d\n",u,v);
	if (inEdges.find(v) == inEdges.end())
 	{
//		printf("add edge list for node: %d\n",v);
		inEdges[v] = new std::set<int>();

    	inEdges[v]->insert(u);
	}
    else
	{
    	inEdges[v]->insert(u);
	}

  
  	if(readAsUndirected)
	{
		if (inEdges.find(u) == inEdges.end())
		{
			inEdges[u] = new std::set<int>();
			inEdges[u]->insert(v);
		}
		else
		{
			if(inEdges[u]->find(v)==inEdges[u]->end()){inEdges[u]->insert(v);}
		}

	}

  
  
  }
  fclose(file);
  printf("nv = %d\n", nv);
 
  
//Trip-----

  for (int i = 0; i< nv; i++)
  {
	  printf("neighbors of %d\n",i);
	  if (inEdges.find(i) != inEdges.end()) {
      std::set<int>::const_iterator it1, it2,it3;
      std::set<int>::const_iterator first = inEdges[i]->begin(), last = inEdges[i]->end();
      for (it1 = first; it1 != last; it1 ++)
        for (it2 = first; it2 != it1; it2 ++) {
			for (it3 = first; it3!=it2; it3++){

			  u = *it3;
			  v = *it2;
			  w = *it1;
//			  printf("add triple %d %d %d\n",u,v,w);			  
			  assert(u < v);
				assert (v< w);
			  if (counter_t.find(std::make_tuple(u, v,w)) == counter_t.end())
				counter_t[std::make_tuple(u, v,w)] = 1;
			  else
				counter_t[std::make_tuple(u, v,w)] ++;
			}
        }
      if (i % 1000 == 0) printf("i = %d\n", i);
    }
  }



  // initialize heap

  //PC---
/*  std::unordered_map<std::pair<int, int>, int>::const_iterator it;
  for (it = counter.begin(); it != counter.end(); it++) {
    PairCount pc(it->first.first, it->first.second, it->second);
	printf("u: %d v: %d count: %d\n",it->first.first,it->first.second,it->second);
	printf("pc: %d %d %d\n",pc.fst,pc.snd,pc.cnt);
	heap.insert(pc);
  }
*/

//Trip----
  map_t::const_iterator it_t;
  for (it_t = counter_t.begin(); it_t != counter_t.end(); it_t++) {
	  N::key_t keyval = it_t->first;
	  TripCount tc(std::get<0>(keyval),std::get<1>(keyval),std::get<2>(keyval),it_t->second);
	  if(it_t->second>1){

	  	//printf("u: %d v: %d w: %d count: %d\n",std::get<0>(keyval),std::get<1>(keyval),std::get<2>(keyval),it_t->second);
	  }
		printf("tc: %d %d %d %d\n",tc.fst,tc.snd,tc.thr,tc.cnt);

		heap_t.insert(tc);
  }
 



auto start = std::chrono::steady_clock::now();

//pc---------
/*  int saved = 0;
  int scoreSaved=  0;
  std::map<int, int> depths;
  for (int i = 0;i<nsteps; i++) {
	  PairCount pc = *heap.begin();
		int preDepth = 0;
		if (depths.find(pc.fst) != depths.end())
		  preDepth = std::max(preDepth, depths[pc.fst]);
		if (depths.find(pc.snd) != depths.end())
		  preDepth = std::max(preDepth, depths[pc.snd]);
		depths[nv] = preDepth + 1;
		saved += pc.cnt;

		printf("pc[%d]: fst(%d) snd(%d) depth(%d) cnt(%d) acc_save(%d)\n", i, pc.fst, pc.snd, preDepth + 1, pc.cnt, saved);
		scoreSaved += (pc.cnt-1);//account for 
		printf("cumulativesave: %d\n",scoreSaved);
	  
		heap.erase(heap.begin());
		for (int j = 0; j < nv; j++)
			if (inEdges.find(j) != inEdges.end()) {
			std::set<int>* list = inEdges[j];
			if ((list->find(pc.fst) != list->end())
			&&  (list->find(pc.snd) != list->end())) {
				  list->erase(pc.fst);
				  list->erase(pc.snd);
					printf("j: %d\n",j);
				  // update counters
				//in addition to removing the pair u,v we reduce counts for any other pairs into j that include u or v
			  std::set<int>::const_iterator it;
			  for (it = list->begin(); it != list->end(); it++) {
				printf("sub pair counts: %d,%d and %d,%d\n",*it,pc.fst,*it,pc.snd);
				 sub_pair_count(*it, pc.fst, counter, heap);
				sub_pair_count(*it, pc.snd, counter, heap);
		//		add_pair_count(*it, nv, counter, heap);
			  }
			  list->insert(nv);
			}
		  }
		nv ++;
  }
*/
  //tc------
  //
 
  int savedt = 0;
  int scoreSavedt=  0;
  std::map<int, int> depthst;
  for (int i = 0;i<nsteps; i++) {
	  TripCount tc = *heap_t.begin();
 		
	    if (tc.cnt < 2)  
        {
			printf("no useful cliques\n");
            break;
        }
	  
 
 
 		int preDepth = 0;
		if (depthst.find(tc.fst) != depthst.end())
		  preDepth = std::max(preDepth, depthst[tc.fst]);
		if (depthst.find(tc.snd) != depthst.end())
		  preDepth = std::max(preDepth, depthst[tc.snd]);
		depthst[nv] = preDepth + 1;
		savedt += tc.cnt;

		printf("tc[%d]: fst(%d) snd(%d) thr(%d) depth(%d) cnt(%d) acc_save(%d)\n", i, tc.fst, tc.snd,tc.thr, preDepth + 1, tc.cnt, savedt);
		scoreSavedt += 2*(tc.cnt-1);//account for 
		printf("cumulativesave: %d\n",scoreSavedt);
	  
		heap_t.erase(heap_t.begin());
		for (int j = 0; j < nv; j++)
			if (inEdges.find(j) != inEdges.end()) {
			std::set<int>* list = inEdges[j];
			if ((list->find(tc.fst) != list->end())
			&&  (list->find(tc.snd) != list->end())
			&& (list->find(tc.thr)  != list->end())) {
				  list->erase(tc.fst);
				  list->erase(tc.snd);
				  list->erase(tc.thr);
					printf("j: %d\n",j);
				  // update counters
				//in addition to removing the pair u,v we reduce counts for any other pairs into j that include u or v
			  std::set<int>::const_iterator it;
			  int u = tc.fst;
			  int v  = tc.snd;
			  int w = tc.thr;
			  for (it = list->begin(); it != list->end(); it++) {
						 if(*it != u && *it !=v && *it !=w)
						 {
						 	printf("remove with u,v,w pairs: %d\n",*it);
							sub_trip_count(*it,u,v,counter_t,heap_t);
							sub_trip_count(*it,u,w,counter_t,heap_t);
							sub_trip_count(*it,v,w,counter_t,heap_t);
						 }
			  }
			  std::set<int>::const_iterator it1,it2;
			  for(it1 = list->begin();it1!=list->end();it1++)
			  {
				  for(it2 =list->begin();it2!=it1;it2++)
				  {
					  printf("pair to remove u,v,w with: %d %d\n",*it1,*it2);
					  sub_trip_count(*it1,*it2,u,counter_t,heap_t);
					  sub_trip_count(*it1,*it2,v,counter_t,heap_t);
					  sub_trip_count(*it1,*it2,w,counter_t,heap_t);
				  }
			  }
			  list->insert(nv);
			}
		  }
		nv ++;
  }

 	auto end = std::chrono::steady_clock::now();
 int mstime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        printf("runtime: %d\n",mstime);
		printf("scoreSaved: %d ratio:  %f\n",scoreSavedt,(float)scoreSavedt/(float)ne);
}

