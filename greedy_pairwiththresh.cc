#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <unordered_map>
#include <set>
#include <assert.h>
#include <chrono>
#include <iostream>
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

struct PairCount {
  PairCount(int u, int v, int _cnt) : fst(u), snd(v), cnt(_cnt) {}
  int fst, snd, cnt;
};

struct pair_count_compare {
  bool operator()(const PairCount& lhs, const PairCount& rhs) const {
    if (lhs.cnt != rhs.cnt) return (lhs.cnt > rhs.cnt);
    if (lhs.fst != rhs.fst) return (lhs.fst < rhs.fst);
    return (lhs.snd < rhs.snd);
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

int main(int argc, char** argv)
{
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

	
	int u, v;
  int nv = 0;
  int ne = 0;
  std::map<int, std::set<int>* > inEdges;
  std::unordered_map<std::pair<int, int> , int> counter;
  std::set<PairCount, pair_count_compare> heap;
 bool readAsUndirected = true;

 
 std::map<int, std::set<int>* > outEdges;
  std::unordered_map<int,int> degmap;
  std::unordered_map<int,int> degmapout;
 
 
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
		degmap[v]=1;
	}
    else
	{
    	inEdges[v]->insert(u);
		degmap[v]++;
	}

 if(outEdges.find(u) == outEdges.end())
    {
        outEdges[u] = new std::set<int>();
        outEdges[u]->insert(v);
        degmapout[u]=  1;
    }
    else
    {
        if(outEdges[u]->find(v)==outEdges[u]->end()){
            outEdges[u] ->insert(v);

            degmapout[u] ++;
        }

    }


	//so do v,u as well
	if(readAsUndirected)
	{
		if (inEdges.find(u) == inEdges.end())
		{
			inEdges[u] = new std::set<int>();
			inEdges[u]->insert(v);
		}

		else
        {
            if(inEdges[u]->find(v)==inEdges[u]->end()){
                inEdges[u]->insert(v);

                degmap[u] ++;
            }

        }


        if(outEdges.find(v) == outEdges.end())
        {
            outEdges[v] = new std::set<int>();
            outEdges[v]->insert(u);
            degmapout[v]=  1;
        }
        else
        {
            if(outEdges[v]->find(u)==outEdges[v]->end()){
                outEdges[v] ->insert(u);

                degmapout[v] ++;
            }

        }

	}
  

  
  
  }
  fclose(file);
  printf("nv = %d\n", nv);

/*    for(int i =1; i<nv;i++)
  {
  		printf("vertex i: %d\n",i);
		std::set<int>::const_iterator init;
		std::set<int>::const_iterator outit;
		for (init = inEdges[i]->begin(); init!=inEdges[i]->end();init ++)
		{
			printf("in: %d\n",*init);
		}
}*/	
	int degthresh = 1;
	int usecount = 0;
	int ignorecount = 0;
	for (int i = 0; i< nv; i++)
	  if (inEdges.find(i) != inEdges.end()) {
      std::set<int>::const_iterator it1, it2;
      std::set<int>::const_iterator first = inEdges[i]->begin(), last = inEdges[i]->end();
      for (it1 = first; it1 != last; it1 ++)
        for (it2 = first; it2 != it1; it2 ++) {
          u = *it2;
          v = *it1;

		  assert(u < v);
		  if(degmapout[u]>degthresh && degmapout[v]>degthresh)
		  {
				usecount++;
		  	//printf("out degs: %d %d\n",degmapout[u],degmapout[v]);
			  if (counter.find(std::make_pair(u, v)) == counter.end())
				counter[std::make_pair(u, v)] = 1;
			  else
				counter[std::make_pair(u, v)] ++;
		  }
		  else{
			  ignorecount++;
		  }
        }
      if (i % 1000 == 0) printf("i = %d\n", i);
    }
printf("used pairs: %d ignored: %d\n",usecount,ignorecount);
  // initialize heap
  std::unordered_map<std::pair<int, int>, int>::const_iterator it;
  for (it = counter.begin(); it != counter.end(); it++) {
    PairCount pc(it->first.first, it->first.second, it->second);
	//printf("u: %d v: %d count: %d\n",it->first.first,it->first.second,it->second);
	//printf("pc: %d %d %d\n",pc.fst,pc.snd,pc.cnt);
	heap.insert(pc);
  }
 
auto start = std::chrono::steady_clock::now();


  int saved = 0;
  int scoreSaved=  0;
  std::map<int, int> depths;
  for (int i = 0;i<nsteps; i++) {
		if(heap.empty())
		{
			printf("no more heap elements!\n");
			break;
		}
	  PairCount pc = *heap.begin();
	   	
		
		if (pc.cnt < 2) 
		{
			printf("no more good heap elements\n");
			break;
		}


		
		int preDepth = 0;
		if (depths.find(pc.fst) != depths.end())
		  preDepth = std::max(preDepth, depths[pc.fst]);
		if (depths.find(pc.snd) != depths.end())
		  preDepth = std::max(preDepth, depths[pc.snd]);
		depths[nv] = preDepth + 1;
		saved += pc.cnt;

		printf("pc[%d]: fst(%d) snd(%d) depth(%d) cnt(%d) acc_save(%d) udeg(%d) vdeg(%d)\n", i, pc.fst, pc.snd, preDepth + 1, pc.cnt, saved,degmapout[pc.fst],degmapout[pc.snd]);
		printf("topdegs: %d %d\n",degmapout[pc.fst],degmapout[pc.snd]);
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
 	auto end = std::chrono::steady_clock::now();
 int mstime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        printf("runtime: %d\n",mstime);
	if(ne==0)
		ne =1;
		printf("scoreSaved: %d ratio: %f\n",scoreSaved, (float)scoreSaved/(float)ne);
}

