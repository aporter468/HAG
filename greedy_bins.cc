#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <unordered_map>
#include <set>
#include <assert.h>
#include <algorithm>    // std::make_heap, std::pop_heap, std::push_heap, std::sort_heap
#include <iterator>     // std::front_inserter
#include <vector>
#include <chrono>
#include <iostream>
#include <list>
#include "matching.h"


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

using namespace std;
struct NborSet {
  NborSet(int v, int _cnt) : rec(v), cnt(_cnt) {}
  int rec, cnt;
};


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

//higher cnt better, otherwise use the node id rank
struct nbor_count_compare {
  bool operator()(const NborSet& lhs, const NborSet& rhs) const {
    if (lhs.cnt != rhs.cnt) return (lhs.cnt > rhs.cnt);
    return (lhs.rec < rhs.rec);
  }
};

//done for everyone up front  1x
void add_nbor_count( int v,
                    std::unordered_map<int,int>& degmap,
                    std::set<NborSet, nbor_count_compare>& heap)
{
	
	int countval = degmap[v];
	NborSet ns(v,countval);
	heap.insert(ns);

}

void printItset(int k, std::vector<std::unordered_map<std::pair<int,int> , int>::iterator>& itset)
{
	std::cout<<"itset: ";
	for (int i =0; i<k; i++)
		{
				std::cout<<itset[i]->first.first<<","<<itset[i]->first.second<<" ";
		}
	std::cout<<endl;

}


int makeAndScoreGraph(int k, 
        std::vector<std::unordered_map<std::pair<int,int> , int>::iterator>& itset, 
        std::map<std::pair<int,int>, std::set<int>* >& sharedOuts)
 {	
	std::unordered_map<pair<int,int>, int> endsToVert;
	set<pair<int,int> > edges;
	set<int> verts;
	unordered_map<int, pair<int,int> > reverseVertMap;
	
	printItset(k,itset);
	int numnodes = 0;
	for (int i =0;i<k; i++)
    {
        std::cout<<"make end set for i="<<i<<std::endl;
        int first = itset[i]->first.first;
        int second = itset[i]->first.second;
		std::cout<<"i="<<i<<" first: "<<first<<" second :"<<second<<std::endl;
        int count = itset[i]->second;
		std::cout<<"count: "<<count<<endl;

		std::pair<int,int> mypair = std::make_pair(first,second);
		
		cout<<"pair: "<<mypair.first<<" "<<mypair.second<<endl;
		if(sharedOuts.count(std::make_pair(first,second))==0)
		{
			cout<<"key pair deosn't exist"<<endl;
		}
		else
		{
			cout<<"found key pair"<<endl;
		}
       // std::set<int> outset2 =  *(sharedOuts[std::make_pair(first,second)]);
        //cout<<"got regular outset."<<endl;
		std::set<int> *outset = sharedOuts[std::make_pair(first,second)];
		std::cout<<"outset size: "<<outset->size()<<" count: "<<count<<std::endl;
	
		
		set<int>::const_iterator  tempIter =outset->begin();
		while(tempIter!=outset->end())
		{
			cout<<"temp iter step: "<<*tempIter<<endl;
			tempIter++;
		}
		
		
		set<int>::const_iterator outsetIter = outset->begin();
       while(outsetIter!=outset->end())
		{
			int recnode = *outsetIter;
			cout<<"recnode: "<<recnode<<endl;
			if(endsToVert.find(make_pair(first,recnode))==endsToVert.end())
			{
				endsToVert[make_pair(first,recnode)] = numnodes;
				verts.insert(numnodes);	
				reverseVertMap[numnodes] = make_pair(first,recnode);
				numnodes++;
			
			}

			if(endsToVert.find(make_pair(second,recnode))==endsToVert.end())
			{
				endsToVert[make_pair(second,recnode)] = numnodes;
				verts.insert(numnodes);
				reverseVertMap[numnodes] =make_pair(second,recnode);
				numnodes++;
			}
			
			cout<<"numnodes = "<<numnodes<<endl;			
			
			int v1 = endsToVert[make_pair(first,recnode)];
			int v2 = endsToVert[make_pair(second,recnode)];
			edges.insert(make_pair(v1,v2));
			cout<<"edges inserted"<<endl;
			cout<<"iter =end?"<<(outsetIter ==outset->end())<<endl;
			outsetIter++;
			cout<<"next iter"<<endl;
        }
		cout<<"end of outsetiter test"<<endl;
    }
	cout<<"sending "<<verts.size()<<" verts with "<<edges.size()<<" edges."<<endl;
    int matchSize=	getMaxMatch(verts,edges);
	cout<<"match size: "<<matchSize<<endl;
	return matchSize;
}
int recursiveSearch(
		int iterindex,
		int k,
        std::vector<std::unordered_map<std::pair<int,int> , int>::iterator>& itset,
		std::map<std::pair<int,int>, std::set<int>* >& sharedOuts,
	  std::unordered_map<std::pair<int, int> , int>& counter)
{
	/*
	if(iterindex ==k)//this recursion level didn't set an additional iter at all
	{	

		 int score = makeAndScoreGraph(k,itset,sharedOuts);
		 return score;//base case always returns value, nothing to max over...
	}*/
	//else
	//{
		int maxscore = 0;
		while(itset[iterindex]!=counter.end())
		{

			
			
			//---------
			//
				/*cout<<"recursive search for iter: "<<iterindex<<endl;

				cout<<"base case printout test:"<<endl;
				int u =5;// itset[i]->first.first;
				int v =8;// itset[i]->first.second;
				cout<<"u,v: "<<u<<" "<<v<<endl;
				std::set<int> *outset = sharedOuts[std::make_pair(u,v)];
				std::cout<<"outset size: "<<outset->size()<<std::endl;


				 set<int>::const_iterator  tempIter =outset->begin();
				while(tempIter!=outset->end())
				{
					cout<<"temp iter step: "<<*tempIter<<endl;
					tempIter++;
				}*/	
//---------------------

			


		int subscore = 0;
			if(iterindex==k-1)
			{
				subscore = makeAndScoreGraph(k,itset,sharedOuts);
			}
			else
			{

				int nextIterLoc = iterindex+1;
				itset[nextIterLoc] = itset[iterindex];
				itset[nextIterLoc]++;
			
				subscore = recursiveSearch(iterindex+1,k,itset,sharedOuts,counter);
			}
			if (subscore >maxscore)
			{
				maxscore = subscore;
			}
			itset[iterindex]++;
		}
		cout<<"max at depth "<<iterindex<<" is "<<maxscore<<endl;
		return maxscore;
	//}



}

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
  heap.erase(pc);
  counter[std::make_pair(u, v)] = oldVal - 1;
  pc.cnt = oldVal - 1;
  heap.insert(pc);
}


void sub_nbor_count(int v, std::unordered_map<int,int>& degmap,
						std::set<NborSet,nbor_count_compare>& heap,int newval)
{
	int oldval = degmap[v];
	NborSet ins(v,oldval);
	heap.erase(ins);
	ins.cnt = newval;
	heap.insert(ins);
}
int main(int argc, char** argv)
{
  //FILE* file = fopen("IMDB-MULTI/IMDB-MULTI_A.txt", "r");
 //FILE* file = fopen("REDDITBINARY/REDDITBINARY_A.txt", "r");
//FILE* file = fopen("PROTEINS/PROTEINS_A.txt", "r");
	//FILE* file = fopen("COLLAB/COLLAB_A.txt", "r");
  //FILE* file = fopen("BZR_MD/BZR_MD_A.txt", "r");
//FILE* file = fopen("amazon0302.txt","r");
//FILE* file = fopen("smalltest.txt","r");
// FILE* file = fopen("../syntheticGraphs/testset1/graph_1.txt","r");
FILE* file = fopen(argv[1],"r");	
int nsteps = atoi(argv[2]);
int readundr = atoi(argv[3]);
	
	
 	bool readAsUndirected = false;
	if(readundr==1)
			readAsUndirected = true;





int u, v;
  int nv = 0;
  int cumulativeScore=0;
  std::map<int, std::set<int>* > inEdges;
  
  std::map<int, std::set<int>* > outEdges;
  std::unordered_map<std::pair<int, int> , int> counter;
  std::map<std::pair<int,int>,std::set<int>* > sharedOuts;
  std::unordered_map<int,int> degmap;
  
  std::unordered_map<int,int> degmapout;
  std::set<NborSet, nbor_count_compare> degheap;
  std::set<PairCount, pair_count_compare> heap;

  while (fscanf(file, "%d %d", &u, &v) != EOF) {
    if (std::max(u, v) >= nv)
      nv = std::max(u, v) + 1;
	printf("edge: %d %d\n",u,v);
	if (inEdges.find(v) == inEdges.end())
 	{
		printf("add inedges for v=%d\n",v);
		inEdges[v] = new std::set<int>();
		inEdges[v]->insert(u);
		degmap[v] = 1;//rely on inEdges to match degmap on index sets
	}
    else
	{
		if(inEdges[v]->find(u)==inEdges[v]->end()){
    	inEdges[v]->insert(u);
		degmap[v] ++;
		}
	}


	if(outEdges.find(u) == outEdges.end())
	{
		printf("add outedges for u=%d\n",u);
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
 	
	if(readAsUndirected)
	{
		if (inEdges.find(u) == inEdges.end())
		{

		printf("add inedges for u=%d\n",u);
			inEdges[u] = new std::set<int>();
			inEdges[u]->insert(v);
		}
		else
		{
			if(inEdges[u]->find(v)==inEdges[u]->end()){
				inEdges[u]->insert(v);
			}

		}
		if(outEdges.find(v) == outEdges.end())
		{
		
		printf("add outedges for v=%d\n",v);
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
	printf("after edge add: uin: %d uout %d vin %d vout %d\n",inEdges[u]->size(),outEdges[u]->size(),inEdges[v]->size(),outEdges[v]->size());
  



  }

  fclose(file);
  printf("nv = %d\n", nv);
  
  
  for (int i = 0; i< nv; i++)
	  if (inEdges.find(i) != inEdges.end()) {
      std::set<int>::const_iterator it1, it2;
      std::set<int>::const_iterator first = inEdges[i]->begin(), last = inEdges[i]->end();
	  printf("in-pairs for i=%d set size %d\n",i,inEdges[i]->size());
	  for (it1 = first; it1 != last; it1 ++)
        for (it2 = first; it2 != it1; it2 ++) {
          u = *it2;
          v = *it1;
		  printf("checking pair (u,v)=(%d,%d)\n",u,v);
          assert(u < v);
          if (counter.find(std::make_pair(u, v)) == counter.end())
		  {
			  counter[std::make_pair(u, v)] = 1;
		  }
          else
		  {
            counter[std::make_pair(u, v)] ++;
		  }

		  if(sharedOuts.find(std::make_pair(u,v)) == sharedOuts.end())
		  {
//			  cout<<"adding out-set for "<<u<<","<<v<<endl;
			   sharedOuts[std::make_pair(u,v)] = new std::set<int>();
		  }
		    cout<<i<<" in outset for "<<u<<","<<v<<endl;
			sharedOuts[std::make_pair(u,v)]->insert(i);
		   	cout<<"check outset:"<<endl;


        }
      if (i % 1000 == 0) printf("i = %d\n", i);
    }

  // initialize heap
  std::unordered_map<std::pair<int, int>, int>::const_iterator it;
 
  for (it = counter.begin(); it != counter.end(); it++) {
    PairCount pc(it->first.first, it->first.second, it->second);
    heap.insert(pc);
  }

  std::unordered_map<int,int>::const_iterator  it2;
  for (it2 = degmap.begin();it2!=degmap.end();it2++)
  {
	  NborSet ins(it2->first,it2->second);
	  degheap.insert(ins);
  }
  
  auto start = std::chrono::steady_clock::now();
  
  //loop over all sets of k pairs - use "counter"

const int k = nsteps;
std::cout<<"k="<<nsteps<<std::endl;
//std::unordered_map<std::pair<int, int> , int>::iterator pcsetItset = new std::unordered_map<std::pair<int, int> , int>::iterator[k];

std::vector<std::unordered_map<std::pair<int,int> , int>::iterator > itset(k);
for (int i =0; i<k; i++)
{
	std::unordered_map<std::pair<int,int> , int>::iterator startiter = counter.begin();
	itset[i] = startiter;
}


int matchingscore=  recursiveSearch(0,k,itset,sharedOuts,counter);
printf("best matching: %d\n",matchingscore);

int costsavedbymatch = matchingscore - k;
if(matchingscore==0 && k >1)
{
	printf("score of 0 for k=%d, trying k=1:\n");
	for (int i =0; i<k; i++)
	{
		std::unordered_map<std::pair<int,int> , int>::iterator startiter = counter.begin();
		itset[i] = startiter;
	}
	
	matchingscore=  recursiveSearch(0,1,itset,sharedOuts,counter);
	costsavedbymatch = matchingscore - 1;
	printf("new score: %d saved: %d\n",matchingscore,costsavedbymatch);
}

	if (costsavedbymatch<0)
		{
			costsavedbymatch  =0;
		}

	printf("cost saved by match: %d t\n",costsavedbymatch);
	auto end = std::chrono::steady_clock::now();
    int mstime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("runtime: %d\n",mstime);
//delete [] pcsetItset;
}

