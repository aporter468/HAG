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
	printf("read undirected? %d\n",readAsUndirected);


  auto start = std::chrono::steady_clock::now();
int u, v;
  int nv = 0;
  int cumulativeScore=0;
  std::map<int, std::set<int>* > inEdges;
  
  std::map<int, std::set<int>* > outEdges;
  std::unordered_map<std::pair<int, int> , int> counter;
  std::unordered_map<int,int> degmap;
  std::unordered_map<int,int> degmapout;
  std::set<NborSet, nbor_count_compare> degheap;
  std::set<PairCount, pair_count_compare> heap;
  int nedges = 0;
  int inedgecount = 0;
  while (fscanf(file, "%d %d", &u, &v) != EOF) {
	  nedges++;
    if (std::max(u, v) >= nv)
      nv = std::max(u, v) + 1;
//	printf("edge: %d %d\n",u,v);
	if (inEdges.find(v) == inEdges.end())
 	{
//		printf("add edge list for node: %d\n",v);
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
			inedgecount++;
			degmap[u] = 1;//rely on inEdges to match degmap on index sets
		}
		else
		{
			if(inEdges[u]->find(v)==inEdges[u]->end()){
				inEdges[u]->insert(v);
				inedgecount++;
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
	printf("nedges: %d inedgecount: %d\n",nedges,inedgecount);
  fclose(file);
/*  printf("nvert = %d\n", nv);
  for(int i =1; i<nv;i++)
  {
  		printf("vertex i: %d\n",i);
		std::set<int>::const_iterator init;
		std::set<int>::const_iterator outit;
		for (init = inEdges[i]->begin(); init!=inEdges[i]->end();init ++)
		{
			printf("in: %d\n",*init);
		}
		for(outit=outEdges[i]->begin(); outit!=outEdges[i]->end();outit ++)
		{
			printf("out: %d\n",*outit);
		}
  }*/
  
/*  for (int i = 0; i< nv; i++)
	  if (inEdges.find(i) != inEdges.end()) {
      std::set<int>::const_iterator it1, it2;
      std::set<int>::const_iterator first = inEdges[i]->begin(), last = inEdges[i]->end();
      for (it1 = first; it1 != last; it1 ++)
        for (it2 = first; it2 != it1; it2 ++) {
          u = *it2;
          v = *it1;
          assert(u < v);
          if (counter.find(std::make_pair(u, v)) == counter.end())
            counter[std::make_pair(u, v)] = 1;
          else
            counter[std::make_pair(u, v)] ++;
        }
      if (i % 1000 == 0) printf("i = %d\n", i);
    }

  // initialize heap
  std::unordered_map<std::pair<int, int>, int>::const_iterator it;
 
  for (it = counter.begin(); it != counter.end(); it++) {
    PairCount pc(it->first.first, it->first.second, it->second);
    heap.insert(pc);
  }
*/
  std::unordered_map<int,int>::const_iterator  it2;
  for (it2 = degmap.begin();it2!=degmap.end();it2++)
  {
	  NborSet ins(it2->first,it2->second);
	  degheap.insert(ins);
  }
  
//  auto start = std::chrono::steady_clock::now();


  int saved = 0;
  std::map<int, int> depths;
 // for (int i = 0; i<nsteps; i++) {
  int i =0;
  int vertsused = 0;
  while(!degheap.empty() && vertsused<nsteps) 
  {
	
	  NborSet topins =  *degheap.begin();
		printf("i=%d top node: %d with deg: %d\n",i, topins.rec,topins.cnt);
		if (topins.cnt<2)
		{
			vertsused  =nsteps; 
			break;
		}
	    sub_nbor_count(topins.rec,degmap,degheap,0);//put at bottom of heap but deg isn't actually 0
		int v = topins.rec;
		int vdeg = topins.cnt;
		std::set<int>::const_iterator edgeiter;
		std::set<int>::const_iterator  first  = inEdges[v]->begin(),last = inEdges[v]->end();
				
  		std::set<NborSet, nbor_count_compare> outdegheap;
		//get the out-neighborhoods of topNode's in-neighbors
		for(edgeiter = first; edgeiter!=last;edgeiter++)
		{
			u = *edgeiter;
			int udeg = degmapout[u];
			NborSet ins(u,udeg);
		  	outdegheap.insert(ins);
		//	printf("v: %d u: %d u deg: %d\n",v,u,udeg);
		
		}

		//maximize the overlap-setsize function

		std::set<int>* sharedOuts;
		sharedOuts = new std::set<int>();
		std::set<int>* prevSharedOuts;
		prevSharedOuts = new std::set<int>();
		int prevscore =  0;
		std::vector<int>* usedSenders;
		usedSenders = new std::vector<int>();
		int senderCount = 0;
		for (int j = 0; j<vdeg; j++)
		{
			NborSet topsend = *outdegheap.begin();
			printf("top sender: %d with  %d\n",topsend.rec,topsend.cnt);
			if(topsend.cnt==0)
			{
				j = vdeg+1;
			}
			else
			{
				topsend.cnt = 0;
				sub_nbor_count(topsend.rec,degmapout,outdegheap,0);
				if(j==0)
				{
				//	sharedOuts = outEdges[topsend.rec];
					std::copy(outEdges[topsend.rec]->begin(),outEdges[topsend.rec]->end(),std::inserter(*sharedOuts,sharedOuts->begin()));
				printf("shareOuts start:  %d\n",sharedOuts->size());
						for (auto it = sharedOuts->begin(); it!=sharedOuts->end(); ++it)
						{
							printf("%d ",*it);

						}
					printf("\n"); 
				}
				else
				{
					*prevSharedOuts = *sharedOuts;
					std::set<int>* newset = outEdges[topsend.rec];
					std::set<int> newintersect;// = new std::set<int>();

					std::set_intersection(sharedOuts->begin(),sharedOuts->end(),newset->begin(),newset->end(),std::inserter(newintersect,newintersect.begin()));

	//				printf("new receiver set: size = %d\n",newintersect.size());
				/*    for (auto it=newintersect.begin(); it != newintersect.end(); ++it)
					{
						printf("%d ",*it);
					}
					printf("\n");*/

					*sharedOuts = newintersect;
					printf("new shared outs size: %d\n",sharedOuts ->size());

				}

				int score = (j+1)*sharedOuts->size() - (j+1+sharedOuts->size()-1);
				
				printf("score of %d from %d receiving from %d senders\n",score,sharedOuts->size(),j+1);
			
			if(score<prevscore)
				{
					//use previous setup and break
					printf("score of %d was better with %d outs.\n",prevscore,prevSharedOuts->size());
					cumulativeScore+=prevscore;
		//			printf("cumulativescore: %d\n",cumulativeScore);
					j=vdeg+1;
				}
				else if(score==0 && j>0)//first will always be 0 but if doesn't go up it can't go down to be caught by previous if
				{
					senderCount = 0;
					j = vdeg+1;
				}

				else if (j==(vdeg-1))//last iter so use score
				{
					printf("use last score of %d\n",score);
					cumulativeScore+=score;
					prevscore=score;
					usedSenders->push_back(topsend.rec);
					senderCount++;
					*prevSharedOuts  = *sharedOuts;
				}
				else
				{
				
					prevscore = score;
					usedSenders->push_back(topsend.rec);
					senderCount++;
					printf("finished updates to keep score as option\n");				
				}
			  }
		}


  //if it hasn't changed still need a printout

	   
		if(senderCount>0)
		{
	
	  	printf("cumulativescore: %d\n",cumulativeScore);
			vertsused++;
			printf("best sender set: %d send to %d for score %d\n",senderCount,prevSharedOuts->size(),prevscore);
			printf("edgessaved: %d\n",prevscore);
				//leave outdeg heap alone b/c its reset at the outer loop anyway, but remove edges and figure out how degmap+degheap should be changed in the end for the receivers actually used

			for (auto it = prevSharedOuts->begin(); it != prevSharedOuts->end(); ++it)
			{
				int receiver = *it;
				printf("deleting for receiver: %d\n",receiver);
				for(auto it2 = usedSenders->begin(); it2!=usedSenders->end();++it2)
				{
					int sender = *it2;
					printf("deleting  for sender: %d\n",sender);
					//out side of the edge
					std::set<int>::iterator deleteIter;
					printf("deleting outedge\n");
					deleteIter = outEdges[sender]->find(receiver);
					printf("found outedge iter to %d\n",receiver);
					if(deleteIter == outEdges[sender]->end())
					{
						printf("outedge missing node!\n");
					}
					outEdges[sender]->erase(deleteIter);
					printf("remove from degmap\n");
					degmapout[sender]--;
					//in side of the edge
					printf("deleting in edge\n");
					deleteIter = inEdges[receiver]->find(sender);
					inEdges[receiver]->erase(deleteIter);
				}
				printf("deleting for %d done\n",receiver);
			int newdeg = degmap[receiver] - senderCount;
			
		    sub_nbor_count(receiver,degmap,degheap,newdeg);	
			degmap[receiver] = newdeg;

		}
	}
	
	/*PairCount pc = *heap.begin();
    int preDepth = 0;
    if (depths.find(pc.fst) != depths.end())
      preDepth = std::max(preDepth, depths[pc.fst]);
    if (depths.find(pc.snd) != depths.end())
      preDepth = std::max(preDepth, depths[pc.snd]);
    depths[nv] = preDepth + 1;
    saved += pc.cnt;
    printf("pc[%d]: fst(%d) snd(%d) depth(%d) cnt(%d) acc_save(%d)\n", i, pc.fst, pc.snd, preDepth + 1, pc.cnt, saved);
    if (pc.cnt < 3) break;
    heap.erase(heap.begin());
    for (int j = 0; j < nv; j++)
      if (inEdges.find(j) != inEdges.end()) {
        std::set<int>* list = inEdges[j];
        if ((list->find(pc.fst) != list->end())
        &&  (list->find(pc.snd) != list->end())) {
          list->erase(pc.fst);
          list->erase(pc.snd);
          // update counters
          std::set<int>::const_iterator it;
          for (it = list->begin(); it != list->end(); it++) {
            sub_pair_count(*it, pc.fst, counter, heap);
            sub_pair_count(*it, pc.snd, counter, heap);
            add_pair_count(*it, nv, counter, heap);
          }
          list->insert(nv);
        }
      }
    nv ++;
 */
		i++;
  }
  printf("finalscore: %d\n",cumulativeScore);
  if (nedges == 0) 
	  nedges = 1;
  printf("finalratio: %f\n",(float)cumulativeScore/(float)nedges);
	auto end = std::chrono::steady_clock::now();
    int mstime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("runtime: %d\n",mstime);

}

