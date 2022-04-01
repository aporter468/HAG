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
	printf("heap erase :%d %d\n",ins.rec,ins.cnt);
	heap.erase(ins);
	ins.cnt = newval;
	printf("insert: %d %d\n",ins.rec,ins.cnt);
	heap.insert(ins);
}

void print_degree_distrib(std::unordered_map<int,int>& degmap,std::string printtag,  int nvert)
{
	for(int i =0; i<nvert;i++)
	{
		if(degmap.find(i)!=degmap.end())
		{
			printf("%s %d\n",printtag.c_str(),degmap[i]);
		}
	}
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
int nsteps2 = atoi(argv[3]);
	int u, v;
  int nv = 0;
  int cumulativeScore=0;
  std::map<int, std::set<int>* > inEdges;
  
  std::map<int, std::set<int>* > outEdges;

  std::unordered_map<int,int> degmap;
  
  std::unordered_map<int,int> degmapout;
  std::set<NborSet, nbor_count_compare> degheap;


  while (fscanf(file, "%d %d", &u, &v) != EOF) {
    if (std::max(u, v) >= nv)
      nv = std::max(u, v) + 1;
	printf("edge: %d %d\n",u,v);
	if (inEdges.find(v) == inEdges.end())
 	{
		printf("add edge list for node: %d\n",v);
		inEdges[v] = new std::set<int>();
		inEdges[v]->insert(u);
		degmap[v] = 1;//rely on inEdges to match degmap on index sets
	}
    else
	{
    	inEdges[v]->insert(u);
		degmap[v] ++;
	}


	if(outEdges.find(u) == outEdges.end())
	{
		outEdges[u] = new std::set<int>();
		outEdges[u]->insert(v);
		degmapout[u]=  1;
	}
	else
	{
		outEdges[u] ->insert(v);
		degmapout[u] ++;
	}
  
  }

  fclose(file);
  printf("nv = %d\n", nv);
 
  std::unordered_map<int,int>::const_iterator  it2;
  for (it2 = degmap.begin();it2!=degmap.end();it2++)
  {
	  NborSet ins(it2->first,it2->second);
	  degheap.insert(ins);
  }
  
  auto start = std::chrono::steady_clock::now();
	print_degree_distrib(degmap,"deg1",nv);

  for (int i = 0; i<nsteps; i++) {
		NborSet topins =  *degheap.begin();
		printf("top node: %d with deg: %d\n", topins.rec,topins.cnt);
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
			else if (j==(vdeg-1))//last iter so use score
            {
                printf("use last score of %d\n",score);
                cumulativeScore+=score;
                prevscore=score;
                usedSenders->push_back(topsend.rec);
                senderCount++;
            }
			else
			{
			
				prevscore = score;
				usedSenders->push_back(topsend.rec);
				senderCount++;
				
			}
		}

  //if it hasn't changed still need a printout
	  printf("cumulativescore: %d\n",cumulativeScore);
	   
		if(senderCount>0)
		{
			printf("best sender set: %d send to %d for score %d\n",senderCount,prevSharedOuts->size(),prevscore);
			printf("edgessaved: %d\n",prevscore);
				//leave outdeg heap alone b/c its reset at the outer loop anyway, but remove edges and figure out how degmap+degheap should be changed in the end for the receivers actually used

			for (auto it = prevSharedOuts->begin(); it != prevSharedOuts->end(); ++it)
			{
				int receiver = *it;
			//	printf("deleting for receiver: %d\n",receiver);
				for(auto it2 = usedSenders->begin(); it2!=usedSenders->end();++it2)
				{
					int sender = *it2;
					//printf("deleting  for sender: %d\n",sender);
					//out side of the edge
					std::set<int>::iterator deleteIter;
					deleteIter = outEdges[sender]->find(receiver);
					outEdges[sender]->erase(deleteIter);
					//in side of the edge
					deleteIter = inEdges[receiver]->find(sender);
					inEdges[receiver]->erase(deleteIter);
				}

			int newdeg = degmap[receiver] - senderCount;
			
		    sub_nbor_count(receiver,degmap,degheap,newdeg);	
			degmap[receiver] = newdeg;

		}
	}
	
  }

	printf("cumulativesave: %d\n",cumulativeScore);
	auto end = std::chrono::steady_clock::now();
    int mstime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("runtime: %d\n",mstime);




//experimental step: cut off high deg nodes - search by out for now

  std::set<NborSet, nbor_count_compare> degheap2;
  std::unordered_map<int,int>::const_iterator  it3;
	for (it3 = degmapout.begin();it3!=degmapout.end();it3++)
  {
	  printf("insert %d with degree %d\n",it3->first,it3->second);
	  NborSet ins(it3->first,it3->second);
	  degheap2.insert(ins);
  }
	printf("made new heap\n");
	int topdelete = 10;
   for (int i = 0; i<topdelete; i++) {
	   printf("deletion #%d\n",i);
		NborSet topins =  *degheap2.begin();
		int clearnode = topins.rec;
		printf("top node: %d with deg: %d\n", topins.rec,topins.cnt);
	    sub_nbor_count(clearnode,degmapout,degheap2,0);
  		NborSet topnew = *degheap2.begin();
		printf("new top top node: %d %d\n",topnew.rec,topnew.cnt);
		
		degmapout[clearnode] = 0;
		//remove this node as an in-edge as needed, and decrease degmap
		printf("degmap cleared,heap updated\n");

		std::set<int>::const_iterator  first  = outEdges[clearnode]->begin(),last =outEdges[clearnode]->end();
				
		std::set<int>::const_iterator edgeiter;
		for(edgeiter = first; edgeiter!=last;edgeiter++)
		{
			u = *edgeiter;
			printf("receiver u: %d found? %d\n",u,(inEdges[u]->find(clearnode)==inEdges[u]->end()));
			if(inEdges[u]->find(clearnode)!=inEdges[u]->end())
			{

				inEdges[u]->erase(inEdges[u]->find(clearnode));
				degmap[u]--;
			}
		
		}		
		printf("clearing neighbors of %d\n",clearnode);
		outEdges[clearnode]->clear();	
		inEdges[clearnode]->clear();
		printf("remaining out edges? %d in? %d degmapout: %d\n",(outEdges[clearnode]->begin()==outEdges[clearnode]->end()),(inEdges[clearnode]->begin()==inEdges[clearnode]->end()),degmapout[clearnode]);
	   
   }	 

   printf("call print degdist\n");
	print_degree_distrib(degmap,"deg2",nv);
   
//step 2:
//build heap post-step 1
//
//
 
  std::unordered_map<std::pair<int, int> , int> counter;
  std::set<PairCount, pair_count_compare> heap;
  int counterinsert = 0;
//nv hasn't been reduced, but we can ignore zero-degree nodes
  for (int i = 0; i< nv; i++)
	  if (inEdges.find(i) != inEdges.end()) {
      std::set<int>::const_iterator it1, it2;
      std::set<int>::const_iterator first = inEdges[i]->begin(), last = inEdges[i]->end();
      for (it1 = first; it1 != last; it1 ++)
        for (it2 = first; it2 != it1; it2 ++) {
          u = *it2;
          v = *it1;
          assert(u < v);
          if (counter.find(std::make_pair(u, v)) == counter.end())
		  {
            counter[std::make_pair(u, v)] = 1;
			counterinsert++;
		}
          else {
            counter[std::make_pair(u, v)] ++;
		  }
        }
      if (i % 1000 == 0) printf("i = %d\n", i);
    }

  // initialize heap
  std::unordered_map<std::pair<int, int>, int>::const_iterator it;
 int newheapsize = 0; 
  for (it = counter.begin(); it != counter.end(); it++) {
    PairCount pc(it->first.first, it->first.second, it->second);
    heap.insert(pc);
	newheapsize++;
  }
  printf("new heap size: %d %d\n",newheapsize,counterinsert);

  auto start2 = std::chrono::steady_clock::now();

 int saved = 0;
int scoreSaved = 0;
 std::map<int, int> depths;
printf("starting orig algo\n");
  for(int i =0; i<nsteps2;i++)
  {
	PairCount pc = *heap.begin();
    int preDepth = 0;
    if (depths.find(pc.fst) != depths.end())
      preDepth = std::max(preDepth, depths[pc.fst]);
    if (depths.find(pc.snd) != depths.end())
      preDepth = std::max(preDepth, depths[pc.snd]);
    depths[nv] = preDepth + 1;
    saved += pc.cnt;
    printf("pc[%d]: fst(%d) snd(%d) depth(%d) cnt(%d) acc_save(%d)\n", i, pc.fst, pc.snd, preDepth + 1, pc.cnt, saved);
    if (pc.cnt < 3) break;
	scoreSaved += (pc.cnt-1);
	cumulativeScore+=(pc.cnt-1);
	printf("cumulativesave: %d\n",cumulativeScore);
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
 

  }
	auto end2 = std::chrono::steady_clock::now();
    int mstime2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
		printf("step2 runtime: %d\n",mstime2);


  int totalscore = cumulativeScore + scoreSaved;
  printf("total score: %d\n",totalscore);
  printf("total time: %d\n",mstime+mstime2);
}

