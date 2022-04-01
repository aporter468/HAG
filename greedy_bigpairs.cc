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



void useSenderPair(int i1,int i2,std::vector<int> &rankedSenders,
		std::unordered_map<int,int>& degmapout)
{
	int v1 = rankedSenders.at(i1);
	int v2 = rankedSenders.at(i2);
	printf("senders: %d %d degs: %d %d product: %d\n",v1,v2,degmapout[v1],degmapout[v2],degmapout[v1]*degmapout[v2]);
}

/*void del_nbor_count(int v, std::unordered_map<int,int>& degmap, std::set<NborSet,nbor_count_compare>& heap)
{

	int oldval = degmap[v];
    printf("delete nbor %d  with%d\n",v,oldval);
     NborSet topout =  *heap.begin();
     printf("top before sub: %d\n",topout.rec);
    NborSet ins(v,oldval);
    heap.erase(topout);
		  NborSet topoutafter =  *heap.begin();
     printf("top after sub: %d\n",topoutafter.rec);

}*/
void sub_nbor_count(int v, std::unordered_map<int,int>& degmap,
						std::set<NborSet,nbor_count_compare>& heap,int newval)
{
	int oldval = degmap[v];
//	printf("sub nbor %d  from %d to %d\n",v,oldval,newval);
	 NborSet topout =  *heap.begin();
	NborSet ins(v,oldval);
	heap.erase(ins);
	ins.cnt = newval;
	heap.insert(ins);
		NborSet topoutafter =  *heap.begin();
	
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
  std::set<NborSet, nbor_count_compare> outdegheap;
	 
  
  int nedges = 0;
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
			degmap[u] = 1;//rely on inEdges to match degmap on index sets
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
  printf("nedges: %d nedges: %d\n",nedges,nedges);
  fclose(file);
/* printf("nvert = %d\n", nv);
  for(int i =1; i<nv;i++)
  {
  		printf("vertex i: %d\n",i);
		std::set<int>::const_iterator init;
		std::set<int>::const_iterator outit;
		if(inEdges.find(i)!=inEdges.end()){
		for (init = inEdges[i]->begin(); init!=inEdges[i]->end();init ++)
		{
			printf("in: %d\n",*init);
		}}

		if(outEdges.find(i)!=outEdges.end()){
		for(outit=outEdges[i]->begin(); outit!=outEdges[i]->end();outit ++)
		{
			printf("out: %d\n",*outit);
		}
		}
  }
  */
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
  

  std::unordered_map<int,int>::const_iterator  it2,it3;
  for (it2 = degmap.begin();it2!=degmap.end();it2++)
  {
	  NborSet ins(it2->first,it2->second);
	  degheap.insert(ins);
  	

	}
  for (it3 = degmapout.begin();it3!=degmapout.end();it3++)
  {
	  NborSet ins(it3->first,it3->second);
		
		outdegheap.insert(ins);

  }

	std::vector<int>* rankedSenders;
	rankedSenders = new std::vector<int>();
	while(!outdegheap.empty())
	{
		NborSet topout = *outdegheap.begin();
		rankedSenders ->push_back(topout.rec);
		printf("next sender: %d  with deg %d\n",topout.rec,degmapout[topout.rec]);
		outdegheap.erase(topout);

	}


//  auto start = std::chrono::steady_clock::now();


  int saved = 0;
  std::map<int, int> depths;
 // for (int i = 0; i<nsteps; i++) {
  int i1 =0;
  int i2=  1;
  int tempProduct = degmapout[rankedSenders->at(i1+1)]*degmapout[rankedSenders->at(i1+2)];
  int vertsused = 0;
  printf("start main loop: initial degs: %d %d temp: %d\n",degmapout[rankedSenders->at(i1)],degmapout[rankedSenders->at(i2)],tempProduct);
  while( vertsused<nsteps && i1<(rankedSenders->size()-1) 
		  &&!((i1>=(rankedSenders->size()-2))&&(i2>=(rankedSenders->size()-1))) ) 
  {

	  
	  
	  std::vector<int>* usedSenders;
	  usedSenders = new std::vector<int>();
	
	 
 /* if(i2>=rankedSenders->size())
  	  {
			printf("i2 is too big? %d\n",i2);
		  i1++;
		  i2=i1+1;
		  if((i1+2)<rankedSenders->size()){

			  tempProduct = degmapout[rankedSenders->at(i1+1)]*degmapout[rankedSenders->at(i1+2)];
		}
		  else
		   {
		   		tempProduct = 0;
		   }
	}
	printf("start loop testing i1 %d i2 %d\n",i1,i2);
	  if( degmapout[rankedSenders->at(i1)]*degmapout[rankedSenders->at(i2)] > tempProduct)
	  {
		  printf("test product: %d tempProduct: %d\n",degmapout[rankedSenders->at(i1)]*degmapout[rankedSenders->at(i2)],tempProduct);
		  //TODO: use curent i1,i2
		  useSenderPair(i1,i2,*rankedSenders,degmapout);
		  usedSenders->push_back(rankedSenders->at(i1));
		  usedSenders->push_back(rankedSenders->at(i2));
		  printf("added senders to used senders\n");
		  i2++;
		 
	  }
	  else
	  {
		  printf("moving to next i1: %d to %d and %d to %d\n",i1,i1+1,i2,i1+2);
		  i1 = i1+1;
		  i2 = i1+1;//previous i1 + 2
			//use this i1,i2
		  useSenderPair(i1,i2,*rankedSenders,degmapout);
		  usedSenders->push_back(rankedSenders->at(i1));
		  usedSenders->push_back(rankedSenders->at(i2));
	
		  if((i1+2)<rankedSenders->size())
		  {
			  printf("in range to updated tempProduct\n");
			  tempProduct = degmapout[rankedSenders->at(i1+1)]*degmapout[rankedSenders->at(i1+2)];
			  	printf("using i1+2=%d; tempProduct now %d\n",(i1+2),tempProduct);
		  }
		  else
		  {
			  tempProduct = 0;
		  }
		  i2++;

	  }
		
*/
		//new -------
		usedSenders->push_back(rankedSenders->at(i1));
		usedSenders->push_back(rankedSenders->at(i2));
		printf("using i1,i2=(%d,%d) = (%d,%d)\n",i1,i2,rankedSenders->at(i1),rankedSenders->at(i2));

	  i1++;
	  i2=i1+1;
	
		//end new------
	    int v = usedSenders->at(0);
		std::set<int>* sharedOuts;
		sharedOuts = new std::set<int>();
		std::copy(outEdges[v]->begin(),outEdges[v]->end(),std::inserter(*sharedOuts,sharedOuts->begin()));		
		
		printf("first node outdeg: %d\n",sharedOuts->size());
		
		int u = usedSenders->at(1);
		std::set<int>* newset = outEdges[u];
		std::set<int> newintersect;// = new std::set<int>();
		
		printf("Second node outdeg: %d\n",newset->size());
		//intersection
		std::set_intersection(sharedOuts->begin(),sharedOuts->end(),newset->begin(),newset->end(),std::inserter(newintersect,newintersect.begin()));

		printf("intersection size: %d\n",newintersect.size());
	    int outsize = newintersect.size();
                    *sharedOuts = newintersect;

	if(outsize>1){	  		cumulativeScore += (outsize-1);}
			printf("cumulativescore: %d\n",cumulativeScore);
	
		printf("stepscore %d\n",outsize-1);
  //if it hasn't changed still need a printout
	if(outsize>1)
		{
	
		printf("best sender set: %d send to %d for score %d\n",2,outsize,outsize-1);

			for (auto it = sharedOuts->begin(); it != sharedOuts->end(); ++it)
			{
				int receiver = *it;
//				printf("deleting for receiver: %d\n",receiver);
				for(auto it2 = usedSenders->begin(); it2!=usedSenders->end();++it2)
				{
					int sender = *it2;
//					printf("deleting  for sender: %d\n",sender);
					//out side of the edge
					std::set<int>::iterator deleteIter;
//					printf("deleting outedge\n");
					deleteIter = outEdges[sender]->find(receiver);
//					printf("found outedge iter to %d\n",receiver);
					if(deleteIter == outEdges[sender]->end())
					{
						printf("outedge missing node!\n");
					}
					outEdges[sender]->erase(deleteIter);
//					printf("remove from degmap\n");
					degmapout[sender]--;
					//in side of the edge
//					printf("deleting in edge\n");
					deleteIter = inEdges[receiver]->find(sender);
					inEdges[receiver]->erase(deleteIter);
				}
//				printf("deleting for %d done\n",receiver);
			int newdeg = degmap[receiver] - 2;
			
		    sub_nbor_count(receiver,degmap,degheap,newdeg);	
			degmap[receiver] = newdeg;

		}
	    for(auto it3 = usedSenders->begin(); it3!=usedSenders->end(); ++it3)
		{
			int newdeg = degmapout[*it3] - outsize;
		}

	}
		
		
		vertsused++;	
		printf("end of loop pass; i1: %d i2: %d out of: %d\n",i1,i2,rankedSenders->size());
  }
  printf("finalscore: %d\n",cumulativeScore);
  if (nedges == 0) 
	  nedges = 1;
  printf("finalratio: %f\n",(float)cumulativeScore/(float)nedges);
	auto end = std::chrono::steady_clock::now();
    int mstime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		printf("runtime: %d\n",mstime);

}

