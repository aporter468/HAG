/* Copyright 2019 Stanford University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "gnn.h"
#include <unordered_map>
#define COUNT_MERGE_THRESHOLD 2

namespace std{
  template<>
  struct hash<std::pair<V_ID, V_ID> >
  {
    size_t operator()(const std::pair<V_ID, V_ID>& p) const
    {
      size_t res = 17;
      res = res * 31 + hash<V_ID>()(p.first);
      res = res * 31 + hash<V_ID>()(p.second);
      return res;
    }
  };
}

struct PairCount {
  PairCount(V_ID u, V_ID v, V_ID _cnt) : fst(u), snd(v), cnt(_cnt) {}
  V_ID fst, snd, cnt;
};

struct pair_count_compare {
  bool operator()(const PairCount& lhs, const PairCount& rhs) const {
    if (lhs.cnt != rhs.cnt) return (lhs.cnt > rhs.cnt);
    if (lhs.fst != rhs.fst) return (lhs.fst < rhs.fst);
    return (lhs.snd < rhs.snd);
  }
};

void add_pair_count(V_ID u, V_ID v,
                    std::unordered_map<std::pair<V_ID, V_ID>, V_ID>& counter,
                    std::set<PairCount, pair_count_compare>& heap)
{
  if (u > v) {V_ID w = u; u = v; v = w;}
  if (counter.find(std::make_pair(u, v)) == counter.end()) {
    counter[std::make_pair(u, v)] = 1;
    PairCount pc(u, v, 1);
    heap.insert(pc);
  } else {
    V_ID oldVal = counter[std::make_pair(u, v)];
    PairCount pc(u, v, oldVal);
    heap.erase(pc);
    counter[std::make_pair(u, v)] = oldVal + 1;
    pc.cnt = oldVal + 1;
    heap.insert(pc);
  }
}

void sub_pair_count(V_ID u, V_ID v,
                    std::unordered_map<std::pair<V_ID, V_ID>, V_ID>& counter,
                    std::set<PairCount, pair_count_compare>& heap)
{
  if (u > v) {V_ID w = u; u = v; v = w;}
  V_ID oldVal = counter[std::make_pair(u, v)];
  PairCount pc(u, v, oldVal);
  heap.erase(pc);
  counter[std::make_pair(u, v)] = oldVal - 1;
  pc.cnt = oldVal - 1;
  heap.insert(pc);
}

void transfer_graph(std::map<V_ID, std::set<V_ID>*>& orgList,
                    std::map<V_ID, std::set<V_ID>*>& optList,
                    std::vector<std::pair<V_ID, V_ID> >& ranges,
                    V_ID nv, E_ID ne, V_ID maxDepth, V_ID width, V_ID& newNv)
{
  printf("CP#1: nv = %d\n", nv);
  std::unordered_map<std::pair<V_ID, V_ID>, V_ID> counter;
  std::set<PairCount, pair_count_compare> heap;
  for (V_ID i = 0; i < nv; i++)
    if (orgList.find(i) != orgList.end()) {
      std::set<V_ID>::const_iterator it1, it2,
          first = orgList[i]->begin(), last = orgList[i]->end();
      for (it1 = first; it1 != last; it1 ++)
        for (it2 = first; it2 != it1; it2 ++) {
          V_ID u = *it2, v = *it1;
          assert(u < v);
          if (counter.find(std::make_pair(u, v)) == counter.end())
            counter[std::make_pair(u, v)] = 1;
          else
            counter[std::make_pair(u, v)] ++;
        }
      if (i % 1000 == 0) printf("v = %d\n", i);
    }
  std::unordered_map<std::pair<V_ID, V_ID>, V_ID>::const_iterator it;
  for (it = counter.begin(); it != counter.end(); it++) {
    PairCount pc(it->first.first, it->first.second, it->second);
    heap.insert(pc);
  }
  V_ID* depths = (V_ID*) malloc(width * sizeof(V_ID));
  V_ID v = nv;
  int saved = 0;
  while (v < nv + width) {
    if (heap.empty())
      break;
    PairCount pc = *heap.begin();
    heap.erase(heap.begin());
    if (pc.cnt <= COUNT_MERGE_THRESHOLD) break;
    saved += pc.cnt - 1;
    printf("[%d] fst = %d snd = %d cnt = %d acc_save(%d)\n",
           v, pc.fst, pc.snd, pc.cnt, saved);
    V_ID preDepth = (pc.fst < nv) ? 0 : depths[pc.fst - nv];
    if (pc.snd >= nv)
      preDepth = std::max(preDepth, depths[pc.snd - nv]);
    if (preDepth >= maxDepth) continue;
    depths[v - nv] = preDepth + 1;
    for (V_ID j = 0; j < v; j++)
      if (orgList.find(j) != orgList.end()) {
        std::set<V_ID>* list = orgList[j];
        if ((list->find(pc.fst) != list->end())
        &&  (list->find(pc.snd) != list->end())) {
          list->erase(pc.fst);
          list->erase(pc.snd);
          // update counters
          std::set<V_ID>::const_iterator it;
          for (it = list->begin(); it != list->end(); it++) {
            sub_pair_count(*it, pc.fst, counter, heap);
            sub_pair_count(*it, pc.snd, counter, heap);
            add_pair_count(*it, v, counter, heap);
          }
          list->insert(v);
        }
      }
    orgList[v] = new std::set<int>();
    orgList[v]->insert(pc.fst);
    orgList[v]->insert(pc.snd);
    v = v + 1;
  }
  newNv = v;
  // Reorder vertices by their depths
  // newIdxs[i] means vertex i is assigned with new ID newIdxs[i]
  V_ID* newIdxs = (V_ID*) malloc(newNv * sizeof(V_ID));
  for (V_ID i = 0; i < nv; i++)
    newIdxs[i] = i;
  V_ID nextIdx = nv;
  for (V_ID d = 1; d <= maxDepth; d++) {
    V_ID rangeLeft = nextIdx;
    for (V_ID i = nv; i < newNv; i++)
      if (depths[i - nv] == d)
        newIdxs[i] = nextIdx++;
    V_ID rangeRight = nextIdx - 1;
    if (rangeRight >= rangeLeft)
      ranges.push_back(std::make_pair(rangeLeft, rangeRight));
  }
  assert(nextIdx == newNv);
  // Construct optList
  for (V_ID i = 0; i < newNv; i++)
    if (orgList.find(i) != orgList.end()) {
      std::set<V_ID>::const_iterator it, first = orgList[i]->begin(),
                                     last = orgList[i]->end();
      for (it = first; it != last; it ++) {
        if (optList.find(newIdxs[i]) == optList.end())
          optList[newIdxs[i]] = new std::set<V_ID>();
        optList[newIdxs[i]]->insert(newIdxs[*it]);
      }
    }
  //for (int i = 0; i < ranges.size(); i++)
  //  printf("[%d] left = %d right = %d\n", i, ranges[i].first, ranges[i].second);
  //for (V_ID i = 0; i < newNv; i++)
  //  if (optList.find(i) != optList.end()) {
  //    std::set<V_ID>::const_iterator it, first = orgList[i]->begin(),
  //                                   last = orgList[i]->end();
  //    for (it = first; it != last; it++)
  //      printf("%d -> %d\n", *it, i);
  //  }
  free(depths);
  free(newIdxs);
}

