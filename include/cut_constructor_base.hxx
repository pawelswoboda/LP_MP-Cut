#ifndef LP_MP_CUT_CONSTRUCTOR_BASE_HXX
#define LP_MP_CUT_CONSTRUCTOR_BASE_HXX

#include "LP_MP.h"
#include "lifted_factors_messages.hxx"

#ifdef LP_MP_PARALLEL
#include <omp.h>
#endif

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <list>
#include <future> // or async?

#include "union_find.hxx"
#include "graph.hxx"

#include "external/max_flow.hxx"

// basic methods for multicut and max-cut constructor. The corresponding {multicut|max-cut} constructors shall inherit from the base constructors defined here.

// to do: 
// - rename AddUnaryFactor to add_edge
// - rename odd_3_wheel to quadruplet
// - rename odd_3_bicycle_wheel to quintuplet

namespace LP_MP {

template<
  typename CUT_CONSTRUCTOR, // CRTP
  typename FACTOR_MESSAGE_CONNECTION, 
  typename EDGE_FACTOR, typename TRIPLET_FACTOR,
  typename UNARY_TRIPLET_MESSAGE_0, typename UNARY_TRIPLET_MESSAGE_1, typename UNARY_TRIPLET_MESSAGE_2
  >
class cut_constructor_base {
public:
   using cut_constructor_base_type = cut_constructor_base<CUT_CONSTRUCTOR, FACTOR_MESSAGE_CONNECTION, EDGE_FACTOR, TRIPLET_FACTOR, UNARY_TRIPLET_MESSAGE_0, UNARY_TRIPLET_MESSAGE_1, UNARY_TRIPLET_MESSAGE_2>;
   using FMC = FACTOR_MESSAGE_CONNECTION;

   using edge_factor = EDGE_FACTOR;
   using triplet_factor = TRIPLET_FACTOR;
   using edge_triplet_message_0 = UNARY_TRIPLET_MESSAGE_0;
   using edge_triplet_message_1 = UNARY_TRIPLET_MESSAGE_1;
   using edge_triplet_message_2 = UNARY_TRIPLET_MESSAGE_2;

   template<typename SOLVER>
   cut_constructor_base(SOLVER& pd)
   : lp_(&pd.GetLP())
   //unaryFactors_(100,hash::array2),
   //tripletFactors_(100,hash::array3)
   {}

   cut_constructor_base(const cut_constructor_base& o)
      : unaryFactors_(o.unaryFactors_),
      unaryFactorsVector_(o.unaryFactorsVector_),
      tripletFactors_(o.tripletFactors_),
      noNodes_(o.noNodes_),
      lp_(o.lp_) 
   {}

   bool get_edge_label(const INDEX i0, const INDEX i1) const
   {
      assert(i0 < i1 && i1 < noNodes_);
      assert(HasUnaryFactor(i0,i1));
      auto* f = GetUnaryFactor(i0,i1);
      return f->GetFactor()->primal()[0];
   }

   virtual edge_factor* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost) // declared virtual so that derived class notices when unary factor is added
   {
      assert(i1 < i2);
      noNodes_ = std::max(noNodes_,i2+1);
      assert(!HasUnaryFactor(i1,i2));
      
      auto* u = new edge_factor();
      (*u->GetFactor())[0] = cost;
      lp_->AddFactor(u);
      auto it = unaryFactors_.insert(std::make_pair(std::array<INDEX,2>{i1,i2}, u)).first;
      unaryFactorsVector_.push_back(std::make_pair(std::array<INDEX,2>{i1,i2}, u));

      if(it != unaryFactors_.begin()) {
         auto prevIt = it;
         --prevIt;
         assert(prevIt->second != u);
         lp_->AddFactorRelation(prevIt->second, u);
      }
      auto nextIt = it;
      ++nextIt;
      if(nextIt != unaryFactors_.end()) {
         assert(nextIt->second != u);
         lp_->AddFactorRelation(u, nextIt->second);
      }

      return u;
   }
   edge_factor* GetUnaryFactor(const INDEX i1, const INDEX i2) const {
      assert(HasUnaryFactor(i1,i2));
      return unaryFactors_.find(std::array<INDEX,2>{i1,i2})->second;
   }
   INDEX number_of_edges() const { return unaryFactors_.size(); }
   INDEX number_of_triplets() const { return tripletFactors_.size(); }

   template<typename MESSAGE_CONTAINER>
   MESSAGE_CONTAINER* LinkUnaryTriplet(edge_factor* u, triplet_factor* t) 
   {
      auto* m = new MESSAGE_CONTAINER(u, t);
      lp_->AddMessage(m);
      return m;
   }

   
   virtual triplet_factor* AddTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) // declared virtual so that derived constructor notices when triplet factor is added
   {
      assert(i1 < i2 && i2 < i3);
      assert(!HasTripletFactor(i1,i2,i3));
      if(!HasUnaryFactor(i1,i2)) {
         AddUnaryFactor(i1,i2,0.0);
      }
      if(!HasUnaryFactor(i1,i3)) {
         AddUnaryFactor(i1,i3,0.0);
      }
      if(!HasUnaryFactor(i2,i3)) {
         AddUnaryFactor(i2,i3,0.0);
      }
      assert(HasUnaryFactor(i1,i2) && HasUnaryFactor(i1,i3) && HasUnaryFactor(i2,i3));
      auto* t = new triplet_factor();
      lp_->AddFactor(t);
      tripletFactors_.insert(std::make_pair( std::array<INDEX,3>{i1,i2,i3}, t ));
      // use following ordering of unary and triplet factors: triplet comes after edge factor (i1,i2) and before (i2,i3)
      auto* before = GetUnaryFactor(i1,i2);
      lp_->AddFactorRelation(before,t);
      auto* middle = GetUnaryFactor(i1,i3);
      lp_->AddFactorRelation(middle,t);
      auto* after = GetUnaryFactor(i2,i3);
      lp_->AddFactorRelation(t,after);
      // get immediate predeccessor and successor and place new triplet in between
      //auto succ = tripletFactors_.upper_bound(std::make_tuple(i1,i2,i3));
      //if(succ != tripletFactors_.end()) {
      //   assert(t != succ->second);
      //   lp_->AddFactorRelation(t,succ->second);
      //}
      // link with all three unary factors
      LinkUnaryTriplet<edge_triplet_message_0>(before, t);
      LinkUnaryTriplet<edge_triplet_message_1>(middle, t);
      LinkUnaryTriplet<edge_triplet_message_2>(after, t);
      return t;
   }
   bool HasUnaryFactor(const std::tuple<INDEX,INDEX> e) const 
   {
      return HasUnaryFactor(std::get<0>(e), std::get<1>(e));
      assert(std::get<0>(e) < std::get<1>(e));
      return unaryFactors_.find(std::array<INDEX,2>{std::get<0>(e),std::get<1>(e)}) != unaryFactors_.end();
   }
   bool HasUnaryFactor(const INDEX i1, const INDEX i2) const 
   {
      assert(i1 < i2 && i2 < noNodes_);
      return (unaryFactors_.find(std::array<INDEX,2>{i1,i2}) != unaryFactors_.end());
   }
   bool HasTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) const 
   {
      assert(i1 < i2 && i2 < i3 && i3 < noNodes_);
      return (tripletFactors_.find(std::array<INDEX,3>{i1,i2,i3}) != tripletFactors_.end());
   }

   triplet_factor* GetTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3) const 
   {
      assert(HasTripletFactor(i1,i2,i3));
      return tripletFactors_.find(std::array<INDEX,3>{i1,i2,i3})->second;
   }


   std::tuple<INDEX,INDEX> GetEdge(const INDEX i1, const INDEX i2) const
   {
      return std::make_tuple(std::min(i1,i2), std::max(i1,i2));
   }

   REAL get_edge_cost(const INDEX i1, const INDEX i2) const
   {
      assert(HasUnaryFactor(i1,i2));
      return *(unaryFactors_.find(std::array<INDEX,2>{i1,i2})->second->GetFactor());
   }

   INDEX NumOfUnaryFactors() const {
      return unaryFactorsVector_.size();
   }

   std::array<INDEX,2> GetEdge(const INDEX e) const {
      return unaryFactorsVector_[e].first;
   }

   REAL GetEdgeCost(const INDEX e) const {
      return *(unaryFactorsVector_[e].second->GetFactor());
   }

   struct triplet_candidate {
      std::array<INDEX,3> nodes;
      REAL cost;
      bool operator<(const triplet_candidate& o) const {
         return this->cost > o.cost;
      }
   };

   // search for violated triplets, e.g. triplets with one negative edge and two positive ones.
   INDEX find_violated_triplets(const INDEX max_triplets_to_add)
   {
      std::vector<INDEX> adjacency_list_count(noNodes_,0);
      // first determine size for adjacency_list
      for(auto& it : unaryFactorsVector_) {
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         adjacency_list_count[i]++;
         adjacency_list_count[j]++; 
      }
      two_dim_variable_array<std::tuple<INDEX,REAL>> adjacency_list(adjacency_list_count);
      std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);
      for(auto& it : unaryFactorsVector_) {
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         const REAL cost_ij = (*it.second->GetFactor())[0];
         assert(i<j);
         adjacency_list[i][adjacency_list_count[i]] = std::make_tuple(j,cost_ij);
         adjacency_list_count[i]++;
         adjacency_list[j][adjacency_list_count[j]] = std::make_tuple(i,cost_ij);
         adjacency_list_count[j]++;
      }

      // Sort the adjacency list, for fast intersections later
      auto adj_sort = [](const auto a, const auto b) { return std::get<0>(a) < std::get<0>(b); };

#pragma omp parallel for schedule(guided)
      for(int i=0; i < adjacency_list.size(); i++) {
         std::sort(adjacency_list[i].begin(), adjacency_list[i].end(), adj_sort);
      }

      // Iterate over all of the edge intersection sets
      // do zrobienia: parallelize
      // we will intersect two adjacency list by head node, but we want to preserve the costs of either edge pointing to head node
      using intersection_type = std::tuple<INDEX,REAL,REAL>;
      auto merge = [](const auto a, const auto b) -> intersection_type { 
         assert(std::get<0>(a) == std::get<0>(b));
         return std::make_tuple(std::get<0>(a), std::get<1>(a), std::get<1>(b)); 
      };
      std::vector<triplet_candidate> triplet_candidates;
#pragma omp parallel
      {
         std::vector<intersection_type> commonNodes(noNodes_);
         std::vector<triplet_candidate> triplet_candidates_per_thread;
#pragma omp for  schedule(guided) nowait
         for(INDEX c=0; c<unaryFactorsVector_.size(); ++c) {
            const REAL cost_ij = (*unaryFactorsVector_[c].second->GetFactor())[0];
            const INDEX i = std::get<0>(unaryFactorsVector_[c].first);
            const INDEX j = std::get<1>(unaryFactorsVector_[c].first);

            // Now find all neighbors of both i and j to see where the triangles are
            // TEMP TEMP -- fails at i=0, j=1, on i==3.
            auto intersects_iter_end = set_intersection_merge(adjacency_list[i].begin(), adjacency_list[i].end(), adjacency_list[j].begin(), adjacency_list[j].end(), commonNodes.begin(), adj_sort, merge);

            for(auto n=commonNodes.begin(); n != intersects_iter_end; ++n) {
               const INDEX k = std::get<0>(*n);

               // Since a triplet shows up three times as an edge plus
               // a node, we only consider it for the case when i<j<k 
               if(!(j<k))
                  continue;
               const REAL cost_ik = std::get<1>(*n);
               const REAL cost_jk = std::get<2>(*n);

               const REAL lb = std::min(0.0, cost_ij) + std::min(0.0, cost_ik) + std::min(0.0, cost_jk);
               const REAL best_labeling = static_cast<CUT_CONSTRUCTOR*>(this)->triplet_cost(cost_ij, cost_ik, cost_jk);
               assert(lb <= best_labeling + eps);
               const REAL guaranteed_dual_increase = best_labeling - lb;
               if(guaranteed_dual_increase > eps) {
                 triplet_candidates_per_thread.push_back({std::array<INDEX,3>({i,j,k}), guaranteed_dual_increase});
               } 
            }
         }
#pragma omp critical
         {
            triplet_candidates.insert(triplet_candidates.end(), triplet_candidates_per_thread.begin(), triplet_candidates_per_thread.end()); 
         }
      }
      std::sort(triplet_candidates.begin(), triplet_candidates.end());

      if(triplet_candidates.size() > 0 && diagnostics()) {
         std::cout << "best triplet candidate in triplet search has guaranteed dual improvement " << triplet_candidates[0].cost << "\n";
      }

      INDEX triplets_added = 0;
      for(const auto& triplet_candidate : triplet_candidates) {
         const INDEX i = triplet_candidate.nodes[0];
         const INDEX j = triplet_candidate.nodes[1];
         const INDEX k = triplet_candidate.nodes[2];
         if(!HasTripletFactor(i,j,k)) {
            AddTripletFactor(i,j,k);
            triplets_added++;
            if(triplets_added > max_triplets_to_add) {
               break;
            } 
         }
      }

      return triplets_added;
   }

   template<typename ITERATOR>
   void cycle_normal_form(ITERATOR cycle_begin, ITERATOR cycle_end) const
   {
      assert(std::distance(cycle_begin, cycle_end) >= 3);
      //assert(std::distance(cycle_begin, cycle_end)%2 == 1);
      // first search for smallest entry and make it first
      std::rotate(cycle_begin, std::min_element(cycle_begin, cycle_end), cycle_end);
      // now two choices left: we can traverse cycle in forward or backward direction. Choose direction such that second entry is smaller than in reverse directoin.
      if(*(cycle_begin+1) > *(cycle_end - 1)) {
         std::reverse(cycle_begin+1, cycle_end);
      }
   } 

   template<typename ITERATOR>
   void triangulate_cycle(const REAL cost, ITERATOR path_begin, ITERATOR path_end, std::vector<triplet_candidate>& candidates)
   {
      assert(std::distance(path_begin, path_end) >= 3);
      cycle_normal_form(path_begin, path_end);
      const INDEX first_node = *path_begin;
      for(auto it=path_begin+1; it+1!=path_end; ++it) {
         if(first_node != *it && first_node != *(it+1)) {
            std::array<INDEX,3> nodes({first_node, *it, *(it+1)});
            std::sort(nodes.begin(), nodes.end());
            assert(HasUniqueValues(nodes));
            candidates.push_back({nodes, cost});
         }
      } 
   }

   void write_labeling_into_factors(const std::vector<char>& labeling) {
      assert(labeling.size() <= unaryFactorsVector_.size());
      for(INDEX c=0; c<labeling.size(); ++c) {
         auto* f = unaryFactorsVector_[c].second;
         f->GetFactor()->primal()[0] = labeling[c];
         f->propagate_primal_through_messages();
      }

      // possibly, additional edges have been added because of tightening. infer labeling of those from union find datastructure
      if(labeling.size() < unaryFactorsVector_.size()) {
         UnionFind uf(noNodes_);
         for(INDEX c=0; c<labeling.size(); ++c) {
            edge_factor* f = unaryFactorsVector_[c].second; 
            if(f->GetFactor()->primal()[0] == false) {
               // connect components 
               const INDEX i = std::get<0>(unaryFactorsVector_[c].first);
               const INDEX j = std::get<1>(unaryFactorsVector_[c].first);
               uf.merge(i,j);
            }
         }
         if(debug()) {
            std::cout << "built union find structure, propagate information now\n";
         }
         for(INDEX c=labeling.size(); c<unaryFactorsVector_.size(); ++c) {
            edge_factor* f = unaryFactorsVector_[c].second; 
            const INDEX i = std::get<0>(unaryFactorsVector_[c].first);
            const INDEX j = std::get<1>(unaryFactorsVector_[c].first);
            if(uf.connected(i,j)) {
               f->GetFactor()->primal()[0] = false;
            } else {
               f->GetFactor()->primal()[0] = true;
            }
            f->propagate_primal_through_messages();
         }
      }
   }

   template<typename STREAM>
   void WritePrimal(STREAM& s)
   {
      assert(no_original_edges_ <= unaryFactorsVector_.size());
      for(INDEX e=0; e<std::min(no_original_edges_, INDEX(unaryFactorsVector_.size())); ++e) {
         const INDEX i = unaryFactorsVector_[e].first[0];
         const INDEX j = unaryFactorsVector_[e].first[1];
         auto* t = unaryFactorsVector_[e].second->GetFactor();
         const bool cut = t->primal()[0];
         s << i << " " << j << " " << cut << "\n";
      } 
   }

   INDEX Tighten(const INDEX max_factors_to_add)
   {
      no_original_edges_ = std::min(no_original_edges_, INDEX(unaryFactorsVector_.size()));

      if(number_of_edges() > 2) {
         if(diagnostics()) {
            std::cout << "Search for violated triplet constraints\n";
         }
         std::size_t triplets_added = find_violated_triplets(max_factors_to_add);
         if(diagnostics()) {
            std::cout << "Added " << triplets_added << " triplet(s) out of " <<  max_factors_to_add << " by searching for triplets\n"; 
         }
         if(triplets_added < 0.6*max_factors_to_add) {
            if(diagnostics()) {
               std::cout << "Additionally search via shortest paths for violated constraints\n";
            }
            triplets_added += static_cast<CUT_CONSTRUCTOR*>(this)->find_violated_cycles(max_factors_to_add - triplets_added);
            if(diagnostics()) {
               std::cout << "Added " << triplets_added << " triplet(s) out of " <<  max_factors_to_add << " in total\n";
            }
         }
         return triplets_added;
      } else {
         return 0;
      }
   }

   void End() 
   {
      // wait for the primal rounding to finish.
      if(debug()) {
         std::cout << "wait for primal computation to end\n";
      }
      if(primal_handle_.valid()) {
         auto labeling = primal_handle_.get();
         write_labeling_into_factors(labeling);
      }
   }


   struct edge {
     const INDEX i,j;
     const REAL cost;
   };

   void round()
   {
      if(diagnostics()) {
         std::cout << "compute primal cut\n";
      }
      std::vector<edge> edges;
      edges.reserve(unaryFactorsVector_.size());

      for(const auto& e : unaryFactorsVector_) {
        const INDEX i = e.first[0];
        const INDEX j = e.first[1];
        const REAL cost = e.second->GetFactor()->operator[](0);
        edges.push_back({i,j,cost});
      }

      primal_handle_ = std::async(std::launch::async, call_round, this, std::move(edges)); 
   }

   static std::vector<char> call_round(cut_constructor_base_type* class_ptr, std::vector<edge> edges)
   {
	   return static_cast<CUT_CONSTRUCTOR*>(class_ptr)->round(std::move(edges));
   }

   // skeleton for asynchronous calling of primal rounding
   void ComputePrimal()
   {
     // in the first round, just start rounding without reading out any solution.
     if(!primal_handle_.valid()) { 
       round();
       return;
     }

     const auto primal_state = primal_handle_.wait_for(std::chrono::seconds(0));

     if(primal_state == std::future_status::deferred) {
       assert(false); // this should not happen, we launch immediately!
       throw std::runtime_error("asynchronuous primal rounding was deferred, but this should not happen");
     } else if(primal_state == std::future_status::ready) {

       if(debug()) {
         std::cout << "collect rounding result\n";
       }
       auto labeling = primal_handle_.get();
       write_labeling_into_factors(labeling);

       if(debug()) {
         std::cout << "restart primal rounding\n";
       }
       round();

     } else {
       if(debug()) {
         std::cout << "rounding is currently running.\n";
       }
     }
   }


protected:
   std::map<std::array<INDEX,2>, edge_factor*> unaryFactors_; // actually unary factors in multicut are defined on edges. assume first index < second one
   INDEX no_original_edges_ = std::numeric_limits<INDEX>::max();
   //std::unordered_map<std::array<INDEX,2>, edge_factor*> unaryFactors_; // actually unary factors in multicut are defined on edges. assume first index < second one
   std::vector<std::pair<std::array<INDEX,2>, edge_factor*>> unaryFactorsVector_; // we store a second copy of unary factors for faster iterating
   // sort triplet factors as follows: Let indices be i=(i1,i2,i3) and j=(j1,j2,j3). Then i<j iff i1+i2+i3 < j1+j2+j3 or for ties sort lexicographically
   struct tripletComp {
      bool operator()(const std::tuple<INDEX,INDEX,INDEX> i, const std::tuple<INDEX,INDEX,INDEX> j) const
      {
         const INDEX iSum = std::get<0>(i) + std::get<1>(i) + std::get<2>(i);
         const INDEX jSum = std::get<0>(j) + std::get<1>(j) + std::get<2>(j);
         if(iSum < jSum) return true;
         else if(iSum > jSum) return false;
         else return i<j; // lexicographic comparison
      }
   };
   std::unordered_map<std::array<INDEX,3>, triplet_factor*> tripletFactors_; // triplet factors are defined on cycles of length three
   INDEX noNodes_ = 0;

   LP* lp_;

   decltype(std::async(std::launch::async, call_round, nullptr, std::vector<edge>{})) primal_handle_;
};

template<
   class CUT_CONSTRUCTOR, // CRTP
   class BASE_CONSTRUCTOR, 
   typename ODD_3_WHEEL_FACTOR,
   typename TRIPLET_ODD_3_WHEEL_MESSAGE_012, typename TRIPLET_ODD_3_WHEEL_MESSAGE_013, typename TRIPLET_ODD_3_WHEEL_MESSAGE_023, typename TRIPLET_ODD_3_WHEEL_MESSAGE_123
>
class odd_wheel_constructor_base : public BASE_CONSTRUCTOR {
public:
   using FMC = typename BASE_CONSTRUCTOR::FMC;
   using base_constructor = BASE_CONSTRUCTOR;

   using odd_3_wheel_factor_container = ODD_3_WHEEL_FACTOR;
   using triplet_odd_3_wheel_message_012_container = TRIPLET_ODD_3_WHEEL_MESSAGE_012;
   using triplet_odd_3_wheel_message_013_container = TRIPLET_ODD_3_WHEEL_MESSAGE_013;
   using triplet_odd_3_wheel_message_023_container = TRIPLET_ODD_3_WHEEL_MESSAGE_023;
   using triplet_odd_3_wheel_message_123_container = TRIPLET_ODD_3_WHEEL_MESSAGE_123;

   using base_constructor::base_constructor;

   // add triplet indices additionally to tripletIndices_
   virtual typename base_constructor::edge_factor* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost)
   {
      typename base_constructor::edge_factor* u = base_constructor::AddUnaryFactor(i1,i2,cost);   
      return u;
   }

   // add triplet indices additionally to tripletIndices_
   virtual typename base_constructor::triplet_factor* AddTripletFactor(const INDEX i1, const INDEX i2, const INDEX i3)
   {
      assert(i1 < i2 && i2 < i3);
      assert(i3 < base_constructor::noNodes_);
      auto* t = base_constructor::AddTripletFactor(i1,i2,i3);
      triplet_vector_.push_back(std::make_tuple(i1,i2,i3,t));
      return t;
   }

   struct triplet_connection {
      std::array<INDEX,2> nodes;
      typename base_constructor::triplet_factor* f;

      bool operator<(const triplet_connection& b) const
      {
         if(this->nodes[0] != b.nodes[0]) {
            return this->nodes[0] < b.nodes[0];
         } else {
            return this->nodes[1] < b.nodes[1];
         }
      }
   };

   using triplet_connections = two_dim_variable_array<triplet_connection>;
   triplet_connections compute_connected_triplets() const
   {
#ifdef LP_MP_PARALLEL
      std::vector<std::atomic<INDEX>> no_triplets_per_node(this->noNodes_, 0);
#else
      std::vector<INDEX> no_triplets_per_node(this->noNodes_, 0);
#endif

#pragma omp parallel for schedule(guided)
      for(INDEX t=0; t<triplet_vector_.size(); ++t) {
         const INDEX i = std::get<0>(triplet_vector_[t]);
         const INDEX j = std::get<1>(triplet_vector_[t]);
         const INDEX k = std::get<2>(triplet_vector_[t]);
         no_triplets_per_node[i]++;
         no_triplets_per_node[j]++;
         no_triplets_per_node[k]++; 
      }

      triplet_connections connected_triplets(no_triplets_per_node);
      std::fill(no_triplets_per_node.begin(), no_triplets_per_node.end(), 0);
      for(INDEX t=0; t<triplet_vector_.size(); ++t) {
         const INDEX i = std::get<0>(triplet_vector_[t]);
         const INDEX j = std::get<1>(triplet_vector_[t]);
         const INDEX k = std::get<2>(triplet_vector_[t]);
         auto* f = std::get<3>(triplet_vector_[t]);

         connected_triplets[i][ no_triplets_per_node[i] ] = triplet_connection({j,k, f});
         no_triplets_per_node[i]++;
         connected_triplets[j][ no_triplets_per_node[j] ] = triplet_connection({i,k, f});
         no_triplets_per_node[j]++;
         connected_triplets[k][ no_triplets_per_node[k] ] = triplet_connection({i,j, f});
         no_triplets_per_node[k]++;
      }

#pragma omp parallel for schedule(guided)
      for(INDEX i=0; i<connected_triplets.size(); ++i) {
         std::sort(connected_triplets[i].begin(), connected_triplets[i].end());
      }
      return std::move(connected_triplets);
   }
   
   odd_3_wheel_factor_container* add_odd_3_wheel_factor(const INDEX i0, const INDEX i1, const INDEX i2, const INDEX i3)
   {
      assert(!has_odd_3_wheel_factor(i0,i1,i2,i3));

      auto* f = new odd_3_wheel_factor_container();
      this->lp_->AddFactor(f);
      odd_3_wheel_factors_.insert(std::make_pair(std::array<INDEX,4>({i0,i1,i2,i3}), f));

      if(!this->HasTripletFactor(i0,i1,i2)) {
         this->AddTripletFactor(i0,i1,i2);
      }
      auto* t012 = this->GetTripletFactor(i0,i1,i2);
      auto* m012 = new triplet_odd_3_wheel_message_012_container(t012, f);
      this->lp_->AddMessage(m012);
      this->lp_->AddFactorRelation(t012, f);

      if(!this->HasTripletFactor(i0,i1,i3)) {
         this->AddTripletFactor(i0,i1,i3);
      }
      auto* t013 = this->GetTripletFactor(i0,i1,i3);
      auto* m013 = new triplet_odd_3_wheel_message_013_container(t013, f);
      this->lp_->AddMessage(m013);

      if(!this->HasTripletFactor(i0,i2,i3)) {
         this->AddTripletFactor(i0,i2,i3);
      }
      auto* t023 = this->GetTripletFactor(i0,i2,i3);
      auto* m023 = new triplet_odd_3_wheel_message_023_container(t023, f);
      this->lp_->AddMessage(m023);

      if(!this->HasTripletFactor(i1,i2,i3)) {
         this->AddTripletFactor(i1,i2,i3);
      }
      auto* t123 = this->GetTripletFactor(i1,i2,i3);
      auto* m123 = new triplet_odd_3_wheel_message_123_container(t123, f);
      this->lp_->AddMessage(m123);
      this->lp_->AddFactorRelation(f, t012);

      return f;
   }
   bool has_odd_3_wheel_factor(const INDEX i0, const INDEX i1, const INDEX i2, const INDEX i3) const
   {
      assert(i0 < i1 && i1 < i2 && i2 < i3);
      return odd_3_wheel_factors_.find(std::array<INDEX,4>({i0,i1,i2,i3})) != odd_3_wheel_factors_.end();
   }
   odd_3_wheel_factor_container* get_odd_3_wheel_factor(const INDEX i0, const INDEX i1, const INDEX i2, const INDEX i3) const
   {
      assert(has_odd_3_wheel_factor(i0,i1,i2,i3));
      return odd_3_wheel_factors_.find(std::array<INDEX,4>({i0,i1,i2,i3}))->second; 
   }

   // what is the lowest threshold so that an odd cycle exists
   REAL compute_odd_cycle_threshold(
         std::unordered_map<INDEX,INDEX>& origToCompressedNode,
         std::vector<INDEX>& compressedToOrigNode,
         std::vector<std::tuple<INDEX,INDEX,REAL>>& compressedEdges 
         )
   {
      std::sort(compressedEdges.begin(), compressedEdges.end(), [](auto a, auto b) { return std::get<2>(a) > std::get<2>(b); });

      const INDEX noCompressedNodes = origToCompressedNode.size();
      const INDEX noBipartiteCompressedNodes = 2*noCompressedNodes;
      UnionFind uf(noBipartiteCompressedNodes); 
      // construct bipartite graph based on triangles
      for(auto& e : compressedEdges) {
         const INDEX jc = std::get<0>(e);
         const INDEX kc = std::get<1>(e);
         uf.merge(jc,noCompressedNodes + kc);
         uf.merge(noCompressedNodes + jc,kc);
         if(uf.connected(jc, noCompressedNodes + jc)) {
            assert(uf.connected(kc, noCompressedNodes + kc));
            return std::get<2>(e);
         }
      }
      return -std::numeric_limits<REAL>::infinity(); // no constraint found
   }

   struct odd_3_wheel_candidate {
      std::array<INDEX,4> nodes;
      REAL cost;
   };

   template<typename ITERATOR>
   void triangulate_odd_wheel(const INDEX i, const REAL cost, ITERATOR path_begin, ITERATOR path_end, std::vector<odd_3_wheel_candidate>& candidates)
   {
      assert(std::distance(path_begin, path_end) >= 3);
      this->cycle_normal_form(path_begin, path_end);
      const INDEX first_node = *path_begin;
      for(auto it=path_begin+1; it+1!=path_end; ++it) {
         std::array<INDEX,4> nodes({i,first_node, *it, *(it+1)});
         std::sort(nodes.begin(), nodes.end());
         assert(HasUniqueValues(nodes));
         candidates.push_back({nodes, cost});
      } 
   }

   // explicitly enumerate all 3-wheels (complete 4-graphs) that are present in the graph and check whether adding them would guarantee dual bound increase
   INDEX find_odd_3_wheels(const INDEX max_factors_to_add)
   {
      std::vector<INDEX> adjacency_list_count(this->noNodes_,0);
      // first determine size for adjacency_list
      for(auto& it : this->unaryFactorsVector_) {
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         adjacency_list_count[i]++;
         adjacency_list_count[j]++; 
      }
      two_dim_variable_array<INDEX> adjacency_list(adjacency_list_count);
      std::fill(adjacency_list_count.begin(), adjacency_list_count.end(), 0);
      for(auto& it : this->unaryFactorsVector_) {
         const INDEX i = std::get<0>(it.first);
         const INDEX j = std::get<1>(it.first);
         assert(i<j);
         adjacency_list[i][adjacency_list_count[i]] = j;
         adjacency_list_count[i]++;
         adjacency_list[j][adjacency_list_count[j]] = i;
         adjacency_list_count[j]++;
      }

      // Sort the adjacency list, for fast intersections later
#pragma omp parallel for  schedule(guided)
      for(int i=0; i < adjacency_list.size(); i++) {
         std::sort(adjacency_list[i].begin(), adjacency_list[i].end());
      } 

      std::vector<odd_3_wheel_candidate> odd_3_wheel_candidates;

#pragma omp parallel
      {
         std::vector<odd_3_wheel_candidate> odd_3_wheel_candidates_local;
         typename odd_3_wheel_factor_container::FactorType test_odd_3_wheel_factor;

         std::vector<INDEX> commonNodes(this->noNodes_);
#pragma omp for schedule(guided) nowait
         for(INDEX i=0; i<triplet_vector_.size(); ++i) {
            const INDEX i1 = std::get<0>(triplet_vector_[i]);
            const INDEX i2 = std::get<1>(triplet_vector_[i]);
            const INDEX i3 = std::get<2>(triplet_vector_[i]);
            auto* f = std::get<3>(triplet_vector_[i])->GetFactor();
            // search for node node k such that there exists an edge from i1, i2 and i3 to it.
            // For this, intersect adjacency lists coming from i1,i2 and i3
            // to do: do not intersect twice but do it at once
            auto intersects_iter_end = std::set_intersection(adjacency_list[i1].begin(), adjacency_list[i1].end(), adjacency_list[i2].begin(), adjacency_list[i2].end(), commonNodes.begin());
            intersects_iter_end = std::set_intersection(commonNodes.begin(), intersects_iter_end, adjacency_list[i3].begin(), adjacency_list[i3].end(), commonNodes.begin());
            for(auto it=commonNodes.begin(); it!=intersects_iter_end; ++it) {
               std::fill(test_odd_3_wheel_factor.begin(), test_odd_3_wheel_factor.end(), 0.0);
               const INDEX k = *it; 
               std::array<INDEX,4> nodes({i1,i2,i3,k});
               std::sort(nodes.begin(), nodes.end());
               if(!has_odd_3_wheel_factor(nodes[0], nodes[1], nodes[2], nodes[3])
                     && this->HasTripletFactor(nodes[0], nodes[1], nodes[2])
                     && this->HasTripletFactor(nodes[0], nodes[1], nodes[3])
                     && this->HasTripletFactor(nodes[0], nodes[2], nodes[3])
                     && this->HasTripletFactor(nodes[1], nodes[2], nodes[3]) 
                     && k == nodes[3] // otherwise we will iterate over the same set of nodes 4 times
                 ) { 
                  // compute guaranteed dual increase from adding odd 3 wheel on nodes.
                  // For this purpose, go over all trplet factors acting on nodes and compute lb over each of them separately as well as the best labeling.
                  // The latter is obtained by reparametrizing everything into a test odd 3 wheel factor and then computing its lower bound.
                  // Question: Can adding an odd 3 wheel be advantageous, even if not all sub-triplets are present or is then necessarily some sub-triplet already violated? If the former is true, we need to take edge factors into account as well. but it might get messay then.

                  //auto* f_i1_k = this->GetUnaryFactor(std::min(i1,k), std::max(i1,k))->GetFactor();
                  //auto* f_i2_k = this->GetUnaryFactor(std::min(i2,k), std::max(i2,k))->GetFactor();
                  //auto* f_i3_k = this->GetUnaryFactor(std::min(i3,k), std::max(i3,k))->GetFactor();

                  //lb += f_i1_k->LowerBound();
                  //lb += f_i2_k->LowerBound();
                  //lb += f_i3_k->LowerBound();

                  REAL lb = 0.0;
                  // go over all 3-subset of nodes containing k
                  {
                     auto* sub_f = this->GetTripletFactor(nodes[0], nodes[1], nodes[2])->GetFactor();
                     lb += sub_f->LowerBound();
                     typename triplet_odd_3_wheel_message_012_container::MessageType m;
                     m.RepamRight(test_odd_3_wheel_factor, *sub_f); 
                  }
                  {
                     auto* sub_f = this->GetTripletFactor(nodes[0], nodes[1], nodes[3])->GetFactor();
                     lb += sub_f->LowerBound();
                     typename triplet_odd_3_wheel_message_013_container::MessageType m;
                     m.RepamRight(test_odd_3_wheel_factor, *sub_f); 
                  }
                  { 
                     auto* sub_f = this->GetTripletFactor(nodes[0], nodes[2], nodes[3])->GetFactor();
                     lb += sub_f->LowerBound();
                     typename triplet_odd_3_wheel_message_023_container::MessageType m;
                     m.RepamRight(test_odd_3_wheel_factor, *sub_f); 
                  }
                  { 
                     auto* sub_f = this->GetTripletFactor(nodes[1], nodes[2], nodes[3])->GetFactor();
                     lb += sub_f->LowerBound();
                     typename triplet_odd_3_wheel_message_123_container::MessageType m;
                     m.RepamRight(test_odd_3_wheel_factor, *sub_f); 
                  } 
                  const REAL best_labeling_cost = test_odd_3_wheel_factor.LowerBound();
                  assert(lb <= best_labeling_cost + eps);
                  if(lb < best_labeling_cost - eps) {
                     odd_3_wheel_candidates_local.push_back({nodes, best_labeling_cost - lb});
                  }
               }
            }
         }
#pragma omp critical
         odd_3_wheel_candidates.insert(odd_3_wheel_candidates.end(), odd_3_wheel_candidates_local.begin(), odd_3_wheel_candidates_local.end());
      }

      std::sort(odd_3_wheel_candidates.begin(), odd_3_wheel_candidates.end(), [](const auto& a, const auto& b) { return a.cost > b.cost; });
      if(odd_3_wheel_candidates.size() > 2) {
         for(INDEX i=0; i<odd_3_wheel_candidates.size()-1; ++i) {
            assert(odd_3_wheel_candidates[i].cost > odd_3_wheel_candidates[i+1].cost);
         }
      }

      INDEX factors_added = 0;
      for(INDEX i=0; i<odd_3_wheel_candidates.size(); ++i) {
         auto& nodes = odd_3_wheel_candidates[i].nodes;
         if(!has_odd_3_wheel_factor(nodes[0], nodes[1], nodes[2], nodes[3])) {
            add_odd_3_wheel_factor(nodes[0], nodes[1], nodes[2], nodes[3]);
            ++factors_added;
         }
         if(factors_added > max_factors_to_add) {
            break;
         }
      }

      if(odd_3_wheel_candidates.size() > 0 && diagnostics()) {
         std::cout << "added " << factors_added << " by local odd 3 wheel search with guaranteed dual improvement " << odd_3_wheel_candidates[0].cost << "\n";
      }

      return factors_added;
   }

   INDEX Tighten(const INDEX max_factors_to_add)
   {
      const INDEX triplets_added = BASE_CONSTRUCTOR::Tighten(max_factors_to_add);
      if(this->number_of_edges() < 2) { return 0; }
      if(triplets_added > 0.1*max_factors_to_add) {
        return triplets_added;
      } else {
        const INDEX odd_3_wheels_added = find_odd_3_wheels(max_factors_to_add - triplets_added);
        if(odd_3_wheels_added > 0.1*max_factors_to_add) {
          return odd_3_wheels_added + triplets_added;
        } else {
          const INDEX odd_wheels_added = static_cast<CUT_CONSTRUCTOR*>(this)->find_odd_wheels(max_factors_to_add - triplets_added - odd_3_wheels_added);
          if(diagnostics()) {
            std::cout << "Added " << odd_wheels_added << " factors for odd wheel constraints\n";
          }
          return odd_wheels_added + odd_3_wheels_added + triplets_added;
        }
      }
   }

protected:
   std::unordered_map<std::array<INDEX,4>, odd_3_wheel_factor_container*> odd_3_wheel_factors_;
   std::vector<std::tuple<INDEX,INDEX,INDEX,typename base_constructor::triplet_factor*>> triplet_vector_;
};

template<
   typename CUT_CONSTRUCTOR, // CRTP
   typename ODD_WHEEL_CONSTRUCTOR, 
   typename ODD_BICYCLE_3_WHEEL_FACTOR,
   typename ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0123,
   typename ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0124,
   typename ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0134,
   typename ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0234,
   typename ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_1234
   >
class odd_bicycle_wheel_constructor_base : public ODD_WHEEL_CONSTRUCTOR {
public:
   using base_constructor = ODD_WHEEL_CONSTRUCTOR;
   using FMC = typename base_constructor::FMC;
   using odd_bicycle_3_wheel_factor_container = ODD_BICYCLE_3_WHEEL_FACTOR;
   using odd_3_wheel_odd_bicycle_3_wheel_message_0123_container = ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0123;
   using odd_3_wheel_odd_bicycle_3_wheel_message_0124_container = ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0124;
   using odd_3_wheel_odd_bicycle_3_wheel_message_0134_container = ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0134;
   using odd_3_wheel_odd_bicycle_3_wheel_message_0234_container = ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_0234;
   using odd_3_wheel_odd_bicycle_3_wheel_message_1234_container = ODD_WHEEL_ODD_BICYCLE_3_WHEEL_MESSAGE_1234;

   using base_constructor::base_constructor;

   template<typename MSG_TYPE>
   auto* connect_odd_3_wheel_factor(odd_bicycle_3_wheel_factor_container* f, const std::array<INDEX,4> idx)
   {
      assert(idx[0] < idx[1] && idx[1] < idx[2] && idx[2] < idx[3]);
      if(!this->has_odd_3_wheel_factor(idx[0], idx[1], idx[2], idx[3])) {
         this->add_odd_3_wheel_factor(idx[0], idx[1], idx[2], idx[3]);
      }
      auto* o = this->get_odd_3_wheel_factor(idx[0], idx[1], idx[2], idx[3]);
      auto* m = new MSG_TYPE(o, f);
      this->lp_->AddMessage(m);
      return o;
   }

   auto* add_odd_bicycle_3_wheel(const INDEX i0, const INDEX i1, const INDEX i2, const INDEX i3, const INDEX i4)
   {
      std::array<INDEX,5> idx({i0,i1,i2,i3,i4});
      return add_odd_bicycle_3_wheel(idx); 
   }

   auto* add_odd_bicycle_3_wheel(const std::array<INDEX,5> idx) {
      assert(!has_odd_bicycle_3_wheel(idx));
      auto* f = new odd_bicycle_3_wheel_factor_container();
      this->lp_->AddFactor(f);
      odd_bicycle_3_wheel_factors_.insert(std::make_pair(idx, f));

      auto* o_0123 = connect_odd_3_wheel_factor<odd_3_wheel_odd_bicycle_3_wheel_message_0123_container>(f, {idx[0], idx[1], idx[2], idx[3]});
      this->lp_->AddFactorRelation(o_0123, f);
      connect_odd_3_wheel_factor<odd_3_wheel_odd_bicycle_3_wheel_message_0124_container>(f, {idx[0], idx[1], idx[2], idx[4]});
      connect_odd_3_wheel_factor<odd_3_wheel_odd_bicycle_3_wheel_message_0134_container>(f, {idx[0], idx[1], idx[3], idx[4]});
      connect_odd_3_wheel_factor<odd_3_wheel_odd_bicycle_3_wheel_message_0234_container>(f, {idx[0], idx[2], idx[3], idx[4]});
      auto* o_1234 = connect_odd_3_wheel_factor<odd_3_wheel_odd_bicycle_3_wheel_message_1234_container>(f, {idx[1], idx[2], idx[3], idx[4]});
      this->lp_->AddFactorRelation(o_1234, f);

      return f;
   }

   bool has_odd_bicycle_3_wheel(const std::array<INDEX,5>& idx) const
   {
      assert(idx[0] < idx[1] && idx[1] < idx[2] && idx[2] < idx[3] && idx[3] < idx[4]);
      return odd_bicycle_3_wheel_factors_.find(idx) != odd_bicycle_3_wheel_factors_.end(); 
   }

/*
   // ij is the axle, uv is the wheel edge
   // labelings with edge ij and uv cut and (i) iu and jv or (ii) iv and ju not cut must have smaller cost than other configurations
   REAL compute_edge_cost(const INDEX i, const INDEX j, const INDEX u, const INDEX v)
   {
      std::array<INDEX,4> idx{i,j,u,v};
      std::sort(idx.begin(), idx.end());
      // can it happen that odd bicycle wheel inequalities can tighten the polytope but odd wheel inequalities are not also violated? Only in this case the below construction makes sense
      if(!this->has_odd_3_wheel_factor(idx[0], idx[1], idx[2], idx[3])) { // create a fake odd 3 wheel factor and reparamaetrize all underlying triplets into it

        typename odd_3_wheel_factor = base_constructor::odd_3_wheel_factor_container::FactorType;
        typename triplet_odd_3_wheel_message_012 = base_constructor::triplet_odd_3_wheel_message_012::MessageType;
        typename triplet_odd_3_wheel_message_013 = base_constructor::triplet_odd_3_wheel_message_013::MessageType;
        typename triplet_odd_3_wheel_message_023 = base_constructor::triplet_odd_3_wheel_message_023::MessageType;
        typename triplet_odd_3_wheel_message_123 = base_constructor::triplet_odd_3_wheel_message_123::MessageType;

        odd_3_wheel_factor f;

        // possibly do not create new messages, but use static functions in these messages directly
        if(this->HasTripletFactor(idx[0],idx[1],idx[2])) {
          const auto& t = *(this->GetTripletFactor(idx[0], idx[1], idx[2])->GetFactor());
          triplet_odd_3_wheel_message_012 m;
          m.RepamRight(f,t);
        }

        if(this->HasTripletFactor(idx[0],idx[1],idx[3])) {
          const auto& t = *(this->GetTripletFactor(idx[0], idx[1], idx[3])->GetFactor());
          triplet_odd_3_wheel_message_013 m;
          m.RepamRight(f,t);
        }

        if(this->HasTripletFactor(idx[0],idx[2],idx[3])) {
          const auto& t = *(this->GetTripletFactor(idx[0], idx[2], idx[3])->GetFactor());
          triplet_odd_3_wheel_message_023 m;
          m.RepamRight(f,t);
        }

        if(this->HasTripletFactor(idx[1],idx[2],idx[3])) {
          const auto& t = *(this->GetTripletFactor(idx[1], idx[2], idx[3])->GetFactor());
          triplet_odd_3_wheel_message_123 m;
          m.RepamRight(f,t);
        }
        return edge_cost_from_3_wheel(i,j,u,v, f);

      } else {
        auto& f = *(this->get_odd_3_wheel_factor(idx[0], idx[1], idx[2], idx[3])->GetFactor());
        return compute_edge_cost_from_odd_3_wheel(i,j,u,v, f);
      }
   }

   // ij is the axle. compute minimum over labelings where edge ij and uv is on - minimum over labelings where above condition is not true
   template<typename ODD_3_WHEEL_FACTOR>
   REAL compute_edge_cost_from_odd_3_wheel(const INDEX i, const INDEX j, const INDEX u, const INDEX v, const ODD_3_WHEEL_FACTOR& f)
   {
      assert(i < j && u < v);
      std::array<INDEX,4> idx{i,j,u,v};
      std::sort(idx.begin(), idx.end());
      static_cast<CUT_CONSTRUCTOR*>(this)->compute_edge_cost_from_odd_3_wheel(i,j,u,v,f);
      REAL min_participating_labelings = std::numeric_limits<REAL>::infinity();
      REAL min_non_participating_labelings;
      // non participating labelings have always != 4 cut edges
      if(CUT_TYPE == cut_type::multicut) {
         min_non_participating_labelings = std::min({0.0, f[0], f[1], f[2], f[3], f[7], f[8], f[9], f[10], f[11], f[12], f[13]});
      } else if(CUT_TYPE == cut_type::maxcut) {`
         min_non_participating_labelings = std::min({0.0, f[0], f[1], f[2], f[3]}); 
      } else {
         assert(false);
      }
      // for participating edges, the axle and wheel edge must be cut
      if(idx[0] == i && idx[1] == j) { // first and sixth edge
         assert(idx[2] == u && idx[3] == v);
         min_participating_labelings = std::min(f[5], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[4]);
      } else if(idx[0] == i && idx[2] == j) { // second and fifth
         assert(idx[1] == u && idx[3] == v);
         min_participating_labelings = std::min(f[4], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[5]);
      } else if(idx[0] == i && idx[3] == j) { // third and fourth
         assert(idx[1] == u && idx[2] == v);
         min_participating_labelings = std::min(f[4], f[5]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[6]);
      } else if(idx[1] == i && idx[2] == j) { // fourth and third
         assert(idx[0] == u && idx[3] == v);
         min_participating_labelings = std::min(f[4], f[5]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[6]);
      } else if(idx[1] == i && idx[3] == j) { // fifth and second
         assert(idx[0] == u && idx[2] == v);
         min_participating_labelings = std::min(f[4], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[5]);
      } else if(idx[2] == i && idx[3] == j) { // sixth and first
         assert(idx[0] == u && idx[1] == v);
         min_participating_labelings = std::min(f[5], f[6]);
         min_non_participating_labelings = std::min(min_non_participating_labelings,f[4]);
      } else {
         assert(false);
      }
      //assert(false);

      return min_participating_labelings - min_non_participating_labelings; 
   }
*/

   using triangle_intersection_type = std::tuple<INDEX,INDEX, typename ODD_WHEEL_CONSTRUCTOR::triplet_factor*, typename ODD_WHEEL_CONSTRUCTOR::triplet_factor*>;

   void compute_triangles( // triangle are pyramids, though, search for better name
         const INDEX i, const INDEX j, const REAL minTh, 
         const typename ODD_WHEEL_CONSTRUCTOR::triplet_connections& connected_triplets,
         std::vector<triangle_intersection_type>& common_edges,
         std::unordered_map<INDEX,INDEX>& origToCompressedNode, 
         std::vector<INDEX>& compressedToOrigNode, 
         std::vector<std::tuple<INDEX,INDEX,REAL>>& compressedEdges)
   {
      assert(i < j);
      origToCompressedNode.clear();
      compressedToOrigNode.clear();
      compressedEdges.clear();

      auto merge = [](const auto a, const auto b) -> triangle_intersection_type { 
         assert(a.nodes[0] == b.nodes[0] && a.nodes[1] == b.nodes[1]);
         return std::make_tuple(a.nodes[0], a.nodes[1], a.f, b.f); 
      };

      // find all triangles ijk

      // find all edges uv such that there exist edge triplets iuv and juv. 
      // this is done by sorting all triplets which have node i and node j, and intersecting the set
      auto intersects_iter_end = set_intersection_merge(
            connected_triplets[i].begin(), connected_triplets[i].end(),
            connected_triplets[j].begin(), connected_triplets[j].end(),
            common_edges.begin(), [](const auto& a, const auto& b) { return a.operator<(b); }, merge);

      for(auto n=common_edges.begin(); n != intersects_iter_end; ++n) {
         const INDEX u = std::get<0>(*n);
         const INDEX v = std::get<1>(*n);
         assert(u < v);
         const auto& iuv = std::get<2>(*n)->GetFactor();
         const auto& juv = std::get<3>(*n)->GetFactor();

         const REAL dual_increase = this->compute_edge_cost(i,j,u,v);

         if(dual_increase >= minTh) { // add edge uv to bipartite graph

            if(origToCompressedNode.find(u) == origToCompressedNode.end()) {
               origToCompressedNode.insert(std::make_pair(u, origToCompressedNode.size()));
               compressedToOrigNode.push_back(u);
            }
            if(origToCompressedNode.find(v) == origToCompressedNode.end()) {
               origToCompressedNode.insert(std::make_pair(v, origToCompressedNode.size()));
               compressedToOrigNode.push_back(v);
            }
            const INDEX uc = origToCompressedNode[u];
            const INDEX vc = origToCompressedNode[v];
            assert(uc != vc);
            compressedEdges.push_back(std::make_tuple(uc,vc, dual_increase));
         }
      } 
   }

   struct bicycle_candidate {
      std::array<INDEX,5> idx;
      REAL cost;
   };

   template<typename ITERATOR>
   void triangulate_odd_bicycle_wheel(const INDEX i, const INDEX j, const REAL cost, ITERATOR path_begin, ITERATOR path_end, std::vector<bicycle_candidate>& candidates)
   {
      assert(i < j);
      assert(std::distance(path_begin, path_end) >= 3);
      this->cycle_normal_form(path_begin, path_end);
      const INDEX first_node = *path_begin;
      for(auto it=path_begin+1; it+1!=path_end; ++it) {
         std::array<INDEX,5> nodes({i,j,first_node, *it, *(it+1)});
         std::sort(nodes.begin(), nodes.end());
         assert(HasUniqueValues(nodes));
         candidates.push_back({nodes, cost});
      }
   }


   INDEX find_violated_odd_bicycle_wheels(const INDEX max_factors_to_add)
   {
      if(this->number_of_edges() > 2) {
         // preprocessing: sort triplets for fast intersection later
         auto connected_triplets = this->compute_connected_triplets();

         std::vector<std::tuple<INDEX,REAL>> threshold(this->unaryFactorsVector_.size()); // edge number and threshold
         std::vector<bicycle_candidate> odd_bicycle_candidates;
         // given a cut axle edge (negative cost), find cut wheel edges such that among the four spokes exactly two are cut and two are connected.
     
#pragma omp parallel
         {
            using intersection_type = std::tuple<INDEX,INDEX, typename ODD_WHEEL_CONSTRUCTOR::TripletFactorContainer*, typename ODD_WHEEL_CONSTRUCTOR::TripletFactorContainer*>;
            std::vector<intersection_type> common_edges(this->number_of_edges()); // possibly this is a bit large!
            auto merge = [](const auto a, const auto b) -> intersection_type { 
               assert(std::get<0>(a) == std::get<0>(b) && std::get<1>(a) == std::get<1>(b));
               return std::make_tuple(std::get<0>(a), std::get<1>(a), std::get<2>(a), std::get<2>(b)); 
            };
            std::vector<bicycle_candidate> odd_bicycle_candidates_local;

            std::unordered_map<INDEX,INDEX> origToCompressedNode;
            std::vector<INDEX> compressedToOrigNode;
            std::vector<std::tuple<INDEX,INDEX,REAL>> compressedEdges;

#pragma omp for schedule(guided) nowait
            for(INDEX e=0; e<this->unaryFactorsVector_.size(); ++e) {

               // edge ij will be treated as axle of odd bicycle wheel
               const INDEX i = std::get<0>(this->unaryFactorsVector_[e])[0];
               const INDEX j = std::get<0>(this->unaryFactorsVector_[e])[1];
               const REAL cost_ij = std::get<1>(this->unaryFactorsVector_[e])->GetFactor()->operator[](0); 
               if(cost_ij < -eps) {
                  //origToCompressedNode.clear();
                  //compressedToOrigNode.clear();
                  //compressedEdges.clear(); 

                  compute_triangles(i, j, eps, connected_triplets, common_edges, origToCompressedNode, compressedToOrigNode, compressedEdges); 
                  const REAL th = this->compute_odd_cycle_threshold(origToCompressedNode, compressedToOrigNode, compressedEdges);

                  if(th > eps) {
                     auto path = this->compute_path_in_bipartite_graph(origToCompressedNode, compressedToOrigNode, compressedEdges, th);
                     triangulate_odd_bicycle_wheel(i,j, th, path.begin(), path.end(), odd_bicycle_candidates_local);
                  }
               } 
            }
#pragma omp critical
            odd_bicycle_candidates.insert(odd_bicycle_candidates.end(), odd_bicycle_candidates_local.begin(), odd_bicycle_candidates_local.end());
         }

         std::sort(odd_bicycle_candidates.begin(), odd_bicycle_candidates.end(), [](const auto& a, const auto& b) { return a.cost > b.cost; });
         INDEX no_factors_added = 0;
         for(INDEX i=0; i<odd_bicycle_candidates.size(); ++i) {
            if(!has_odd_bicycle_3_wheel( odd_bicycle_candidates[i].idx )) {
               add_odd_bicycle_3_wheel( odd_bicycle_candidates[i].idx );
               ++no_factors_added;
               if(no_factors_added >= max_factors_to_add) {
                  break;
               }
            } 
         }
         return no_factors_added;
      } else {
         return 0;
      }
   }

   INDEX Tighten(const INDEX max_factors_to_add)
   {
      const INDEX factors_added = ODD_WHEEL_CONSTRUCTOR::Tighten(max_factors_to_add);
      if(factors_added < 0.1*max_factors_to_add) {
         if(diagnostics()) {
            std::cout << "search for odd bicycle wheels\n";
         }
         const INDEX odd_bicycles_added = static_cast<CUT_CONSTRUCTOR*>(this)->find_violated_odd_bicycle_wheels(max_factors_to_add - factors_added);
         if(diagnostics()) {
            std::cout << "added " << odd_bicycles_added << " odd bicycle wheels\n";
         }
         return odd_bicycles_added + factors_added;
      } else {
         return factors_added;
      }
   }
private:

   std::unordered_map<std::array<INDEX,5>, odd_bicycle_3_wheel_factor_container*> odd_bicycle_3_wheel_factors_;
};


template< class BASE_CONSTRUCTOR, // CRTP
	class CUT_CONSTRUCTOR,
	typename LIFTED_CUT_FACTOR,
	typename CUT_EDGE_LIFTED_FACTOR_MSG, typename LIFTED_EDGE_LIFTED_FACTOR_MSG
      >
      class lifted_constructor : public CUT_CONSTRUCTOR {
         public:
            using base_constructor = BASE_CONSTRUCTOR;
	    using class_type = lifted_constructor<BASE_CONSTRUCTOR, CUT_CONSTRUCTOR, LIFTED_CUT_FACTOR, CUT_EDGE_LIFTED_FACTOR_MSG, LIFTED_EDGE_LIFTED_FACTOR_MSG>;
            using LiftedCutFactorContainer = LIFTED_CUT_FACTOR;
            using CutEdgeLiftedFactorMessageContainer = CUT_EDGE_LIFTED_FACTOR_MSG;
            using CutEdgeLiftedFactorMessage = typename CutEdgeLiftedFactorMessageContainer::MessageType;
            using LiftedEdgeLiftedFactorMessageContainer = LIFTED_EDGE_LIFTED_FACTOR_MSG;
            using LiftedEdgeLiftedFactorMessage = typename LiftedEdgeLiftedFactorMessageContainer::MessageType;

            // do zrobienia: use this everywhere instead of std::array<INDEX,2>
            struct Edge : public std::array<INDEX,2> {
               Edge(const INDEX i, const INDEX j) : std::array<INDEX,2>({std::min(i,j), std::max(i,j)}) {}
            };
            using CutId = std::vector<Edge>;

            template<typename SOLVER>
               lifted_constructor(SOLVER& pd) : CUT_CONSTRUCTOR(pd) {}

            virtual typename CUT_CONSTRUCTOR::edge_factor* AddUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost)
            {
               assert(i1<i2);
               auto* f = CUT_CONSTRUCTOR::AddUnaryFactor(i1,i2,cost);
               if(!addingTighteningEdges) {
                  baseEdges_.push_back({i1, i2, f});
               } else { // all auxiliary edges may be regarded as lifted ones
                  AddLiftedUnaryFactor(i1, i2, cost);
               }
               return f;
            }
            typename CUT_CONSTRUCTOR::edge_factor* AddLiftedUnaryFactor(const INDEX i1, const INDEX i2, const REAL cost)
            {
               auto* f = CUT_CONSTRUCTOR::AddUnaryFactor(i1, i2, cost);
               liftedEdges_.push_back({i1,i2,f}); 
               return f;
            }

            bool HasCutFactor(const CutId& cut) 
            {
               assert(std::is_sorted(cut.begin(), cut.end()));
               return liftedFactors_.find(cut) != liftedFactors_.end();
            }

            bool HasLiftedEdgeInCutFactor(const CutId& cut, const INDEX i1, const INDEX i2)
            {
               assert(i1<i2);
               assert(HasCutFactor(cut));
               const auto& edgeList = liftedFactors_[cut].second;
               // speed this up by holding edge list sorted
               return std::find(edgeList.begin(),edgeList.end(),Edge({i1,i2})) != edgeList.end();
            }

            // do zrobienia: provide AddCutFactor(const CutId& cut, const INDEX i1, const INDEX i2) as well
            LiftedCutFactorContainer* AddCutFactor(const CutId& cut)
            {
               assert(!HasCutFactor(cut));
               //std::cout << "Add cut with edges ";
               //for(auto i : cut) { std::cout << "(" << std::get<0>(i) << "," << std::get<1>(i) << ");"; } std::cout << "\n";
               auto* f = new LiftedCutFactorContainer(cut.size());
               CUT_CONSTRUCTOR::lp_->AddFactor(f);
               // connect the cut edges
               for(INDEX e=0; e<cut.size(); ++e) {
                  auto* unaryFactor = CUT_CONSTRUCTOR::GetUnaryFactor(cut[e][0],cut[e][1]);
                  auto* m = new CutEdgeLiftedFactorMessageContainer(CutEdgeLiftedFactorMessage(e),unaryFactor,f);
                  CUT_CONSTRUCTOR::lp_->AddMessage(m);
               }
               liftedFactors_.insert(std::make_pair(cut,std::make_pair(f,std::vector<Edge>())));
               return f;
            }

            LiftedCutFactorContainer* GetCutFactor(const CutId& cut)
            {
               assert(HasCutFactor(cut));
               return liftedFactors_[cut].first;
            }

            void AddLiftedEdge(const CutId& cut, const INDEX i1, const INDEX i2)
            {
               assert(!HasLiftedEdgeInCutFactor(cut,i1,i2));
               assert(HasCutFactor(cut));
               assert(i1<i2);
               //std::cout << "Add lifted edge (" << i1 << "," << i2 << ") to cut\n";
               auto& c = liftedFactors_[cut];
               auto* f = c.first;
               auto& edgeList = c.second;
               auto* unaryFactor = CUT_CONSTRUCTOR::GetUnaryFactor(i1,i2);
               f->GetFactor()->IncreaseLifted();
               auto* m = new LiftedEdgeLiftedFactorMessageContainer(LiftedEdgeLiftedFactorMessage(edgeList.size() + cut.size()), unaryFactor, f);
               CUT_CONSTRUCTOR::lp_->AddMessage(m);
               c.second.push_back(Edge({i1,i2}));
            }



            INDEX Tighten(const INDEX max_factors_to_add)
            {
               const bool prevMode = addingTighteningEdges;
               addingTighteningEdges = true;
               assert(max_factors_to_add > 5); //otherwise the below arrangement makes little sense.
               const INDEX noBaseConstraints = CUT_CONSTRUCTOR::Tighten(std::ceil(0.8*max_factors_to_add));
               //return noBaseConstraints;
               INDEX noLiftingConstraints = 0;
               if(debug()) {
                  std::cout << "number of cut constraints = " << liftedFactors_.size() << "\n";
               }
               if(noBaseConstraints < max_factors_to_add) {
                  REAL th = FindViolatedCutsThreshold(max_factors_to_add - noBaseConstraints);
                  if(th >= 0.0) {
                     noLiftingConstraints = FindViolatedCuts(th, max_factors_to_add - noBaseConstraints);
                     if(diagnostics()) {
                        std::cout << "added " << noLiftingConstraints << " lifted cut factors.\n";
                     }
                  }
               }
               addingTighteningEdges = prevMode;
               return noBaseConstraints + noLiftingConstraints;
            }

            REAL FindViolatedCutsThreshold(const INDEX max_triplets_to_add)
            {
               // make one function to reuse allocated datastructures.
               UnionFind uf(this->noNodes_);
               std::vector<std::tuple<INDEX,INDEX,REAL>> edges;
               edges.reserve(baseEdges_.size());
               for(const auto& e : baseEdges_) {
                  const INDEX i = e.i;
                  const INDEX j = e.j;
                  const REAL weight = e.weight();
                  if(weight > 0) {
                     uf.merge(i,j);
                  } else {
                     edges.push_back(std::make_tuple(i,j,weight));
                  }
               }
               std::sort(edges.begin(),edges.end(), [] (auto& a, auto& b)->bool { return std::get<2>(a) > std::get<2>(b); });
               std::vector<std::list<std::tuple<INDEX,REAL>>> liftedEdges(this->noNodes_);
               for(auto& e : liftedEdges_) {
                  const INDEX i = e.i;
                  const INDEX j = e.j;
                  const REAL weight = e.weight();
                  if(weight > 0.0) {
                     liftedEdges[i].push_front(std::make_tuple(j,weight));
                     liftedEdges[j].push_front(std::make_tuple(i,weight));
                  }
               }
               auto edge_sort = [] (const std::tuple<INDEX,REAL>& a, const std::tuple<INDEX,REAL>& b)->bool { return std::get<0>(a) < std::get<0>(b); };
               for(INDEX i=0; i<liftedEdges.size(); ++i) {
                  liftedEdges[i].sort(edge_sort);
               }
               REAL prevWeight = 0.0;
               REAL maxTh = -std::numeric_limits<REAL>::infinity();
               for(INDEX e=0; e<edges.size(); ++e) {
                  const INDEX i = std::get<0>(edges[e]);
                  const INDEX j = std::get<1>(edges[e]);
                  const REAL weight = std::get<2>(edges[e]);
                  if(uf.find(i) != uf.find(j)) {
                     uf.merge(i,j);
                     const INDEX c = uf.find(i);
                     // (i) merge edges emanating from same connected component
                     if(i != c) {
                        liftedEdges[c].merge(liftedEdges[i],edge_sort);
                     } 
                     if(j != c) {
                        liftedEdges[c].merge(liftedEdges[j],edge_sort);
                     }
                     std::transform(liftedEdges[c].begin(), liftedEdges[c].end(), liftedEdges[c].begin(), 
                           [&uf,&maxTh,weight] (std::tuple<INDEX,REAL> a)->std::tuple<INDEX,REAL> {
                           const INDEX cc = uf.find(std::get<0>(a));
                           if(cc != std::get<0>(a)) {
                           const REAL th = std::min(-weight, std::get<1>(a));
                           maxTh = std::max(th,maxTh); // do zrobienia: correct?
                           }
                           std::get<0>(a) = cc;
                           return a; 
                           });
                     // do zrobienia: unnecessary sort here: has been sorted more or less after merge, only parallel edges need to be sorted
                     liftedEdges[c].sort([] (const std::tuple<INDEX,REAL>& a, const std::tuple<INDEX,REAL>& b)->bool { 
                           if(std::get<0>(a) != std::get<0>(b)) {
                           return std::get<0>(a) < std::get<0>(b); 
                           }
                           return std::get<1>(a) > std::get<1>(b); // this ensures that remove removes copies with smaller weight. Thus, the largest one only remains.
                           });
                     // (ii) remove lifted edges that have been collapsed. do zrobienia: record weight of those
                     liftedEdges[c].remove_if([&uf,c,maxTh] (const auto& a)->bool { 
                           if(std::get<1>(a) < maxTh) return true; // we need not take this edge into account anymore, as it would lead to smaller minimum dual improvement.
                           // first check whether lifted edge belonged to different connected components before the last merge operation. If yes, update threshold
                           const INDEX cc = uf.find(std::get<0>(a));
                           return uf.find(std::get<0>(a)) == c; 
                           });
                     // (iii) take maximal representative of edges that are parallel now
                     liftedEdges[c].unique([] (const auto& a, const auto& b)->bool { return std::get<0>(a) == std::get<0>(b); }); 
                  }

               }

               //if(maxTh != -std::numeric_limits<REAL>::infinity()) {
               //   std::cout << "\nnon-trivial cut factor with weight = " << maxTh << "\n\n";
               //}
               return 0.1*maxTh;
            }

            INDEX FindViolatedCuts(const INDEX minDualIncrease, const INDEX noConstraints)
            {
               // update weight of base edges
               //for(auto& e : baseEdges_) {
               //   e.w = CUT_CONSTRUCTOR::GetUnaryFactor(e.i,e.j)->operator[](0);
               //}
               //std::sort(baseEdges_.begin(), baseEdges_.end(), [](const Edge& e1, const Edge& e2) { return e1.weight() < e2.weight(); });
               UnionFind uf(CUT_CONSTRUCTOR::noNodes_);
               for(const auto& e : baseEdges_) {
                  if(e.weight() >= -minDualIncrease) {
                     uf.merge(e.i,e.j);
                  }
               }

               // build reduced graph with connected components as nodes and edge weights as number of edges with weight < -minDualIncrease

               // union find indices are not contiguous. Make them so, to use them as identifiers for connected components
               std::map<INDEX,INDEX> ufIndexToContiguous;
               for(INDEX i=0; i<CUT_CONSTRUCTOR::noNodes_; ++i) {
                  const INDEX ufIndex = uf.find(i);
                  if(ufIndexToContiguous.find(ufIndex) == ufIndexToContiguous.end()) {
                     ufIndexToContiguous.insert(std::make_pair(ufIndex,ufIndexToContiguous.size()));
                  }
               }
               const INDEX ccNodes = ufIndexToContiguous.size();

               std::map<INDEX,INDEX> origToCompressedNode; // union find index to compressed node indices // do zrobienia: use hash map
               for(INDEX i=0; i<CUT_CONSTRUCTOR::noNodes_; ++i) {
                  const INDEX ufIndex = uf.find(i);
                  const INDEX collapsedIndex = ufIndexToContiguous[ufIndex];
                  origToCompressedNode.insert(std::make_pair(i,collapsedIndex));
               }

               INDEX ccEdges = 0;
               std::map<Edge,std::vector<Edge>> ccToBaseEdges;
               for(const auto& e : baseEdges_) {
                  if(!uf.connected(e.i,e.j)) {
                     const INDEX i = ufIndexToContiguous[uf.find(e.i)];
                     const INDEX j = ufIndexToContiguous[uf.find(e.j)];
                     //if(ccEdgeCap.find({i,j}) == cc.EdgeCap.end()) {
                     //   ccEdgeCap.insert(std::make_pair(std::array<INDEX,2>(i,j),1));
                     //}
                     //ccEdgeCap[std::array<INDEX,2>(i,j)]++;
                     if(ccToBaseEdges.find(Edge(i,j)) == ccToBaseEdges.end()) {
                        ++ccEdges;
                        ccToBaseEdges.insert(std::make_pair(Edge(i,j),std::vector<Edge>()));
                     }
                     ccToBaseEdges[Edge(i,j)].push_back(Edge(e.i,e.j));
                  }
               }
               BKMaxFlow::Graph<int,int,int> maxFlow(ccNodes, ccEdges); 
               maxFlow.add_node(ccNodes);
               for(auto& e : ccToBaseEdges) {
                  const INDEX i = e.first.operator[](0);
                  const INDEX j = e.first.operator[](1);
                  const INDEX cap = e.second.size();
                  maxFlow.add_edge(i,j,cap,cap);
               }

               // note: this can possibly be made faster by caching the weight
               std::sort(liftedEdges_.begin(), liftedEdges_.end(), [](const weighted_edge& e1, const weighted_edge& e2) { return e1.weight() > e2.weight(); });
               INDEX factorsAdded = 0;

               // note that currently possibly multiple max flow computations are performed, when lifted edges come from the same connected components. This is superfluous and searches could be remembered and reused.
               const int capacityMax = baseEdges_.size()+1;
               for(const auto& liftedEdge : liftedEdges_) {
                  if(factorsAdded >= noConstraints) { 
		     if(debug()) {
                       std::cout << "maximal number of constraints to add reached";
		     }
                     break; 
                  }
                  if(liftedEdge.weight() > minDualIncrease) {
                     const INDEX i = origToCompressedNode[liftedEdge.i];
                     const INDEX j = origToCompressedNode[liftedEdge.j];
                     if(!uf.connected(i,j)) {
                        // find minimum cut in unweighted graph containing only base edges with weight < eps
                        maxFlow.add_tweights(i,capacityMax,0);
                        maxFlow.add_tweights(j,0,capacityMax);
                        const INDEX noCutEdges = maxFlow.maxflow();
                        assert(noCutEdges > 0 && noCutEdges < baseEdges_.size()); // otherwise there is no path from i to j or all paths were collapsed
                        std::vector<Edge> minCut;
                        minCut.reserve(noCutEdges);
                        // now do a dfs from i on those vertices which are in the same segment as i. The edges from those to vertices in segment of j form a minimum cut.
                        std::stack<INDEX> q;
                        std::vector<bool> visited(ccNodes,false);
                        q.push(i);
                        while(!q.empty()) {
                           const INDEX v = q.top();
                           q.pop();
                           if(visited[v]) {
                              continue;
                           }
                           visited[v] = true;
                           auto* a = maxFlow.get_first_arc(v);
                           while(a != nullptr) {
                              int v_test, w; // do zrobienia: use proper type
                              maxFlow.get_arc_ends(a, v_test, w);
                              assert(v != w);
                              assert(v_test == v);
                              if(maxFlow.what_segment(INDEX(v)) != maxFlow.what_segment(INDEX(w))) {
                                 // expand all edges that were collapsed into (v,w) in the original graph
                                 //spdlog::get("logger")->info() << "edge in mincut: " << v << "," << w;
                                 for(const auto& e : ccToBaseEdges[Edge(INDEX(v),INDEX(w))]) {
                                    //spdlog::get("logger")->info() << " expanded edge : " << e[0] << "," << e[1];
                                    minCut.push_back(Edge(e[0],e[1]));
                                 }
                              } else if(!visited[INDEX(w)]) {
                                 q.push(INDEX(w));
                              }
                              a = a->next;
                           }
                        }
                        assert(minCut.size() == noCutEdges);
                        std::sort(minCut.begin(),minCut.end()); // unique form of cut

                        // check if minimum cut is already present
                        if(!HasCutFactor(minCut)) {
                           auto* f = AddCutFactor(minCut);
                           AddLiftedEdge(minCut,liftedEdge.i,liftedEdge.j);
                           ++factorsAdded;
                        } else {
                           auto* MinCutFactor = GetCutFactor(minCut);
                           if(!HasLiftedEdgeInCutFactor(minCut,liftedEdge.i,liftedEdge.j)) {
                              AddLiftedEdge(minCut,liftedEdge.i,liftedEdge.j);
                              ++factorsAdded;
                           }
                        }

                        // restore original terminal weights
                        maxFlow.add_tweights(i,-capacityMax,0);
                        maxFlow.add_tweights(j,0,-capacityMax);
                     }
                  }
               }

               return factorsAdded;
            }

            // check if all lifted edges are primally consistent by asserting that a path of zero values exists in the ground graph whenever lifted edge is zero
            bool CheckPrimalConsistency() const
            {
               const bool multicutConsistent = CUT_CONSTRUCTOR::CheckPrimalConsistency();
               if(!multicutConsistent) {
                  return false;
               }

               //collect connectivity information with union find w.r.t. base edges
               UnionFind uf(CUT_CONSTRUCTOR::noNodes_);
               for(const auto& e : baseEdges_) {
                  if(e.f->GetFactor()->primal()[0] == false) {
                     uf.merge(e.i,e.j);
                  }
               }
               for(const auto& e : liftedEdges_) {
                  if(e.f->GetFactor()->primal()[0] == false) {
                     if(!uf.connected(e.i,e.j)) {
                        return false;
                     }
                  }
               }
               return true;
            }

   static std::vector<char> call_round(class_type* class_ptr, std::vector<Edge> base_edges, std::vector<Edge> lifted_edges, std::vector<REAL> edge_values)
   {
	   return static_cast<BASE_CONSTRUCTOR*>(class_ptr)->round(base_edges, lifted_edges, edge_values);
   }

   void round()
   {
      if(diagnostics()) {
         std::cout << "compute lifted multicut primal with GAEC + KLj\n";
      }
      std::vector<REAL> edgeValues;
      edgeValues.reserve(this->baseEdges_.size() + this->liftedEdges_.size());

      if(debug()) {
         std::cout << "# base edges = " << baseEdges_.size() << ", # lifted edges = " << liftedEdges_.size() << "\n";
      }

      // do zrobienia: initalize the graph structures only once
      std::vector<Edge> base_edges;
      std::vector<Edge> lifted_edges;
      std::vector<REAL> edge_values;
      for(const auto& e : baseEdges_) {
	      base_edges.push_back({e.i, e.j});
	      lifted_edges.push_back({e.i,e.j});
	      edge_values.push_back(e.f->GetFactor()->operator[](0));
      }
      for(const auto& e : liftedEdges_) {
	      lifted_edges.push_back({e.i,e.j});
	      edge_values.push_back(e.f->GetFactor()->operator[](0));
      }

      primal_handle_ = std::async(std::launch::async, call_round, this, std::move(base_edges), std::move(lifted_edges), std::move(edge_values));
   }


   void ComputePrimal()
   {
      if(!primal_handle_.valid()) { 
         round();
         return;
      }

      const auto primal_state = primal_handle_.wait_for(std::chrono::seconds(0));

      if(primal_state == std::future_status::deferred) {
         assert(false); // this should not happen, we launch immediately!
         throw std::runtime_error("asynchronuous lifted rounding was deferred, but this should not happen");
      } else if(primal_state == std::future_status::ready) {

         if(debug()) {
            std::cout << "collect lifted rounding result\n";
         }
         auto labeling = primal_handle_.get();
         this->write_labeling_into_factors(labeling);

         if(debug()) {
            std::cout << "restart primal rounding\n";
         }
         round();

      } else {
         if(debug()) {
            std::cout << "lifted rounding is currently running.\n";
         }
      }
   } 

            
         private:
   decltype(std::async(std::launch::async, call_round, nullptr, std::vector<Edge>{}, std::vector<Edge>{}, std::vector<REAL>{})) primal_handle_;

            //struct Edge {INDEX i; INDEX j; REAL w;}; // replace by WeightedEdge
            struct weighted_edge {
               INDEX i; 
               INDEX j; 
               typename CUT_CONSTRUCTOR::edge_factor* f;
               REAL weight() const { return (*f->GetFactor())[0]; }
            };
            bool addingTighteningEdges = false; // controls whether edges are added to baseEdges_
            std::vector<weighted_edge> baseEdges_;
            std::vector<weighted_edge> liftedEdges_;

            std::vector<std::vector<INDEX>> cutEdgesLiftedFactors_;
            std::vector<std::vector<INDEX>> liftedEdgesLiftedFactors_;

            std::map<CutId,std::pair<LiftedCutFactorContainer*,std::vector<Edge>>> liftedFactors_;
      };



} // end namespace LP_MP

#endif // LP_MP_CUT_CONSTRUCTOR_BASE_HXX


