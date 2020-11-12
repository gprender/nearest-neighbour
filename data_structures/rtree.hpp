// rtree.hpp

#include <memory>
#include <vector>

#include "spatial.hpp" 

#pragma once

namespace spatial {

    template<typename T>
    class Rtree {
        private:
            struct Entry;
            struct LeafEntry;
            struct InternalEntry;

            class Node{
                public:
                    index_t load;  // # of POINTS in the node + all subnodes
                    int m;  // # of ENTRIES in the current node
                    std::vector<std::shared_ptr<Entry>> entries;

                    // Is there a better way to access Rtree class members?
                    Rtree& rtree;

                    Node(Rtree& rt);
                    // ~Node();
                    bool insert(index_t point_idx);
                    void split(int const branch_idx);
                    int choose_branch(Point const p) const;
                    bool is_leaf() const;
            };

            // Polymorphism to deal with the two types of tree entries
            struct Entry { 
                Rectangle bounding_box; 
                virtual ~Entry() { }
            };
            struct LeafEntry : Entry { index_t idx; };
            struct InternalEntry : Entry { std::shared_ptr<Node> node; };

            std::shared_ptr<InternalEntry> root_entry;
            std::vector<Datum<T>> data;

            void split_root();

        public:
            Rtree();
            // ~Rtree();
            void build(std::vector<T> const raw_data);
            void insert(Datum<T> new_datum);
            std::vector<T> query_knn(
                unsigned const k, coord_t const x, coord_t const y
            ) const;
            index_t get_load() const;
    };
}
