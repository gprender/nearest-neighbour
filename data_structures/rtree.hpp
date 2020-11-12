// rtree.hpp

#include <memory>
#include <vector>
#include <variant>

#include "spatial.hpp" 

#pragma once

namespace spatial {

    template<typename T>
    class Rtree {
        private:
            class Entry;

            class Node {
                public:
                    index_t load;
                    std::vector<Entry> entries;

                    Node();
                    ~Node();
                    bool insert(Datum<T> const& datum);
                    void split(int const branch_idx);
                    int choose_branch(Point const p) const;
                    void pick_seeds(std::vector<Entry> const& entry_choices);
                    void distribute(std::vector<Entry>& leftovers);
                    int pick_next(std::vector<Entry> const& leftovers) const;
                    bool is_leaf() const;
            };

            /**
             * For the sake of abstraction from the std::variant polymorphism,
             * we're keeping the member variables private, at the cost of using
             * shared_ptr rather than unique_ptr for the std::variant contents.
             */
            class Entry {
                private:
                    Rectangle _bounding_box;
                    std::variant<
                        std::shared_ptr<Datum<T>>, 
                        std::shared_ptr<Node>
                    > _contents;

                public:
                    Entry(Rectangle rect, std::shared_ptr<Datum<T>> datum): 
                        _bounding_box(rect),
                        _contents(datum)
                    { }

                    Entry(Rectangle rect, std::shared_ptr<Node> node): 
                        _bounding_box(rect),
                        _contents(node)
                    { }

                    Rectangle get_mbb() const { return _bounding_box; }

                    void set_mbb(Rectangle const rect) { 
                        _bounding_box = rect; 
                    }

                    std::shared_ptr<Datum<T>> get_datum() const { 
                        return std::get<std::shared_ptr<Datum<T>>>(
                            _contents
                        ); 
                    }

                    std::shared_ptr<Node> get_node() const { 
                        return std::get<std::shared_ptr<Node>>(
                            _contents
                        ); 
                    }
            };

            std::unique_ptr<Entry> root_entry;
            std::vector<Datum<T>> data;

            void split_root();

        public:
            Rtree();
            ~Rtree();
            void build(std::vector<T> const& raw_data);
            void insert(Datum<T> const& new_datum);
            std::vector<T> query_knn(
                unsigned const k, coord_t const x, coord_t const y
            ) const;
            index_t get_load() const;
    };
}
