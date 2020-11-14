// rtree.hpp

#include <queue>
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

                    bool is_leaf_entry() const {
                        return std::holds_alternative
                            <std::shared_ptr<Datum<T>>>(
                                _contents
                            );
                    }
            };

            /**
             * A thin wrapper around a std::priority_queue for Entries.
             * Helps to simplify the distance browsing code.
             */
            class EntryPQ {
                private:
                    struct EntryPQE {
                        Entry entry;
                        coord_t dist;
                    };

                    struct Closer {
                        bool operator()(EntryPQE const a, EntryPQE const b) {
                            return (a.dist > b.dist);
                        }
                    };

                    std::priority_queue<
                        EntryPQE, 
                        std::vector<EntryPQE>, 
                        Closer
                    > pq;

                    Point query_point;

                public:
                    EntryPQ(Point p):
                        query_point(p)
                    { }

                    void push(Entry e) {
                        pq.push(
                            (EntryPQE){e, distance(query_point, e.get_mbb())}
                        );
                    }

                    EntryPQE pop() {
                        auto const pqe = pq.top();
                        pq.pop();
                        return pqe;
                    }

                    EntryPQE peek() { return pq.top(); }

                    /**
                     * Push all of an entry's children onto the priority queue.
                     */
                    void expand(Entry e) {
                        for (auto const& child : e.get_node()->entries) {
                            push(child);
                        }
                    }
            };

            /**
             * Similar to above class, but with Datum elements. 
             */
            class DatumPQ {
                private:
                    struct DatumPQE {
                        Datum<T> datum;
                        coord_t dist;
                    };

                    struct Farther {
                        bool operator()(DatumPQE const a, DatumPQE const b) {
                            return (a.dist < b.dist);
                        }
                    };

                    std::priority_queue<
                            DatumPQE, 
                            std::vector<DatumPQE>, 
                            Farther
                    > pq;

                    Point query_point;

                public:
                    DatumPQ(Point p):
                        query_point(p)
                    { }
                
                    void push(Datum<T> d) {
                        pq.push(
                            (DatumPQE){d, distance(query_point, d.point)}
                        );
                    }

                    DatumPQE pop() {
                        auto const pqe = pq.top();
                        pq.pop();
                        return pqe;
                    }

                    DatumPQE peek() { return pq.top(); }

                    /**
                     * Conditionally push a datum onto the priority queue,
                     * if it's closer than the top (furthest) element.
                     */
                    void choose(Datum<T> d) {
                        coord_t const new_dist = distance(query_point, d.point);
                        if (pq.top().dist > new_dist) {
                            pq.pop();
                            pq.push((DatumPQE){d, new_dist});
                        }
                    }

                    unsigned size() { return pq.size(); }

                    bool empty() { return (pq.size() == 0); }
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
