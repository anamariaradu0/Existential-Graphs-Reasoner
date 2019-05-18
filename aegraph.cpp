// Copyright 2019 Luca Istrate, Danut Matei
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <utility>
#include <cassert>
#include "./aegraph.h"

std::string strip(std::string s) {
    // deletes whitespace from the beginning and end of the string
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::pair<std::string, std::string> split_first(std::string s,
    char delimiter = ',') {
    // returns a pair that contains: <first_cut, rest_of_graph>

    int numOpen = 0;

    int len = s.size();
    for (int i = 0; i < len; i++) {
        char c = s[i];
        if (c == delimiter && numOpen == 0)
            return std::make_pair(
                    strip(s.substr(0, i)), strip(s.substr(i+1, len)));
        if (c == '[')
            numOpen += 1;
        if (c == ']')
            numOpen -= 1;
    }

    return {strip(s), std::string()};
}


std::vector<std::string> split_level(std::string s, char delimiter = ',') {
    // splits 's' into separate entities (atoms, subgraphs)

    std::vector<std::string> result;
    auto snd = s;
    while (true) {
        auto p = split_first(snd, delimiter);
        auto fst = p.first;
        snd = p.second;

        result.push_back(fst);

        if (snd.empty())
            return result;
    }
}


int AEGraph::num_subgraphs() const {
    return subgraphs.size();
}


int AEGraph::num_atoms() const {
    return atoms.size();
}


int AEGraph::size() const {
    return num_atoms() + num_subgraphs();
}


bool AEGraph::operator<(const AEGraph& other) const {
    return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
    return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
    return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
    // offers an easier way of accessing the nested graphs
    if (index < num_subgraphs()) {
        return subgraphs[index];
    }

    if (index < num_subgraphs() + num_atoms()) {
        return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
    }

    return AEGraph("()");
}

std::ostream& operator<<(std::ostream &out, const AEGraph &g) {
    out << g.repr();
    return out;
}

AEGraph::AEGraph(std::string representation) {
    // constructor that creates an AEGraph structure from a
    // serialized representation
    char left_sep = representation[0];
    char right_sep = representation[representation.size() - 1];

    assert((left_sep == '(' && right_sep == ')')
        || (left_sep == '[' && right_sep == ']'));

    // if the left separator is '(' then the AEGraph is the entire
    // sheet of assertion
    if (left_sep == '(') {
        is_SA = true;
    } else {
        is_SA = false;
    }

    // eliminate the first pair of [] or ()
    representation = representation.substr(1, representation.size() - 2);

    // split the graph into separate elements
    auto v = split_level(representation);
    // add them to the corresponding vector
    for (auto s : v) {
        if (s[0] != '[') {
            atoms.push_back(s);
        } else {
            subgraphs.push_back(AEGraph(s));
        }
    }

    // also internally sort the new graph
    this->sort();
}

std::string AEGraph::repr() const {
    // returns the serialized representation of the AEGraph

    std::string left, right;
    if (is_SA) {
        left = '(';
        right = ')';
    } else {
        left = '[';
        right = ']';
    }

    std::string result;
    for (auto subgraph : subgraphs) {
        result += subgraph.repr() + ", ";
    }

    int len = atoms.size();
    if (len != 0) {
        for (int i = 0; i < len - 1; i++) {
            result += atoms[i] + ", ";
        }
        result += atoms[len - 1];
    } else {
        if (subgraphs.size() != 0)
            return left + result.substr(0, result.size() - 2) + right;
    }

    return left + result + right;
}


void AEGraph::sort() {
    std::sort(atoms.begin(), atoms.end());

    for (auto& sg : subgraphs) {
        sg.sort();
    }

    std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
    // checks if an atom is in a graph
    if (find(atoms.begin(), atoms.end(), other) != atoms.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

bool AEGraph::contains(const AEGraph& other) const {
    // checks if a subgraph is in a graph
    if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const std::string other)
    const {
    // returns all paths in the tree that lead to an atom like <other>
    std::vector<std::vector<int>> paths;

    int len_atoms = num_atoms();
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_atoms; i++) {
        if (atoms[i] == other && size() > 1) {
            paths.push_back({i + len_subgraphs});
        }
    }

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i].contains(other)) {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const AEGraph& other)
    const {
    // returns all paths in the tree that lead to a subgraph like <other>
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i] == other && size() > 1) {
            paths.push_back({i});
        } else {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::possible_double_cuts() const {
    // returns all paths in the tree that lead to a possible double cut
    std::vector<std::vector<int>> paths_to_cuts;
    int len_subgraphs = num_subgraphs();
    
    // checks if a subgraph has exactly one child, that being a subgraph
    // as well
    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i].num_subgraphs() == 1 &&
                                            subgraphs[i].num_atoms() == 0) {
            paths_to_cuts.push_back({i});
        }
        auto r = subgraphs[i].possible_double_cuts();
        for (auto& v : r) {
            v.insert(v.begin(), i);
        }
        copy(r.begin(), r.end(), back_inserter(paths_to_cuts));
    }

    return paths_to_cuts;
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
    // erases a double cut from the tree
    AEGraph updated_graph(repr());
    int len_path = where.size();
    
    // recursively reconstructs the new subgraphs and atoms
    // vectors
    if (len_path == 1) {
        AEGraph temp = subgraphs[where[0]].subgraphs[0];

        int len_atoms = temp.num_atoms();
        for (int i = 0; i < len_atoms; i++) {
            updated_graph.atoms.push_back(temp.atoms[i]);
        }

        int len_subgraphs = temp.num_subgraphs();
        for (int i = 0; i < len_subgraphs; i++) {
            updated_graph.subgraphs.push_back(temp.subgraphs[i]);
        }

        updated_graph.subgraphs.erase(updated_graph.subgraphs.begin() +
                                    where[0]);
    } else {
        int index = where[0];
        where.erase(where.begin());
        updated_graph.subgraphs[index] =
                            updated_graph.subgraphs[index].double_cut(where);
    }

    return updated_graph;
}


std::vector<std::vector<int>> AEGraph::possible_erasures(int level) const {
    // returns all paths in the tree that lead to a possible erasure
    std::vector<std::vector<int>> paths_to_erasures;
    
    // checks if the level is even
    if (level % 2 && (size() > 1 || is_SA)) {
        for (int i = 0; i < size(); i++) {
            paths_to_erasures.push_back({i});
        }
    }

    // recursively checks every next level for the condition
    int len_subgraphs = num_subgraphs();
    for (int i = 0; i < len_subgraphs; i++) {
        auto r = subgraphs[i].possible_erasures(level + 1);
        for (auto& v : r) {
            v.insert(v.begin(), i);
        }
        copy(r.begin(), r.end(), back_inserter(paths_to_erasures));
    }

    return paths_to_erasures;
}


AEGraph AEGraph::erase(std::vector<int> where) const {
    // erases a subgraph/atom from the tree
    AEGraph updated_graph(repr());
    int len_path = where.size();
    int index = where[0];
    
    // recursively deletes the subgraph/atom from the respective position
    // by modifying the corresponding vectors
    if (len_path == 1) {
        int len_subgraphs = num_subgraphs();
        int len_atoms = num_atoms();

        if (index < len_subgraphs) {
            updated_graph.subgraphs.erase(updated_graph.subgraphs.begin() +
                                        index);
        } else if (index < len_subgraphs + len_atoms) {
            updated_graph.atoms.erase(updated_graph.atoms.begin() + index -
                                    len_subgraphs);
        }
    } else {
        where.erase(where.begin());
        updated_graph.subgraphs[index] =
                                    updated_graph.subgraphs[index].erase(where);
    }

    return updated_graph;
}

std::vector<std::vector<int>> AEGraph::possible_deiterations() const {
    // returns a vector of the paths to all possible deiterations
    std::vector<std::vector<int>> paths_to_deiters;
    int len_subgraphs = num_subgraphs();
    int len_atoms = num_atoms();
    
    // for each subgraph and atom, we chech whether it is found in one
    // of its children
    for (int i = 0; i < len_subgraphs + len_atoms; i++) {
        for (int j = 0; j < len_subgraphs; j++) {
            if (i != j) {
                if (i < len_subgraphs) {
                    // checks whether there are two identical subgraphs
                    // with the same parent
                    if (subgraphs[j] == subgraphs[i]) {
                        paths_to_deiters.push_back({j});
                    }
                    
                    // checks whether a subgraph contains another
                    // identical to it
                    if (subgraphs[j].contains(subgraphs[i])) {
                        auto r = subgraphs[j].get_paths_to(subgraphs[i]);
                        for (auto& v : r) {
                            v.insert(v.begin(), j);
                        }
                        copy(r.begin(), r.end(),
                                            back_inserter(paths_to_deiters));
                    }
                } else {
                    // checks for atoms
                    if (subgraphs[j].contains(atoms[i - len_subgraphs])) {
                        auto r = subgraphs[j].get_paths_to(atoms[i -
                                                            len_subgraphs]);
                        for (auto& v : r) {
                            v.insert(v.begin(), j);
                        }
                        copy(r.begin(), r.end(),
                                            back_inserter(paths_to_deiters));
                    }
                }
            }
        }
    
        // reccursive call of function
        if (i < len_subgraphs) {
            auto r = subgraphs[i].possible_deiterations();
            for (auto& v : r) {
                v.insert(v.begin(), i);
            }
            copy(r.begin(), r.end(), back_inserter(paths_to_deiters));
        }
    }
    
    // checks whether there are two identical atoms with the same parent
    for (int i = 0; i < len_atoms; i++) {
        for (int j = 0; j < len_atoms; j++) {
            if (i != j && atoms[i] == atoms[j]) {
                paths_to_deiters.push_back({j});
            }
        }
    }
    
    // eliminating duplicates
    std::sort(paths_to_deiters.begin(), paths_to_deiters.end());
    int no_paths = paths_to_deiters.size();
    for (int i = 2; i < no_paths; i++) {
        if (paths_to_deiters[i] == paths_to_deiters[i - 1]) {
            paths_to_deiters.erase(paths_to_deiters.begin() + i);
        }
    }

    return paths_to_deiters;
}

AEGraph AEGraph::deiterate(std::vector<int> where) const {
    // elimination deiterations using the erase method
    return erase(where);
}

