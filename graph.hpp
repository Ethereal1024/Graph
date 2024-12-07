#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <bits/stdc++.h>

const int INF = 1e9;

namespace COLOR {
    const int WHITE = -1;
    const int GRAY = 0;
    const int RED = 1;
    const int ORANGE = 2;
    const int YELLOW = 3;
    const int GREEN = 4;
    const int CYAN = 5;
    const int BLUE = 6;
    const int PURPLE = 7;
    const int BLACK = 8;
}  // namespace Color

struct Edge;
enum EdgeType { EE,
                OE,
                OO };

typedef struct Node {
    std::string label_;
    int color_ = -1;
    std::vector<Edge*> connections_;

    int degree();

    void print_node();

    std::vector<Node*> get_neighbors();
} Node;

typedef struct Edge {
    int length;
    int color_ = -1;
    std::pair<Node*, Node*> endpoints_;

    Node* another(const std::string& label);

    EdgeType edge_type();

    std::vector<Edge*> get_neighbors();

    void print_edge();
} Edge;

typedef struct Path {
    int dist_;
    std::vector<std::string> track_;

    void print_path();
} Path;

typedef struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const;
} pair_hash;

typedef std::pair<std::string, std::string> StrPair;

class Stainer;

class Graph {
public:
    Graph() {}

    Graph(const std::vector<Node>& nodes, const std::vector<Edge>& edges);

    Graph(const std::vector<Node>& nodes);

    Graph(const int& num);

    Graph(const std::vector<std::string>& labels);

    void add_node(const std::string& label);

    void add_node(const std::vector<std::string>& labels);

    void delete_node(const std::string& label);

    void connect(const std::string& label1, const std::string& label2, const int& length);

    void disconnect(const std::string& label1, const std::string& label2);

    void display_nodes();

    void display_edges();

    std::unordered_map<std::string, Node> get_nodes() const;

    std::unordered_map<StrPair, Edge, pair_hash> get_edges() const;

    std::vector<std::string> get_nodes_labels() const;

    Edge* edge(const std::string& label1, const std::string& label2);

    void operator=(const Graph& srcGraph);

private:
    friend class Optimizer;
    std::unordered_map<std::string, Node> nodes_;
    std::unordered_map<StrPair, Edge, pair_hash> edges_;

    void init_edges(const std::vector<Edge>& edges);
    void init_nodes(const std::vector<Node>& nodes);
};

class Optimizer {
public:
    explicit Optimizer(Graph& graph) : graph_(&graph) {};

    std::unordered_map<std::string, Path> dijkstra(std::string label);

    std::vector<Node*> odd_nodes();

    std::vector<int> get_epoch_dist_logger();

    bool is_euler();

    void bleach_node(std::string label);

    void bleach_all_nodes();  

    void DFS(Node* srcNode, std::vector<Path>& paths, const int& barrierColor, const int& targetColor);

    std::vector<Path> odd_node_connections(const std::string& label);

    std::string find_center(const std::vector<std::string>& labelNodes, int& distLog);

    void graph_k_means(std::vector<std::string>& startCenters, const int& epoch);

private:
    Graph* graph_;

    std::vector<int> epoch_dist_logger_;

    void DFS_step(Node* currentNode, int currentDistance, std::vector<std::string>& currentTrack, std::vector<Path>& paths, const int& barrierColor, const int& targetColor);

    void graph_k_means_step(std::vector<std::string>& centers, std::vector<int>* distLogs);
};

#endif