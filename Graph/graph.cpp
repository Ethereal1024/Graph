#include "graph.hpp"

#include <bits/stdc++.h>

int Node::degree() {
    return connections_.size();
}

void Node::print_node() {
    std::cout << "label: " << label_ << ", degree: " << connections_.size() << ", color: " << color_ << std::endl;
}

std::vector<Node*> Node::get_neighbors() {
    std::vector<Node*> neighbors;
    for (auto connection : connections_) {
        neighbors.push_back(connection->another(label_));
    }
    return neighbors;
}

Node* Edge::another(const std::string& label) {
    if (endpoints_.first->label_ == label) {
        return endpoints_.second;
    } else if (endpoints_.second->label_ == label) {
        return endpoints_.first;
    } else {
        std::cerr << "Error: This edge doesn't have such endpoint." << std::endl;
        return nullptr;
    }
}

EdgeType Edge::edge_type() {
    return static_cast<EdgeType>(
        (endpoints_.first->degree() & 1) + (endpoints_.second->degree() & 1));
}

std::vector<Edge*> Edge::get_neighbors() {
    std::vector<Edge*> neighbors;
    for (auto connection : endpoints_.first->connections_) {
        neighbors.push_back(connection);
    }

    for (auto connection : endpoints_.second->connections_) {
        neighbors.push_back(connection);
    }
    return neighbors;
}

void Edge::print_edge() {
    const int side_l = 5;
    int str_l = std::to_string(length).size() + endpoints_.first->label_.size() + endpoints_.second->label_.size();
    std::cout << endpoints_.first->label_;
    int left_l = str_l >> 1;
    int right_l = str_l - left_l;
    for (int i = 0; i < side_l - left_l; i++) {
        std::cout << "-";
    }
    std::cout << length;
    for (int i = 0; i < side_l - right_l + 1; i++) {
        std::cout << "-";
    }
    std::cout << endpoints_.second->label_;
    std::cout << "  color: " << color_ << std::endl;
}

void Path::print_path() {
    std::cout << dist_ << " (";
    auto& track = track_;
    for (int i = 0; i < track.size(); i++) {
        std::cout << track[i];
        if (i != track.size() - 1) {
            std::cout << "->";
        }
    }
    std::cout << ")";
    std::cout << "\n";
}

template <class T1, class T2>
std::size_t pair_hash::operator()(const std::pair<T1, T2>& p) const {
    auto hash1 = std::hash<T1>{}(p.first);
    auto hash2 = std::hash<T2>{}(p.second);
    return hash1 ^ hash2;
}

Graph::Graph(const std::vector<Node>& nodes, const std::vector<Edge>& edges) {
    init_nodes(nodes);
    init_edges(edges);
}

Graph::Graph(const std::vector<Node>& nodes) {
    init_nodes(nodes);
}

Graph::Graph(const int& num) {
    if (num <= 26) {
        for (int i = 0; i < num; i++) {
            char tmp_label = 'A' + i;
            std::string label(1, tmp_label);
            nodes_[label] = {label, -1, {}};
        }
    } else {
        for (int i = 0; i < num; i++) {
            std::string tmp_label = "N" + std::to_string(i);
            nodes_[tmp_label] = {tmp_label, -1, {}};
        }
    }
}

Graph::Graph(const std::vector<std::string>& labels) {
    add_node(labels);
}

void Graph::add_node(const std::string& label) {
    if (nodes_.find(label) == nodes_.end()) {
        nodes_[label] = {label, -1, {}};
    } else {
        return;
    }
}

void Graph::add_node(const std::vector<std::string>& labels) {
    for (auto label : labels) {
        if (nodes_.find(label) == nodes_.end()) {
            nodes_[label] = {label, -1, {}};
        } else {
            std::cerr << "Error : node can not be created repeatedly." << std::endl;
        }
    }
}

void Graph::delete_node(const std::string& label) {
    if (nodes_.find(label) != nodes_.end()) {
        for (auto neighbor : nodes_[label].get_neighbors()) {
            disconnect(label, neighbor->label_);
        }
        nodes_.erase(label);
    } else {
        std::cerr << "Error : node does not exist." << std::endl;
    }
}

void Graph::connect(const std::string& label1, const std::string& label2, const int& length) {
    if (nodes_.find(label1) != nodes_.end() && nodes_.find(label2) != nodes_.end()) {
        StrPair label_key({label1, label2});
        if (edges_.find(label_key) == edges_.end()) {
            Node *node1 = &nodes_[label1], *node2 = &nodes_[label2];
            edges_[label_key] = Edge({length, -1, {node1, node2}});
            node1->connections_.push_back(&edges_[label_key]);
            node2->connections_.push_back(&edges_[label_key]);

            auto connection0 = node1->connections_[0];
        } else {
            return;
            // std::cerr << "Error: connection already exists." << std::endl;
        }
    } else {
        std::cerr << "Error: label not found." << std::endl;
    }
}

void Graph::disconnect(const std::string& label1, const std::string& label2) {
    if (nodes_.find(label1) != nodes_.end() && nodes_.find(label2) != nodes_.end()) {
        auto label_key = edges_.find(StrPair({label1, label2}));
        if (label_key == edges_.end()) {
            label_key = edges_.find(StrPair({label2, label1}));
        }

        if (label_key != edges_.end()) {
            Node *node1 = &nodes_[label1], *node2 = &nodes_[label2];

            auto remove_connection = [](Node* node, Node* other_node) {
                for (auto it = node->connections_.begin(); it != node->connections_.end();) {
                    if ((*it)->endpoints_.first == other_node || (*it)->endpoints_.second == other_node) {
                        it = node->connections_.erase(it);  // 使用erase并更新迭代器
                    } else {
                        ++it;
                    }
                }
            };

            remove_connection(node1, node2);
            remove_connection(node2, node1);
            edges_.erase(label_key);
        } else {
            std::cerr << "Error: connection never created." << std::endl;
        }
    } else {
        std::cerr << "Error: label not found." << std::endl;
    }
}

void Graph::display_nodes() {
    for (auto it : nodes_) {
        Node& node = it.second;
        node.print_node();
    }
}

void Graph::display_edges() {
    for (auto it : edges_) {
        Edge& edge = it.second;
        edge.print_edge();
    }
}

std::unordered_map<std::string, Node> Graph::get_nodes() const {
    return nodes_;
}

std::unordered_map<StrPair, Edge, pair_hash> Graph::get_edges() const {
    return edges_;
}

std::vector<std::string> Graph::get_nodes_labels() const {
    std::vector<std::string> result;
    for (auto elem : nodes_) {
        result.push_back(elem.second.label_);
    }
    return result;
}

Edge* Graph::edge(const std::string& label1, const std::string& label2) {
    if (nodes_.find(label1) != nodes_.end() && nodes_.find(label2) != nodes_.end()) {
        auto label_key = edges_.find(StrPair({label1, label2}));
        if (label_key == edges_.end()) {
            label_key = edges_.find(StrPair({label2, label1}));
        }

        if (label_key != edges_.end()) {
            return &label_key->second;
        } else {
            // std::cerr << "Error: Edge not found." << std::endl;
            return nullptr;
        }
    } else {
        // std::cerr << "Error: label not found." << std::endl;
        return nullptr;
    }
}

Edge* Graph::edge(const StrPair& labelPair) {
    return edge(labelPair.first, labelPair.second);
}

Edge* Graph::edge(const std::pair<Node*, Node*>& nodePair) {
    return edge(nodePair.first->label_, nodePair.second->label_);
}

void Graph::operator=(const Graph& srcGraph) {
    nodes_ = srcGraph.get_nodes();
    edges_ = srcGraph.get_edges();
}

void Graph::init_edges(const std::vector<Edge>& edges) {
    for (auto edge : edges) {
        std::pair<std::string, std::string> endpoints = {
            edge.endpoints_.first->label_, edge.endpoints_.second->label_};
        edges_[endpoints] = edge;
    }
}

void Graph::init_nodes(const std::vector<Node>& nodes) {
    for (auto node : nodes) {
        add_node(node.label_);
        nodes_[node.label_].connections_ = node.connections_;
    }
}

std::unordered_map<std::string, Path> Optimizer::dijkstra(std::string label) {
    auto it = graph_->nodes_.find(label);
    std::unordered_map<std::string, Path> min_distances;
    if (it != graph_->nodes_.end()) {
        std::unordered_map<std::string, Path> distances;
        for (auto nodeKV : graph_->nodes_) {
            distances[nodeKV.second.label_] = {INF, {}};
            min_distances[nodeKV.second.label_] = {INF, {label}};
        }

        Node* currentNode = &(it->second);
        int currentWeight = 0;
        std::vector<std::string> currentTrack = {label};
        distances.erase(label);
        min_distances[label] = {0, {label}};

        while (!distances.empty()) {
            for (auto neighbor : currentNode->get_neighbors()) {
                auto it = distances.find(neighbor->label_);
                int updateLen = graph_->edge(currentNode->label_, neighbor->label_)->length + currentWeight;
                std::vector<std::string> updateTrack = currentTrack;
                updateTrack.push_back(neighbor->label_);
                if (it != distances.end() && it->second.dist_ > updateLen) {
                    it->second.dist_ = updateLen;
                    it->second.track_ = updateTrack;
                }
            }

            int minLen = INF;
            std::string minLabel = "";
            for (auto element : distances) {
                if (element.second.dist_ < minLen) {
                    minLen = element.second.dist_;
                    minLabel = element.first;
                }
            }

            currentWeight = minLen;
            currentNode = &graph_->nodes_[minLabel];
            currentTrack = distances[minLabel].track_;

            if (minLen < min_distances[minLabel].dist_) {
                min_distances[minLabel].dist_ = minLen;
                min_distances[minLabel].track_ = distances[minLabel].track_;
            }
            distances.erase(minLabel);
        }

    } else {
        std::cerr << "Error : node does not exist." << std::endl;
    }
    return min_distances;
}

std::vector<Node*> Optimizer::odd_nodes() {
    std::vector<Node*> oddNodes;
    for (auto it : graph_->nodes_) {
        Node* node = &(graph_->nodes_[it.first]);
        if (node->degree() & 1) {
            oddNodes.push_back(node);
        }
    }
    return oddNodes;
}

std::vector<int> Optimizer::get_epoch_dist_logger() {
    return epoch_dist_logger_;
}

bool Optimizer::is_euler() {
    int n_odd = odd_nodes().size();
    return (n_odd == 0 || n_odd == 2);
}

void Optimizer::bleach_node(std::string label) {
    graph_->nodes_[label].color_ = -1;
}

void Optimizer::bleach_all_nodes() {
    for (auto element : graph_->nodes_) {
        element.second.color_ = -1;
    }
}

void Optimizer::bleach_all_edges() {
    for (auto element : graph_->edges_) {
        element.second.color_ = -1;
    }
}

void Optimizer::DFS_step(Node* currentNode, int currentDistance, std::vector<std::string>& currentTrack, std::vector<Path>& paths, const int& barrierColor, const int& targetColor) {
    if (graph_->nodes_[currentNode->label_].color_ == barrierColor) {
        return;
    }
    if (graph_->nodes_[currentNode->label_].color_ == targetColor) {
        paths.push_back({currentDistance, currentTrack});
        graph_->nodes_[currentNode->label_].color_ = barrierColor;
        return;
    }

    graph_->nodes_[currentNode->label_].color_ = barrierColor;

    auto neighbors = currentNode->get_neighbors();
    for (auto neighbor : neighbors) {
        int newEdgeLen = graph_->edge(currentNode->label_, neighbor->label_)->length;
        currentTrack.push_back(neighbor->label_);
        DFS_step(&graph_->nodes_[neighbor->label_], currentDistance + newEdgeLen, currentTrack, paths, barrierColor, targetColor);
        currentTrack.pop_back();
    }
}

void Optimizer::DFS(Node* srcNode, std::vector<Path>& paths, const int& barrierColor, const int& targetColor) {
    std::vector<std::string> track = {srcNode->label_};
    DFS_step(srcNode, 0, track, paths, barrierColor, targetColor);
}

std::vector<Path> Optimizer::odd_node_connections(const std::string& label) {
    std::vector<Node*> oddNodes = odd_nodes();
    auto it = oddNodes.begin();
    while (it != oddNodes.end()) {
        if ((*it)->label_ == label) {
            break;
        }
        ++it;
    }
    if (it == oddNodes.end()) {
        std::cerr << "Error: this node is not odd node." << std::endl;
        return {};
    }

    for (auto node : oddNodes) {
        node->color_ = COLOR::RED;
    }
    bleach_node(label);

    std::vector<Path> result;
    DFS(&graph_->nodes_[label], result, COLOR::GRAY, COLOR::RED);
    bleach_all_nodes();

    return result;
}

std::string Optimizer::find_center(const std::vector<std::string>& labelNodes, int& distLog) {
    std::unordered_map<std::string, int> sumDists;
    for (std::string label : labelNodes) {
        std::unordered_map<std::string, Path> nodeDij = dijkstra(label);
        int sumDist = 0;
        for (std::string other : labelNodes) {
            sumDist += nodeDij[other].dist_;
        }
        sumDists[label] = sumDist;
    }

    int minSumDist = INF;
    std::string result;
    for (auto element : sumDists) {
        if (element.second < minSumDist) {
            minSumDist = element.second;
            result = element.first;
        }
    }

    distLog = minSumDist;
    return result;
}

void Optimizer::graph_k_means_step(std::vector<std::string>& centers, std::vector<int>* distLogs) {
    std::vector<std::vector<std::string>> groups(centers.size());
    // 为每个节点重新分配质心
    for (auto& element : graph_->nodes_) {
        Node& node = element.second;
        std::unordered_map<std::string, Path> nodeDij = dijkstra(node.label_);
        int minDist = INF;
        int master_index;
        for (int i = 0; i < centers.size(); i++) {
            int dist = nodeDij[centers[i]].dist_;
            if (dist < minDist) {
                minDist = dist;
                master_index = i;
            }
        }
        node.color_ = master_index;
        groups[master_index].push_back(node.label_);
    }

    distLogs->clear();
    distLogs->resize(centers.size(), 0);
    for (int i = 0; i < centers.size(); i++) {
        int distLog;
        centers[i] = find_center(groups[i], distLog);
        (*distLogs)[i] = distLog;
    }
    // for (auto dist : *distLogs) {
    //     std::cout << dist << " ";
    // }
    // std::cout << "\n";
}

void Optimizer::graph_k_means(std::vector<std::string>& startCenters, const int& epochs) {
    epoch_dist_logger_.clear();
    for (int i = 0; i < epochs; i++) {
        std::vector<int> distLogs;
        graph_k_means_step(startCenters, &distLogs);
        int sumDist = 0;
        for (auto dist : distLogs) {
            sumDist += dist;
        }
        epoch_dist_logger_.push_back(sumDist);
    }
}

std::unordered_map<int, Graph> Optimizer::separate_graph_by_colors() {
    std::unordered_map<int, Graph> separatedGraphs;
    for (auto elem : graph_->nodes_) {
        auto& node = elem.second;
        if (separatedGraphs.find(node.color_) == separatedGraphs.end()) {
            separatedGraphs[node.color_] = Graph();
        }

        separatedGraphs[node.color_].add_node(node.label_);
        for (const auto& neighbor : node.get_neighbors()) {
            if (neighbor->color_ == node.color_) {
                separatedGraphs[node.color_].add_node(neighbor->label_);
                int length = graph_->edge(node.label_, neighbor->label_)->length;
                separatedGraphs[node.color_].connect(node.label_, neighbor->label_, length);
            }
        }
    }
    return separatedGraphs;
}

void Optimizer::node_enclosure(const std::string& srcNodeLabel, int numArea) {
    bleach_all_nodes();
    Node& srcNode = graph_->nodes_[srcNodeLabel];
    graph_->nodes_[srcNodeLabel].color_ = numArea;
    std::vector<std::shared_ptr<std::queue<Node*>>> targetQueues;
    std::vector<Node*> neighbors = srcNode.get_neighbors();
    int size = numArea < neighbors.size() ? numArea : neighbors.size();
    for (int i = 0; i < size; ++i) {
        auto currentQueue = std::make_shared<std::queue<Node*>>();
        currentQueue->push(neighbors[i]);
        std::cout << currentQueue->front()->label_ << std::endl;
        targetQueues.push_back(currentQueue);
    }
    std::vector<int> totalLengths(targetQueues.size(), 0);

    bool end = 0;
    while (!end) {
        end = 1;
        for (int i = 0; i < targetQueues.size(); ++i) {
            auto& currentQueue = targetQueues[i];
            if (!currentQueue->empty()) {
                end = 0;
                Node* currentNode = currentQueue->front();
                graph_->nodes_[currentNode->label_].color_ = i;
                for (auto& neighbor : currentNode->get_neighbors()) {
                    int postColor = graph_->nodes_[neighbor->label_].color_;
                    if (postColor == numArea)
                        continue;
                    if (postColor == -1 || totalLengths[postColor] > totalLengths[i]) {
                        currentQueue->push(neighbor);
                        totalLengths[i] += graph_->edge(currentNode->label_, neighbor->label_)->length;
                    }
                }
                currentQueue->pop();
            }
        }
    }
}

bool Optimizer::are_connected(const std::string& src, const std::string& dist, int edgeColor) {
    graph_->nodes_[dist].color_ = -2;
    std::vector<Path> paths;
    DFS(&graph_->nodes_[src], paths, -3, -2);
    if (paths.size() < 1) {
        return false;
    }
    bool connectionExist = 0;
    for (auto path : paths) {
        const std::vector<std::string>& track = path.track_;
        bool allColorCorrect = 1;
        for (int i = 0; i < path.track_.size() - 1; ++i) {
            if (graph_->edge(track[i], track[i + 1])->color_ != edgeColor) {
                allColorCorrect = 0;
            }
        }
        if (allColorCorrect) {
            connectionExist = 1;
        }
    }
    return connectionExist;
}

void Optimizer::edge_enclosure(const std::string& srcNodeLabel, int numArea) {
    bleach_all_edges();
    Node& srcNode = graph_->nodes_[srcNodeLabel];
    graph_->nodes_[srcNodeLabel].color_ = numArea;
    std::vector<std::shared_ptr<std::queue<Edge*>>> targetQueues;
    std::vector<Node*> neighbors = srcNode.get_neighbors();
    int size = numArea < neighbors.size() ? numArea : neighbors.size();
    std::vector<int> totalLengths(size, 0);
    // std::cout << "E" << std::endl;
    for (int i = 0; i < size; ++i) {
        auto currentQueue = std::make_shared<std::queue<Edge*>>();
        Edge* neighborEdge = graph_->edge(srcNodeLabel, neighbors[i]->label_);
        currentQueue->push(neighborEdge);
        targetQueues.push_back(currentQueue);
        totalLengths[i] += neighborEdge->length;
    }

    // std::cout << "A" << std::endl;
    bool end = 0;
    while (!end) {
        end = 1;
        for (int i = 0; i < targetQueues.size(); ++i) {
            auto& currentQueue = targetQueues[i];
            if (!currentQueue->empty()) {
                end = 0;
                Edge* currentEdge = currentQueue->front();
                int originColor = graph_->edge(currentEdge->endpoints_)->color_;
                graph_->edge(currentEdge->endpoints_)->color_ = i;

                if (originColor != -1) {
                    auto nodes = currentEdge->endpoints_;
                    std::string label1 = nodes.first->label_;
                    std::string label2 = nodes.second->label_;
                    bool condition1 = are_connected(label1, srcNodeLabel, i);
                    bool condition2 = are_connected(label1, srcNodeLabel, originColor);
                    bool condition3 = are_connected(label2, srcNodeLabel, i);
                    bool condition4 = are_connected(label2, srcNodeLabel, originColor);
                    bool condition = (condition1 || condition2) && (condition3 || condition4);
                    // bool condition = true;
                    if (!condition) {
                        graph_->edge(currentEdge->endpoints_)->color_ = originColor;
                        currentQueue->pop();
                        continue;
                    }
                }
                // std::cout << "B" << std::endl;
                for (auto& neighbor : currentEdge->get_neighbors()) {
                    int postColor = graph_->edge(neighbor->endpoints_)->color_;
                    if (postColor == numArea)
                        continue;
                    if (postColor == -1 || totalLengths[postColor] > totalLengths[i]) {
                        currentQueue->push(neighbor);
                        totalLengths[i] += neighbor->length;
                    }
                }
                currentQueue->pop();
            }
        }
    }
}

void Optimizer::flood(const std::string& srcNodeLabel, int numBranch) {
    Graph pipeGraph = *graph_;
    Node& srcNode = pipeGraph.nodes_[srcNodeLabel];
    auto edgeCmp = [](const Edge* a, const Edge* b) {
        return a->length > b->length;
    };

    bleach_all_edges();
    std::vector<Node*> neighbors = srcNode.get_neighbors();
    int size = numBranch < neighbors.size() ? numBranch : neighbors.size();
    std::vector<int> totalLengths(size, 0);
    std::vector<Edge*> floodingEdge;
    for (int i = 0; i < size; ++i) {
        Edge* neighborEdge = graph_->edge(srcNodeLabel, neighbors[i]->label_);
        neighborEdge->color_ = i;
        floodingEdge.push_back(neighborEdge);
    }
    std::sort(floodingEdge.begin(), floodingEdge.end(), edgeCmp);

    while (!floodingEdge.empty()) {
        for (auto& edge : floodingEdge) {
            edge->print_edge();
        }
        system("pause");
        int backLen = floodingEdge.back()->length;
        for (auto& edge : floodingEdge) {
            edge->length -= backLen;
        }
        while (floodingEdge.back()->length == 0) {
            Edge* edge = *floodingEdge.begin();
            for (auto& neighbor : edge->get_neighbors()) {
                if (neighbor->color_ == -1) {
                    floodingEdge.insert(floodingEdge.begin(), neighbor);
                    neighbor->color_ = edge->color_;
                }
            }
            floodingEdge.pop_back();
        }
        std::sort(floodingEdge.begin(), floodingEdge.end(), edgeCmp);
    }
}