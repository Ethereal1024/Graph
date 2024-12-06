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
    std::cout << endpoints_.second->label_ << std::endl;
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
            std::cerr << "Error: connection already exists." << std::endl;
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

void Optimizer::DFS_step(Node* currentNode, int currentDistance, std::vector<std::string>& currentTrack, std::vector<Path>& paths, const int& barrierColor, const int& targetColor) {
    if (currentNode->color_ == barrierColor) {
        return;
    }
    if (currentNode->color_ == targetColor) {
        paths.push_back({currentDistance, currentTrack});
        currentNode->color_ = barrierColor;
        return;
    }

    currentNode->color_ = barrierColor;

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

/*
std::vector<Path> Optimizer::separate_euler_path() {
    std::vector<Node*> oddNodes = odd_nodes();
    std::vector<std::vector<Path>> oddNodePaths;

    auto pathCmp = [](const Path& a, const Path& b) {
        return a.dist_ < b.dist_;
    };

    auto vectorSizeCmp = [](const std::vector<Path>& a, const std::vector<Path>& b) {
        return a.size() > b.size();
    };

    for (Node* oddNode : oddNodes) {
        std::vector<Path> currentOddNodePath = odd_node_connections(oddNode->label_);
        std::sort(currentOddNodePath.begin(), currentOddNodePath.end(), pathCmp);
        oddNodePaths.push_back(currentOddNodePath);
    }
    std::sort(oddNodePaths.begin(), oddNodePaths.end(), vectorSizeCmp);

    std::vector<Path> result;
    while(!oddNodePaths.empty()) {
        result.push_back(oddNodePaths.back().back());
    }
}
*/

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

    delete distLogs;
    distLogs = new std::vector<int>(centers.size());
    for (int i = 0; i < centers.size(); i++) {
        int distLog;
        centers[i] = find_center(groups[i], distLog);
        (*distLogs)[i] = distLog;
    }
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

float Drawer::distance_on_canvas(const Point& p1, const Point& p2) {
    if (p1.x == INF || p2.x == INF) {
        return INF;
    }
    int deltaX = p1.x - p2.x;
    int deltaY = p1.y - p2.y;
    return sqrt(deltaX * deltaX + deltaY * deltaY);
}

float Drawer::distance_on_graph(const std::string& label1, const std::string& label2) {
    Edge* edge = graph_->edge(label1, label2);
    if (edge == nullptr) {
        return INF;
    } else {
        return edge->length;
    }
}

std::pair<Drawer::Point, Drawer::Point> Drawer::circle_intersection(const Point& p1, const Point& p2, const float& r1, const float& r2) {
    float x1 = p1.x;
    float x2 = p2.x;
    float y1 = p1.y;
    float y2 = p2.y;

    // std::cout << "\np1: (" << p1.x << ", " << p1.y << ")\n";
    // std::cout << "p2: (" << p2.x << ", " << p2.y << ")\n";
    // std::cout << "r1: " << r1 << "\n";
    // std::cout << "r2: " << r2 << "\n";

    float originDistance = distance_on_canvas(p1, p2);
    // std::cout << originDistance << "\n";
    if (originDistance >= (r1 + r2)) {

        // std::cout << originDistance << std::endl;

        float x = x1 + (x2 - x1) * (r1 / r2);
        float y = y1 + (y2 - y1) * (r1 / r2);
        return std::pair<Point, Point>({{x, y}, {x, y}});
    }

    if (x1 == x2) {
        float x = x1;
        float yUp, yDown, rUp, rDown;
        if (y1 < y2) {
            yUp = y2;
            yDown = y1;
            rUp = r2;
            rDown = r1;
        } else {
            yUp = y1;
            yDown = y2;
            rUp = r1;
            rDown = r2;
        }
        float h = yUp - yDown;
        float t = abs((rUp * rUp - rDown * rDown + h * h) / (2 * h));
        float d = sqrt(rUp * rUp - t * t);
        float y = yUp - t;
        Point result1 = {x - d, y};
        Point result2 = {x + d, y};
        return std::pair<Point, Point>({result1, result2});
    }

    else if (y1 == y2) {
        float y = y1;
        float xLeft, xRight, rLeft, rRight;
        if (x1 < x2) {
            xLeft = x1;
            xRight = x2;
            rLeft = r1;
            rRight = r2;
        } else {
            xLeft = x2;
            xRight = x1;
            rLeft = r2;
            rRight = r1;
        }
        float d = xRight - xLeft;
        float t = abs((rRight * rRight - rLeft * rLeft + d * d) / (2 * d));
        float h = sqrt(rRight * rRight - t * t);
        float x = xRight - t;
        Point result1 = {x, y - h};
        Point result2 = {x, y + h};
        return std::pair<Point, Point>({result1, result2});
    }

    else {
        float k1 = (y1 - y2) / (x1 - x2);
        float k2 = -(x1 - x2) / (y1 - y2);
        float b1 = (x1 * y2 - x2 * y1) / (x1 - x2);
        float b2 = -(r1 * r1 - x1 * x1 - y1 * y1 - r2 * r2 + x2 * x2 + y2 * y2) / (2 * (y1 - y2));

        // std::cout << k1 << " " << b1 << "\n" << k2 << " " << b2 << "\n";

        float xMid = -(b1 - b2) / (k1 - k2);
        float xLeft, xRight;
        if (x1 < x2) {
            xLeft = x1 - r1;
            xRight = x2 + r2;
        } else {
            xLeft = x2 - r2;
            xRight = x1 + r1;
        }

        float begin1 = xLeft, end1 = xMid, begin2 = xRight, end2 = xMid;
        float test1 = (begin1 + end1) / 2;
        float test2 = (begin2 + end2) / 2;

        auto search = [&](float& test, float& begin, float& end, const float& r, const Point& p) {
            while (abs(end - begin) > 0.01) {
                // std::cout << begin << " " << test << " " << end << "\n";
                // system("pause");
                Point testPoint = {test, k2 * test + b2};
                if (distance_on_canvas(testPoint, p) < r) {
                    end = test;
                } else if (distance_on_canvas(testPoint, p) > r) {
                    begin = test;
                } else {
                    break;
                }
                test = (begin + end) / 2;
            }
        };

        search(test1, begin1, end1, r1, p1);
        search(test2, begin2, end2, r2, p2);

        Point result1 = {test1, k2 * test1 + b2};
        Point result2 = {test2, k2 * test2 + b2};
        return std::pair<Point, Point>({result1, result2});
    }
}

void Drawer::fix_node_step(const std::string& nodeLabel) {
    if (coordinates_[nodeLabel].fixed()) {
        return;
    }

    std::vector<std::string> fixedNeighbors;
    for (auto it = coordinates_.begin(); fixedNeighbors.size() < 3 && it != coordinates_.end(); ++it) {
        if (it->second.fixed() && graph_->edge(it->first, nodeLabel) != nullptr) {
            fixedNeighbors.push_back(it->first);
        }
    }

    if (fixedNeighbors.size() == 0) {
        coordinates_[nodeLabel] = {0, 0};
    }

    else if (fixedNeighbors.size() == 1) {
        std::string neighborLabel = fixedNeighbors[0];
        std::vector<Node*> neighborNeighbors = graph_->nodes_[neighborLabel].get_neighbors();
        int index = 0;
        for (int i = 0; i < neighborNeighbors.size(); i++) {
            if (neighborNeighbors[i]->label_ == nodeLabel) {
                index = i;
            }
        }
        float len = graph_->edge(neighborLabel, nodeLabel)->length;
        Point pointDelta;
        if (index % 4 == 0) {
            pointDelta = Point({len, 0});
        } else if (index % 4 == 1) {
            pointDelta = Point({0, len});
        } else if (index % 4 == 2) {
            pointDelta = Point({-len, 0});
        } else {
            pointDelta = Point({0, -len});
        }
        coordinates_[nodeLabel] = coordinates_[neighborLabel] + pointDelta;
    }

    else {
        std::string label1 = fixedNeighbors[0];
        std::string label2 = fixedNeighbors[1];

        const Point& point1 = coordinates_[label1];
        const Point& point2 = coordinates_[label2];

        float len1 = graph_->edge(label1, nodeLabel)->length;
        float len2 = graph_->edge(label2, nodeLabel)->length;

        std::pair<Point, Point> choise = circle_intersection(point1, point2, len1, len2);

        Point point3;
        std::string label3;
        if (fixedNeighbors.size() == 2) {
            std::vector<Node*> node1Neighbors = graph_->nodes_[label1].get_neighbors();
            for (auto neighbor : node1Neighbors) {
                Point& candidate = coordinates_[neighbor->label_];
                if (candidate.fixed()) {
                    point3 = candidate;
                    label3 = neighbor->label_;
                    break;
                }
            }
        } else {
            point3 = coordinates_[fixedNeighbors[2]];
        }

        float graphDistance = distance_on_graph(nodeLabel, label3);
        float cost1 = abs(graphDistance - distance_on_canvas(choise.first, point3));
        float cost2 = abs(graphDistance - distance_on_canvas(choise.second, point3));

        coordinates_[nodeLabel] = cost1 < cost2 ? choise.first : choise.second;
    }
    // std::cout << "(" << coordinates_[nodeLabel].x << "," << coordinates_[nodeLabel].y << ")" << std::endl;
}

std::pair<float, float> Drawer::get_coordinate(const std::string& nodeLabel) {
    Point coordinate = coordinates_[nodeLabel];
    if (!coordinate.fixed()) {
        std::cerr << "Error: node not fixed." << std::endl;
    }
    float x = std::round(coordinate.x * 100.0f) / 100.0f;
    float y = std::round(coordinate.y * 100.0f) / 100.0f;
    return std::pair<float, float>({x, y});
}

void Drawer::fix_node(const std::string srcLabel) {
    std::queue<std::string> fixQueue;
    fixQueue.push(srcLabel);
    graph_->nodes_[srcLabel].color_ = COLOR::RED;
    while (!fixQueue.empty()) {
        std::string& label = fixQueue.front();
        // std::cout << label << " ";
        
        fix_node_step(label);

        // auto coordinate = coordinates_[label];
        // std::cout << "\\node[draw, circle] (" << label << ") at (" << coordinate.x / 4 << ", " << coordinate.y / 4;
        // if (label.size() > 1) {
        //     label = "\\" + label;
        // }
        // std::cout << ") {$" << label << "$};\n";
        // system("pause");

        std::vector<Node*> neighbors = graph_->nodes_[label].get_neighbors();
        for (Node* neighbor : neighbors) {
            const std::string& neighborLabel = neighbor->label_;
            if (graph_->nodes_[neighborLabel].color_ != COLOR::RED) {
                fixQueue.push(neighborLabel);
                graph_->nodes_[neighborLabel].color_ = COLOR::RED;
            }
        }

        fixQueue.pop();
    }
    // std::cout << "\n";
}
