#include "graph.hpp"

typedef struct EdgeCmd {
    std::string p1;
    std::string p2;
    int length;
} EdgeCommand;

std::vector<EdgeCmd> edgeCmds = {
    {"Q", "H", 8},
    {"Q", "S", 6},
    {"Q", "T", 6},
    {"T", "R", 3},
    {"T", "U", 4},
    {"H", "Gamma", 6},
    {"H", "V", 6},
    {"S", "V", 4},
    {"S", "U", 3},
    {"U", "P", 2},
    {"R", "P", 5},
    {"R", "s", 6},
    {"P", "z", 6},
    {"s", "y", 6},
    {"y", "z", 4},
    {"V", "Omega", 4},
    {"Omega", "W", 3},
    {"W", "O", 3},
    {"z", "A", 8},
    {"O", "A", 8},
    {"Gamma", "gamma", 2},
    {"gamma", "alpha", 5},
    {"alpha", "beta", 5},
    {"alpha", "Z", 2},
    {"Omega", "omega", 2},
    {"omega", "X", 4},
    {"omega", "Z", 1},
    {"Z", "Y", 5},
    {"beta", "Y", 2},
    {"gamma", "eta", 3},
    {"eta", "beta", 3},
    {"eta", "Pi", 4},
    {"Gamma", "Pi", 8},
    {"Pi", "G", 4},
    {"beta", "phi", 6},
    {"G", "phi", 5},
    {"phi", "Phi", 4},
    {"Y", "Phi", 6},
    {"Y", "X", 1},
    {"X", "N", 5},
    {"O", "N", 5},
    {"N", "M", 5},
    {"Phi", "M", 8},
    {"G", "F", 9},
    {"M", "L", 10},
    {"F", "L", 16},
    {"F", "I", 10},
    {"A", "B", 14},
    {"B", "K", 10},
    {"L", "K", 14},
    {"K", "J", 6},
    {"I", "J", 16},
    {"I", "a", 12},
    {"J", "f", 12},
    {"a", "f", 20},
    {"a", "b", 14},
    {"b", "c", 20},
    {"c", "d", 8},
    {"d", "e", 8},
    {"e", "f", 16},
    {"c", "g", 8},
    {"d", "h", 8},
    {"g", "h", 6},
    {"g", "theta", 3},
    {"theta", "h", 3},
    {"h", "i", 16},
    {"i", "C", 20},
    {"C", "D", 12},
    {"D", "f", 10},
    {"D", "E", 8},
    {"C", "E", 7},
    {"i", "o", 12},
    {"o", "n", 8},
    {"n", "m", 5},
    {"o", "j", 8},
    {"j", "m", 8},
    {"m", "k", 10},
    {"j", "k", 5},
    {"k", "l", 3},
    {"l", "j", 5},
    {"l", "p", 10},
    {"p", "q", 16},
    {"p", "r", 20},
    {"r", "t", 4},
    {"t", "v", 10},
    {"t", "u", 6},
    {"u", "v", 6},
    {"u", "x", 5},
    {"v", "w", 4},
    {"x", "w", 5},
    {"x", "s", 6},
    {"w", "y", 6},
    {"E", "e", 16}};

void test() {
    Graph graph(6);
    graph.connect("A", "C", 8);
    graph.connect("A", "E", 5);
    graph.connect("A", "B", 11);
    graph.connect("B", "D", 8);
    graph.connect("C", "D", 2);
    graph.connect("B", "E", 5);
    graph.connect("E", "D", 3);
    graph.connect("E", "F", 4);
    graph.connect("B", "F", 2);
    graph.display_edges();
    graph.display_nodes();
    std::cout
        << "\n";

    Optimizer optim(graph);

    std::cout << graph.edge("E", "F")->length << std::endl;
    std::unordered_map<std::string, Path> min_distances = optim.dijkstra("B");
    for (auto element : min_distances) {
        std::cout << element.first << ": " << element.second.dist_ << std::endl;
        std::cout << "(";
        auto& track = element.second.track_;
        for (int i = 0; i < track.size(); i++) {
            std::cout << track[i];
            if (i != track.size() - 1) {
                std::cout << "->";
            }
        }
        std::cout << ")";
        std::cout << "\n";
    }
}

Graph* createGraphByCommands() {
    Graph* graph = new Graph;
    for (auto command : edgeCmds) {
        graph->add_node(command.p1);
        graph->add_node(command.p2);
        graph->connect(command.p1, command.p2, command.length);
    }
    return graph;
}

int main() {
    Graph graph = *createGraphByCommands();
    // Optimizer optim(graph);
    // std::vector<std::string> startCenters = {"z", "j", "d", "i", "r", "beta"};
    // optim.graph_k_means(startCenters, 30);
    // // graph.display_nodes();
    // for(auto it : graph.get_nodes()) {
    //     std::cout << it.second.label_ << " " << it.second.color_ << "\n";
    // }
    // for(auto it : startCenters) {
    //     std::cout << it << " ";
    // }
    // std::cout << "\n";
    Drawer drawer(graph);
    auto nodeLabels = graph.get_nodes_labels();
    drawer.fix_node("r");
    for (auto label : nodeLabels) {
        std::pair<float, float> coordinate = drawer.get_coordinate(label);
        // if(label == "Z" || label == "c" || label == "b" || label == "I" || label == "F" || label == "a" || label == "alpha"){
        //     continue;
        // }
        // std::cout << label << ": (" << coordinate.first << ", " << coordinate.second << ")" << std::endl;
        
        // std::cout << "\\node[draw, circle] (" << label << ") at (" << coordinate.first / 10 << ", " << coordinate.second / 10;
        // if (label.size() > 1) {
        //     label = "\\" + label;
        // }
        // std::cout << ") {$" << label << "$};\n";
    }
    return 0;
}