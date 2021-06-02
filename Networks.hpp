#pragma once

#include <algorithm>
#include <limits>
#include <cmath>
#include <deque>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <utility> //* pair
#include <vector>

#include "CSV.hpp"
#include "linearAlgebra.hpp"
#include "pcg_random.hpp"

template <typename T>
struct Node {
  public:
    //* Member variables
    T index;
    std::set<T> neighbors;

    //* Constructor
  public:
    Node() {}
    Node(const T& t_index) : index(t_index) {}

    //* Comparison operators
    bool operator==(const Node& t_node) const {
        return this->index == t_node.index;
    }
    bool operator<(const Node& t_node) const {
        return this->index < t_node.index;
    }
};

template <typename T>
struct Network {
    //* Member variables
  public:
    std::string type;
    T size{0};
    unsigned long long linkSize{0};
    double meanDegree{0.0};
    std::vector<std::set<T>> adjacency;

    //* Constructor
  public:
    Network() {}
    Network(const T& t_size) : size(t_size) {
        adjacency.assign(t_size, std::set<T>());
    }

    //* Reset the whole network withiout size information
    void clear();

    //* Whether link between two input nodes is already exists
    bool linkExists(const T&, const T&) const;

    //* Add single link betweeen two input nodes
    //* If there is already a link, this function does nothing
    void addLink(const T&, const T&);

    //* Delete single link between two input nodes
    //* If there is no link, this function does nothing
    void deleteLink(const T&, const T&);

    //* Get number of degree
    const std::map<T, T> getDegreeDist() const;

    //* Get distance between all nodes
    const std::vector<std::vector<T>> getFullDistance() const;

    //* Print total information of network
    void print(const std::string& t_fileName = "") const;

    //* Print adjecency matrix
    void printAdjacency(const std::string& t_fileName = "", const bool t_append = false) const;

    //* Print degree distribution
    void printDegreeDist(const std::string& t_fileName = "", const std::string& t_seperator = ",", const std::string& t_secondSeperator = "\n") const;

    //* Load network from adjacency matrix
    void loadAdjacency(const std::string& t_fileName = "");

}; //* End of struct Network

template <typename T>
void Network<T>::clear() {
    //* Clear information about number of links
    linkSize = 0;
    meanDegree = 0;

    //* Clear adjacency information
    for (auto& neighbors : adjacency) {
        neighbors.clear();
    }
}

template <typename T>
bool Network<T>::linkExists(const T& t_index1, const T& t_index2) const {
    if (adjacency[t_index1].size() <= adjacency[t_index2].size()) {
        const auto candidates = adjacency[t_index1];
        auto it = candidates.find(t_index2);
        return it != candidates.end() ? true : false;
    } else {
        const auto candidates = adjacency[t_index2];
        auto it = candidates.find(t_index1);
        return it != candidates.end() ? true : false;
    }
}

template <typename T>
void Network<T>::addLink(const T& t_index1, const T& t_index2) {
    if (linkExists(t_index1, t_index2)) {
        return;
    }

    //* Update linksize information
    ++linkSize;

    //* Update adjacency information
    adjacency[t_index1].emplace(t_index2);
    adjacency[t_index2].emplace(t_index1);
}

template <typename T>
void Network<T>::deleteLink(const T& t_index1, const T& t_index2) {
    if (!linkExists(t_index1, t_index2)) {
        return;
    }

    //* Update linksize information
    --linkSize;

    //* Update adjacency information
    adjacency[t_index1].erase(t_index2);
    adjacency[t_index2].erase(t_index1);
}

template <typename T>
const std::map<T, T> Network<T>::getDegreeDist() const {
    std::map<T, T> degreeDist;
    for (T index = 0; index < size; ++index) {
        ++degreeDist[adjacency[index].size()];
    }
    return degreeDist;
}

template <typename T>
const std::vector<std::vector<T>> Network<T>::getFullDistance() const {
    //* Initialize distance as infinity (maximum value)
    std::vector<std::vector<T>> fullDistance(size, std::vector<T>(size, std::numeric_limits<T>::max()));

    //* Get shortest path using BFS algorithm
    for (unsigned index = 0; index < size; ++index) {
        //* Assing queue and visited for single index
        std::deque<unsigned> queue;
        std::vector<bool> visited(size, false);
        std::vector<unsigned> distance(size, 0);

        //* Start with initial node
        queue.push_back(index);
        visited[index] = true;

        //* Iterate until reaching every possible nodes
        while (!queue.empty()) {
            const unsigned next = queue.front();
            queue.pop_front();
            fullDistance[index][next] = distance[next];

            //* Check first visited neighbor of 'next'
            for (const unsigned& nextNeighbor : adjacency[next]) {
                if (!visited[nextNeighbor]) {
                    distance[nextNeighbor] = distance[next] + 1;
                    queue.push_back(nextNeighbor);
                    visited[nextNeighbor] = true;
                }
            }
        }
    }
    return fullDistance;
}

template <typename T>
void Network<T>::print(const std::string& t_fileName) const {
    //* Specify out stream: given file or terminal
    std::streambuf* buf;
    std::ofstream file;
    if (t_fileName.size()) {
        file.open(t_fileName);
        buf = file.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    //* Print total data
    std::ostream out(buf);
    out << "Network Type: " << type << "\n";
    out << "Network unsigned: " << size << "\n";
    out << "Total Link unsigned: " << linkSize << "\n";
    out << "Adjacency Matrix\n";
    for (T i = 0; i < size; ++i) {
        out << i << ": ";
        for (const T& neighbor : adjacency[i]) {
            out << neighbor << ", ";
        }
        out << "\n";
    }
}

template <typename T>
void Network<T>::printAdjacency(const std::string& t_fileName, const bool t_append) const {
    //* Specify out stream: given file or terminal
    std::streambuf* buf;
    std::ofstream file;
    if (t_fileName.size()) {
        t_append ? file.open(t_fileName, std::ios_base::app) : file.open(t_fileName);
        buf = file.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    //* Print adjacency matrix of network
    std::ostream out(buf);
    for (T index = 0; index < size; ++index) {
        for (T neighbor = 0; neighbor < size - 1; ++neighbor) {
            adjacency[index].find(neighbor) != adjacency[index].end() ? out << 1 << "," : out << 0 << ",";
        }
        adjacency[index].find(size - 1) != adjacency[index].end() ? out << 1 << "\n" : out << 0 << "\n";
    }
}

template <typename T>
void Network<T>::loadAdjacency(const std::string& t_fileName) {
    //* Read raw adjacency matrix
    std::vector<std::vector<int>> rawAdj;
    CSV::read(t_fileName, rawAdj);

    //* Get information
    size = rawAdj.size();
    adjacency.assign(size, std::set<T>());
    linkSize = 0;
    for (T index=0; index<size; ++index){
        for (T neighbor = 0; neighbor < size - 1; ++neighbor){
            if (rawAdj[index][neighbor]){
                adjacency[index].emplace_hint(adjacency[index].end(), neighbor);
                ++linkSize;
            }
        }
    }
    meanDegree = (double)linkSize / size;
    linkSize /= 2;
}

template <typename T>
void Network<T>::printDegreeDist(const std::string& t_fileName, const std::string& t_seperator, const std::string& t_secondSeperator) const {
    const std::map<T, T> degreeDist = getDegreeDist();

    //* Specify out stream: given file or terminal
    std::streambuf* buf;
    std::ofstream file;
    if (t_fileName.size()) {
        file.open(t_fileName);
        buf = file.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    //* Print degree distribution of network
    std::ostream out(buf);
    for (auto it = degreeDist.begin(); it != degreeDist.end(); ++it) {
        out << it->first << t_seperator << it->second << t_secondSeperator;
    }
}

//* Erdos-Renyi network
namespace ER {
template <typename T>
Network<T> generate(const T& t_size, const double& t_probability, pcg32& t_randomEngine) {
    Network<T> ER(t_size);
    std::uniform_real_distribution<double> probabilityDistribution(0, 1);
    for (T index1 = 0; index1 < t_size; ++index1) {
        for (T index2 = index1 + 1; index2 < t_size; ++index2) {
            if (probabilityDistribution(t_randomEngine)) {
                ER.addLink(index1, index2);
            }
        }
    }
    ER.type = "ER";
    ER.meanDegree = 2.0 * ER.linkSize / t_size;
    return ER;
}

template <typename T>
Network<T> generate(const T& t_size, const unsigned long long& t_linkSize, pcg32& t_randomEngine) {
    Network<T> ER(t_size);
    std::uniform_int_distribution<T> indexDistribution(0, t_size - 1);
    while (ER.linkSize < t_linkSize) {
        T index1, index2;
        do {
            index1 = indexDistribution(t_randomEngine);
            index2 = indexDistribution(t_randomEngine);
        } while (index1 == index2);
        ER.addLink(index1, index2);
    }
    ER.type = "ER";
    ER.meanDegree = 2.0 * t_linkSize / t_size;
    return ER;
}
} // namespace ER

//* Random Regular network
namespace RR {
template <typename T>
Network<T> generate(const T& t_size, const T& t_degree, pcg32& t_randomEngine) {
    Network<T> RR(t_size);
    std::deque<T> stubs;
    for (T index = 0; index < t_size; ++index) {
        for (T degree = 0; degree < t_degree; ++degree) {
            stubs.emplace_back(index);
        }
    }
    std::shuffle(stubs.begin(), stubs.end(), t_randomEngine);
    while (!stubs.empty()) {
        const T index1 = stubs.front();
        stubs.pop_front();
        const T index2 = stubs.front();
        stubs.pop_front();
        RR.addLink(index1, index2);
    }
    RR.type = "RR";
    RR.meanDegree = (double)t_degree;
    return RR;
}
} // namespace RR

//* Scale free network
namespace SF {
//* power law distribution with input exponent excluding t_lower and t_upper
template <typename T>
T randomPowerLawDistribution(const int& t_lower, const T& t_upper, const double& t_exponent, const double& t_prob) {
    return std::pow((std::pow(t_upper + 0.5, t_exponent + 1) - std::pow(t_lower - 0.5, t_exponent + 1)) * t_prob + std::pow(t_lower - 0.5, t_exponent + 1), 1.0 / (t_exponent + 1)) + 0.5;
}

template <typename T>
Network<T> generate(const T& t_size, const unsigned long long& t_linkSize, const double& t_degreeExponent, pcg32& t_randomEngine) {
    Network<T> SF(t_size);
    const double weightExponent = 1.0 / (t_degreeExponent - 1.0);
    std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
    while (SF.linkSize < t_linkSize) {
        T index1, index2;
        do {
            index1 = randomPowerLawDistribution(1, t_size, -weightExponent, probabilityDistribution(t_randomEngine)) - 1;
            index2 = randomPowerLawDistribution(1, t_size, -weightExponent, probabilityDistribution(t_randomEngine)) - 1;
        } while (index1 == index2);
        SF.addLink(index1, index2);
    }
    SF.type = "SF";
    SF.meanDegree = 2.0 * t_linkSize / t_size;
    return SF;
}
} // namespace SF

//* Scale free network by Chung-Lu
namespace CL {
//* weight distribution with input weight
template <typename T>
T weightSampling(const std::vector<double>& t_weight, const double& t_prob, const T& t_dummy) {
    return (T)(std::lower_bound(t_weight.begin(), t_weight.end(), t_prob * t_weight.back()) - t_weight.begin());
}

template <typename T>
Network<T> generate(const T& t_size, const unsigned long long& t_linkSize, const double& t_degreeExponent, pcg32& t_randomEngine) {
    Network<T> CL(t_size);
    const double weightExponent = 1.0 / (t_degreeExponent - 1.0);
    const double correction = weightExponent < 0.5 ? 1.0 : std::pow(10.0 * std::sqrt(2.0) * (1.0 - weightExponent), 1.0 / weightExponent) * std::pow(t_size, 1.0 - 1.0 / (2.0 - weightExponent));

    std::vector<double> weight(t_size, 0.0);
    weight[0] = std::pow(correction, -1.0 * weightExponent);
    for (T i = 1; i < t_size; ++i) {
        weight[i] = weight[i - 1] + std::pow(i + correction, -1.0 * weightExponent);
    }

    std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
    while (CL.linkSize < t_linkSize) {
        T index1, index2;
        do {
            index1 = weightSampling(weight, probabilityDistribution(t_randomEngine), t_size);
            index2 = weightSampling(weight, probabilityDistribution(t_randomEngine), t_size);
        } while (index1 == index2);
        CL.addLink(index1, index2);
    }
    CL.type = "CL";
    CL.meanDegree = 2.0 * t_linkSize / t_size;
    return CL;
}
} // namespace CL

//* Weighted node
template <typename T>
struct WNode {
  public:
    //* Member variables
    T index;
    std::map<T, double> neighbors;

    //* Constructor
    WNode() {}
    WNode(const T& t_index) : index(t_index) {}

    //* Comparison operator
    bool operator==(const WNode& t_wnode) const {
        return this->index == t_wnode.index;
    }
    bool operator<(const WNode& t_wnode) const {
        return this->index < t_wnode.index;
    }
};

//* Weighted network
template <typename T>
struct WNetwork {
  public:
    //* Member variables
    std::string type;
    T size{0};
    unsigned long long linkSize{0};
    double meanDegree{0.0};
    std::vector<std::map<T, double>> wadjacency;

  public:
    //* Constructor
    WNetwork() {}
    WNetwork(const T& t_size) : size(t_size) {
        wadjacency.assign(t_size, std::map<T, double>());
        for (T index = 0; index < t_size; ++index) {
            wadjacency[index][index] = 0.0;
        }
    }

    //* Reset the whole network
    void clear();

    //* Return weight if link exists or return 0
    const double getWeight(const T&, const T&) const;

    //* Add single link
    void addLink(const T&, const T&, const double& t_weight = 1.0);

    //* Delete single link
    void deleteLink(const T&, const T&);

    //* Get normal degree distribution
    const std::map<T, T> getDegreeDist() const;

    //* Get node weight distribution
    const std::map<double, T> getNodeWeightDist() const;

    //* Print Adjacency matrix with bool (except node weight)
    void printAdjacency(const std::string& t_fileName = "", const bool t_append = false) const;

    //* Print Adjacency matrix with weight
    void printWAdjacency(const std::string& t_fileName = "", const bool t_append = false) const;

    //* Prind degree distribution
    void printDegreeDist(const std::string& t_fileName = "", const std::string& t_seperator = ",", const std::string& t_secondSeperator = "\n") const;

    //* Prind node weight distribution
    void printNodeWeightDist(const std::string& t_fileName = "", const std::string& t_seperator = ",", const std::string& t_secondSeperator = "\n") const;

}; //* End of struct WNetwork

template <typename T>
void WNetwork<T>::clear() {
    //* Clear link size information
    linkSize = 0;
    meanDegree = 0.0;

    //* Clear adjacency information
    for (T index = 0; index < size; ++index) {
        wadjacency[index].clear();
        wadjacency[index][index] = 0.0;
    }
}

template <typename T>
const double WNetwork<T>::getWeight(const T& t_index1, const T& t_index2) const {
    if (wadjacency[t_index1].size() <= wadjacency[t_index2].size()) {
        const std::map<T, double> candidates = wadjacency[t_index1];
        const auto it = candidates.find(t_index2);
        return it != candidates.end() ? it->second : 0.0;
    } else {
        const std::map<T, double> candidates = wadjacency[t_index2];
        const auto it = candidates.find(t_index1);
        return it != candidates.end() ? it->second : 0.0;
    }
}

template <typename T>
void WNetwork<T>::addLink(const T& t_index1, const T& t_index2, const double& t_weight) {
    if (getWeight(t_index1, t_index2)) {
        return;
    }
    //* Update linksize information
    ++linkSize;

    //* Update weighted adjacency information
    wadjacency[t_index1][t_index2] += t_weight;
    wadjacency[t_index2][t_index1] += t_weight;
}

template <typename T>
void WNetwork<T>::deleteLink(const T& t_index1, const T& t_index2) {
    if (!getWeight(t_index1, t_index2)) {
        return;
    }
    //* Update linksize information
    --linkSize;

    //* Update weighted adjacency information
    wadjacency[t_index1].erase(t_index2);
    wadjacency[t_index2].erase(t_index1);
}

template <typename T>
const std::map<T, T> WNetwork<T>::getDegreeDist() const {
    std::map<T, T> degreeDist;
    for (T index = 0; index < size; ++index) {
        ++degreeDist[wadjacency[index].size() - 1];
    }
    return degreeDist;
}

template <typename T>
const std::map<double, T> WNetwork<T>::getNodeWeightDist() const {
    std::map<double, T> nodeWeightDist;
    for (T index = 0; index < size; ++index) {
        ++nodeWeightDist[wadjacency[index].at(index)];
    }
    return nodeWeightDist;
}

template <typename T>
void WNetwork<T>::printAdjacency(const std::string& t_fileName, const bool t_append) const {
    //* Specify out stream: given file or terminal
    std::streambuf* buf;
    std::ofstream file;
    if (t_fileName.size()) {
        t_append ? file.open(t_fileName, std::ios_base::app) : file.open(t_fileName);
        buf = file.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    //* Print adjacency matrix of network
    std::ostream out(buf);
    for (T index = 0; index < size; ++index) {
        for (T neighbor = 0; neighbor < size - 1; ++neighbor) {
            (wadjacency[index].find(neighbor) != wadjacency[index].end() && neighbor != index) ? out << 1 << "," : out << 0 << ",";
        }
        (wadjacency[index].find(size - 1) != wadjacency[index].end() && size - 1 != index) ? out << 1 << "\n" : out << 0 << "\n";
    }
}

template <typename T>
void WNetwork<T>::printWAdjacency(const std::string& t_fileName, const bool t_append) const {
    //* Specify out stream: given file or terminal
    std::streambuf* buf;
    std::ofstream file;
    if (t_fileName.size()) {
        t_append ? file.open(t_fileName, std::ios_base::app) : file.open(t_fileName);
        buf = file.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    //* Print weighted adjacency matrix of network
    std::ostream out(buf);
    for (T index = 0; index < size; ++index) {
        for (T neighbor = 0; neighbor < size - 1; ++neighbor) {
            neighbor != index ? out << getWeight(index, neighbor) << "," : out << 0 << ",";
        }
        index != size - 1 ? out << getWeight(index, size - 1) << "\n" : out << 0 << "\n";
    }
}

template <typename T>
void WNetwork<T>::printDegreeDist(const std::string& t_fileName, const std::string& t_seperator, const std::string& t_secondSeperator) const {
    const std::map<T, T> degreeDist = getDegreeDist();

    //* Specify out stream: given file or terminal
    std::streambuf* buf;
    std::ofstream file;
    if (t_fileName.size()) {
        file.open(t_fileName);
        buf = file.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    //* Print degree distribution of network
    std::ostream out(buf);
    for (auto it = degreeDist.begin(); it != degreeDist.end(); ++it) {
        out << it->first << t_seperator << it->second << t_secondSeperator;
    }
}

template <typename T>
void WNetwork<T>::printNodeWeightDist(const std::string& t_fileName, const std::string& t_seperator, const std::string& t_secondSeperator) const {
    const std::map<double, T> nodeWeightDist = getNodeWeightDist();

    //* Specify out stream: given file or terminal
    std::streambuf* buf;
    std::ofstream file;
    if (t_fileName.size()) {
        file.open(t_fileName);
        buf = file.rdbuf();
    } else {
        buf = std::cout.rdbuf();
    }

    //* Print node distribution of network
    std::ostream out(buf);
    for (auto it = nodeWeightDist.begin(); it != nodeWeightDist.end(); ++it) {
        out << it->first << t_seperator << it->second << t_secondSeperator;
    }
}

//* Weighted scale free network by Chung-Lu
//* Degree distribution and node weight distribution follows same power law
namespace WCL {
template <typename T>
T weightSampling(const std::vector<double>& t_weight, const double& t_prob, const T& t_dummy) {
    return (unsigned)(std::lower_bound(t_weight.begin(), t_weight.end(), t_prob * t_weight.back()) - t_weight.begin());
}

template <typename T>
WNetwork<T> generate(const T& t_size, const double& t_meanPopulation, const double& t_meanDegree, const double& t_degreeExponent, const double& t_linkWeightExponent, pcg32& t_randomEngine) {
    WNetwork<T> WCL(t_size);
    //* Genearting conventional scale free network by Chung-Lu
    const double weightExponent = 1.0 / (t_degreeExponent - 1.0);
    const double correction = weightExponent < 0.5 ? 1.0 : std::pow(10.0 * std::sqrt(2.0) * (1.0 - weightExponent), 1.0 / weightExponent) * std::pow(t_size, 1.0 - 1.0 / (2.0 - weightExponent));

    std::vector<double> weight(t_size, 0.0);
    weight[0] = std::pow(correction, -1.0 * weightExponent);
    for (T i = 1; i < t_size; ++i) {
        weight[i] = weight[i - 1] + std::pow(i + correction, -1.0 * weightExponent);
    }

    std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
    while (WCL.linkSize < t_meanDegree * t_size / 2.0) {
        T index1, index2;
        do {
            index1 = weightSampling(weight, probabilityDistribution(t_randomEngine), t_size);
            index2 = weightSampling(weight, probabilityDistribution(t_randomEngine), t_size);
        } while (index1 == index2);
        WCL.addLink(index1, index2);
    }

    WCL.type = "CL";
    WCL.meanDegree = t_meanDegree;

    //* Add weight to node and link
    const unsigned long long population = t_size * t_meanPopulation;
    for (T index = 0; index < t_size; ++index) {
        //* Get weight of node and link
        const T degree = WCL.wadjacency[index].size() - 1;
        const double nodeWeight = degree / (2.0 * WCL.linkSize);
        std::map<T, double> linkWeight;
        for (auto it = WCL.wadjacency[index].begin(); it != WCL.wadjacency[index].end(); ++it) {
            const T neighbor = it->first;
            if (neighbor != index) {
                const T neighborDegree = WCL.wadjacency[neighbor].size() - 1;
                linkWeight[neighbor] = std::pow((double)degree * neighborDegree, t_linkWeightExponent);
            }
        }
        const double totalLinkWeight = linearAlgebra::accumulate(linkWeight);

        //* Add weight to weighted adjacency matrix
        for (auto it = WCL.wadjacency[index].begin(); it != WCL.wadjacency[index].end(); ++it) {
            if (it->first != index) {
                WCL.wadjacency[index][it->first] = linkWeight[it->first] / totalLinkWeight;
            } else {
                WCL.wadjacency[index][it->first] = std::floor(population * nodeWeight + 0.5);
            }
        }
    }
    return WCL;
}
} // namespace WCL