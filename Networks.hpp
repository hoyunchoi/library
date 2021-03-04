#pragma once

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <utility>  //* pair
#include <vector>

#include "linearAlgebra.hpp"
#include "pcg_random.hpp"

struct Node {
    //* Member variables
    unsigned m_index;
    std::set<unsigned> m_neighbors;

    //* Constructor
    Node() {}
    Node(const unsigned& t_index) : m_index(t_index) {}

    //* Comparison operators
    bool operator==(const Node& t_node) const {
        return this->m_index == t_node.m_index;
    }
    bool operator<(const Node& t_node) const {
        return this->m_index < t_node.m_index;
    }
};

struct Network {
    //* Member variables
    std::string m_type;
    unsigned m_size{0};
    unsigned m_linkSize{0};
    double m_meanDegree{0.0};
    std::vector<std::set<unsigned>> m_adjacency;

    //* Constructor
    Network() {}
    Network(const unsigned& t_size) : m_size(t_size) {
        m_adjacency.assign(t_size, std::set<unsigned>());
    }

    //* Reset the whole network withiout size information
    void clear() {
        //* Clear linksize information
        m_linkSize = 0;

        //* Clear adjacency information
        for (auto& neighbors : m_adjacency) {
            neighbors.clear();
        }
    }

    //* Whether link is already exists
    bool linkExists(const unsigned& t_index1, const unsigned& t_index2) const {
        if (m_adjacency[t_index1].size() <= m_adjacency[t_index2].size()) {
            const auto candidates = m_adjacency[t_index1];
            auto it = candidates.find(t_index2);
            return it != candidates.end() ? true : false;
        } else {
            const auto candidates = m_adjacency[t_index2];
            auto it = candidates.find(t_index1);
            return it != candidates.end() ? true : false;
        }
    }

    //* Add single link
    void addLink(const unsigned& t_index1, const unsigned& t_index2) {
        //* Update linksize information
        if (!linkExists(t_index1, t_index2)) {
            ++m_linkSize;
        }

        //* Update adjacency information
        m_adjacency[t_index1].emplace(t_index2);
        m_adjacency[t_index2].emplace(t_index1);
    }

    //* Delete single link
    void deleteLink(const unsigned& t_index1, const unsigned& t_index2) {
        //* Update linksize information
        if (linkExists(t_index1, t_index2)) {
            --m_linkSize;
        }

        //* Update adjacency information
        m_adjacency[t_index1].erase(t_index2);
        m_adjacency[t_index2].erase(t_index1);
    }

    //* Print total information of network
    void print() {
        std::cout << "Network Type: " << m_type << "\n";
        std::cout << "Network unsigned: " << m_size << "\n";
        std::cout << "Total Link unsigned: " << m_linkSize << "\n";
        std::cout << "Adjacency Matrix\n";
        for (unsigned i = 0; i < m_size; ++i) {
            std::cout << i << ": ";
            for (const unsigned& neighbor : m_adjacency[i]) {
                std::cout << neighbor << ", ";
            }
            std::cout << "\n";
        }
    }
};  //* End of struct Network

//* Erdos-Renyi network
namespace ER {
Network generate(const unsigned& t_size, const double& t_probability, pcg32& t_randomEngine) {
    Network ER(t_size);
    std::uniform_real_distribution<double> probabilityDistribution(0, 1);
    for (unsigned index1 = 0; index1 < t_size; ++index1) {
        for (unsigned index2 = index1 + 1; index2 < t_size; ++index2) {
            if (probabilityDistribution(t_randomEngine)) {
                ER.addLink(index1, index2);
            }
        }
    }
    ER.m_type = "ER";
    ER.m_meanDegree = 2.0 * ER.m_linkSize / t_size;
    return ER;
}

Network generate(const unsigned& t_size, const unsigned& t_linkSize, pcg32& t_randomEngine) {
    Network ER(t_size);
    std::uniform_int_distribution<int> indexDistribution(0, t_size - 1);
    while (ER.m_linkSize < t_linkSize) {
        unsigned index1, index2;
        do {
            index1 = indexDistribution(t_randomEngine);
            index2 = indexDistribution(t_randomEngine);
        } while (index1 == index2);
        ER.addLink(index1, index2);
    }
    ER.m_type = "ER";
    ER.m_meanDegree = 2.0 * t_linkSize / t_size;
    return ER;
}
}  // namespace ER

//* Random Regular network
namespace RR {
Network generate(const unsigned& t_size, const unsigned& t_degree, pcg32& t_randomEngine) {
    Network RR(t_size);
    std::deque<unsigned> stubs;
    for (unsigned index = 0; index < t_size; ++index) {
        for (unsigned degree = 0; degree < t_degree; ++degree) {
            stubs.emplace_back(index);
        }
    }
    std::shuffle(stubs.begin(), stubs.end(), t_randomEngine);
    while (!stubs.empty()) {
        unsigned index1 = stubs.front();
        stubs.pop_front();
        unsigned index2 = stubs.front();
        stubs.pop_front();
        RR.addLink(index1, index2);
    }
    RR.m_type = "RR";
    RR.m_meanDegree = (double)t_degree;
    return RR;
}
}  // namespace RR

//* Scale free network
namespace SF {
unsigned randomPowerLawDistribution(const int& t_lower, const int& t_upper, const double& t_exponent, const double& t_prob) {
    std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
    return std::pow((std::pow(t_upper + 0.5, t_exponent + 1) - std::pow(t_lower - 0.5, t_exponent + 1)) * t_prob + std::pow(t_lower - 0.5, t_exponent + 1), 1.0 / (t_exponent + 1)) + 0.5;
}

Network generate(const unsigned& t_size, const unsigned& t_linkSize, const double& t_degreeExponent, pcg32& t_randomEngine) {
    Network SF(t_size);
    const double weightExponent = 1.0 / (t_degreeExponent - 1.0);
    std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
    while (SF.m_linkSize < t_linkSize) {
        unsigned index1, index2;
        do {
            index1 = randomPowerLawDistribution(1, t_size, -weightExponent, probabilityDistribution(t_randomEngine)) - 1;
            index2 = randomPowerLawDistribution(1, t_size, -weightExponent, probabilityDistribution(t_randomEngine)) - 1;
        } while (index1 == index2);
        SF.addLink(index1, index2);
    }
    SF.m_type = "SF";
    SF.m_meanDegree = 2.0 * t_linkSize / t_size;
    return SF;
}
}  // namespace SF

//* Scale free network by Chung-Lu
namespace CL {
unsigned weightSampling(const std::vector<double>& t_weight, const double& t_prob) {
    return (unsigned)(std::lower_bound(t_weight.begin(), t_weight.end(), t_prob * t_weight.back()) - t_weight.begin());
}

Network generate(const unsigned& t_size, const unsigned& t_linkSize, const double& t_degreeExponent, pcg32& t_randomEngine) {
    Network CL(t_size);
    const double weightExponent = 1.0 / (t_degreeExponent - 1.0);
    const double correction = weightExponent < 0.5 ? 1.0 : std::pow(10.0 * std::sqrt(2.0) * (1.0 - weightExponent), 1.0 / weightExponent) * std::pow(t_size, 1.0 - 1.0 / (2.0 - weightExponent));

    std::vector<double> weight(t_size, 0.0);
    weight[0] = std::pow(correction, -1.0 * weightExponent);
    for (unsigned i = 1; i < t_size; ++i) {
        weight[i] = weight[i - 1] + std::pow(i + correction, -1.0 * weightExponent);
    }

    std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
    while (CL.m_linkSize < t_linkSize) {
        unsigned index1, index2;
        do {
            index1 = weightSampling(weight, probabilityDistribution(t_randomEngine));
            index2 = weightSampling(weight, probabilityDistribution(t_randomEngine));
        } while (index1 == index2);
        CL.addLink(index1, index2);
    }
    CL.m_type = "CL";
    CL.m_meanDegree = 2.0 * t_linkSize / t_size;
    return CL;
}
}  // namespace CL

//* Weighted node
struct WNode {
    //* Member variables
    unsigned m_index;
    std::map<unsigned, double> m_neighbors;

    //* Constructor
    WNode() {}
    WNode(const unsigned& t_index) : m_index(t_index) {}

    //* Comparison operator
    bool operator<(const WNode& t_wnode) const {
        return this->m_index < t_wnode.m_index;
    }
};

//* Weighted network
struct WNetwork {
    //* Member variables
    std::string m_type;
    unsigned m_size{0};
    unsigned m_linkSize{0};
    double m_meanDegree{0.0};
    std::vector<std::map<unsigned, double>> m_wadjacency;

    //* Constructor
    WNetwork() {}
    WNetwork(const unsigned& t_size) : m_size(t_size) {
        m_wadjacency.assign(t_size, std::map<unsigned, double>());
        for (unsigned index = 0; index < t_size; ++index) {
            m_wadjacency[index][index] = 0.0;
        }
    }

    //* Reset the whole network
    void clear() {
        //* Clear link size information
        m_linkSize = 0;

        //* Clear adjacency information
        for (unsigned index = 0; index < m_size; ++index) {
            m_wadjacency[index].clear();
            m_wadjacency[index][index] = 0.0;
        }
    }

    //* Whether link is already exists
    double getWeight(const unsigned& t_index1, const unsigned& t_index2) const {
        if (m_wadjacency[t_index1].size() <= m_wadjacency[t_index2].size()) {
            const std::map<unsigned, double> candidates = m_wadjacency[t_index1];
            const auto it = candidates.find(t_index2);
            return it != candidates.end() ? it->second : 0.0;
        } else {
            const std::map<unsigned, double> candidates = m_wadjacency[t_index2];
            const auto it = candidates.find(t_index1);
            return it != candidates.end() ? it->second : 0.0;
        }
    }

    //* Add single link
    void addLink(const unsigned& t_index1, const unsigned& t_index2, const double& t_weight = 1.0) {
        //* Update linksize information
        if (!getWeight(t_index1, t_index2)) {
            ++m_linkSize;
        }

        //* Update weighted adjacency information
        m_wadjacency[t_index1].emplace(t_index2, t_weight);
        m_wadjacency[t_index2].emplace(t_index1, t_weight);
    }

    //* Delete single link
    void deleteLink(const unsigned& t_index1, const unsigned& t_index2) {
        //* Get weight information
        const double weight = getWeight(t_index1, t_index2);

        //* Update linksize information
        if (weight) {
            --m_linkSize;
        }

        //* Update weighted adjacency information
        m_wadjacency[t_index1].erase(t_index2);
        m_wadjacency[t_index2].erase(t_index1);
    }

    //* Print Adjacency matrix with weight
    void printWAdjacency(const std::string& t_fileName = "", const bool t_append = false) const {
        //* File name is given: print into the file
        if (t_fileName.size()) {
            std::ofstream file;
            t_append ? file.open(t_fileName, std::ios_base::app) : file.open(t_fileName);
            for (unsigned index = 0; index < m_size; ++index) {
                for (unsigned neighbor = 0; neighbor < m_size - 1; ++neighbor) {
                    neighbor != index ? file << getWeight(index, neighbor) << "," : file << 0 << ",";
                }
                index != m_size - 1 ? file << getWeight(index, m_size - 1) << "\n" : file << 0 << "\n";
            }
            file.close();
        }
        //* File name is not given: print to terminal
        else {
            for (unsigned index = 0; index < m_size; ++index) {
                for (unsigned neighbor = 0; neighbor < m_size - 1; ++neighbor) {
                    neighbor != index ? std::cout << getWeight(index, neighbor) << "," : std::cout << 0 << ",";
                }
                index != m_size - 1 ? std::cout << getWeight(index, m_size - 1) << "\n" : std::cout << 0 << "\n";
            }
        }
    }

    //* Print Adjacency matrix with bool (except node weight)
    void printAdjacency(const std::string& t_fileName = "", const bool t_append = false) const {
        //* File name is given: print into the file
        if (t_fileName.size()) {
            std::ofstream file;
            t_append ? file.open(t_fileName, std::ios_base::app) : file.open(t_fileName);
            for (unsigned index = 0; index < m_size; ++index) {
                for (unsigned neighbor = 0; neighbor < m_size - 1; ++neighbor) {
                    (m_wadjacency[index].find(neighbor) != m_wadjacency[index].end() && neighbor != index) ? file << 1 << "," : file << 0 << ",";
                }
                (m_wadjacency[index].find(m_size - 1) != m_wadjacency[index].end() && m_size - 1 != index) ? file << 1 << "\n" : file << 0 << "\n";
            }
            file.close();
        }
        //* File name is not given: print to terminal
        else {
            for (unsigned index = 0; index < m_size; ++index) {
                for (unsigned neighbor = 0; neighbor < m_size - 1; ++neighbor) {
                    (m_wadjacency[index].find(neighbor) != m_wadjacency[index].end() && neighbor != index) ? std::cout << 1 << "," : std::cout << 0 << ",";
                }
                (m_wadjacency[index].find(m_size - 1) != m_wadjacency[index].end() && m_size - 1 != index) ? std::cout << 1 << "\n" : std::cout << 0 << "\n";
            }
        }
    }

    //* Prind degree distribution
    void printDegreeDist(const std::string& t_fileName = "", const std::string& t_seperator = ",", const std::string& t_secondSeperator = "\n") const {
        using namespace linearAlgebra;

        //* Generate degree distribution
        std::map<int, double> degreeDist;
        for (unsigned index = 0; index < m_size; ++index) {
            ++degreeDist[m_wadjacency[index].size() - 1];
        }
        degreeDist /= (double)linearAlgebra::accumulate(degreeDist);

        //* File name is given: print into the file
        if (t_fileName.size()) {
            std::ofstream file(t_fileName);
            for (auto it = degreeDist.begin(); it != degreeDist.end(); ++it) {
                file << it->first << t_seperator << it->second << t_secondSeperator;
            }
            file.close();
        }
        //* File name is not given: print to terminal
        else {
            for (auto it = degreeDist.begin(); it != degreeDist.end(); ++it) {
                std::cout << it->first << t_seperator << it->second << t_secondSeperator;
            }
        }
    }

    //* Prind node weight distribution
    void printNodeWeightDist(const std::string& t_fileName = "", const std::string& t_seperator = ",", const std::string& t_secondSeperator = "\n") const {
        using namespace linearAlgebra;

        //* Generate node weight distribution
        std::map<double, double> nodeWeightDist;
        for (unsigned index = 0; index < m_size; ++index) {
            ++nodeWeightDist[m_wadjacency[index].at(index)];
        }
        nodeWeightDist /= linearAlgebra::accumulate(nodeWeightDist);

        //* File name is given: print into the file
        if (t_fileName.size()) {
            std::ofstream file(t_fileName);
            for (auto it = nodeWeightDist.begin(); it != nodeWeightDist.end(); ++it) {
                file << it->first << t_seperator << it->second << t_secondSeperator;
            }
            file.close();
        }
        //* File name is not given: print to terminal
        else {
            for (auto it = nodeWeightDist.begin(); it != nodeWeightDist.end(); ++it) {
                std::cout << it->first << t_seperator << it->second << t_secondSeperator;
            }
        }
    }
};  //* End of struct WNetwork

//* Weighted scale free network by Chung-Lu. Degree distribution and node weight distribution follows same power law
namespace WCL {
unsigned weightSampling(const std::vector<double>& t_weight, const double& t_prob) {
    return (unsigned)(std::lower_bound(t_weight.begin(), t_weight.end(), t_prob * t_weight.back()) - t_weight.begin());
}

WNetwork generate(const unsigned& t_size, const unsigned& t_meanPopulation, const double& t_meanDegree, const double& t_degreeExponent, const double& t_linkWeightExponent, pcg32& t_randomEngine) {
    WNetwork CL(t_size);
    //* Genearting conventional scale free network by Chung-Lu
    const double weightExponent = 1.0 / (t_degreeExponent - 1.0);
    const double correction = weightExponent < 0.5 ? 1.0 : std::pow(10.0 * std::sqrt(2.0) * (1.0 - weightExponent), 1.0 / weightExponent) * std::pow(t_size, 1.0 - 1.0 / (2.0 - weightExponent));

    std::vector<double> weight(t_size, 0.0);
    weight[0] = std::pow(correction, -1.0 * weightExponent);
    for (unsigned i = 1; i < t_size; ++i) {
        weight[i] = weight[i - 1] + std::pow(i + correction, -1.0 * weightExponent);
    }

    std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
    while (CL.m_linkSize < t_meanDegree * t_size / 2.0) {
        unsigned index1, index2;
        do {
            index1 = weightSampling(weight, probabilityDistribution(t_randomEngine));
            index2 = weightSampling(weight, probabilityDistribution(t_randomEngine));
        } while (index1 == index2);
        CL.addLink(index1, index2);
    }

    CL.m_type = "CL";
    CL.m_meanDegree = t_meanDegree;

    //* Add weight to node and link
    const unsigned long long population = t_size * t_meanPopulation;
    for (unsigned index = 0; index < t_size; ++index) {
        //* Add node itself to weighted adjacency matrix
        CL.m_wadjacency[index][index] = 0.0;

        //* Get weight of node and link
        const unsigned degree = CL.m_wadjacency[index].size() - 1;
        const double nodeWeight = degree / (2.0 * CL.m_linkSize);
        std::map<unsigned, double> linkWeight;
        for (auto it = CL.m_wadjacency[index].begin(); it != CL.m_wadjacency[index].end(); ++it) {
            const unsigned neighbor = it->first;
            if (neighbor != index) {
                const int neighborDegree = CL.m_wadjacency[neighbor].size() - 1;
                linkWeight[neighbor] = std::pow(degree * neighborDegree, t_linkWeightExponent);
            }
        }
        const double totalLinkWeight = linearAlgebra::accumulate(linkWeight);

        //* Add weight to weighted adjacency matrix
        for (auto it = CL.m_wadjacency[index].begin(); it != CL.m_wadjacency[index].end(); ++it) {
            if (it->first != index) {
                CL.m_wadjacency[index][it->first] = linkWeight[it->first] / totalLinkWeight;
            } else {
                CL.m_wadjacency[index][it->first] = std::floor(population * nodeWeight + 0.5);
            }
        }
    }
    return CL;
}
}  // namespace WCL