#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <deque>
#include <cmath>
#include <set>
#include <utility> //* pair
#include <algorithm>
#include <random>
#include <map>
#include <fstream>

#include "pcg_random.hpp"

typedef int Node;
typedef std::pair<Node, double> WNode; //* weighted node

struct Network
{
    //* Member variables
    int m_size{0};
    int m_linkSize{0};
    std::vector<std::set<Node>> m_adjacentMatrix;
    pcg32 m_randomEngine;

    //* constructor
    Network() {}

    Network(const int &t_size)
    : m_size(t_size)
    {
        m_adjacentMatrix.resize(t_size);
    }

    //* Simple Get Functions
    const int size() const { return m_size; }
    const int linkSize() const { return m_linkSize; }
    const std::set<Node> adjacent(const Node &t_node) const { return m_adjacentMatrix[t_node]; }
    std::vector<std::set<Node>> adjacent() const { return m_adjacentMatrix; }

    //* Show full information of network
    void show(const int& debugMode = 0, const std::string& t_outFileName = "temp.txt")
    {
        if (debugMode == 1 && t_outFileName == "temp.txt")
        {
            for (Node node = 0; node < m_size; ++node)
            {
                for (Node neighbor : m_adjacentMatrix[node])
                {
                    if (neighbor > node)
                    {
                        std::cout << node << "," << neighbor << "\n";
                    }
                }
            }
        }
        else if (debugMode == 1){
            std::ofstream outFile(t_outFileName);
            for (Node node = 0; node < m_size; ++node)
            {
                for (Node neighbor : m_adjacentMatrix[node])
                {
                    if (neighbor > node)
                    {
                        outFile << node << "," << neighbor << "\n";
                    }
                }
            }
        }
        else
        {
            std::cout << "Total number of nodes : " << m_size << ", links : " << m_linkSize << "\n";
            for (Node node = 0; node < m_size; ++node)
            {
                for (Node neighbor : m_adjacentMatrix[node])
                {
                    std::cout << "(" << node << "," << neighbor << ") ";
                }
                m_adjacentMatrix[node].size() ? std::cout << "\n" : std::cout << "()\n";
            }
        }
    }

    //* Show degree distribution
    void showDegree(const std::string& t_outFileName = "temp.txt")
    {
        std::map<int, int> degreeDistribution;
        for (int node = 0; node < m_size; ++node)
        {
            ++degreeDistribution[m_adjacentMatrix[node].size()];
        }
        if (t_outFileName == "temp.txt"){
            for (auto it = degreeDistribution.begin(); it != degreeDistribution.end(); ++it)
            {
                std::cout << it->first << "," << it->second << "\n";
            }
        }
        else
        {
            std::ofstream outFile(t_outFileName);
            for (auto it = degreeDistribution.begin(); it!= degreeDistribution.end(); ++it)
            {
                outFile << it->first << "," << it->second <<"\n";
            }
        }


    }

    //* Whether link is already exists
    bool linkExists(const Node &t_node1, const Node &t_node2)
    {
        std::set<Node> candidates;
        if (m_adjacentMatrix[t_node1].size() < m_adjacentMatrix[t_node2].size())
        {
            candidates = m_adjacentMatrix[t_node1];
            auto it = std::find(candidates.begin(), candidates.end(), t_node2);
            return it != candidates.end() ? true : false;
        }
        else
        {
            candidates = m_adjacentMatrix[t_node2];
            auto it = std::find(candidates.begin(), candidates.end(), t_node1);
            return it != candidates.end() ? true : false;
        }
    }

    //* Add single link
    void addLink(const Node &t_node1, const Node &t_node2)
    {
        //* WNode already exists
        if (linkExists(t_node1, t_node2))
        {
            return;
        }

        //* Update node information
        m_adjacentMatrix[t_node1].insert(t_node2);
        m_adjacentMatrix[t_node2].insert(t_node1);

        //* Update network information
        ++m_linkSize;
    }

    //* Delete single link
    void deleteLink(const Node &t_node1, const Node &t_node2)
    {
        //* WNode does not exists
        if (!linkExists(t_node1, t_node2))
        {
            return;
        }

        //* Update node information
        m_adjacentMatrix[t_node1].erase(t_node2);
        m_adjacentMatrix[t_node2].erase(t_node1);

        //* Update network information
        --m_linkSize;
    }

    //* Reset network
    void clear()
    {
        m_adjacentMatrix.clear();
        m_linkSize = 0;
    }

    //* Rewire single link of one node
    void rewire(const Node &t_node)
    {
        //* Delete link to old neighbor
        auto iter = m_adjacentMatrix[t_node].begin();
        std::uniform_int_distribution<int> neighborDistribution(0, m_adjacentMatrix[t_node].size() - 1);
        std::advance(iter, neighborDistribution(m_randomEngine));
        deleteLink(t_node, *iter);

        //* Add link to new neighbor
        std::uniform_int_distribution<int> nodeDistribution(0, m_size - 1);
        Node newNeighbor;
        do
        {
            newNeighbor = nodeDistribution(m_randomEngine);
        } while (linkExists(t_node, newNeighbor));
        addLink(t_node, newNeighbor);
    }
}; //! End of Network Struct

//* Erdos-Renyi network
struct ER_Network : public Network
{
public:
    //* constructor
    ER_Network() {}

    ER_Network(const int &t_size) : Network(t_size) {}

    ER_Network(const int &t_size, const double &t_probability)
        : Network(t_size)
    {
        m_randomEngine.seed((std::random_device())());
        generate(t_probability);
    }

    ER_Network(const int &t_size, const double &t_probability, const pcg32 &t_randomEngine)
        : Network(t_size)
    {
        m_randomEngine = t_randomEngine;
        generate(t_probability);
    }

    ER_Network(const int &t_size, const int &t_linkSize)
        : Network(t_size)
    {
        m_randomEngine.seed((std::random_device())());
        generate(t_linkSize);
    }

    ER_Network(const int &t_size, const int &t_linkSize, const pcg32 &t_randomEngine)
        : Network(t_size)
    {
        m_randomEngine = t_randomEngine;
        generate(t_linkSize);
    }

    //* Generator
    void generate(const double &t_probability)
    {
        std::uniform_real_distribution<double> probabilityDistribution(0, 1);
        for (Node node1 = 0; node1 < m_size; ++node1)
        {
            for (Node node2 = node1 + 1; node2 < m_size; ++node2)
            {
                if (probabilityDistribution(m_randomEngine) < t_probability)
                {
                    addLink(node1, node2);
                }
            }
        }
    }
    void generate(const int &t_linkSize)
    {
        std::uniform_int_distribution<int> nodeDistribution(0, m_size - 1);
        while (m_linkSize < t_linkSize)
        {
            Node node1, node2;
            do
            {
                node1 = nodeDistribution(m_randomEngine);
                node2 = nodeDistribution(m_randomEngine);
            } while (node1 == node2);
            addLink(node1, node2);
        }
    }
}; //* End of struct ER

//* Random Regular network
struct RR_Network : public Network
{
private:
    //* Member variables
    double m_meanDegree;

public:
    //* Constructor
    RR_Network(const int &t_size, const int &t_meanDegree)
        : Network(t_size), m_meanDegree(t_meanDegree)
    {
        m_randomEngine.seed((std::random_device())());
        generate();
    }

    RR_Network(const int &t_size, const int &t_meanDegree, const pcg32 &t_randomEngine)
        : Network(t_size), m_meanDegree(t_meanDegree)
    {
        m_randomEngine = t_randomEngine;
        generate();
    }

    //* Generator
    void generate()
    {
        std::deque<int> stubs;
        for (Node node = 0; node < m_size; ++node)
        {
            for (int degree = 0; degree < m_meanDegree; ++degree)
            {
                stubs.emplace_back(node);
            }
        }
        std::shuffle(stubs.begin(), stubs.end(), m_randomEngine);
        while (!stubs.empty())
        {
            Node node1 = stubs.front();
            stubs.pop_front();
            Node node2 = stubs.front();
            stubs.pop_front();
            m_adjacentMatrix[node1].insert(node2);
            m_adjacentMatrix[node2].insert(node1);
        }

        m_linkSize = m_size*m_meanDegree/2;
    }
}; //* End of struct RR_Network

//* Static scale free network with Robin hood Algorithm
struct SF_Network : public Network{
private:
    //* Member variables
    double m_degreeExponent;

public:
    //* Constructor
    SF_Network(const int &t_size, const int &t_linkSize, const double &t_degreeExponent)
        : Network(t_size), m_degreeExponent(t_degreeExponent)
    {
        m_randomEngine.seed((std::random_device())());
        generate(t_linkSize);
    }

    SF_Network(const int &t_size, const int &t_linkSize, const double &t_degreeExponent, const pcg32 &t_randomEngine)
        : Network(t_size), m_degreeExponent(t_degreeExponent)
    {
        m_randomEngine = t_randomEngine;
        generate(t_linkSize);
    }

    //* Generator
    void generate(const int &t_linkSize)
    {
        const double weightExponent = 1.0 / (m_degreeExponent - 1.0);
        while (m_linkSize < t_linkSize)
        {
            Node node1, node2;
            do
            {
                node1 = randomPowerLawDistribution(1, m_size, -weightExponent) - 1;
                node2 = randomPowerLawDistribution(1, m_size, -weightExponent) - 1;
            } while (node1 == node2);
            addLink(node1, node2);
        }
    }

    int randomPowerLawDistribution(const int &t_lower, const int &t_upper, const double &t_exponent)
    {
        std::uniform_real_distribution<double> realDistribution(0, 1);
        return std::pow((std::pow(t_upper + 0.5, t_exponent + 1) - std::pow(t_lower - 0.5, t_exponent + 1)) * realDistribution(m_randomEngine) + std::pow(t_lower - 0.5, t_exponent + 1), 1.0 / (t_exponent + 1)) + 0.5;
    }
}; //* End of struct SF_Network

//* Scale free network with Chung Lu Algorithm
struct CL_Network : public Network{
private:
    double m_alpha;
    std::vector<double> m_weight;

public:
    //* Constructor
    CL_Network(const int& t_size, const int& t_linkSize, const double& t_degreeExponent)
        : Network(t_size)
    {
        m_randomEngine.seed((std::random_device())());
        m_alpha = 1.0/(t_degreeExponent-1.0);
        generate(t_linkSize);

    }

    CL_Network(const int& t_size, const int& t_linkSize, const double& t_degreeExponent, const pcg32& t_randomEngine)
        : Network(t_size)
    {
        m_randomEngine = t_randomEngine;
        m_alpha = 1.0/(t_degreeExponent-1.0);
        generate(t_linkSize);
    }

    void generate(const int& t_linkSize){
        const double correction = m_alpha < 0.5 ? 1.0 : std::pow(10.0*std::sqrt(2.0)*(1.0-m_alpha), 1.0/m_alpha) * std::pow(m_size, 1.0-1.0/(2.0-m_alpha));
        m_weight.assign(m_size, 0.0);
        m_weight[0] = std::pow(correction, -1.0*m_alpha);

        for (int i=1; i<m_size; ++i){
            m_weight[i] = m_weight[i-1]+std::pow(i+correction, -1.0*m_alpha);
        }
        std::uniform_real_distribution<double> realDistribution(0, 1);
        while (m_linkSize < t_linkSize){
            Node node1, node2;
            do {
                node1 = weightSampling(realDistribution(m_randomEngine));
                node2 = weightSampling(realDistribution(m_randomEngine));
            } while (node1 == node2);
            addLink(node1, node2);
        }
    }

    Node weightSampling(const double& t_prob){
        return (Node) (std::lower_bound(m_weight.begin(), m_weight.end(), t_prob*m_weight.back())-m_weight.begin());
    }
};//* End of struct CL_Network

//* Network class with merging cluster by Newman-Ziff algorithm
struct NZ_Network{
private:
    //* Size of nodes and links
    int m_size{0};
    int m_linkSize{0};

    //* Maximum Cluster and Second Maximum Cluster
    int m_maximumClusterSize{1};
    int m_secondMaximumClusterSize{1};
    int m_deltaMaximumClusterSize{0};

    //* m_parent[node] : parent of each 'node'
    std::vector<Node> m_parent;

    //* m_sortedCluster[size] : number of cluster of 'size'
    std::map<int,int> m_sortedCluster;

    //* m_birth[root] : birth time of each 'root'
    //* changedAge[root] : {age, size} of cluster with 'root'
    std::vector<Node> m_birth;
    std::vector<std::pair<int, int>> m_changedAge;

public:
    //* Constructor
    NZ_Network() {}

    NZ_Network(const int &t_size)
    : m_size(t_size)
    {
        //! Make every node to root node with size 1
        m_parent.resize(t_size,-1);

        //! Initialize Sorted Cluster
        m_sortedCluster[1] = t_size;

        //! Initialize birth time with 0 and changedAge
        m_birth.resize(t_size);
        m_changedAge.resize(2,std::pair<int, int> {0,0});
    }

    //* Simple get functions
    int getMaximumClusterSize() const {return m_maximumClusterSize;}
    int getSecondMaximumClusterSize() const {return m_secondMaximumClusterSize;}
    int getClusterSize(const Node &t_root) const {return -m_parent[t_root];}
    int getDeltaMaximumClusterSize() const {return m_deltaMaximumClusterSize;}
    std::vector<std::pair<int,int>> getChangedAge() const {return m_changedAge;}
    std::map<int,int> getSortedCluster(const int &excludeNum=1) const
    {
        std::map<int,int> result=m_sortedCluster;

        //! exclude maximum cluster
        --result[m_maximumClusterSize];

        //! exclude second giant
        if (excludeNum-1){
            --result[m_secondMaximumClusterSize];
        }
        return result;
    }

    //* get the root of input node
    Node getRoot(const Node &t_node){
        //* t_node is root
        if (m_parent[t_node] < 0){
            return t_node;
        }
        //* recursively find node
        return m_parent[t_node] = getRoot(m_parent[t_node]);
    }

    //* Merge two clusters
    void merge(const Node &t_root1, const Node &t_root2){
        //! update link size
        m_linkSize++;

        //! Save age and reset birth time
        m_changedAge[0] = std::pair<int,int>{m_linkSize-m_birth[t_root1], -m_parent[t_root1]};
        m_changedAge[1] = std::pair<int,int>{m_linkSize-m_birth[t_root2], -m_parent[t_root2]};
        m_birth[t_root1] = m_linkSize;
        m_birth[t_root2] = m_linkSize;

        //! Get the size of each clusters
        const int size1 = -m_parent[t_root1];
        const int size2 = -m_parent[t_root2];
        const int newSize = size1+size2;

        //! Update Parent
        m_parent[t_root1] -= size2;
        m_parent[t_root2] = t_root1;

        //! Update Sorted Cluster
        --m_sortedCluster[size1];
        if (m_sortedCluster[size1]==0){
            m_sortedCluster.erase(size1);
        }
        --m_sortedCluster[size2];
        if (m_sortedCluster[size2]==0){
            m_sortedCluster.erase(size2);
        }
        ++m_sortedCluster[newSize];

        //! find maximum cluster
        if (m_maximumClusterSize < newSize){
            m_deltaMaximumClusterSize = newSize-m_maximumClusterSize;
            m_maximumClusterSize = newSize;
        }
        else{
            m_deltaMaximumClusterSize = 0;
        }

        //! find second giant cluster
        for (auto it=m_sortedCluster.rbegin(); it!=m_sortedCluster.rend(); ++it){
            if (it->first!=m_maximumClusterSize || it->second>1){
                m_secondMaximumClusterSize=it->first;
                break;
            }
        }
    }

    //* Mean cluster size
    double getMeanClusterSize() const{
        const int firstMoment = m_size-m_maximumClusterSize;
        double secondMoment = 0;
        for (auto it=m_sortedCluster.begin(); it!=m_sortedCluster.end(); ++it){
            secondMoment += pow(it->first,2)*it->second;
        }
        //! exclude infinite size cluster
        secondMoment -= pow(m_maximumClusterSize,2);

        return secondMoment/firstMoment;
    }
};



//* Weighted Network
struct WNetwork
{
    //* Member variables
    int m_size{0};
    int m_linkSize{0};
    std::vector<std::set<WNode>> m_adjacentMatrix;
    pcg32 m_randomEngine;

    //* constructor
    WNetwork() {}

    WNetwork(const int &t_size)
        : m_size(t_size)
    {
        m_adjacentMatrix.resize(t_size);
    }

    //* Simple get function
    const int size() const { return m_size; }
    const int linkSize() const { return m_linkSize; }
    const std::set<WNode> adjacent(const Node &t_node) const { return m_adjacentMatrix[t_node];}
    const std::vector<std::set<WNode>> adjacent() const { return m_adjacentMatrix; }

    //* Show full information of network
    void show(const int debugMode = 0)
    {
        if (debugMode == 1)
        {
            for (Node node = 0; node < m_size; ++node)
            {
                for (WNode link : m_adjacentMatrix[node])
                {
                    if (link.first > node)
                    {
                        std::cout << node << "," << link.first << "," << link.second << "\n";
                    }
                }
            }
        }
        else
        {
            std::cout << "Total number of nodes : " << m_size << ", links : " << m_linkSize << "\n";
            for (Node node = 0; node < m_size; ++node)
            {
                for (WNode link : m_adjacentMatrix[node])
                {
                    std::cout << "(" << node << "," << link.first << ")" << link.second << ", ";
                }
                m_adjacentMatrix[node].size() ? std::cout << "\n" : std::cout << "()\n";
            }
        }
    }

    //* Show degree distribution
    void showDegree()
    {
        std::map<int, int> degreeDistribution;
        for (int node = 0; node < m_size; ++node)
        {
            ++degreeDistribution[m_adjacentMatrix[node].size()];
        }
        for (auto it = degreeDistribution.begin(); it != degreeDistribution.end(); ++it)
        {
            // std::cout<<"("<<it->first<<","<<it->second<<"),";
            std::cout << it->first << "," << it->second << "\n";
        }
        std::cout << "\n";
    }

    //* Whether link is already exists
    bool linkExists(const Node &t_node1, const Node &t_node2)
    {
        std::set<Node> candidates;
        if (m_adjacentMatrix[t_node1].size() < m_adjacentMatrix[t_node2].size())
        {
            for (WNode link : m_adjacentMatrix[t_node1])
            {
                candidates.insert(link.first);
            }
            auto it = std::find(candidates.begin(), candidates.end(), t_node2);
            return it != candidates.end() ? true : false;
        }
        else
        {
            for (WNode link : m_adjacentMatrix[t_node2])
            {
                candidates.insert(link.first);
            }
            auto it = std::find(candidates.begin(), candidates.end(), t_node1);
            return it != candidates.end() ? true : false;
        }
    }

    //* Add single link
    void addLink(const Node &t_node1, const Node &t_node2, const double &t_weight)
    {
        //* link already exists
        if (linkExists(t_node1, t_node2))
        {
            return;
        }

        //* Update node information
        m_adjacentMatrix[t_node1].insert({t_node2, t_weight});
        m_adjacentMatrix[t_node2].insert({t_node1, t_weight});

        //* Update network information
        ++m_linkSize;
    }

    //* Delete single link
    void deleteLink(const Node &t_node1, const Node &t_node2)
    {
        //* WNode does not exists
        if (!linkExists(t_node1, t_node2))
        {
            return;
        }

        //* Update node information
        const std::set<WNode> neighbors = m_adjacentMatrix[t_node1];
        double weight;
        for (WNode neighbor : neighbors){
            if (neighbor.first == t_node2){
                weight = neighbor.second;
                break;
            }
        }
        m_adjacentMatrix[t_node1].erase({t_node2, weight});
        m_adjacentMatrix[t_node2].erase({t_node1, weight});

        //* Update network information
        --m_linkSize;
    }

    //* Reset network
    void clear()
    {
        m_adjacentMatrix.clear();
        m_linkSize = 0;
    }

    //* Rewire single link of one node with preserving weight
    void rewire(const Node &t_node)
    {
        //* Delete link to old neighbor
        auto iter = m_adjacentMatrix[t_node].begin();
        std::uniform_int_distribution<int> neighborDistribution(0, m_adjacentMatrix[t_node].size() - 1);
        std::advance(iter, neighborDistribution(m_randomEngine));
        const WNode deletedLink = *iter;
        deleteLink(t_node, deletedLink.first);

        //* Add link to new neighbor
        std::uniform_int_distribution<int> nodeDistribution(0, m_size - 1);
        Node newNeighbor;
        do
        {
            newNeighbor = nodeDistribution(m_randomEngine);
        } while (linkExists(t_node, newNeighbor));
        addLink(t_node, newNeighbor, deletedLink.second);
    }
}; //* End of struct WNetwork

struct ER_WNetwork : public WNetwork
{
private:
    //* Member variables
    std::normal_distribution<double> m_weightDistribution;

public:
    //* constructor
    ER_WNetwork() {}

    ER_WNetwork(const int &t_size) : WNetwork(t_size) {}

    ER_WNetwork(const int &t_size, const double &t_probability, const double &t_meanWeight, const double &t_stdWeight)
        : WNetwork(t_size)
    {
        m_randomEngine.seed((std::random_device())());
        generate(t_probability, t_meanWeight, t_stdWeight);
    }

    ER_WNetwork(const int &t_size, const double &t_probability, const double &t_meanWeight, const double &t_stdWeight, const pcg32 &t_randomEngine)
        : WNetwork(t_size)
    {
        m_randomEngine = t_randomEngine;
        generate(t_probability, t_meanWeight, t_stdWeight);
    }

    ER_WNetwork(const int &t_size, const int &t_linkSize, const double &t_meanWeight, const double &t_stdWeight)
        : WNetwork(t_size)
    {
        m_randomEngine.seed((std::random_device())());
        generate(t_linkSize, t_meanWeight, t_stdWeight);
    }

    ER_WNetwork(const int &t_size, const int &t_linkSize, const double &t_meanWeight, const double &t_stdWeight, const pcg32 &t_randomEngine)
        : WNetwork(t_size)
    {
        m_randomEngine = t_randomEngine;
        generate(t_linkSize, t_meanWeight, t_stdWeight);
    }

    //* Generator
    void generate(const double &t_probability, const double &t_meanWeight, const double &t_stdWeight)
    {
        m_weightDistribution.param(std::normal_distribution<double>::param_type(t_meanWeight, t_stdWeight));
        std::uniform_real_distribution<double> probabilityDistribution(0, 1);
        for (Node node1 = 0; node1 < m_size; ++node1)
        {
            for (Node node2 = node1 + 1; node2 < m_size; ++node2)
            {
                if (probabilityDistribution(m_randomEngine) < t_probability)
                {
                    addLink(node1, node2, m_weightDistribution(m_randomEngine));
                }
            }
        }
    }
    void generate(const int &t_linkSize, const double &t_meanWeight, const double &t_stdWeight)
    {
        m_weightDistribution.param(std::normal_distribution<double>::param_type(t_meanWeight, t_stdWeight));
        std::uniform_int_distribution<int> nodeDistribution(0, m_size - 1);
        while (m_linkSize < t_linkSize)
        {
            Node node1, node2;
            do
            {
                node1 = nodeDistribution(m_randomEngine);
                node2 = nodeDistribution(m_randomEngine);
            } while (node1 == node2);
            addLink(node1, node2, m_weightDistribution(m_randomEngine));
        }
    }
}; //* End of struct ER_WNetwork

struct SF_WNetwork : public WNetwork
{
private:
    //* Member variables
    double m_degreeExponent;
    std::normal_distribution<double> m_weightDistribution;

public:
    //* Constructor
    SF_WNetwork(const int &t_size, const int &t_linkSize, const double &t_degreeExponent, const double &t_meanWeight, const double &t_stdWeight)
        : WNetwork(t_size), m_degreeExponent(t_degreeExponent)
    {
        m_randomEngine.seed((std::random_device())());
        generate(t_linkSize, t_degreeExponent, t_meanWeight, t_stdWeight);
    }

    SF_WNetwork(const int &t_size, const int &t_linkSize, const double &t_degreeExponent, const double &t_meanWeight, const double &t_stdWeight, const pcg32 &t_randomEngine)
        : WNetwork(t_size), m_degreeExponent(t_degreeExponent)
    {
        m_randomEngine = t_randomEngine;
        generate(t_linkSize, t_degreeExponent, t_meanWeight, t_stdWeight);
    }

    //* Generator
    void generate(const int &t_linkSize, const double &t_degreeExponent, const double &t_meanWeight, const double &t_stdWeight)
    {
        m_weightDistribution.param(std::normal_distribution<double>::param_type(t_meanWeight, t_stdWeight));
        const double weightExponent = 1.0 / (m_degreeExponent - 1.0);
        while (m_linkSize < t_linkSize)
        {
            Node node1, node2;
            do
            {
                node1 = randomPowerLawDistribution(1, m_size, weightExponent) - 1;
                node2 = randomPowerLawDistribution(1, m_size, weightExponent) - 1;
            } while (node1 == node2);
            addLink(node1, node2, m_weightDistribution(m_randomEngine));
        }
    }

    int randomPowerLawDistribution(const int &t_lower, const int &t_upper, const double &t_exponent)
    {
        std::uniform_real_distribution<double> realDistribution(0, 1);
        return std::pow((std::pow(t_upper + 0.5, t_exponent + 1) - std::pow(t_lower - 0.5, t_exponent + 1)) * realDistribution(m_randomEngine) + std::pow(t_lower - 0.5, t_exponent + 1), 1.0 / (t_exponent + 1)) + 0.5;
    }
};//* end of struct SF_WNetwork