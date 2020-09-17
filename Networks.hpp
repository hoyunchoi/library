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
typedef std::pair<Node, double> WNode; //* weighted link with node,weight

struct Network
{
    //* Member variables
    int m_size{0};
    int m_linkSize{0};
    std::vector<std::set<Node>> m_adjacentMatrix;

    //* constructor
    Network() {}

    Network(const int &t_size)
    : m_size(t_size){
        m_adjacentMatrix.resize(t_size);
    }

    //* Simple Get Functions
    const int size() const { return m_size;}
    const int linkSize() const { return m_linkSize;}
    const std::set<Node> adjacent(const Node &t_node) const { return m_adjacentMatrix[t_node];}
    std::vector<std::set<Node>> adjacent() const { return m_adjacentMatrix;}

    //* Show full information of network
    void printAdjacent(const std::string& t_outFileName = "default"){
        if (t_outFileName=="default"){
            std::cout << "Total number of nodes: " << m_size << ", links: "<< m_linkSize << "\n";
            for (Node node=0; node<m_size; ++node){
                std::cout << node << ":";
                for (Node neighbor : m_adjacentMatrix[node]){
                    std::cout << neighbor << ", ";
                }
                std::cout << "\n";
            }
        }
        else{
            std::ofstream outFile(t_outFileName);
            outFile << "Total number of nodes: " << m_size << ", links: "<< m_linkSize << "\n";
            for (Node node=0; node<m_size; ++node){
                outFile << node << ":";
                for (Node neighbor : m_adjacentMatrix[node]){
                    outFile << neighbor << ", ";
                }
                outFile << "\n";
            }
        }
    }

    //* Show degree distribution
    void printDegree(const std::string& t_outFileName = "default"){
        std::map<int, int> degreeDistribution;
        for (int node = 0; node < m_size; ++node){
            ++degreeDistribution[m_adjacentMatrix[node].size()];
        }
        if (t_outFileName == "default"){
            for (auto it = degreeDistribution.begin(); it != degreeDistribution.end(); ++it){
                std::cout << it->first << "," << it->second << "\n";
            }
        }
        else{
            std::ofstream outFile(t_outFileName);
            for (auto it = degreeDistribution.begin(); it!= degreeDistribution.end(); ++it){
                outFile << it->first << "," << it->second <<"\n";
            }
        }


    }

    //* Whether link is already exists
    bool linkExists(const Node &t_node1, const Node &t_node2){
        std::set<Node> candidates;
        if (m_adjacentMatrix[t_node1].size() < m_adjacentMatrix[t_node2].size()){
            candidates = m_adjacentMatrix[t_node1];
            auto it = std::find(candidates.begin(), candidates.end(), t_node2);
            return it != candidates.end() ? true : false;
        }
        else{
            candidates = m_adjacentMatrix[t_node2];
            auto it = std::find(candidates.begin(), candidates.end(), t_node1);
            return it != candidates.end() ? true : false;
        }
    }

    //* Add single link
    void addLink(const Node &t_node1, const Node &t_node2){
        //* Update node information
        m_adjacentMatrix[t_node1].insert(t_node2);
        m_adjacentMatrix[t_node2].insert(t_node1);

        //* Update network information
        ++m_linkSize;
    }

    //* Delete single link
    void deleteLink(const Node &t_node1, const Node &t_node2){
        //* Update node information
        m_adjacentMatrix[t_node1].erase(t_node2);
        m_adjacentMatrix[t_node2].erase(t_node1);

        //* Update network information
        --m_linkSize;
    }

    //* Reset network
    void clear(){
        m_adjacentMatrix.clear();
        m_linkSize = 0;
    }

    //* Rewire single link of one node
    void rewire(const Node& t_node, const Node& t_oldNeighbor, const Node& t_newNeighbor){
        //* Delete link to old neighbor
        deleteLink(t_node, t_oldNeighbor);

        //* Add link to new neighbor
        addLink(t_node, t_newNeighbor);
    }
}; //* End of Network Struct

//* Erdos-Renyi network
namespace ER{
    Network generate(const int& t_size, const double& t_probability, pcg32& t_randomEngine){
        Network ER(t_size);
        std::uniform_real_distribution<double> probabilityDistribution(0, 1);
        for (Node node1 = 0; node1 < ER.size(); ++node1){
            for (Node node2 = node1+1; node2 < ER.size(); ++node2){
                if (probabilityDistribution(t_randomEngine)){
                    ER.addLink(node1, node2);
                }
            }
        }
        return ER;
    }

    Network generate(const int& t_size, const int& t_linkSize, pcg32& t_randomEngine){
        Network ER(t_size);
        std::uniform_int_distribution<int> nodeDistribution(0, ER.size() - 1);
        while (ER.linkSize() < t_linkSize){
            Node node1, node2;
            do{
                node1 = nodeDistribution(t_randomEngine);
                node2 = nodeDistribution(t_randomEngine);
            } while(node1 == node2);
            ER.addLink(node1, node2);
        }
        return ER;
    }
}//* End of namespace ER

namespace RR{
    Network generate(const int& t_size, const int& t_meanDegree, pcg32& t_randomEngine){
        Network RR(t_size);
        std::deque<int> stubs;
        for (Node node=0; node<RR.size(); ++node){
            for (int degree=0; degree<t_meanDegree; ++degree){
                stubs.emplace_back(node);
            }
        }
        std::shuffle(stubs.begin(), stubs.end(), t_randomEngine);
        while (!stubs.empty()){
            Node node1 = stubs.front();
            stubs.pop_front();
            Node node2 = stubs.front();
            stubs.pop_front();
            RR.addLink(node1, node2);
        }
        return RR;
    }

}//* End of namespace RR

namespace SF{
    Node randomPowerLawDistribution(const int& t_lower, const int& t_upper, const double& t_exponent, const double& t_prob){
        std::uniform_real_distribution<double> realDistribution(0.0, 1.0);
        return std::pow((std::pow(t_upper + 0.5, t_exponent + 1) - std::pow(t_lower - 0.5, t_exponent + 1)) * t_prob + std::pow(t_lower - 0.5, t_exponent + 1), 1.0 / (t_exponent + 1)) + 0.5;
    }

    Network generate(const int& t_size, const int& t_linkSize, const double& t_degreeExponent, pcg32& t_randomEngine){
        Network SF(t_size);
        const double weightExponent = 1.0/(t_degreeExponent-1.0);
        std::uniform_real_distribution<double> realDistribution(0.0, 1.0);
        while(SF.linkSize() < t_linkSize){
            Node node1, node2;
            do{
                node1 = randomPowerLawDistribution(1, SF.size(), -weightExponent, realDistribution(t_randomEngine))-1;
                node2 = randomPowerLawDistribution(1, SF.size(), -weightExponent, realDistribution(t_randomEngine))-1;
            } while(node1 == node2);
            SF.addLink(node1, node2);
        }
        return SF;
    }
}//* End of namespace SF

namespace CL{
    Node weightSampling (const std::vector<double>& t_weight, const double& t_prob){
        return (Node) (std::lower_bound(t_weight.begin(), t_weight.end(), t_prob*t_weight.back())-t_weight.begin());
    }

    Network generate(const int& t_size, const int& t_linkSize, const double& t_degreeExponent, pcg32& t_randomEngine){
        Network CL(t_size);
        const double weightExponent = 1.0/(t_degreeExponent-1.0);
        const double correction = weightExponent < 0.5 ? 1.0 : std::pow(10.0*std::sqrt(2.0)*(1.0-weightExponent), 1.0/weightExponent) * std::pow(CL.size(), 1.0-1.0/(2.0-weightExponent));

        std::vector<double> weight(CL.size(), 0.0);
        weight[0] = std::pow(correction, -1.0*weightExponent);
        for (int i=1; i<CL.size(); ++i){
            weight[i] = weight[i-1]+std::pow(i+correction, -1.0*weightExponent);
        }

        std::uniform_real_distribution<double> realDistribution(0.0, 1.0);
        while (CL.linkSize() < t_linkSize){
            Node node1, node2;
            do{
                node1 = weightSampling(weight, realDistribution(t_randomEngine));
                node2 = weightSampling(weight, realDistribution(t_randomEngine));
            } while(node1 == node2);
        }
        return CL;
    }
}//* End of namespace CL


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
    std::map<int, int> m_sortedCluster;

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
        //* Make every node to root node with size 1
        m_parent.resize(t_size,-1);

        //* Initialize Sorted Cluster
        m_sortedCluster[1] = t_size;

        //* Initialize birth time with 0 and changedAge
        m_birth.resize(t_size);
        m_changedAge.resize(2,std::pair<int, int> {0,0});
    }

    //* Simple get functions
    int getMaximumClusterSize() const {return m_maximumClusterSize;}
    int getClusterSize(const Node &t_root) const {return -m_parent[t_root];}
    int getDeltaMaximumClusterSize() const {return m_deltaMaximumClusterSize;}
    std::vector<std::pair<int,int>> getChangedAge() const {return m_changedAge;}


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

    }

    //* Calculate second Maximum Cluster Size
    void processSecondMaximumClusterSize(){
        for (auto it=m_sortedCluster.rbegin(); it != m_sortedCluster.rend(); ++it){
            if (it->first != m_maximumClusterSize || it->second>1){
                m_secondMaximumClusterSize = it->first;
                break;
            }
        }
    }

    //* Get second Maximum Cluster Size
    int getSecondMaximumClusterSize(){
        processSecondMaximumClusterSize();
        return m_secondMaximumClusterSize;
    }

    //* Get sorted Cluster
    std::map<int,int> getSortedCluster(const int &excludeNum=1){
        std::map<int,int> result=m_sortedCluster;

        //! exclude maximum cluster
        --result[m_maximumClusterSize];

        //! exclude second giant
        if (excludeNum-1){
            processSecondMaximumClusterSize();
            --result[m_secondMaximumClusterSize];
        }
        return result;
    }

    //* Get Mean cluster size
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
struct WNetwork{
    //* Member variables
    int m_size{0};
    int m_linkSize{0};
    std::vector<std::set<WNode>> m_adjacentMatrix;

    //* Constructor
    WNetwork() {}

    WNetwork(const int& t_size)
    : m_size(t_size){
        m_adjacentMatrix.resize(t_size);
    }

    //* Simple Get functions
    const int size() const {return m_size;}
    const int linkSize() const {return m_linkSize;}
    const std::set<WNode> adjacent(const Node& t_node) const {return m_adjacentMatrix[t_node];}
    std::vector<std::set<WNode>> adjacent() const {return m_adjacentMatrix;}

    //* Get functions
    double weight(const Node& t_node1, const Node& t_node2){
        double weight = 0.0;
        if (m_adjacentMatrix[t_node1].size() < m_adjacentMatrix[t_node2].size()){
            for (WNode neighbor : m_adjacentMatrix[t_node1]){
                if (neighbor.first == t_node2){
                    weight = neighbor.second;
                    break;
                }
            }
        }
        else{
            for (WNode neighbor : m_adjacentMatrix[t_node2]){
                if (neighbor.first == t_node1){
                    weight = neighbor.second;
                    break;
                }
            }
        }
        return weight;
    }

    //* Whether link is already exists
    bool linkExists(const Node& t_node1, const Node& t_node2, const double& t_weight){
        if (m_adjacentMatrix[t_node1].size() < m_adjacentMatrix[t_node2].size()){
            const std::set<WNode> candidates = m_adjacentMatrix[t_node1];
            auto it = std::find(candidates.begin(), candidates.end(), std::make_pair(t_node2, t_weight));
            return it != candidates.end() ? true : false;
        }
        else{
            const std::set<WNode> candidates = m_adjacentMatrix[t_node2];
            auto it = std::find(candidates.begin(), candidates.end(), std::make_pair(t_node1, t_weight));
            return it != candidates.end() ? true : false;
        }
    }

    //* Add single link with given weight
    void addLink(const Node& t_node1, const Node& t_node2, const double& t_weight){
        //* update adjacent information
        m_adjacentMatrix[t_node1].insert(std::make_pair(t_node2, t_weight));
        m_adjacentMatrix[t_node2].insert(std::make_pair(t_node1, t_weight));

        //* update link size
        ++m_linkSize;
    }

    //* Delete single link with given weight
    void deleteLink(const Node& t_node1, const Node& t_node2, const double& t_weight){
        //* Update adjacent information
        m_adjacentMatrix[t_node1].erase(std::make_pair(t_node2, t_weight));
        m_adjacentMatrix[t_node2].erase(std::make_pair(t_node1, t_weight));

        //* Update link size
        --m_linkSize;
    }

    void deleteLink(const Node& t_node1, const Node& t_node2){
        //* Find weight of input link
        const double weight = this->weight(t_node1, t_node2);
        deleteLink(t_node1, t_node2, weight);
    }

    //* Rewire single link of one node with preserving weight
    void rewire(const Node& t_node, const Node& t_oldNeighbor, const Node& t_newNeighbor){
        //* Delete link to old neighbor
        const double weight = this->weight(t_node, t_oldNeighbor);
        deleteLink(t_node, t_oldNeighbor, weight);

        //* Add link to new neighbor
        addLink(t_node, t_newNeighbor, weight);
    }

    //* Reset Network
    void clear(){
        m_adjacentMatrix.clear();
        m_linkSize = 0;
    }


    //* Show full information of network
    void printAdjacent(const std::string& t_outFileName = "default"){
        if (t_outFileName=="default"){
            std::cout << "Total number of nodes: " << m_size << ", links: "<< m_linkSize << "\n";
            for (Node node=0; node<m_size; ++node){
                std::cout << node << ":";
                for (WNode neighbor : m_adjacentMatrix[node]){
                    std::cout << "(" << neighbor.first << "," << neighbor.second << "), ";
                }
                std::cout << "\n";
            }
        }
        else{
            std::ofstream outFile(t_outFileName);
            outFile << "Total number of nodes: " << m_size << ", links: "<< m_linkSize << "\n";
            for (Node node=0; node<m_size; ++node){
                outFile << node << ":";
                for (WNode neighbor : m_adjacentMatrix[node]){
                    outFile << "(" << neighbor.first << "," << neighbor.second << "), ";
                }
                outFile << "\n";
            }
        }
    }

    void printDegree(const std::string& t_outFileName = "default"){
        std::map<int, int> degreeDistribution;
        for (int node=0; node<m_size; ++node){
            ++degreeDistribution[m_adjacentMatrix[node].size()];
        }
        if (t_outFileName == "default"){
            for (auto it=degreeDistribution.begin(); it!=degreeDistribution.end(); ++it){
                std::cout << it->first << "," << it->second << "\n";
            }
        }
        else{
            std::ofstream outFile(t_outFileName);
            for (auto it=degreeDistribution.begin(); it!=degreeDistribution.end(); ++it){
                outFile << it->first << "," << it->second << "\n";
            }
        }
    }
};//* End of struct WNetwork