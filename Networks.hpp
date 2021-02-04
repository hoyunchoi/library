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
#include "linearAlgebra.hpp"

using Size = unsigned;
using Adjacency = std::vector<std::set<Size>>;
using WAdjacency = std::vector<std::map<Size, double>>;

struct Node
{
    //* Member variables
    Size m_index;
    std::set<Size> m_neighbors;

    //* Constructor
    Node(){}
    Node(const Size& t_index) : m_index(t_index){}

    //* Ordering operator for set
    bool operator< (const Node& t_node) const {
        return this->m_index < t_node.m_index;
    }

};

struct Network
{
    //* Member variables
    std::string m_type;
    Size m_size{0};
    Size m_linkSize{0};
    double m_meanDegree{0.0};
    Adjacency m_adjacency;

    //* Constructor
    Network(){}
    Network(const Size& t_size) : m_size(t_size){
        m_adjacency.assign(t_size, std::set<Size>());
    }

    //* Reset the whole network withiout size information
    void clear(){
        //* Clear linksize information
        m_linkSize = 0;

        //* Clear adjacency information
        for (auto& neighbors : m_adjacency){
            neighbors.clear();
        }
    }

    //* Whether link is already exists
    bool linkExists(const Size& t_index1, const Size& t_index2) const {
        if (m_adjacency[t_index1].size() <= m_adjacency[t_index2].size()){
            const auto candidates = m_adjacency[t_index1];
            auto it = candidates.find(t_index2);
            return it != candidates.end() ? true : false;
        }
        else{
            const auto candidates = m_adjacency[t_index2];
            auto it = candidates.find(t_index1);
            return it != candidates.end() ? true : false;
        }
    }

    //* Add single link
    void addLink(const Size& t_index1, const Size& t_index2){
        //* Update linksize information
        if (!linkExists(t_index1, t_index2)){
            ++m_linkSize;
        }

        //* Update adjacency information
        m_adjacency[t_index1].emplace(t_index2);
        m_adjacency[t_index2].emplace(t_index1);
    }

    //* Delete single link
    void deleteLink(const Size& t_index1, const Size& t_index2){
        //* Update linksize information
        if (linkExists(t_index1, t_index2)){
            --m_linkSize;
        }

        //* Update adjacency information
        m_adjacency[t_index1].erase(t_index2);
        m_adjacency[t_index2].erase(t_index1);
    }

    //* Print total information of network
    void print(){
        std::cout << "Network Type: " << m_type << "\n";
        std::cout << "Network Size: " << m_size << "\n";
        std::cout << "Total Link Size: " << m_linkSize << "\n";
        std::cout << "Adjacency Matrix\n";
        for (Size i=0; i<m_size; ++i){
            std::cout << i << ": ";
            for (const Size& neighbor : m_adjacency[i]){
                std::cout << neighbor << ", ";
            }
            std::cout << "\n";
        }
    }
};//* End of struct Network


//* Erdos-Renyi network
namespace ER{
    Network generate(const Size& t_size, const double& t_probability, pcg32& t_randomEngine){
        Network ER(t_size);
        std::uniform_real_distribution<double> probabilityDistribution(0, 1);
        for (Size index1 = 0; index1 < t_size; ++index1){
            for (Size index2 = index1+1; index2 < t_size; ++index2){
                if (probabilityDistribution(t_randomEngine)){
                    ER.addLink(index1, index2);
                }
            }
        }
        ER.m_type = "ER";
        ER.m_meanDegree = 2.0*ER.m_linkSize/t_size;
        return ER;
    }

    Network generate(const Size& t_size, const Size& t_linkSize, pcg32& t_randomEngine){
        Network ER(t_size);
        std::uniform_int_distribution<int> indexDistribution(0, t_size - 1);
        while (ER.m_linkSize < t_linkSize){
            Size index1, index2;
            do{
                index1 = indexDistribution(t_randomEngine);
                index2 = indexDistribution(t_randomEngine);
            } while(index1 == index2);
            ER.addLink(index1, index2);
        }
        ER.m_type = "ER";
        ER.m_meanDegree = 2.0*t_linkSize/t_size;
        return ER;
    }
}//* End of namespace ER


//* Random Regular network
namespace RR{
    Network generate(const Size& t_size, const Size& t_degree, pcg32& t_randomEngine){
        Network RR(t_size);
        std::deque<Size> stubs;
        for (Size index=0; index<t_size; ++index){
            for (Size degree=0; degree<t_degree; ++degree){
                stubs.emplace_back(index);
            }
        }
        std::shuffle(stubs.begin(), stubs.end(), t_randomEngine);
        while (!stubs.empty()){
            Size index1 = stubs.front();
            stubs.pop_front();
            Size index2 = stubs.front();
            stubs.pop_front();
            RR.addLink(index1, index2);
        }
        RR.m_type = "RR";
        RR.m_meanDegree = (double)t_degree;
        return RR;
    }
}//* End of namespace RR



//* Scale free network
namespace SF{
    Size randomPowerLawDistribution(const int& t_lower, const int& t_upper, const double& t_exponent, const double& t_prob){
        std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
        return std::pow((std::pow(t_upper+0.5, t_exponent+1) - std::pow(t_lower-0.5, t_exponent+1)) * t_prob + std::pow(t_lower-0.5, t_exponent+1), 1.0 / (t_exponent+1)) + 0.5;
    }

    Network generate(const Size& t_size, const Size& t_linkSize, const double& t_degreeExponent, pcg32& t_randomEngine){
        Network SF(t_size);
        const double weightExponent = 1.0/(t_degreeExponent-1.0);
        std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
        while(SF.m_linkSize < t_linkSize){
            Size index1, index2;
            do{
                index1 = randomPowerLawDistribution(1, t_size, -weightExponent, probabilityDistribution(t_randomEngine))-1;
                index2 = randomPowerLawDistribution(1, t_size, -weightExponent, probabilityDistribution(t_randomEngine))-1;
            } while(index1 == index2);
            SF.addLink(index1, index2);
        }
        SF.m_type = "SF";
        SF.m_meanDegree = 2.0*t_linkSize/t_size;
        return SF;
    }
}//* End of namespace SF

//* Scale free network by Chung-Lu
namespace CL{
    Size weightSampling (const std::vector<double>& t_weight, const double& t_prob){
        return (Size) (std::lower_bound(t_weight.begin(), t_weight.end(), t_prob*t_weight.back())-t_weight.begin());
    }

    Network generate(const Size& t_size, const Size& t_linkSize, const double& t_degreeExponent, pcg32& t_randomEngine){
        Network CL(t_size);
        const double weightExponent = 1.0/(t_degreeExponent-1.0);
        const double correction = weightExponent < 0.5 ? 1.0 : std::pow(10.0*std::sqrt(2.0)*(1.0-weightExponent), 1.0/weightExponent) * std::pow(t_size, 1.0-1.0/(2.0-weightExponent));

        std::vector<double> weight(t_size, 0.0);
        weight[0] = std::pow(correction, -1.0*weightExponent);
        for (Size i=1; i<t_size; ++i){
            weight[i] = weight[i-1]+std::pow(i+correction, -1.0*weightExponent);
        }

        std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
        while (CL.m_linkSize < t_linkSize){
            Size index1, index2;
            do{
                index1 = weightSampling(weight, probabilityDistribution(t_randomEngine));
                index2 = weightSampling(weight, probabilityDistribution(t_randomEngine));
            } while(index1 == index2);
            CL.addLink(index1, index2);
        }
        CL.m_type = "CL";
        CL.m_meanDegree = 2.0*t_linkSize/t_size;
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
    std::vector<int> m_parent;

    //* m_sortedCluster[size] : number of cluster of 'size'
    std::map<int, int> m_sortedCluster;

    //* m_birth[root] : birth time of each 'root'
    //* changedAge[root] : {age, size} of cluster with 'root'
    std::vector<int> m_birth;
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
    int getClusterSize(const int& t_root) const {return -m_parent[t_root];}
    int getDeltaMaximumClusterSize() const {return m_deltaMaximumClusterSize;}
    std::vector<std::pair<int,int>> getChangedAge() const {return m_changedAge;}


    //* get the root of input node
    int getRoot(const int &t_node){
        //* t_node is root
        if (m_parent[t_node] < 0){
            return t_node;
        }
        //* recursively aaaaa node
        return m_parent[t_node] = getRoot(m_parent[t_node]);
    }

    //* Merge two clusters
    void merge(const int& t_root1, const int& t_root2){
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

        //! aaaaa maximum cluster
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

//* Weighted node
struct WNode
{
    //* Member variables
    Size m_index;
    std::map<Size, double> m_neighbors;

    //* Constructor
    WNode(){}
    WNode(const Size& t_index) : m_index(t_index){}

    //* Ordering operator for set
    bool operator< (const WNode& t_wnode) const {
        return this->m_index < t_wnode.m_index;
    }
};

//* Weighted network
struct WNetwork
{
    //* Member variables
    std::string m_type;
    Size m_size{0};
    Size m_linkSize{0};
    double m_meanDegree{0.0};
    WAdjacency m_wadjacency;

    //* Constructor
    WNetwork(){}
    WNetwork(const Size& t_size) : m_size(t_size){
        m_wadjacency.assign(t_size, std::map<Size, double>());
        for (Size index=0; index<t_size; ++index){
            m_wadjacency[index][index] = 0.0;
        }
    }

    //* Reset the whole network
    void clear(){
        //* Clear link size information
        m_linkSize = 0;

        //* Clear adjacency information
        for (auto& wneighbors : m_wadjacency){
            wneighbors.clear();
        }
    }

    //* Whether link is already exists
    double linkExists(const Size& t_index1, const Size& t_index2) const {
        if (m_wadjacency[t_index1].size() <= m_wadjacency[t_index2].size()){
            const auto candidates = m_wadjacency[t_index1];
            auto it = candidates.find(t_index2);
            return it != candidates.end() ? it->second : 0.0;
        }
        else{
            const auto candidates = m_wadjacency[t_index2];
            auto it = candidates.find(t_index1);
            return it != candidates.end() ? it->second : 0.0;
        }
    }

    //* Add single link
    void addLink(const Size& t_index1, const Size& t_index2, const double& t_weight = 1.0){
        //* Update linksize information
        if (!linkExists(t_index1, t_index2)){
            ++m_linkSize;
        }

        //* Update weighted adjacency information
        m_wadjacency[t_index1].emplace(t_index2, t_weight);
        m_wadjacency[t_index2].emplace(t_index1, t_weight);
    }

    //* Delete single link
    void deleteLink(const Size& t_index1, const Size& t_index2){
        //* Get weight information
        const double weight = linkExists(t_index1, t_index2);

        //* Update linksize information
        if (weight){
            --m_linkSize;
        }

        //* Update weighted adjacency information
        m_wadjacency[t_index1].erase(t_index2);
        m_wadjacency[t_index2].erase(t_index1);
    }

    //* Print Adjacency matrix with weight
    void printWAdjacency(const std::string& t_fileName = "", const std::string& t_seperator = ",", const std::string& t_secondSeperator = "\n", const bool t_append = false) const {
        //* File name is given: print into the file
        if (t_fileName.size()){
            std::ofstream file;
            t_append ? file.open(t_fileName, std::ios_base::app) : file.open(t_fileName);
            for (Size index=0; index<m_size; ++index){
                for (Size neighbor=0; neighbor<m_size; ++neighbor){
                    //* neighbor is connected with index
                    if (m_wadjacency[index].find(neighbor) != m_wadjacency[index].end()){
                        file << m_wadjacency[index].at(neighbor) << t_seperator;
                    }
                    //* neighbor is not connected with index
                    else{
                        file << 0.0 << t_seperator;
                    }
                }
                if (t_secondSeperator == "\n"){
                    file << t_secondSeperator;
                }
            }
            file << "\n";
            file.close();
        }
        //* File name is not given: print to terminal
        else{
            for (Size index=0; index<m_size; ++index){
                for (Size neighbor=0; neighbor<m_size; ++neighbor){
                    //* neighbor is connected with index
                    if (m_wadjacency[index].find(neighbor) != m_wadjacency[index].end()){
                        std::cout << m_wadjacency[index].at(neighbor) << t_seperator;
                    }
                    //* neighbor is not connected with index
                    else{
                        std::cout << 0.0 << t_seperator;
                    }
                }
                if (t_secondSeperator == "\n"){
                    std::cout << t_secondSeperator;
                }
            }
        }
    }

    //* Print Adjacency matrix with bool (except node weight)
    void printAdjacency(const std::string& t_fileName = "", const std::string& t_seperator = ",", const std::string& t_secondSeperator = "\n", const bool t_append = false) const {
        //* File name is given: print into the file
        if (t_fileName.size()){
            std::ofstream file;
            t_append ? file.open(t_fileName, std::ios_base::app) : file.open(t_fileName);
            for (Size index=0; index<m_size; ++index){
                for (Size neighbor=0; neighbor<m_size; ++neighbor){
                    //* neighbor is connected with index
                    if (neighbor != index && m_wadjacency[index].find(neighbor) != m_wadjacency[index].end()){
                        file << 1 << t_seperator;
                    }
                    else{
                        file << 0 << t_seperator;
                    }
                }
                if (t_secondSeperator == "\n"){
                    file << t_secondSeperator;
                }
            }
            file << "\n";
            file.close();
        }
        //* File name is not given: print to terminal
        else{
            for (Size index=0; index<m_size; ++index){
                for (Size neighbor=0; neighbor<m_size; ++neighbor){
                    //* neighbor is connected with index
                    if (neighbor == index){
                        std::cout << m_wadjacency[index].at(index) << t_seperator;
                    }
                    else if (m_wadjacency[index].find(neighbor) != m_wadjacency[index].end()){
                        std::cout << 1 << t_seperator;
                    }
                    else{
                        std::cout << 0 << t_seperator;
                    }
                }
                if (t_secondSeperator == "\n"){
                    std::cout << t_secondSeperator;
                }
            }
        }
    }

    //* Prind degree distribution
    void printDegreeDist (const std::string& t_fileName = "", const std::string& t_seperator = ",", const std::string& t_secondSeperator = "\n") const {
        using namespace linearAlgebra;

        //* Generate degree distribution
        std::map<int, double> degreeDist;
        for (Size index=0; index<m_size; ++index){
            ++degreeDist[m_wadjacency[index].size()-1];
        }
        degreeDist /= (double)linearAlgebra::accumulate(degreeDist);

        //* File name is given: print into the file
        if (t_fileName.size()){
            std::ofstream file(t_fileName);
            for (auto it=degreeDist.begin(); it!=degreeDist.end(); ++it){
                file << it->first << t_seperator << it->second << t_secondSeperator;
            }
            file.close();
        }
        //* File name is not given: print to terminal
        else{
            for (auto it=degreeDist.begin(); it!=degreeDist.end(); ++it){
                std::cout << it->first << t_seperator << it->second << t_secondSeperator;
            }
        }
    }

    //* Prind node weight distribution
    void printNodeWeightDist (const std::string& t_fileName = "", const std::string& t_seperator = ",", const std::string& t_secondSeperator = "\n") const {
        using namespace linearAlgebra;

        //* Generate node weight distribution
        std::map<double, double> nodeWeightDist;
        for (Size index=0; index<m_size; ++index){
            ++nodeWeightDist[m_wadjacency[index].at(index)];
        }
        nodeWeightDist /= linearAlgebra::accumulate(nodeWeightDist);

        //* File name is given: print into the file
        if (t_fileName.size()){
            std::ofstream file(t_fileName);
            for (auto it=nodeWeightDist.begin(); it!=nodeWeightDist.end(); ++it){
                file << it->first << t_seperator << it->second << t_secondSeperator;
            }
            file.close();
        }
        //* File name is not given: print to terminal
        else{
            for (auto it=nodeWeightDist.begin(); it!=nodeWeightDist.end(); ++it){
                std::cout << it->first << t_seperator << it->second << t_secondSeperator;
            }
        }
    }
};//* End of struct WNetwork

//* Weighted scale free network by Chung-Lu. Degree distribution and node weight distribution follows same power law
namespace WCL{
    Size weightSampling(const std::vector<double>& t_weight, const double& t_prob){
        return (Size) (std::lower_bound(t_weight.begin(), t_weight.end(), t_prob*t_weight.back())-t_weight.begin());
    }

    WNetwork generate(const Size& t_size, const Size& t_meanPopulation, const double& t_meanDegree, const double& t_degreeExponent, const double& t_linkWeightExponent, pcg32& t_randomEngine){
        WNetwork CL(t_size);
        //* Genearting conventional scale free network by Chung-Lu
        const double weightExponent = 1.0/(t_degreeExponent-1.0);
        const double correction = weightExponent < 0.5 ? 1.0 : std::pow(10.0*std::sqrt(2.0)*(1.0-weightExponent), 1.0/weightExponent) * std::pow(t_size, 1.0-1.0/(2.0-weightExponent));

        std::vector<double> weight(t_size, 0.0);
        weight[0] = std::pow(correction, -1.0*weightExponent);
        for (Size i=1; i<t_size; ++i){
            weight[i] = weight[i-1] + std::pow(i+correction, -1.0*weightExponent);
        }

        std::uniform_real_distribution<double> probabilityDistribution(0.0, 1.0);
        while(CL.m_linkSize < t_meanDegree * t_size / 2.0){
            Size index1, index2;
            do{
                index1 = weightSampling(weight, probabilityDistribution(t_randomEngine));
                index2 = weightSampling(weight, probabilityDistribution(t_randomEngine));
            } while(index1 == index2);
            CL.addLink(index1, index2);
        }

        CL.m_type = "CL";
        CL.m_meanDegree = t_meanDegree;

        //* Add weight to node and link
        const unsigned long long population = t_size * t_meanPopulation;
        for (Size index=0; index<t_size; ++index){
            //* Add node itself to weighted adjacency matrix
            CL.m_wadjacency[index][index] = 0.0;

            //* Get weight of node and link
            const Size degree = CL.m_wadjacency[index].size() - 1;
            const double nodeWeight = degree / (2.0 * CL.m_linkSize);
            std::map<Size, double> linkWeight;
            for (auto it = CL.m_wadjacency[index].begin(); it != CL.m_wadjacency[index].end(); ++it){
                const Size neighbor = it->first;
                if (neighbor != index){
                    const int neighborDegree = CL.m_wadjacency[neighbor].size() - 1;
                    linkWeight[neighbor] = std::pow(degree * neighborDegree, t_linkWeightExponent);
                }
            }
            const double totalLinkWeight = linearAlgebra::accumulate(linkWeight);

            //* Add weight to weighted adjacency matrix
            for (auto it = CL.m_wadjacency[index].begin(); it != CL.m_wadjacency[index].end(); ++it){
                if (it->first != index){
                    CL.m_wadjacency[index][it->first] = linkWeight[it->first] / totalLinkWeight;
                }
                else{
                    CL.m_wadjacency[index][it->first] = std::floor(population * nodeWeight + 0.5);
                }
            }
        }
        return CL;
    }
}//* End of namespace WCL