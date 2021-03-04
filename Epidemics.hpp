#pragma once

#include <map>
#include <random>  //uniform distribution
#include <set>
#include <string>
#include <vector>

#include "CSV.hpp"
#include "Networks.hpp"
#include "pcg_random.hpp"
#include "stringFormat.hpp"

struct Node_Epidemic : public Node {
    //* Member variables
    std::string m_state;
    double m_transitionRate{0.0};

    //* Generator
    Node_Epidemic() {}
    Node_Epidemic(const unsigned& t_index) : Node(t_index) {}
    Node_Epidemic(const unsigned& t_index, const std::string& t_state) : Node(t_index), m_state(t_state) {}
};

/*
    SIR Model simulation
    S+I -> I+I with rate SI_II
    I -> R with rate I_R
*/
namespace SIR {
//* pre-defined parameters
const std::string rootDirectory = "../data/epidemics/SIR/";
std::map<std::string, int> stateToInt = {{"S", 0}, {"I", 1}, {"R", 2}};
unsigned seedSize;

//* Network parameters
std::string networkType;
unsigned networkSize;
unsigned meanDegree;
std::vector<Node_Epidemic> nodes;

//* SIR rate parameter
double SI_II;
double I_R;

//* Random parameter
int randomEngineSeed;
pcg32 randomEngine;
std::uniform_real_distribution<double> probabilityDistribution(0, 1);

//* File name convention
const std::string fileName() {
    const std::string fileName = "N" + to_stringWithExponent((double)networkSize, 1) + ",M" + std::to_string(meanDegree) + ",SIII" + to_stringWithPrecision(SI_II, 2) + ",IR" + to_stringWithPrecision(I_R, 2);

    return randomEngineSeed == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(randomEngineSeed) + ".txt";
}

//* Set Network
void setNetwork(const Network& t_network) {
    networkType = t_network.m_type;
    networkSize = t_network.m_size;
    meanDegree = t_network.m_meanDegree;

    nodes.clear();
    nodes.reserve(networkSize);
    for (unsigned index = 0; index < networkSize; ++index) {
        Node_Epidemic node(index, "S");
        node.m_neighbors = t_network.m_adjacency[index];
        nodes.emplace_back(node);
    }
}

//* Set rate parameters
void setRate(const double& t_SI_II, const double& t_I_R, const unsigned& t_seedSize = 1) {
    SI_II = t_SI_II;
    I_R = t_I_R;
    seedSize = t_seedSize;
}

//* Simulate SIR model using 4th-order Runge-Kutta method
namespace RK4 {
//* pre-defined parameters
const int recentR_length = 5;
const double err = 1e-7;
const double deltaT = 1e-2;

//* Parameters used for RK4
unsigned iteration;
std::vector<double> recentR;
double currentTime;
double ratioS;
double ratioI;
double ratioR;

//* Initialize
void initialize(const int& t_randomEngineSeed, const pcg32& t_randomEngine) {
    randomEngineSeed = t_randomEngineSeed;
    randomEngine = t_randomEngine;
    iteration = 0;
    currentTime = 0.0;
    recentR.assign(recentR_length, -1);
    ratioS = 1.0 - (double)seedSize / networkSize;
    ratioI = (double)seedSize / networkSize;
    ratioR = 0.0;
}

//* RK-4
const double dotS(const double& t_S, const double& t_I) { return -SI_II * meanDegree * t_S * t_I; }
const double dotI(const double& t_S, const double& t_I) { return SI_II * meanDegree * t_S * t_I - I_R * t_I; }
const double dotR(const double& t_I) { return I_R * t_I; }

//* update one step
void update() {
    const double ratioS1 = dotS(ratioS, ratioI);
    const double ratioI1 = dotI(ratioS, ratioI);
    const double ratioR1 = dotR(ratioI);

    const double ratioS2 = dotS(ratioS + ratioS1 * deltaT / 2, ratioI + ratioI1 * deltaT / 2);
    const double ratioI2 = dotI(ratioS + ratioS1 * deltaT / 2, ratioI + ratioI1 * deltaT / 2);
    const double ratioR2 = dotR(ratioI + ratioI1 * deltaT / 2);

    const double ratioS3 = dotS(ratioS + ratioS2 * deltaT / 2, ratioI + ratioI2 * deltaT / 2);
    const double ratioI3 = dotI(ratioS + ratioS2 * deltaT / 2, ratioI + ratioI2 * deltaT / 2);
    const double ratioR3 = dotR(ratioI + ratioI2 * deltaT / 2);

    const double ratioS4 = dotS(ratioS + ratioS3 * deltaT, ratioI + ratioI3 * deltaT);
    const double ratioI4 = dotI(ratioS + ratioS3 * deltaT, ratioI + ratioI3 * deltaT);
    const double ratioR4 = dotR(ratioI + ratioI3 * deltaT);

    ratioS += (ratioS1 + 2 * ratioS2 + 2 * ratioS3 + ratioS4) * deltaT / 6;
    ratioI += (ratioI1 + 2 * ratioI2 + 2 * ratioI3 + ratioI4) * deltaT / 6;
    ratioR += (ratioR1 + 2 * ratioR2 + 2 * ratioR3 + ratioR4) * deltaT / 6;

    ++iteration;
    currentTime += deltaT;
    recentR[iteration % recentR_length] = ratioR;
}  //* End of function SIR::RK4::update

//* Check if the state reached equilibrium
bool equilibrium() {
    if (recentR[recentR_length - 1] < 0) {
        return false;
    }
    if (ratioI < 1.0 / networkSize || ratioR > 1.0 - 1.0 / networkSize) {
        return true;
    }
    const double average = std::accumulate(recentR.begin(), recentR.end(), 0.0) / recentR_length;
    for (auto& r : recentR) {
        if (fabs(r - average) > err) {
            return false;
        }
    }
    return true;
}  //* End of function SIR::RK4::equilibrium

double run(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/RK4/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Write the result
    writeFile << currentTime << "," << ratioS << "," << ratioI << "," << ratioR << "\n";
    while (!equilibrium()) {
        update();
        writeFile << currentTime << "," << ratioS << "," << ratioI << "," << ratioR << "\n";
    }
    writeFile.close();
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SIR::RK4::run
}  // namespace RK4

//* Simulate SIR model using Gillespi Algorithm
namespace GA {
//* Parameters used for GA
std::set<unsigned> reactingIndex;
double deltaT, currentTime;
int numS;
int numI;
int numR;

//* Update transition rate of single node
void updateTransitionRate(const unsigned& t_index) {
    const int intState = stateToInt[nodes[t_index].m_state];
    switch (intState) {
        //* S process
        case 0: {
            unsigned infectiousNeighbor = 0;
            for (const unsigned& neighbor : nodes[t_index].m_neighbors) {
                if (nodes[neighbor].m_state == "I") {
                    ++infectiousNeighbor;
                }
            }
            nodes[t_index].m_transitionRate = SI_II * infectiousNeighbor;
            break;
        }
        //* I process
        case 1: {
            nodes[t_index].m_transitionRate = I_R;
            break;
        }
        //* R Process
        case 2: {
            nodes[t_index].m_transitionRate = 0.0;
            break;
        }
    }
}

//* Update transition rate of all reacting nodes
void updateTransitionRate() {
    for (const unsigned& index : reactingIndex) {
        updateTransitionRate(index);
    }
}

//* Initiallize
void initialize(const int& t_randomEngineSeed, const pcg32& t_randomEngine) {
    randomEngineSeed = t_randomEngineSeed;
    randomEngine = t_randomEngine;

    currentTime = 0.0;
    numS = networkSize - seedSize;
    numI = seedSize;
    numR = 0;
    reactingIndex.clear();
    for (unsigned index = 0; index < seedSize; ++index) {
        nodes[index].m_state = "I";
        reactingIndex.emplace(index);
        for (const unsigned& neighbor : nodes[index].m_neighbors) {
            reactingIndex.emplace(neighbor);
        }
    }
    updateTransitionRate();
}

//* update one step for every nodes
void syncUpdate() {
    deltaT = 1e-2;

    //* Do reactions according to each transition rate and add I into reacting
    std::set<unsigned> newReactingIndex;
    for (const unsigned& index : reactingIndex) {
        const int intState = stateToInt[nodes[index].m_state];
        const double transitionProb = 1.0 - std::exp(-1.0 * nodes[index].m_transitionRate * deltaT);
        switch (intState) {
            //* S process
            case 0: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "I";
                    --numS;
                    ++numI;
                    newReactingIndex.emplace(index);
                }
                break;
            }
            //* I process
            case 1: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "R";
                    --numI;
                    ++numR;
                } else {
                    newReactingIndex.emplace(index);
                }
                break;
            }
        }
    }

    //* Add neighbor S of I node into reacting nodes
    reactingIndex = newReactingIndex;
    for (const unsigned& index : newReactingIndex) {
        for (const unsigned& neighbor : nodes[index].m_neighbors) {
            if (nodes[neighbor].m_state == "S") {
                reactingIndex.emplace(neighbor);
            }
        }
    }

    //* Update time and transition rate
    currentTime += deltaT;
    updateTransitionRate();
}  //* End of function syncUpdate

//* update one step for single node
void asyncUpdate() {
    //* Calculate total tarnsition rate and delta time
    double totalTransitionRate = 0.0;
    for (const unsigned& index : reactingIndex) {
        totalTransitionRate += nodes[index].m_transitionRate;
    }
    deltaT = std::log(1.0 / probabilityDistribution(randomEngine)) / totalTransitionRate;
    totalTransitionRate *= probabilityDistribution(randomEngine);

    //* Choose target node to be reacted
    std::vector<unsigned> shuffledReactingIndex(reactingIndex.begin(), reactingIndex.end());
    std::shuffle(shuffledReactingIndex.begin(), shuffledReactingIndex.end(), randomEngine);
    double cumulativeTransitionRate = 0.0;
    unsigned target = 0;
    for (unsigned i = 0; i < shuffledReactingIndex.size(); ++i) {
        target = shuffledReactingIndex[i];
        cumulativeTransitionRate += nodes[target].m_transitionRate;
        if (cumulativeTransitionRate > totalTransitionRate) {
            break;
        }
    }

    //* React target node and update transition rate
    const int intState = stateToInt[nodes[target].m_state];
    switch (intState) {
        //* S process
        case 0: {
            nodes[target].m_state = "I";
            --numS;
            ++numI;
            updateTransitionRate(target);
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    reactingIndex.emplace(neighbor);
                    updateTransitionRate(neighbor);
                }
            }
            break;
        }
        //* I process
        case 1: {
            nodes[target].m_state = "R";
            --numI;
            ++numR;
            reactingIndex.erase(target);
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    nodes[neighbor].m_transitionRate -= SI_II;
                    if (!nodes[neighbor].m_transitionRate) {
                        reactingIndex.erase(neighbor);
                    }
                }
            }
            break;
        }
    }

    //* Update time
    currentTime += deltaT;

}  //* End of function SIR::GA::asyncUpdate

double syncRun(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/GA_sync/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    while (numI > 0) {
        syncUpdate();
        writeFile << currentTime << "," << numS / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    }
    writeFile.close();
    const double ratioR = numR / (double)networkSize;
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SIR::GA::syncRun

double asyncRun(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/GA_async/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    while (numI > 0) {
        asyncUpdate();
        writeFile << currentTime << "," << numS / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    }
    writeFile.close();
    const double ratioR = numR / (double)networkSize;
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SIR::GA::asyncRun
}  // namespace GA
}  // namespace SIR

/*
    SWIR Model simulation
    S+I -> W+I with rate SI_WI
    W+I -> I+I with rate WI_II
    I -> R with rate I_R
*/
namespace SWIR {
//* pre-defined parameters
const std::string rootDirectory = "../data/epidemics/SWIR/";
std::map<std::string, int> stateToInt = {{"S", 0}, {"W", 1}, {"I", 2}, {"R", 3}};
unsigned seedSize;

//* Network parameters
std::string networkType;
unsigned networkSize;
unsigned meanDegree;
std::vector<Node_Epidemic> nodes;

//* SWIR rate parameter
double SI_WI, WI_II, I_R;

//* Random parameter
int randomEngineSeed;
pcg32 randomEngine;
std::uniform_real_distribution<double> probabilityDistribution(0, 1);

//* File name convention
const std::string fileName() {
    const std::string fileName = "N" + to_stringWithExponent((double)networkSize, 1) + ",M" + std::to_string(meanDegree) + ",SIWI" + to_stringWithPrecision(SI_WI, 2) + ".WIII" + to_stringWithPrecision(WI_II, 2) + ",IR" + to_stringWithPrecision(I_R, 2);

    return randomEngineSeed == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(randomEngineSeed) + ".txt";
}

//* Set Network
void setNetwork(const Network& t_network) {
    networkType = t_network.m_type;
    networkSize = t_network.m_size;
    meanDegree = t_network.m_meanDegree;

    nodes.clear();
    nodes.reserve(networkSize);
    for (unsigned index = 0; index < networkSize; ++index) {
        Node_Epidemic node(index, "S");
        node.m_neighbors = t_network.m_adjacency[index];
        nodes.emplace_back(node);
    }
}

//* Get model parameters
void setRate(const double& t_SI_WI, const double& t_WI_II, const double& t_I_R, const unsigned& t_seedSize = 1) {
    SI_WI = t_SI_WI;
    WI_II = t_WI_II;
    I_R = t_I_R;
    seedSize = t_seedSize;
}

//* Simulate SWIR model using 4th-order Runge-Kutta method
namespace RK4 {
//* pre-defined parameters
const int recentR_length = 5;
const double err = 1e-7;
const double deltaT = 1e-2;

//* Parameters used for RK4
unsigned iteration;
std::vector<double> recentR;
double currentTime;
double ratioS, ratioW, ratioI, ratioR;

//* Initialize
void initialize(const int& t_randomEngineSeed, const pcg32& t_randomEngine) {
    randomEngineSeed = t_randomEngineSeed;
    randomEngine = t_randomEngineSeed;
    iteration = 0;
    currentTime = 0.0;
    recentR.assign(recentR_length, -1);
    ratioS = 1.0 - (double)seedSize / networkSize;
    ratioW = 0.0;
    ratioI = (double)seedSize / networkSize;
    ratioR = 0.0;
}

//* RK-4
const double dotS(const double& t_S, const double& t_I) { return -SI_WI * meanDegree * t_S * t_I; }
const double dotW(const double& t_S, const double& t_W, const double& t_I) { return SI_WI * meanDegree * t_S * t_I - WI_II * meanDegree * t_W * t_I; }
const double dotI(const double& t_W, const double& t_I) { return WI_II * meanDegree * t_W * t_I - I_R * t_I; }
const double dotR(const double& t_I) { return I_R * t_I; }

//* Update one step
void update() {
    const double ratioS1 = dotS(ratioS, ratioI);
    const double ratioW1 = dotW(ratioS, ratioW, ratioI);
    const double ratioI1 = dotI(ratioW, ratioI);
    const double ratioR1 = dotR(ratioI);

    const double ratioS2 = dotS(ratioS + deltaT / 2 * ratioS1, ratioI + deltaT / 2 * ratioI1);
    const double ratioW2 = dotW(ratioS + deltaT / 2 * ratioS1, ratioW + deltaT / 2 * ratioW1, ratioI + deltaT / 2 * ratioI1);
    const double ratioI2 = dotI(ratioW + deltaT / 2 * ratioW1, ratioI + deltaT / 2 * ratioI1);
    const double ratioR2 = dotR(ratioI + deltaT / 2 * ratioI1);

    const double ratioS3 = dotS(ratioS + deltaT / 2 * ratioS2, ratioI + deltaT / 2 * ratioI2);
    const double ratioW3 = dotW(ratioS + deltaT / 2 * ratioS2, ratioW + deltaT / 2 * ratioW2, ratioI + deltaT / 2 * ratioI2);
    const double ratioI3 = dotI(ratioW + deltaT / 2 * ratioW2, ratioI + deltaT / 2 * ratioI2);
    const double ratioR3 = dotR(ratioI + deltaT / 2 * ratioI2);

    const double ratioS4 = dotS(ratioS + deltaT * ratioS3, ratioI + deltaT * ratioI3);
    const double ratioW4 = dotW(ratioS + deltaT * ratioS3, ratioW + deltaT * ratioW3, ratioI + deltaT * ratioI3);
    const double ratioI4 = dotI(ratioW + deltaT * ratioW3, ratioI + deltaT * ratioI3);
    const double ratioR4 = dotR(ratioI + deltaT * ratioI3);

    ratioS += deltaT / 6 * (ratioS1 + 2 * ratioS2 + 2 * ratioS3 + ratioS4);
    ratioW += deltaT / 6 * (ratioW1 + 2 * ratioW2 + 2 * ratioW3 + ratioW4);
    ratioI += deltaT / 6 * (ratioI1 + 2 * ratioI2 + 2 * ratioI3 + ratioI4);
    ratioR += deltaT / 6 * (ratioR1 + 2 * ratioR2 + 2 * ratioR3 + ratioR4);

    ++iteration;
    currentTime += deltaT;
    recentR[iteration % recentR_length] = ratioR;
}  //* End of function SWIR::RK4::update

//* Check if the state reached equilibrium
bool equilibrium() {
    if (recentR[recentR_length - 1] < 0) {
        return false;
    }
    // if (ratioI < 1.0/networkSize || ratioR > 1.0-1.0/networkSize){
    //     return true;
    // }
    const double average = std::accumulate(recentR.begin(), recentR.end(), 0.0) / recentR_length;
    for (auto& r : recentR) {
        if (fabs(r - average) > err) {
            return false;
        }
    }
    return true;
}  //* End of function SEIR::RK4::equilibrium

double run(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/RK4/" + fileName();
    std::ofstream writeFile(writeFileName);
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Write the result
    writeFile << currentTime << "," << ratioS << "," << ratioW << "," << ratioI << "," << ratioR << "\n";
    while (!equilibrium()) {
        update();
        writeFile << currentTime << "," << ratioS << "," << ratioW << "," << ratioI << "," << ratioR << "\n";
    }
    writeFile.close();
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SEIR::RK4::run
}  // namespace RK4

//* Simulate SWIR model using Gillespi Algorithm
namespace GA {
//* Parameters used for GA
std::set<unsigned> reactingIndex;
double deltaT, currentTime;
int numS, numW, numI, numR;

//* Update transition rate of single node
void updateTransitionRate(const unsigned& t_index) {
    const int intState = stateToInt[nodes[t_index].m_state];
    switch (intState) {
        //* S process
        case 0: {
            unsigned infectiousNeighbor = 0;
            for (const unsigned& neighbor : nodes[t_index].m_neighbors) {
                if (nodes[neighbor].m_state == "I") {
                    ++infectiousNeighbor;
                }
            }
            nodes[t_index].m_transitionRate = SI_WI * infectiousNeighbor;
            break;
        }
        //* W process
        case 1: {
            unsigned infectiousNeighbor = 0;
            for (const unsigned& neighbor : nodes[t_index].m_neighbors) {
                if (nodes[neighbor].m_state == "I") {
                    ++infectiousNeighbor;
                }
            }
            nodes[t_index].m_transitionRate = WI_II * infectiousNeighbor;
            break;
        }
        //* I process
        case 2: {
            nodes[t_index].m_transitionRate = I_R;
            break;
        }
        //* R process
        case 3: {
            nodes[t_index].m_transitionRate = 0.0;
            break;
        }
    }
}

void updateTransitionRate() {
    for (const unsigned& index : reactingIndex) {
        updateTransitionRate(index);
    }
}

//* Initialize
void initialize(const int& t_randomEngineSeed, const pcg32& t_randomEngine) {
    randomEngineSeed = t_randomEngineSeed;
    randomEngine = t_randomEngine;

    currentTime = 0.0;
    numS = networkSize - seedSize;
    numW = 0;
    numI = seedSize;
    numR = 0;
    reactingIndex.clear();
    for (unsigned index = 0; index < seedSize; ++index) {
        nodes[index].m_state = "I";
        reactingIndex.emplace(index);
        for (const unsigned& neighbor : nodes[index].m_neighbors) {
            reactingIndex.emplace(neighbor);
        }
    }
    updateTransitionRate();
}

//* Update one step for every nodes
void syncUpdate() {
    deltaT = 1e-2;

    //* Do reactions according to each transition rate and add I into reacting nodes
    std::set<unsigned> newReactingIndex;
    for (const unsigned& index : reactingIndex) {
        const int intState = stateToInt[nodes[index].m_state];
        const double transitionProb = 1.0 - std::exp(-1.0 * nodes[index].m_transitionRate * deltaT);
        switch (intState) {
            //* S process
            case 0: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "W";
                    --numS;
                    ++numW;
                }
                break;
            }
            //* W process
            case 1: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "I";
                    --numW;
                    ++numI;
                    newReactingIndex.emplace(index);
                }
                break;
            }
            //* I process
            case 2: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "R";
                    --numI;
                    ++numR;
                } else {
                    newReactingIndex.emplace(index);
                }
                break;
            }
        }
    }

    //* Add neighbor of S,W of I node into reacting nodes
    reactingIndex = newReactingIndex;
    for (const unsigned& index : newReactingIndex) {
        for (const unsigned& neighbor : nodes[index].m_neighbors) {
            if (nodes[neighbor].m_state == "S" || nodes[neighbor].m_state == "W") {
                reactingIndex.emplace(neighbor);
            }
        }
    }

    //* Update time and Transition rate of updated reacting nodes
    currentTime += deltaT;
    updateTransitionRate();

}  //* end of function SWIR::GA::syncUpdate

//* Update one step for single node
void asyncUpdate() {
    //* Calculate total transition rate and Delta time
    double totalTransitionRate = 0.0;
    for (const unsigned& index : reactingIndex) {
        totalTransitionRate += nodes[index].m_transitionRate;
    }
    deltaT = std::log(1.0 / probabilityDistribution(randomEngine)) / totalTransitionRate;
    totalTransitionRate *= probabilityDistribution(randomEngine);

    //* Choose target node to be reacted
    std::vector<unsigned> shuffledReactingIndex(reactingIndex.begin(), reactingIndex.end());
    std::shuffle(shuffledReactingIndex.begin(), shuffledReactingIndex.end(), randomEngine);
    double cumulativeTransitionRate = 0.0;
    unsigned target;
    for (unsigned i = 0; i < shuffledReactingIndex.size(); ++i) {
        target = shuffledReactingIndex[i];
        cumulativeTransitionRate += nodes[target].m_transitionRate;
        if (cumulativeTransitionRate > totalTransitionRate) {
            break;
        }
    }

    //* React target node and update transition rate
    const int intState = stateToInt[nodes[target].m_state];
    switch (intState) {
        //* S process
        case 0: {
            nodes[target].m_state = "W";
            --numS;
            ++numW;
            nodes[target].m_transitionRate *= WI_II / SI_WI;
            break;
        }
        //* W process
        case 1: {
            nodes[target].m_state = "I";
            --numW;
            ++numI;
            nodes[target].m_transitionRate = I_R;
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    nodes[neighbor].m_transitionRate += SI_WI;
                    reactingIndex.emplace(neighbor);
                } else if (nodes[neighbor].m_state == "W") {
                    nodes[neighbor].m_transitionRate += WI_II;
                    reactingIndex.emplace(neighbor);
                }
            }
            break;
        }
        //* I process
        case 2: {
            nodes[target].m_state = "R";
            --numI;
            ++numR;
            nodes[target].m_transitionRate = 0.0;
            reactingIndex.erase(target);
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    nodes[neighbor].m_transitionRate -= SI_WI;
                    if (!nodes[neighbor].m_transitionRate) {
                        reactingIndex.erase(neighbor);
                    }
                } else if (nodes[neighbor].m_state == "W") {
                    nodes[neighbor].m_transitionRate -= WI_II;
                    if (!nodes[neighbor].m_transitionRate) {
                        reactingIndex.erase(neighbor);
                    }
                }
            }

            break;
        }
    }

    //* Update time
    currentTime += deltaT;

}  //* end of function SWIR::GA::asyncUpdate

double syncRun(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/GA_sync/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS / (double)networkSize << "," << numW / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    while (numI > 0) {
        syncUpdate();
        writeFile << currentTime << "," << numS / (double)networkSize << "," << numW / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    }
    writeFile.close();
    const double ratioR = numR / (double)networkSize;
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SWIR::GA::syncRun

double asyncRun(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/GA_async/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS / (double)networkSize << "," << numW / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    while (numI > 0) {
        asyncUpdate();
        writeFile << currentTime << "," << numS / (double)networkSize << "," << numW / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    }
    writeFile.close();
    const double ratioR = numR / (double)networkSize;
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SWIR::GA::asyncRun
}  // namespace GA
}  // namespace SWIR

/*
    SEIR Model simulation
    S+I -> E+I with rate SI_EI
    E -> I with rate E_I
    I -> R with rate I_R
*/
namespace SEIR {
//* pre-defined parameters
std::map<std::string, int> stateToInt = {{"S", 0}, {"E", 1}, {"I", 2}, {"R", 3}};
unsigned seedSize;

//* Network parameters
std::string networkType;
unsigned networkSize;
unsigned meanDegree;
std::vector<Node_Epidemic> nodes;

//* SEIR rate parameter
double SI_EI, E_I, I_R;

//* Random parameter
int randomEngineSeed;
pcg32 randomEngine;
std::uniform_real_distribution<double> probabilityDistribution(0, 1);

//* File name convention
std::string rootDirectory = "../data/epidemics/SEIR/";
const std::string fileName() {
    const std::string fileName = "N" + to_stringWithExponent((double)networkSize, 1) + ",M" + std::to_string(meanDegree) + ",SIEI" + to_stringWithPrecision(SI_EI, 2) + ",EI" + to_stringWithPrecision(E_I, 2) + ",IR" + to_stringWithPrecision(I_R, 2);

    return randomEngineSeed == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(randomEngineSeed) + ".txt";
}

//* Set Network
void setNetwork(const Network& t_network) {
    networkType = t_network.m_type;
    networkSize = t_network.m_size;
    meanDegree = t_network.m_meanDegree;

    nodes.clear();
    nodes.reserve(networkSize);
    for (unsigned index = 0; index < networkSize; ++index) {
        Node_Epidemic node(index, "S");
        node.m_neighbors = t_network.m_adjacency[index];
        nodes.emplace_back(node);
    }
}

//* Set rate parameters
void setRate(const double& t_SI_EI, const double& t_E_I, const double& t_I_R, const unsigned& t_seedSize = 1) {
    SI_EI = t_SI_EI;
    E_I = t_E_I;
    I_R = t_I_R;
    seedSize = t_seedSize;
}

//* Simulate SEIR model using 4th-order Runge-Kutta method
namespace RK4 {
//* pre-defined Parameters
const int recentR_length = 5;
const double err = 1e-7;
const double deltaT = 1e-2;

//* Parameters used for RK4
unsigned iteration;
std::vector<double> recentR;
double currentTime;
double ratioS, ratioE, ratioI, ratioR;

//* Initialize
void initialize(const int& t_randomEngineSeed, const pcg32& t_randomEngine) {
    randomEngineSeed = t_randomEngineSeed;
    randomEngine = t_randomEngineSeed;
    iteration = 0;
    currentTime = 0.0;
    recentR.assign(recentR_length, -1);
    ratioS = 1.0 - (double)seedSize / networkSize;
    ratioE = 0.0;
    ratioI = (double)seedSize / networkSize;
    ratioR = 0.0;
}

//* RK-4
const double dotS(const double& t_S, const double& t_I) { return -SI_EI * meanDegree * t_S * t_I; }
const double dotE(const double& t_S, const double& t_E, const double& t_I) { return SI_EI * meanDegree * t_S * t_I - E_I * t_E; }
const double dotI(const double& t_E, const double& t_I) { return E_I * t_E - I_R * t_I; }
const double dotR(const double& t_I) { return I_R * t_I; }

//* Update one step
void update() {
    const double ratioS1 = dotS(ratioS, ratioI);
    const double ratioE1 = dotE(ratioS, ratioE, ratioI);
    const double ratioI1 = dotI(ratioE, ratioI);
    const double ratioR1 = dotR(ratioI);

    const double ratioS2 = dotS(ratioS + deltaT / 2 * ratioS1, ratioI + deltaT / 2 * ratioI1);
    const double ratioE2 = dotE(ratioS + deltaT / 2 * ratioS1, ratioE + deltaT / 2 * ratioE1, ratioI + deltaT / 2 * ratioI1);
    const double ratioI2 = dotI(ratioE + deltaT / 2 * ratioE1, ratioI + deltaT / 2 * ratioI1);
    const double ratioR2 = dotR(ratioI + deltaT / 2 * ratioI1);

    const double ratioS3 = dotS(ratioS + deltaT / 2 * ratioS2, ratioI + deltaT / 2 * ratioI2);
    const double ratioE3 = dotE(ratioS + deltaT / 2 * ratioS2, ratioE + deltaT / 2 * ratioE2, ratioI + deltaT / 2 * ratioI2);
    const double ratioI3 = dotI(ratioE + deltaT / 2 * ratioE2, ratioI + deltaT / 2 * ratioI2);
    const double ratioR3 = dotR(ratioI + deltaT / 2 * ratioI2);

    const double ratioS4 = dotS(ratioS + deltaT * ratioS3, ratioI + deltaT * ratioI3);
    const double ratioE4 = dotE(ratioS + deltaT * ratioS3, ratioE + deltaT * ratioE3, ratioI + deltaT * ratioI3);
    const double ratioI4 = dotI(ratioE + deltaT * ratioE3, ratioI + deltaT * ratioI3);
    const double ratioR4 = dotR(ratioI + deltaT * ratioI3);

    ratioS += deltaT / 6 * (ratioS1 + 2 * ratioS2 + 2 * ratioS3 + ratioS4);
    ratioE += deltaT / 6 * (ratioE1 + 2 * ratioE2 + 2 * ratioE3 + ratioE4);
    ratioI += deltaT / 6 * (ratioI1 + 2 * ratioI2 + 2 * ratioI3 + ratioI4);
    ratioR += deltaT / 6 * (ratioR1 + 2 * ratioR2 + 2 * ratioR3 + ratioR4);

    ++iteration;
    currentTime += deltaT;
    recentR[iteration % recentR_length] = ratioR;
}  //* End of function SEIR::RK4::update

//* Check if the state reached equilibrium
bool equilibrium() {
    if (recentR[recentR_length - 1] < 0) {
        return false;
    }
    if (ratioE + ratioI < 1.0 / networkSize || ratioR > 1.0 - 1.0 / networkSize) {
        return true;
    }
    const double average = std::accumulate(recentR.begin(), recentR.end(), 0.0) / recentR_length;
    for (auto& r : recentR) {
        if (fabs(r - average) > err) {
            return false;
        }
    }
    return true;
}  //* End of function SEIR::RK4::equilibrium

double run(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/RK4/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << ratioS << "," << ratioE << "," << ratioI << "," << ratioR << "\n";
    while (!equilibrium()) {
        update();
        writeFile << currentTime << "," << ratioS << "," << ratioE << "," << ratioI << "," << ratioR << "\n";
    }
    writeFile.close();
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SEIR::RK4::run
}  // namespace RK4

//* Simulate SEIR model using Gillespi Algorithm
namespace GA {
//* Parameters used for GA
std::set<unsigned> reactingIndex;
double deltaT, currentTime;
int numS, numE, numI, numR;

//* Update transition rate of single node
void updateTransitionRate(const unsigned& t_index) {
    const int intState = stateToInt[nodes[t_index].m_state];
    switch (intState) {
        //* S process
        case 0: {
            unsigned infectiousNeighbor = 0;
            for (const unsigned& neighbor : nodes[t_index].m_neighbors) {
                if (nodes[neighbor].m_state == "I") {
                    ++infectiousNeighbor;
                }
            }
            nodes[t_index].m_transitionRate = SI_EI * infectiousNeighbor;
            break;
        }
        //* E process
        case 1: {
            nodes[t_index].m_transitionRate = E_I;
            break;
        }
        //* I process
        case 2: {
            nodes[t_index].m_transitionRate = I_R;
            break;
        }
        //* R process
        case 3: {
            nodes[t_index].m_transitionRate = 0.0;
            break;
        }
    }
}

//* Update transition rate of all reacting nodes
void updateTransitionRate() {
    for (const unsigned& index : reactingIndex) {
        updateTransitionRate(index);
    }
}

//* Initialize
void initialize(const int& t_randomEngineSeed, const pcg32& t_randomEngine) {
    randomEngineSeed = t_randomEngineSeed;
    randomEngine = t_randomEngine;

    currentTime = 0.0;
    numS = networkSize - seedSize;
    numE = 0;
    numI = seedSize;
    numR = 0;
    reactingIndex.clear();
    for (unsigned index = 0; index < seedSize; ++index) {
        nodes[index].m_state = "I";
        reactingIndex.emplace(index);
        for (const unsigned& neighbor : nodes[index].m_neighbors) {
            reactingIndex.emplace(neighbor);
        }
    }
    updateTransitionRate();
}

//* Update one step for every nodes
void syncUpdate() {
    deltaT = 1e-2;

    //* Do reactions according to each transition rate and add E,I into reacting nodes
    std::set<unsigned> newReactingIndex;
    for (const unsigned& index : reactingIndex) {
        const int intState = stateToInt[nodes[index].m_state];
        const double transitionProb = 1.0 - std::exp(-1.0 * nodes[index].m_transitionRate * deltaT);
        switch (intState) {
            //* S process
            case 0: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "E";
                    --numS;
                    ++numE;
                    newReactingIndex.emplace(index);
                }
                break;
            }
            //* E process
            case 1: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "I";
                    --numE;
                    ++numI;
                }
                newReactingIndex.emplace(index);
                break;
            }
            //* I Process
            case 2: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "R";
                    --numI;
                    ++numR;
                } else {
                    newReactingIndex.emplace(index);
                }
                break;
            }
        }
    }

    //* Add neighbor S of I node into reacting nodes
    reactingIndex = newReactingIndex;
    for (const unsigned& index : newReactingIndex) {
        if (nodes[index].m_state == "I") {
            for (const unsigned& neighbor : nodes[index].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    reactingIndex.emplace(neighbor);
                }
            }
        }
    }

    //* Update time and Transition rate
    currentTime += deltaT;
    updateTransitionRate();

}  //* End of function SEIR::GA::syncUpdate

//* Update one step for single node
void asyncUpdate() {
    //* Calculate total transition rate and Delta time
    double totalTransitionRate = 0.0;
    for (const unsigned& index : reactingIndex) {
        totalTransitionRate += nodes[index].m_transitionRate;
    }
    deltaT = std::log(1.0 / probabilityDistribution(randomEngine)) / totalTransitionRate;
    totalTransitionRate *= probabilityDistribution(randomEngine);

    //* Choose target node to be reacted
    std::vector<unsigned> shuffledReactingIndex(reactingIndex.begin(), reactingIndex.end());
    std::shuffle(shuffledReactingIndex.begin(), shuffledReactingIndex.end(), randomEngine);
    double cumulativeTransitionRate = 0.0;
    unsigned target;
    for (unsigned i = 0; i < shuffledReactingIndex.size(); ++i) {
        target = shuffledReactingIndex[i];
        cumulativeTransitionRate += nodes[target].m_transitionRate;
        if (cumulativeTransitionRate > totalTransitionRate) {
            break;
        }
    }

    //* React target node and update transition rate
    const int intState = stateToInt[nodes[target].m_state];
    switch (intState) {
        //* S process
        case 0: {
            nodes[target].m_state = "E";
            --numS;
            ++numE;
            nodes[target].m_transitionRate = E_I;
            break;
        }
        //* E process
        case 1: {
            nodes[target].m_state = "I";
            --numE;
            ++numI;
            nodes[target].m_transitionRate = I_R;
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    reactingIndex.emplace(neighbor);
                    nodes[neighbor].m_transitionRate += SI_EI;
                }
            }
            break;
        }
        //* I Process
        case 2: {
            nodes[target].m_state = "R";
            --numI;
            ++numR;
            nodes[target].m_transitionRate = 0.0;
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    nodes[neighbor].m_transitionRate -= SI_EI;
                    if (!nodes[neighbor].m_transitionRate) {
                        reactingIndex.erase(neighbor);
                    }
                }
            }
            break;
        }
    }

    //* Update Time
    currentTime += deltaT;

}  //* End of functon SEIR::GA::asyncUpdate

double syncRun(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/GA_sync/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS / (double)networkSize << "," << numE / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    while (numI > 0) {
        syncUpdate();
        writeFile << currentTime << "," << numS / (double)networkSize << "," << numE / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    }
    writeFile.close();
    const double ratioR = numR / (double)networkSize;
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SEIR::GA::syncRun

double asyncRun(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/GA_async/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS / (double)networkSize << "," << numE / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    while (numI > 0) {
        asyncUpdate();
        writeFile << currentTime << "," << numS / (double)networkSize << "," << numE / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    }
    writeFile.close();
    const double ratioR = numR / (double)networkSize;
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function SEIR::GA::asyncRun
}  // namespace GA
}  // namespace SEIR

/*
    GSEIR Model simulation
    S+E -> E+E with rate SE_EE
    S+I -> E+I with rate SI_EI
    E -> I with rate E_I
    I -> R with rate I_R
*/
namespace GSEIR {
//* pre-defined Parameters
const std::string rootDirectory = "../data/epidemics/GSEIR/";
std::map<std::string, int> stateToInt = {{"S", 0}, {"E", 1}, {"I", 2}, {"R", 3}};
unsigned seedSize;

//* Network parameters
std::string networkType;
unsigned networkSize;
unsigned meanDegree;
std::vector<Node_Epidemic> nodes;

//* GSEIR rate parameter
double SE_EE, SI_EI, E_I, I_R;

//* Random parameter
int randomEngineSeed;
pcg32 randomEngine;
std::uniform_real_distribution<double> probabilityDistribution(0, 1);

//* File name convention
const std::string fileName() {
    const std::string fileName = "N" + to_stringWithExponent((double)networkSize, 1) + ",M" + std::to_string(meanDegree) + ",SEEE" + to_stringWithPrecision(SE_EE, 2) + ",SIEI" + to_stringWithPrecision(SI_EI, 2) + ",EI" + to_stringWithPrecision(E_I, 2) + ",IR" + to_stringWithPrecision(I_R, 2);

    return randomEngineSeed == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(randomEngineSeed) + ".txt";
}

//* Set Network
void setNetwork(const Network& t_network) {
    networkType = t_network.m_type;
    networkSize = t_network.m_size;
    meanDegree = t_network.m_meanDegree;

    nodes.clear();
    nodes.reserve(networkSize);
    for (unsigned index = 0; index < networkSize; ++index) {
        Node_Epidemic node(index, "S");
        node.m_neighbors = t_network.m_adjacency[index];
        nodes.emplace_back(node);
    }
}

//* Set rate parameters
void setRate(const double& t_SE_EE, const double& t_SI_EI, const double& t_E_I, const double& t_I_R, const unsigned& t_seedSize = 1) {
    SE_EE = t_SE_EE;
    SI_EI = t_SI_EI;
    E_I = t_E_I;
    I_R = t_I_R;
    seedSize = t_seedSize;
}

//* Simulate GSEIR model using 4th-order Runge-Kutta method
namespace RK4 {
//* pre-defined parameters
const int recentR_length = 5;
const double err = 1e-7;
const double deltaT = 1e-2;

//* Parameters used for RK4
unsigned iteration;
std::vector<double> recentR;
double currentTime;
double ratioS;
double ratioE;
double ratioI;
double ratioR;

//* Initialize
void initialize(const int& t_randomEngineSeed, const pcg32& t_randomEngine) {
    randomEngineSeed = t_randomEngineSeed;
    randomEngine = t_randomEngineSeed;
    iteration = 0;
    currentTime = 0.0;
    recentR.assign(recentR_length, -1);
    ratioS = 1.0 - (double)seedSize / networkSize;
    ratioE = 0.0;
    ratioI = (double)seedSize / networkSize;
    ratioR = 0.0;
}

//* RK-4
const double dotS(const double& t_S, const double& t_E, const double& t_I) { return -SE_EE * meanDegree * t_S * t_E - SI_EI * meanDegree * t_S * t_I; }
const double dotE(const double& t_S, const double& t_E, const double& t_I) { return SE_EE * meanDegree * t_S * t_E + SI_EI * meanDegree * t_S * t_I - E_I * t_E; }
const double dotI(const double& t_E, const double& t_I) { return E_I * t_E - I_R * t_I; }
const double dotR(const double& t_I) { return I_R * t_I; }

//* Update one step
void update() {
    const double ratioS1 = dotS(ratioS, ratioE, ratioI);
    const double ratioE1 = dotE(ratioS, ratioE, ratioI);
    const double ratioI1 = dotI(ratioE, ratioI);
    const double ratioR1 = dotR(ratioI);

    const double ratioS2 = dotS(ratioS + deltaT / 2 * ratioS1, ratioE + deltaT / 2 * ratioE1, ratioI + deltaT / 2 * ratioI1);
    const double ratioE2 = dotE(ratioS + deltaT / 2 * ratioS1, ratioE + deltaT / 2 * ratioE1, ratioI + deltaT / 2 * ratioI1);
    const double ratioI2 = dotI(ratioE + deltaT / 2 * ratioE1, ratioI + deltaT / 2 * ratioI1);
    const double ratioR2 = dotR(ratioI + deltaT / 2 * ratioI1);

    const double ratioS3 = dotS(ratioS + deltaT / 2 * ratioS2, ratioE + deltaT / 2 * ratioE2, ratioI + deltaT / 2 * ratioI2);
    const double ratioE3 = dotE(ratioS + deltaT / 2 * ratioS2, ratioE + deltaT / 2 * ratioE2, ratioI + deltaT / 2 * ratioI2);
    const double ratioI3 = dotI(ratioE + deltaT / 2 * ratioE2, ratioI + deltaT / 2 * ratioI2);
    const double ratioR3 = dotR(ratioI + deltaT / 2 * ratioI2);

    const double ratioS4 = dotS(ratioS + deltaT * ratioS3, ratioE + deltaT * ratioE3, ratioI + deltaT * ratioI3);
    const double ratioE4 = dotE(ratioS + deltaT * ratioS3, ratioE + deltaT * ratioE3, ratioI + deltaT * ratioI3);
    const double ratioI4 = dotI(ratioE + deltaT * ratioE3, ratioI + deltaT * ratioI3);
    const double ratioR4 = dotR(ratioI + deltaT * ratioI3);

    ratioS += deltaT / 6 * (ratioS1 + 2 * ratioS2 + 2 * ratioS3 + ratioS4);
    ratioE += deltaT / 6 * (ratioE1 + 2 * ratioE2 + 2 * ratioE3 + ratioE4);
    ratioI += deltaT / 6 * (ratioI1 + 2 * ratioI2 + 2 * ratioI3 + ratioI4);
    ratioR += deltaT / 6 * (ratioR1 + 2 * ratioR2 + 2 * ratioR3 + ratioR4);

    ++iteration;
    currentTime += deltaT;
    recentR[iteration % recentR_length] = ratioR;
}  //* End of function GSEIR::RK4::update

//* Check if the state reached equilibrium
bool equilibrium() {
    if (recentR[recentR_length - 1] < 0) {
        return false;
    }
    if (ratioE + ratioI < 1.0 / networkSize || ratioR > 1.0 - 1.0 / networkSize) {
        return true;
    }
    const double average = std::accumulate(recentR.begin(), recentR.end(), 0.0) / recentR_length;
    for (auto& r : recentR) {
        if (fabs(r - average) > err) {
            return false;
        }
    }
    return true;
}  //* End of function GSEIR::RK4::equilibrium

double run(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/RK4/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Write the result
    writeFile << currentTime << "," << ratioS << "," << ratioE << "," << ratioI << "," << ratioR << "\n";
    while (!equilibrium()) {
        update();
        writeFile << currentTime << "," << ratioS << "," << ratioE << "," << ratioI << "," << ratioR << "\n";
    }
    writeFile.close();
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function GSEIR::RK4::run
}  // namespace RK4

//* Simulate GSEIR model using Gillespi Algorithm
namespace GA {
//* Parameters used for GA
std::set<unsigned> reactingIndex;
double deltaT, currentTime;
unsigned numS, numE, numI, numR;

//* Update transition rate of single node
void updateTransitionRate(const unsigned& t_index) {
    const int intState = stateToInt[nodes[t_index].m_state];
    switch (intState) {
        //* S process
        case 0: {
            unsigned exposedNeighbor = 0;
            unsigned infectiousNeighbor = 0;
            for (const unsigned& neighbor : nodes[t_index].m_neighbors) {
                if (nodes[neighbor].m_state == "E") {
                    ++exposedNeighbor;
                } else if (nodes[neighbor].m_state == "I") {
                    ++infectiousNeighbor;
                }
            }
            nodes[t_index].m_transitionRate = SE_EE * exposedNeighbor + SI_EI * infectiousNeighbor;
            break;
        }

        //* E process
        case 1: {
            nodes[t_index].m_transitionRate = E_I;
            break;
        }
        //* I process
        case 2: {
            nodes[t_index].m_transitionRate = I_R;
            break;
        }
        //* R process
        case 3: {
            nodes[t_index].m_transitionRate = 0.0;
            break;
        }
    }
}

//* Update transition rate of all reacting nodes
void updateTransitionRate() {
    for (const unsigned& index : reactingIndex) {
        updateTransitionRate(index);
    }
}

//* Initialize
void initialize(const int& t_randomEngineSeed, const pcg32& t_randomEngine) {
    randomEngineSeed = t_randomEngineSeed;
    randomEngine = t_randomEngine;

    currentTime = 0.0;
    numS = networkSize - seedSize;
    numE = 0;
    numI = seedSize;
    numR = 0;
    reactingIndex.clear();
    for (unsigned index = 0; index < seedSize; ++index) {
        nodes[index].m_state = "I";
        reactingIndex.emplace(index);
        for (const unsigned& neighbor : nodes[index].m_neighbors) {
            reactingIndex.emplace(neighbor);
        }
    }
    updateTransitionRate();
}

//* Update one step for every nodes
void syncUpdate() {
    deltaT = 1e-1;

    //* Do reactions according to each transition rate and add E,I into reacting nodes
    std::set<unsigned> newReactingIndex;
    for (const unsigned& index : reactingIndex) {
        const int intState = stateToInt[nodes[index].m_state];
        const double transitionProb = 1.0 - std::exp(-1.0 * nodes[index].m_transitionRate * deltaT);
        switch (intState) {
            //* S process
            case 0:
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "E";
                    --numS;
                    ++numE;
                    newReactingIndex.emplace(index);
                }
                break;
            //* E process
            case 1: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "I";
                    --numE;
                    ++numI;
                }
                newReactingIndex.emplace(index);
                break;
            }
            //* I process
            case 2: {
                if (probabilityDistribution(randomEngine) <= transitionProb) {
                    nodes[index].m_state = "R";
                    --numI;
                    ++numR;
                } else {
                    newReactingIndex.emplace(index);
                }
                break;
            }
        }
    }

    //* Add neighbor of S of E,I node into reacting nodes
    reactingIndex = newReactingIndex;
    for (const unsigned& index : newReactingIndex) {
        for (const unsigned& neighbor : nodes[index].m_neighbors) {
            if (nodes[neighbor].m_state == "S") {
                reactingIndex.emplace(neighbor);
            }
        }
    }

    //* Update time and Transition rate of updated reacting nodes
    currentTime += deltaT;
    updateTransitionRate();

}  //* End of function GSEIR::GA::syncUpdate

//* Update one step for single node
void asyncUpdate() {
    //* Calculate total transition rate and Delta time
    double totalTransitionRate = 0.0;
    for (const unsigned& index : reactingIndex) {
        totalTransitionRate += nodes[index].m_transitionRate;
    }
    deltaT = std::log(1.0 / probabilityDistribution(randomEngine)) / totalTransitionRate;
    totalTransitionRate *= probabilityDistribution(randomEngine);

    //* Choose target node to be reacted
    std::vector<unsigned> shuffledReactingIndex(reactingIndex.begin(), reactingIndex.end());
    std::shuffle(shuffledReactingIndex.begin(), shuffledReactingIndex.end(), randomEngine);
    double cumulativeTransitionRate = 0.0;
    unsigned target;
    for (unsigned i = 0; i < shuffledReactingIndex.size(); ++i) {
        target = shuffledReactingIndex[i];
        cumulativeTransitionRate += nodes[target].m_transitionRate;
        if (cumulativeTransitionRate > totalTransitionRate) {
            break;
        }
    }

    //* React target node and update transition rate
    const int intState = stateToInt[nodes[target].m_state];
    switch (intState) {
        //* S process
        case 0: {
            nodes[target].m_state = "E";
            --numS;
            ++numE;
            nodes[target].m_transitionRate = E_I;
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    nodes[neighbor].m_transitionRate += SE_EE;
                    reactingIndex.emplace(neighbor);
                }
            }
            break;
        }
        //* E process
        case 1: {
            nodes[target].m_state = "I";
            --numE;
            ++numI;
            nodes[target].m_transitionRate = I_R;
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    nodes[neighbor].m_transitionRate += SI_EI - SE_EE;
                    reactingIndex.emplace(neighbor);
                }
            }
            break;
        }
        //* I Process
        case 2: {
            nodes[target].m_state = "R";
            --numI;
            ++numR;
            nodes[target].m_transitionRate = 0.0;
            reactingIndex.erase(target);
            for (const unsigned& neighbor : nodes[target].m_neighbors) {
                if (nodes[neighbor].m_state == "S") {
                    nodes[neighbor].m_transitionRate -= SI_EI;
                    if (!nodes[neighbor].m_transitionRate) {
                        reactingIndex.erase(neighbor);
                    }
                }
            }
            break;
        }
    }

    //* Update Time
    currentTime += deltaT;

}  //* End of function GSEIR::GA::asyncUpdate

double syncRun(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/GA_sync/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS / (double)networkSize << "," << numE / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    while (numI > 0) {
        syncUpdate();
        writeFile << currentTime << "," << numS / (double)networkSize << "," << numE / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    }
    writeFile.close();
    const double ratioR = numR / (double)networkSize;
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function GSEIR::GA::syncRun

double asyncRun(const int& t_randomEngineSeed, const pcg32& t_randomEngine, const bool t_deletion = true) {
    //* Initialize model
    initialize(t_randomEngineSeed, t_randomEngine);

    //* Define write file
    const std::string writeFileName = rootDirectory + networkType + "/GA_async/" + fileName();
    std::ofstream writeFile(writeFileName);

    //* Write the result
    writeFile << currentTime << "," << numS / (double)networkSize << "," << numE / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    while (numI > 0) {
        asyncUpdate();
        writeFile << currentTime << "," << numS / (double)networkSize << "," << numE / (double)networkSize << "," << numI / (double)networkSize << "," << numR / (double)networkSize << "\n";
    }
    writeFile.close();
    const double ratioR = numR / (double)networkSize;
    std::cout << writeFileName << ": " << ratioR << "\n";
    if (t_deletion && ratioR < 0.01) {
        CSV::deleteFile(writeFileName);
    }
    return ratioR;
}  //* End of function GSEIR::GA::asyncRun
}  // namespace GA
}  // namespace GSEIR