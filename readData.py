import numpy as np
import pandas as pd
import csv, glob

rootDirectory = "../data/mBFW/"
observables = ["orderParameter", "orderParameter/", "meanClusterSize", "secondGiant" ,"interEventTime" ,"deltaAcceptance" ,"orderParameterDistribution" ,"clusterSizeDistribution" ,"ageDistribution/before" ,"ageDistribution/during", "interEventTimeDistribution/before", "interEventTimeDistribution/during", "deltaUpperBoundDistribution/before", "deltaUpperBoundDistribution/during" ,"deltaAcceptanceDistribution/before","deltaAcceptanceDistribution/during" ,"interEventTime_DeltaAcceptance" ,"upperBound_DeltaAcceptance" ,"deltaUpperBound_DeltaAcceptance"]
relativePath = ["/average/", "logBin/", "/average/", "/average/", "/logBin/", "/logBin/", "/average/", "/logBin/", "/logBin/", "/logBin/", "/logBin/", "/logBin/", "/logBin/", "/logBin/", "/logBin/", "/logBin/", "/logBin/", "/logBin/", "/logBin/" ]

directory = {}
for observable,path in zip(observables, relativePath):
    directory[observable] = rootDirectory+observable+path

#* CSV Reader
def readCSV(t_filename):
    data = pd.read_csv(t_filename, sep=',', header=None)
    data = data.values.transpose()
    if (len(data) == 1):
        return data[0]
    else:
        return tuple([tuple(row) for row in data])

#* File Name Convections
def defaultFileName(t_networkSize, t_acceptanceThreshold):
    return "N{:.1e},G{:.1f}*".format(t_networkSize, t_acceptanceThreshold)

def filename_time(t_networkSize, t_acceptanceThreshold, t_time):
    return "N{:.1e},G{:.1f}*,T{:.4f}*".format(t_networkSize, t_acceptanceThreshold, t_time)

def filename_orderParameter(t_networkSize, t_acceptanceThreshold, t_orderParameter):
    return "N{:.1e},G{:.1f}*,OP{:.4f}*".format(t_networkSize, t_acceptanceThreshold, t_orderParameter)


#* Read Observables
def read(t_observable, t_networkSize, t_acceptanceThreshold, t_checkPoint=None):
    if t_observable == "orderParameterDistribution":
        fileList = np.sort(glob.glob(directory[t_observable] + filename_time(t_networkSize, t_acceptanceThreshold, t_checkPoint)))
    elif t_observable == "clusterSizeDistribution":
        fileList = np.sort(glob.glob(directory[t_observable] + filename_orderParameter(t_networkSize, t_acceptanceThreshold, t_checkPoint)))
    else:
        fileList = np.sort(glob.glob(directory[t_observable] + defaultFileName(t_networkSize, t_acceptanceThreshold)))
    if fileList.size == 0:
        print("There is no file about " + t_observable)
        return
    return readCSV(fileList[-1])

if __name__=="__main__":
    print("This is a module readData.py")