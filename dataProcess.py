import numpy as np

#* Log binning
def intLogBin(raw_x, raw_y):
    assert len(raw_x) == len(raw_y), "At int log binning, need same length of array"

    #* Setup values for double log binning
    min = np.logspace(0, 10, num=101, base=10)
    value = np.sqrt(min[1:] * min[:-1])

    #* Bin the data
    binned = {}
    sampled = {}
    for i in range(len(raw_x)):
        for j in range(len(value)):
            if raw_x[i] < min[j]:
                if value[j] in binned:
                    binned[value[j]] += raw_y[i]
                    sampled[value[j]] += 1
                else:
                    binned[value[j]] = raw_y[i]
                    sampled[value[j]] = 1
                break;
    for value in binned:
        binned[value] /= sampled[value]
    return np.array(list(binned.keys())), np.array(list(binned.values()))

def doubleLogBin(raw_x, raw_y):
    assert len(raw_x) == len(raw_y), "At double log binning, need same length of array"

    #* Setup values for double log binning
    min = np.logspace(-10, 0, num=101, base=10)
    value = np.sqrt(min[1:] * min[:-1])

    #* Bin the data
    binned = {}
    sampled = {}
    for i in range(len(raw_x)):
        for j in range(len(value)):
            if raw_x[i] < min[j]:
                if value[j] in binned:
                    binned[value[j]] += raw_y[i]
                    sampled[value[j]] += 1
                else:
                    binned[value[j]] = raw_y[i]
                    sampled[value[j]] = 1
                break;
    for value in binned:
        binned[value] /= sampled[value]
    return np.array(list(binned.keys())), np.array(list(binned.values()))

#* Linear binning
def doubleLinBin(raw_x, raw_y):
    assert len(raw_x) == len(raw_y), "At double lin binning, need same length of array"

    #* Setup values for double linear binning
    min = np.linspace(0.0, 1.0, 1001)
    value = (min[1:] + min[:-1])/2.0

    #* Bin the data
    binned = {}
    sampled = {}
    for i in range(len(raw_x)):
        for j in range(len(value)):
            if raw_x[i] < min[j]:
                if value[j] in binned:
                    binned[value[j]] += raw_y[i]
                    sampled[value[j]] += 1
                else:
                    binned[value[j]] = raw_y[i]
                    sampled[value[j]] = 1
                break;
    for value in binned:
        binned[value] /= sampled[value]
    return np.array(list(binned.keys())), np.array(list(binned.values()))

