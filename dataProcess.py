import numpy as np

#* Log binning


def distLogBin(raw_x, raw_y, min_exponent=0.0, max_exponent=10.0, delta_exponent=0.1):
    assert len(raw_x) == len(raw_y), "At dist log binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    bins = np.ceil(np.logspace(min_exponent, max_exponent, base=10.0, num=int((max_exponent - min_exponent) / delta_exponent) + 1))
    values = np.sqrt(bins[1:] * bins[:-1])
    binSize = bins[1:] - bins[:-1]
    digitized = np.digitize(raw_x, bins=bins, right=False)
    binned = np.zeros_like(values)
    for i in range(1, len(bins)):
        temp = np.sum(raw_y[digitized==i])
        if temp.size != 0:
            binned[i - 1] = np.sum(temp) / binSize[i - 1]
    values = values[binned != 0]
    binned = binned[binned != 0]
    return values, binned / np.sum(binned)


def avgLogBin(raw_x, raw_y, min_exponent=0.0, max_exponent=10.0, delta_exponent=0.1):
    assert len(raw_x) == len(raw_y), "At avg log binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    bins = np.ceil(np.logspace(min_exponent, max_exponent, base=10.0, num=int((max_exponent - min_exponent) / delta_exponent) + 1))
    values = np.sqrt(bins[1:] * bins[:-1])
    digitized = np.digitize(raw_x, bins=bins, right=False)
    binned = np.zeros_like(values)
    for i in range(1, len(bins)):
        temp = raw_y[digitized == i]
        if temp.size != 0:
            binned[i - 1] = np.mean(raw_y[digitized == i])
    values = values[binned != 0]
    binned = binned[binned != 0]
    return values, binned

#* Linear binning
def distLinBin(raw_x, raw_y, min_val=0.0, max_val=1.0, delta=5e-4):
    assert len(raw_x) == len(raw_y), "At avg lin binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    bins = np.arange(min_val, max_val + delta, delta)
    values = (bins[1:] + bins[:-1]) / 2.0
    digitized = np.digitize(raw_x, bins=bins, right=False)
    binned = np.zeros_like(values)
    for i in range(1, len(bins)):
        temp = raw_y[digitized == i]
        if temp.size != 0:
            binned[i - 1] = np.sum(temp)
    values = values[binned != 0]
    binned = binned[binned != 0]
    return values, binned / np.sum(binned)

def avgLinBin(raw_x, raw_y, min_val=0.0, max_val=1.0, delta=1e-4):
    assert len(raw_x) == len(raw_y), "At avg lin binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    bins = np.arange(min_val, max_val + delta, delta)
    values = (bins[1:] + bins[:-1]) / 2.0
    digitized = np.digitize(raw_x, bins=bins, right=False)
    binned = np.zeros_like(values)
    for i in range(1, len(bins)):
        temp = raw_y[digitized == i]
        if temp.size != 0:
            binned[i - 1] = np.mean(temp)
    return values, binned


if __name__ == "__main__":
    a= np.array([])
    print(a.is_empty())