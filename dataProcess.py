import numpy as np

#* Log binning


def distLogBin(raw_x, raw_y, min_exponent=0.0, max_exponent=10.0, delta_exponent=0.1):
    assert len(raw_x) == len(raw_y), "At dist log binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    bin = np.logspace(min_exponent, max_exponent, num=int((max_exponent - min_exponent) / delta_exponent)+1, base=10.0)
    binned_y, binned_x = np.histogram(raw_x, bins=bin, weights=raw_y, density=True)
    binned_x = np.sqrt(binned_x[1:] * binned_x[:-1])
    return binned_x[binned_y != 0], binned_y[binned_y != 0]


def avgLogBin(raw_x, raw_y, min_exponent=0.0, max_exponent=10.0, delta_exponent=0.1):
    assert len(raw_x) == len(raw_y), "At avg log binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    bin = np.logspace(min_exponent, max_exponent, num=int((max_exponent - min_exponent) / delta_exponent)+1, base=10.0)
    binned_y, binned_x= np.histogram(raw_x, bins=bin, weights=raw_y)
    sampled, _ = np.histogram(raw_x, bins=bin)
    binned_x = np.sqrt(binned_x[1:] * binned_x[:-1])
    return binned_x[sampled!=0], binned_y[sampled!=0]/sampled[sampled!=0]



#* Linear binning
def distLinBin(raw_x, raw_y, min_val=0.0, max_val=1.0, delta=5e-4):
    assert len(raw_x) == len(raw_y), "At avg lin binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    bin = np.arange(min_val, max_val+delta, delta)
    binned_y, binned_x = np.histogram(raw_x, bins=bin, weights=raw_y, density=True)
    binned_x = (binned_x[1:] + binned_x[:-1]) / 2.0
    return binned_x[binned_y != 0], binned_y[binned_y != 0]


def avgLinBin(raw_x, raw_y, min_val=0.0, max_val=1.0, delta=1e-4):
    assert len(raw_x) == len(raw_y), "At avg lin binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    bin = np.arange(min_val, max_val + delta, delta)
    binned_y, binned_x = np.histogram(raw_x, bins=bin, weights=raw_y)
    sampled,_ = np.histogram(raw_x, bins=bin)
    binned_x = (binned_x[1:] + binned_x[:-1]) / 2.0
    return binned_x[sampled!=0], binned_y[sampled!=0]/sampled[sampled!=0]


if __name__ == "__main__":
    a = np.array([123, 2354])
    b = np.array([1, 2])
    bin_a, bin_b = avgLogBin(a,b)
    print(bin_a, bin_b)