import typing
import numpy as np

#* Log binning


def distLogBin(raw_x: np.ndarray,
               raw_y: np.ndarray,
               min_exponent: float = None,
               max_exponent: float = None,
               delta_exponent: float = None) -> typing.Tuple[np.ndarray]:
    """
        Log binning data. Binned data will be normalized
        Args:
            raw_x: x-value of histogram
            raw_y: y-value of histogram
            min_exponent: minimum exponent of resulting data. If not given, becomes minimum exponent of raw_x
            max_exponent: maximum exponent of resulting data. If not given, becomes maximum exponent of raw_x
            delta_exponent: width of resulting histogram. If not given, be the value to make number of bin=100
        Return:
            Binned data with x,y values. Area below the histogram = 1
    """
    assert len(raw_x) == len(raw_y), "At dist log binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"
    if min_exponent is None:
        min_exponent = np.log10(np.min(raw_x))
    if max_exponent is None:
        max_exponent = np.log10(np.max(raw_x))
    if delta_exponent is None:
        delta_exponent = (max_exponent - min_exponent) / 100.0

    bin = np.logspace(min_exponent, max_exponent, num=int((max_exponent - min_exponent) / delta_exponent) + 1, base=10.0)
    binned_y, binned_x = np.histogram(raw_x, bins=bin, weights=raw_y, density=True)
    binned_x = np.sqrt(binned_x[1:] * binned_x[:-1])
    return binned_x[binned_y != 0], binned_y[binned_y != 0]


def avgLogBin(raw_x: np.ndarray,
              raw_y: np.ndarray,
              min_exponent: float = None,
              max_exponent: float = None,
              delta_exponent: float = None) -> typing.Tuple[np.ndarray]:
    """
        Log binning data. Binned data is average of raw_y
        Args:
            raw_x: x-value of histogram
            raw_y: y-value of histogram
            min_exponent: minimum exponent of resulting data. If not given, becomes minimum exponent of raw_x
            max_exponent: maximum exponent of resulting data. If not given, becomes maximum exponent of raw_x
            delta_exponent: width of resulting histogram. If not given, be the value to make number of bin=100
        Return:
            Binned data with x,y values. y values is average of raw_y
    """
    assert len(raw_x) == len(raw_y), "At avg log binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    if min_exponent is None:
        min_exponent = np.log10(np.min(raw_x))
    if max_exponent is None:
        max_exponent = np.log10(np.max(raw_x))
    if delta_exponent is None:
        delta_exponent = (max_exponent - min_exponent) / 100.0

    bin = np.logspace(min_exponent, max_exponent, num=int((max_exponent - min_exponent) / delta_exponent) + 1, base=10.0)
    binned_y, binned_x = np.histogram(raw_x, bins=bin, weights=raw_y)
    sampled, _ = np.histogram(raw_x, bins=bin)
    binned_x = np.sqrt(binned_x[1:] * binned_x[:-1])
    return binned_x[sampled != 0], binned_y[sampled != 0] / sampled[sampled != 0]


#* Linear binning
def distLinBin(raw_x: np.ndarray,
               raw_y: np.ndarray,
               min_val: float = None,
               max_val: float = None,
               delta: float = None,
               keep_zero: bool = False) -> typing.Tuple[np.ndarray]:
    """
        Linear binning data. Binned data will be normalized
        Args:
            raw_x: x-value of histogram
            raw_y: y-value of histogram
            min_val: minimum value of resulting data. If not given, becomes minimum value of raw_x
            max_val: maximum value of resulting data. If not given, becomes maximum value of raw_x
            delta: width of resulting histogram. If not given, be the value to make number of bin=1000
            keep_zero: Whether to keep 0 value bin. Default: False
        Return:
            Binned data with x,y values. Area below the histogram = 1
    """
    assert len(raw_x) == len(raw_y), "At avg lin binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    if min_val is None:
        min_val = np.min(raw_x)
    if max_val is None:
        max_val = np.max(raw_x)
    if delta is None:
        delta = (max_val - min_val) / 1000.0

    bin = np.arange(min_val, max_val + delta, delta)
    binned_y, binned_x = np.histogram(raw_x, bins=bin, weights=raw_y, density=True)
    binned_x = (binned_x[1:] + binned_x[:-1]) / 2.0
    if keep_zero:
        return binned_x, binned_y
    else:
        return binned_x[binned_y != 0], binned_y[binned_y != 0]


def avgLinBin(raw_x: np.ndarray,
              raw_y: np.ndarray,
              min_val: float = None,
              max_val: float = None,
              delta: float = None,
              keep_zero: bool = False) -> typing.Tuple[np.ndarray]:
    """
        Linear binning data. Binned data is average of y-value
        Args:
            raw_x: x-value of histogram
            raw_y: y-value of histogram
            min_val: minimum value of resulting data. If not given, becomes minimum value of raw_x
            max_val: maximum value of resulting data. If not given, becomes maximum value of raw_x
            delta: width of resulting histogram. If not given, be the value to make number of bin=10000
            keep_zero: Whether to keep 0 value bin. Default: False
        Return:
            Binned data with x,y values. y values is average of raw_y
    """
    assert len(raw_x) == len(raw_y), "At avg lin binning, need same length of array"
    assert type(raw_x) is np.ndarray, "At avg lin binning, input need to be numpy array"
    assert type(raw_y) is np.ndarray, "At avg lin binning, input need to be numpy array"

    if min_val is None:
        min_val = np.min(raw_x)
    if max_val is None:
        max_val = np.max(raw_x)
    if delta is None:
        delta = (max_val - min_val) / 10000.0

    bin = np.arange(min_val, max_val + delta, delta)
    binned_y, binned_x = np.histogram(raw_x, bins=bin, weights=raw_y)
    sampled, _ = np.histogram(raw_x, bins=bin)
    binned_x = (binned_x[1:] + binned_x[:-1]) / 2.0

    non_empty_bin = sampled > 0
    binned_y[non_empty_bin] /= sampled[non_empty_bin]
    if keep_zero:
        return binned_x, binned_y
    else:
        return binned_x[non_empty_bin], binned_y[non_empty_bin]


if __name__ == "__main__":
    print("This is module data process")
