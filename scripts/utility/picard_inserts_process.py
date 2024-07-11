import numpy
import sys

from copy import copy


# returns the convolutional or statistical autocorrelation of the passed histogram
# unused/untested currently
def get_autocorr(hist, mode="conv"):
    if mode not in ["conv", "stat"]:
        raise ValueError("'mode' argument must be 'conv'olutional or 'stat'istical")
    x = [0] * (1 + max(map(int, hist.keys())))
    for k, v in hist.iteritems():
        x[int(k)] = v

    if mode == "conv":
        result = numpy.correlate(x, x, mode="same")

    if mode == "stat":
        result = numpy.array(
            [1] + [numpy.corrcoef(x[:-i], x[i:])[0, 1] for i in range(1, len(x))]
        )

    return result


# accepts a histogram as a dictionary of index-value pairs and returns it as an array
# fills in all unseen values as 0
def get_arr(hist):
    x = [0] * (1 + max(map(int, hist.keys())))
    for k, v in hist.items():
        x[int(k)] = v
    return x


# returns a normalized a histogram
# mode=max uses the maximum value (highest peak) as the normalization factor
# mode=sum uses the sum of all values as the normalization factor
def normalize_hist(hist, mode="max"):
    hist = numpy.array(hist)
    if mode not in ["max", "sum"]:
        raise ValueError("'mode' argument must be 'max' or 'sum'")
    if mode == "max":
        result = hist.astype(float) / float(max(hist))
    if mode == "sum":
        result = hist.astype(float) / float(sum(hist))
    return result


# returns the normalized diff of a copy of a passed histogram
# unused/untested currently
def diff_norm(df, idx, hist_name="agg_hist", mode="max", r=(50, 140)):
    try:
        c = copy(get_arr(df[hist_name].iloc[idx])[r[0] : r[1]])
    except Exception:
        c = copy(df[hist_name].iloc[idx][r[0] : r[1]])
    return numpy.diff(normalize_hist(c, mode=mode))


# accepts a histogram and returns a sorted array of the n highest peaks in a given range
# unused/untested currently
def get_peak_loc(hist, r=(100, 106), npeaks=1):
    peaks = []
    try:
        n_highest = sorted(hist[r[0] : r[1]], reverse=True)[:npeaks]
        for n in n_highest:
            peaks.append(hist.index(n))
    except Exception:
        peaks = numpy.nan
    return peaks


def get_large_small_ratio(hist, r=(60, 120)):
    ratio = numpy.nan
    hist = get_arr(hist)
    try:
        rt = (hist[r[1]]) / (hist[r[0]])
        ratio = rt
    except Exception:
        ratio = numpy.nan
    return ratio


# calculate:
# pa = peak at x
# pv = fft value at pa
# vat = value at idx 11
# adp = mean-adjusted peak at x
# adv = mean-adjusted value at idx 11
def f_transform(hist, position=11):
    pa = numpy.nan
    pv = numpy.nan
    vat = numpy.nan
    _adp = numpy.nan
    adv = numpy.nan
    try:
        norm_hist = normalize_hist(get_arr(hist), mode="max")
        mx = list(norm_hist).index(max(norm_hist))
        norm_hist = norm_hist[mx : mx + 105]
        ft = numpy.fft.rfft(numpy.diff(norm_hist))
        m = numpy.argmax(ft)
        adj = numpy.mean(ft)
        pa = m
        pv = ft[pa]
        vat = ft[11]
        _adp = pv - adj
        adv = vat - adj
    except Exception:
        adv = numpy.nan
    return adv


def read_hist(file):
    d = {}
    with open(file) as f:
        for line in f:
            line = line.rstrip()
            try:
                (key, val) = line.split()
                d[int(key)] = int(val)
            except Exception:
                next
    return d


def main():
    inputfile = sys.argv[1]
    hist = read_hist(inputfile)
    fourier_transform_eleven = f_transform(hist)
    large_small_ratio = get_large_small_ratio(hist)

    print("insert-ls-ratio\t%.4f" % (large_small_ratio))
    print("insert-ft-eleven\t%.4f" % (fourier_transform_eleven.real))


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
