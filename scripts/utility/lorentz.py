#!/usr/bin/env python3
# Calculates Lorentz statistics from a frequency table

import sys


def lorentz(seen):
    nums = {float(k): seen[k] for k in seen.keys()}

    total_height = sum(k * seen[k] for k in seen.keys())
    total_width = sum(seen.values())

    slope = total_height / total_width

    robinhood_pos = 0
    robinhood_idx = 0

    width, height, area = 0, 0, 0
    for val in sorted(nums.keys()):
        count = int(nums[val])
        delta_h = val * count
        w = count

        width += w
        height += delta_h
        if delta_h > slope and robinhood_pos == 0:
            robinhood_pos = (width-1) / total_width
            robinhood_idx = (slope*(width-1) - height) / total_height
        area += w * (height - delta_h / 2.0)

    ideal = total_width * total_height / 2.0
    gini = (ideal - area) / ideal

    return {
        "gini-coefficient": gini,
        "robinhood-index": robinhood_idx,
        "robinhood-position": robinhood_pos,
    }


def main():
    seen = {}
    # Read from stdin
    for line in sys.stdin:
        (count, times) = line.split()
        seen[int(count)] = int(times)

    # Print results
    for (key, value) in sorted(lorentz(seen).items()):
        print("%s\t%.4f" % (key, value))

if __name__ == "__main__":
    main()
