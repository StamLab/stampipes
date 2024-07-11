#!/usr/bin/env python3


from collections import namedtuple
import sys
import copy


def process_line(line):
    cols = line.strip().split("\t")
    contains_duplicates = False
    col = ""
    for col in cols:
        if "XA:Z" in col:
            contains_duplicates = True
            break
    sys.stdout.write(line)
    if contains_duplicates:
        aligns = col.split(":")[2].split(";")
        i = 1
        for align in aligns:
            if align:
                i += 1
                new_line = copy.copy(cols)
                fields = align.split(",")
                # for field in fields:
                # print("field", field)
                new_line[0] = cols[0] + "_" + str(i)  # read name
                new_line[1] = int(new_line[1]) | 2 & (~512)  # flag
                new_line[2] = fields[0]  # chrom
                new_line[3] = abs(int(fields[1]))  # pos
                new_line[4] = "30"  # quality
                new_line[5] = fields[2]  # cigar
                new_line[6] = cols[6]  # pair chrom
                new_line[7] = "???"  # pair pos
                new_line[8] = cols[8]  # template length

                # Reverse strand
                if "-" in fields[1]:
                    new_line[1] = new_line[1] | 16 & (~32)
                else:
                    new_line[1] = new_line[1] | 32 & (~16)

                new_line = [str(n) for n in new_line]
                sys.stdout.write("\t".join(new_line) + "\n")


def get_secondary_aligns(read):
    contains_duplicates = False
    cols = read.split("\t")
    for col in cols:
        if "XA:Z" in col:
            contains_duplicates = True
            break
    if not contains_duplicates:
        return []
    aligns = col.split(":")[2].split(";")[:-1]
    # print(aligns)
    return sorted(
        [Alignment(*a.split(",")[:-1]) for a in aligns],
    )


Alignment = namedtuple("Alignment", "chr pos cigar")


# Returns an array of tuples
# Each tuple
def get_aligns(read):
    cols = read.split("\t")
    return [Alignment(cols[2], cols[3], cols[5])] + get_secondary_aligns(read)


def process_pair(r1, r2):
    s1 = get_secondary_aligns(r1)
    s2 = get_secondary_aligns(r2)

    a1 = get_aligns(r1)
    a2 = get_aligns(r2)

    if len(a1) == 1 and len(a2) == 1:
        return

    def valid_pair(i, j):
        # Make sure reads are on same chromosome
        if i.chr != j.chr:
            return False

        p1 = int(i.pos)
        p2 = int(j.pos)

        # Make sure reads are close together
        if abs(abs(p1) - abs(p2)) > 750:
            return False

        # Make sure reads are stranded correctly
        if p1 < 0 and p2 < 0:
            return False
        if p1 > 0 and p2 > 0:
            return False
        return True

    pairs = [(i, j) for i in a1 for j in a2 if valid_pair(i, j)]

    i = 0
    for p in pairs:
        read = p[0]
        mate = p[1]
        i += 1
        c1 = r1.split("\t")
        c2 = r2.split("\t")

        c1[0] += "_" + str(i)
        c1[1] = int(c1[1]) | 2 & (~512)
        c1[2] = read.chr
        c1[3] = abs(int(read.pos))
        c1[4] = 30
        c1[5] = read.cigar
        c1[6] = mate.chr
        c1[7] = abs(int(mate.pos))

        c1[8] = 76  # hack

        c2[0] += "_" + str(i)
        c2[1] = int(c2[1]) | 2 & (~512)
        c2[2] = mate.chr
        c2[3] = abs(int(mate.pos))
        c2[4] = 30
        c2[5] = mate.cigar
        c2[6] = read.chr
        c2[7] = abs(int(read.pos))

        c2[8] = 76  # hack

        # Correct strandedness
        if int(read.pos) > 0:
            c1[1] = c1[1] | 16 & (~32)
            c2[1] = c2[1] | 32 & (~16)
        else:
            c1[1] = c1[1] | 32 & (~16)
            c2[1] = c2[1] | 16 & (~32)

        # Convert to strings
        c1 = [str(c) for c in c1]
        c2 = [str(c) for c in c2]

        sys.stdout.write("\t".join(c1))
        sys.stdout.write("\t".join(c2))

    return
    # c2[0] += "_" + str(i)

    #         new_line[0] = cols[0] + "_" + str(i)        # read name
    #         new_line[1] = int(new_line[1]) | 2 & (~512) # flag
    #         new_line[2] = fields[0]                     # chrom
    #         new_line[3] = abs(int(fields[1]))           # pos
    #         new_line[4] = "30"                          # quality
    #         new_line[5] = fields[2]                     # cigar
    #         new_line[6] = cols[6]                       # pair chrom
    #         new_line[7] = '???'                         # pair pos
    #         new_line[8] = cols[8]                       # template length

    if not s1 and not s2:
        sys.stdout.write(r1)
        sys.stdout.write(r2)
        return

    if s1 and s2:
        # print("r1", r1)
        # print("r2", r2)
        # print("s1", s1)
        # print("s2", s2)
        # combinations = []
        for i in s1:
            for j in s2:
                pass

        return

    # print("r1", r1)
    # print("r2", r2)


def main(args=[]):
    while True:
        try:
            r1 = next(sys.stdin)
            if r1.startswith("@"):
                sys.stdout.write(r1)
                continue
            r2 = next(sys.stdin)
            process_pair(r1, r2)
        except StopIteration:
            break


if __name__ == "__main__":
    main(sys.argv)
