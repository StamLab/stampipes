#!/usr/bin/env python3
def main():
    import argparse
    import random
    import shutil

    import pysam

    parser = argparse.ArgumentParser()
    parser.add_argument("infile")
    parser.add_argument("outfile")
    parser.add_argument("paired_reads_count")
    parser.add_argument("paired_reads_count_to_select")
    args = parser.parse_args()

    # bam file to read from
    infile = args.infile
    # bam file to write to
    outfile = args.outfile
    # number of paired reads to select
    paired_reads_count_to_select = int(args.paired_reads_count_to_select)
    # total number of reads in the input file
    paired_reads_count = int(args.paired_reads_count)

    # if total number of reads less than number of reads to select,
    # then copy input file to output file
    if paired_reads_count_to_select >= paired_reads_count:
        shutil.copyfile(infile, outfile)
        return

    sorted_read_indexes = random.sample(
        range(paired_reads_count), paired_reads_count_to_select
    )
    sorted_read_indexes.sort()

    print("Selecting %d read pairs" % len(sorted_read_indexes))

    if paired_reads_count_to_select > 100:
        print("First 100 indexes to be selected", sorted_read_indexes[:100])
    else:
        print("Indexes to be selected", sorted_read_indexes)

    # input pysam file
    in_alignment_file = pysam.AlignmentFile(infile, "rb")
    # output pysam file
    out_alignment_file = pysam.AlignmentFile(outfile, "wb", template=in_alignment_file)
    # current index of the next paired read in the input file
    current_index_in_file = 0
    # current index of random indexes to select
    current_index_in_indexes = 0
    # dict that store unique read name for paired read that are not to be selected
    # the name is put when one read of pair is retrieved from file
    # and when second read of pair is retrieved the name is removed
    already_read_map = {}
    # the same as 'already_read_map' but stores unique names of pairs to be selected
    already_read_one_read_map = {}
    # boolean flag indicating that we have selected all required reads except the last one
    find_last_mate = False
    # iterating through reads in the input file
    for alignment in in_alignment_file:
        # if read name is in 'already_read_map' means that this read is the second part of paired read
        # so we do not increment 'current_index_in_file', just move on
        if alignment.query_name in already_read_map:
            del already_read_map[alignment.query_name]
            continue
        # if read name is in 'already_read_one_map'
        # means that this read is the second part of paired read that we need to select
        # so we do not increment 'current_index_in_file', just write the read and move on
        if alignment.query_name in already_read_one_read_map:
            out_alignment_file.write(alignment)
            del already_read_one_read_map[alignment.query_name]
            # if we are looking for the last read, then stop iterating
            # otherwise just move on
            if find_last_mate and len(already_read_one_read_map) == 0:
                break
            else:
                continue
        if find_last_mate:
            continue
        # if current read is under next index that we need to select
        if current_index_in_file == sorted_read_indexes[current_index_in_indexes]:
            # store query_name key to find second mate in the pair
            already_read_one_read_map[alignment.query_name] = True
            # write to out file
            out_alignment_file.write(alignment)
            current_index_in_indexes += 1
            if current_index_in_indexes == paired_reads_count_to_select:
                find_last_mate = True
        else:
            already_read_map[alignment.query_name] = True

        current_index_in_file += 1

    out_alignment_file.close()
    in_alignment_file.close()


if __name__ == "__main__":
    main()
