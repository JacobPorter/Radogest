#!/usr/bin/env python
"""
Balance the data by making copies of records until the
taxid distribution is the same.  Also copy the total number of records
until the minimum size of the file is above a threshold.

:Authors:
    Jacob S. Porter <jsporter@virginia.edu>
"""
import argparse
import datetime
import math
import os
import sys
from collections import defaultdict
from multiprocessing import Lock, Pool
from random import shuffle

from SeqIterator import SeqReader, SeqWriter

LOCK = Lock()
INEXACT_MAX = 10000000


def fix_files(dir_name, files, dest_dir, dir_type, min_size=0, inexact=False):
    """Fix the files."""
    fasta_file = [f for f in files if f.endswith("fasta")][0]
    taxid_work = fasta_file.split(".")[0]
    print("Working on: {} {} {}".format(dir_name, files, taxid_work),
          file=sys.stderr)
    sys.stderr.flush()
    fasta_reader = SeqReader(os.path.join(dir_name, fasta_file))
    fasta_records = defaultdict(list)
    for record in fasta_reader:
        taxid = record[0].split(":")[1]
        fasta_records[taxid].append(record)
    fasta_reader.close()
    taxid_counts = {
        taxid: len(fasta_records[taxid])
        for taxid in fasta_records
    }
    taxid_max = max([taxid_counts[taxid] for taxid in taxid_counts])
    new_fasta = defaultdict(list)
    my_size = taxid_max
    taxid_exclude = []
    if inexact:
        exclude_amount = 0
        for taxid in taxid_counts:
            if taxid_counts[taxid] > taxid_max * 0.8:
                taxid_exclude.append(taxid)
                exclude_amount += taxid_counts[taxid]
        my_size = (INEXACT_MAX - exclude_amount) / (len(taxid_counts) -
                                                    len(taxid_exclude))
    else:
        for taxid in fasta_records:
            new_fasta[taxid] = fasta_records[taxid]
            if taxid in taxid_exclude:
                continue
            count = len(new_fasta[taxid])
            pos = 0
            while count < my_size:
                pos %= count
                new_fasta[taxid].append(fasta_records[taxid][0])
                pos += 1
                count += 1
    total_count = sum([len(new_fasta[taxid]) for taxid in new_fasta])
    all_fasta = []
    for taxid in new_fasta:
        for record in new_fasta[taxid]:
            all_fasta.append((record, taxid))
    if total_count < min_size:
        dup_amount = math.ceil(min_size / total_count)
        all_fasta *= dup_amount
    shuffle(all_fasta)
    LOCK.acquire()
    try:
        create_dir = os.path.join(dest_dir, taxid_work, dir_type)
        os.makedirs(create_dir)
    except FileExistsError:
        print("Directory for {} already exists.".format(create_dir),
              file=sys.stderr)
        sys.stderr.flush()
    LOCK.release()
    fasta_writer = SeqWriter(
        open(
            os.path.join(dest_dir, taxid_work, dir_type,
                         taxid_work + ".fasta"), "w"))
    with open(
            os.path.join(dest_dir, taxid_work, dir_type,
                         taxid_work + ".taxid"), "w") as taxid_fd:
        for tup in all_fasta:
            fasta_writer.write(tup[0])
            print(tup[1], file=taxid_fd)
        taxid_fd.flush()
    fasta_writer.flush()
    fasta_writer.close()
    ret_value = (int(taxid_work), dir_type, len(all_fasta))
    print(ret_value, file=sys.stdout)
    sys.stdout.flush()
    return ret_value


def traverse_directory(src_dir,
                       dest_dir,
                       inexact=False,
                       min_size=0,
                       pool_size=20):
    """Traverse a directory and fix the files."""
    with Pool(processes=pool_size) as pool:
        res_list = []
        # lock = Lock()
        for dir_name, subdirs, files in os.walk(src_dir):
            # print(dir_name, subdirs, files, file=sys.stderr)
            dir_type = None
            if dir_name.endswith("train"):
                dir_type = "train"
            if dir_name.endswith("test"):
                dir_type = "test"
            if dir_type and pool_size > 1:
                res = pool.apply_async(
                    fix_files,
                    (dir_name, files, dest_dir, dir_type, min_size, inexact))
                res_list.append(res)
            elif dir_type and pool_size == 1:
                fix_files(dir_name, files, dest_dir, dir_type, min_size,
                          inexact)
        if pool_size > 1:
            for res in res_list:
                res.get()


def main():
    """Parse arguments."""
    tic = datetime.datetime.now()
    parser = argparse.ArgumentParser(description=('Copy data.'))
    parser.add_argument("src_directory",
                        type=str,
                        help=("A directory with fasta files and taxid files."))
    parser.add_argument("dest_directory",
                        type=str,
                        help=("The directory to copy fasta files and "
                              "taxid files."))
    parser.add_argument("--min_size",
                        "-m",
                        type=int,
                        help=("The minimum number of records.  "
                              "If there are fewer records, "
                              "then duplicates will be made."),
                        default=256)
    parser.add_argument("--processes",
                        "-p",
                        type=int,
                        help=('The number of processes to use.'),
                        default=20)
    parser.add_argument("--inexact",
                        "-i",
                        action="store_true",
                        help=("Use an inexact balancing approach."),
                        default=False)
    args = parser.parse_args()
    print(args, file=sys.stderr)
    traverse_directory(args.src_directory, args.dest_directory, args.inexact,
                       args.min_size, args.processes)
    toc = datetime.datetime.now()
    print("The process took time: {}".format(toc - tic), file=sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()


if __name__ == "__main__":
    main()
