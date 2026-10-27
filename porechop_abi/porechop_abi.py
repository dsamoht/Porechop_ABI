"""
Copyright 2017-2022 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop
This module contains the main script for Porechop.
It is executed when a user runs `porechop` after installation or
`porechop-runner.py` directly from the source directory.
This file is part of Porechop. Porechop is free software:
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. Porechop is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
Porechop. If not, see <http://www.gnu.org/licenses/>.
Author: Ryan Wick
Last modified 2022-07-06
"""

import os
import sys
import subprocess
import shutil
import re
from collections import defaultdict
from .misc import load_fasta_or_fastq, print_table, red
from .misc import bold_underline, int_to_str
from .adapters import ADAPTERS, make_full_native_barcode_adapter,\
    make_old_full_rapid_barcode_adapter, make_new_full_rapid_barcode_adapter
from .nanopore_read import NanoporeRead
from .abinitio import launch_ab_initio
from .arg_parser import get_arguments
from .parse_adapter_file import *
from .cpp_function_wrappers import align_regions, find_middle_hits

# Reads per batch of native alignments (bounds the memory held by precomputed results).
NATIVE_CHUNK = 10000


def main():
    args = get_arguments()

    # modifiying global variable ADAPTERS
    global ADAPTERS

    # Clear database if required
    if(args.discard_database):
        ADAPTERS = []

    # Add ab-initio adapter to the database
    if(args.ab_initio):
        if(args.verbosity > 0):
            print('\n' + bold_underline('Ab Initio Phase'),
                  flush=True,
                  file=args.print_dest)
        ADAPTERS += launch_ab_initio(args)
        # If we just want to display inferred adapter, exit here.
        if(args.guess_adapter_only):
            exit()

    # Do we need to add an adapter from a custom adapter file ?
    if(args.custom_adapters):
        ADAPTERS += get_adapters(args.custom_adapters)

    # Letting Porechop do its job.
    reads, check_reads, read_type = load_reads(args.input,
                                               args.verbosity,
                                               args.print_dest,
                                               args.check_reads)
    matching_sets = find_matching_adapter_sets(check_reads,
                                               args.verbosity,
                                               args.end_size,
                                               args.scoring_scheme_vals,
                                               args.print_dest,
                                               args.adapter_threshold,
                                               args.threads)
    matching_sets = fix_up_1d2_sets(matching_sets)

    if args.barcode_dir:
        forward_or_reverse_barcodes = choose_barcoding_kit(matching_sets,
                                                           args.verbosity,
                                                           args.print_dest)
    else:
        forward_or_reverse_barcodes = None

    display_adapter_set_results(matching_sets, args.verbosity, args.print_dest)
    matching_sets = add_full_barcode_adapter_sets(matching_sets)

    if args.verbosity > 0:
        print('\n', file=args.print_dest)

    if matching_sets:
        check_barcodes = (args.barcode_dir is not None)
        find_adapters_at_read_ends(reads, matching_sets, args.verbosity,
                                   args.end_size,
                                   args.extra_end_trim, args.end_threshold,
                                   args.scoring_scheme_vals, args.print_dest,
                                   args.min_trim_size,
                                   args.threads, check_barcodes,
                                   args.barcode_threshold,
                                   args.barcode_diff,
                                   args.require_two_barcodes,
                                   forward_or_reverse_barcodes)
        display_read_end_trimming_summary(reads, args.verbosity,
                                          args.print_dest)

        if not args.no_split:
            find_adapters_in_read_middles(reads, matching_sets, args.verbosity,
                                          args.middle_threshold,
                                          args.extra_middle_trim_good_side,
                                          args.extra_middle_trim_bad_side,
                                          args.scoring_scheme_vals,
                                          args.print_dest,
                                          args.threads,
                                          args.discard_middle)
            display_read_middle_trimming_summary(reads, args.discard_middle,
                                                 args.verbosity,
                                                 args.print_dest)
    elif args.verbosity > 0:
        print('No adapters found - output reads are unchanged from input reads\n',
              file=args.print_dest)

    output_reads(reads, args.format, args.output, read_type, args.verbosity,
                 args.discard_middle, args.min_split_read_size, args.print_dest,
                 args.barcode_dir, args.input, args.untrimmed, args.threads,
                 args.discard_unassigned)


def load_reads(input_file_or_directory, verbosity, print_dest, check_read_count):

    # If the input is a file, just load reads from that file. The check reads will just be the
    # first reads from that file.
    if os.path.isfile(input_file_or_directory):

        if verbosity > 0:
            print('\n' + bold_underline('Loading reads'),
                  flush=True, file=print_dest)
            print(input_file_or_directory, flush=True, file=print_dest)
        reads, read_type = load_fasta_or_fastq(input_file_or_directory)
        if read_type == 'FASTA':
            reads = [NanoporeRead(x[2], x[1], '') for x in reads]
        else:  # FASTQ
            reads = [NanoporeRead(x[4], x[1], x[3]) for x in reads]
        check_reads = reads[:check_read_count]

    # If the input is a directory, assume it's an Albacore directory and search it recursively for
    # fastq files. The check reads will be spread over all of the input files.
    elif os.path.isdir(input_file_or_directory):
        if verbosity > 0:
            print('\n' + bold_underline('Searching for FASTQ files'),
                  flush=True, file=print_dest)
        fastqs = sorted([os.path.join(dir_path, f)
                         for dir_path, _, filenames in os.walk(input_file_or_directory)
                         for f in filenames
                         if f.lower().endswith('.fastq') or f.lower().endswith('.fastq.gz')])
        if not fastqs:
            sys.exit('Error: could not find fastq files in ' +
                     input_file_or_directory)
        reads = []
        read_type = 'FASTQ'
        check_reads = []
        check_reads_per_file = int(round(check_read_count / len(fastqs)))
        for fastq_file in fastqs:
            if verbosity > 0:
                print(fastq_file, flush=True, file=print_dest)
            file_reads, _ = load_fasta_or_fastq(fastq_file)
            file_reads = [NanoporeRead(x[4], x[1], x[3]) for x in file_reads]

            albacore_barcode = get_albacore_barcode_from_path(fastq_file)
            for read in file_reads:
                read.albacore_barcode_call = albacore_barcode
            reads += file_reads
            check_reads += file_reads[:check_reads_per_file]
        if verbosity > 0:
            print('', flush=True, file=print_dest)

    else:
        sys.exit('Error: could not find ' + input_file_or_directory)

    if verbosity > 0:
        print(int_to_str(len(reads)) + ' reads loaded\n\n',
              flush=True, file=print_dest)
    return reads, check_reads, read_type


def get_albacore_barcode_from_path(albacore_path):
    if '/unclassified/' in albacore_path:
        return 'none'
    matches = re.findall('/barcode(\\d\\d)/', albacore_path)
    if matches:
        albacore_barcode_num = matches[-1]
        return 'BC' + albacore_barcode_num
    return None


def find_matching_adapter_sets(check_reads, verbosity, end_size, scoring_scheme_vals, print_dest,
                               adapter_threshold, threads):
    """
    Aligns all of the adapter sets to the start/end of reads to see which (if any) matches best.
    """
    read_count = len(check_reads)
    if verbosity > 0:
        print(bold_underline('Looking for known adapter sets'),
              flush=True, file=print_dest)
        output_progress_line(0, read_count, print_dest)

    search_adapters = [a for a in ADAPTERS if '(full sequence)' not in a.name]
    search_adapter_count = len(search_adapters)

    # All alignments are done natively (in parallel), then applied read by read.
    start_hits = align_regions([read.seq[:end_size] for read in check_reads],
                               [x.start_sequence[1] if x.start_sequence else None
                                for x in search_adapters],
                               scoring_scheme_vals, threads)
    end_hits = align_regions([read.seq[-end_size:] for read in check_reads],
                             [x.end_sequence[1] if x.end_sequence else None
                              for x in search_adapters],
                             scoring_scheme_vals, threads)
    for read_num, read in enumerate(check_reads):
        for i, adapter_set in enumerate(search_adapters):
            read.align_adapter_set(adapter_set, end_size, scoring_scheme_vals,
                                   start_hits[read_num][i], end_hits[read_num][i])
        if verbosity > 0:
            output_progress_line(read_num+1, read_count, print_dest)

    if verbosity > 0:
        output_progress_line(read_count, read_count,
                             print_dest, end_newline=True)

    return [x for x in search_adapters if x.best_start_or_end_score() >= adapter_threshold]


def choose_barcoding_kit(adapter_sets, verbosity, print_dest):
    """
    If the user is sorting reads by barcode bin, choose one barcode configuration (rev comp
    barcodes at the start of the read or at the end of the read) and ignore the other.
    """
    # Tally up scores for forward and reverse barcodes.
    forward_start_or_end, reverse_start_or_end = 0, 0
    forward_start_and_end, reverse_start_and_end = 0, 0
    for adapter_set in adapter_sets:
        if 'barcode' in adapter_set.name.lower():
            if '(forward)' in adapter_set.name.lower():
                forward_start_or_end += adapter_set.best_start_or_end_score()
                forward_start_and_end += adapter_set.best_start_score
                forward_start_and_end += adapter_set.best_end_score
            elif '(reverse)' in adapter_set.name.lower():
                reverse_start_or_end += adapter_set.best_start_or_end_score()
                reverse_start_and_end += adapter_set.best_start_score
                reverse_start_and_end += adapter_set.best_end_score

    if forward_start_or_end == 0 and reverse_start_or_end == 0:
        sys.exit(
            'Error: no barcodes were found, so Porechop cannot perform barcode demultiplexing')

    # If possible, make a decision using each barcode's best start OR end score.
    orientation = None
    if forward_start_or_end > reverse_start_or_end:
        orientation = 'forward'
    elif reverse_start_or_end > forward_start_or_end:
        orientation = 'reverse'

    # If that didn't work (i.e. it's a tie between forward and reverse), then choose based on the
    # sum of both start AND end scores.
    elif forward_start_and_end > reverse_start_and_end:
        orientation = 'forward'
    elif reverse_start_and_end > forward_start_and_end:
        orientation = 'reverse'

    if orientation is None:
        sys.exit('Error: Porechop could not determine barcode orientation')

    if verbosity > 0:
        print('\nBarcodes determined to be in ' +
              orientation + ' orientation', file=print_dest)
    return orientation


def fix_up_1d2_sets(matching_sets):
    """
    The 1D^2 adapters share a lot of common sequence with the old SQK-MAP006_short adapters, so if
    there are strong 1D^2 hits, we can exclude the redundant SQK-MAP006_short adapters.
    """
    matching_set_names = [x.name for x in matching_sets]
    if '1D^2 part 1' in matching_set_names and '1D^2 part 2' in matching_set_names and \
            'SQK-MAP006 Short' in matching_set_names:
        part_1_score = [x for x in matching_sets
                        if x.name == '1D^2 part 1'][0].best_start_or_end_score()
        part_2_score = [x for x in matching_sets
                        if x.name == '1D^2 part 2'][0].best_start_or_end_score()
        sqk_score = [x for x in matching_sets
                     if x.name == 'SQK-MAP006 Short'][0].best_start_or_end_score()
        if part_1_score >= sqk_score and part_2_score >= sqk_score:
            matching_sets = [
                x for x in matching_sets if x.name != 'SQK-MAP006 Short']
    return matching_sets


def display_adapter_set_results(matching_sets, verbosity, print_dest):
    if verbosity < 1:
        return
    table = [['Set', 'Best read start %ID', 'Best read end %ID']]
    row_colours = {}
    matching_set_names = [x.name for x in matching_sets]
    search_adapters = [a for a in ADAPTERS if '(full sequence)' not in a.name]
    for adapter_set in search_adapters:
        start_score = '%.1f' % adapter_set.best_start_score
        end_score = '%.1f' % adapter_set.best_end_score
        table.append([adapter_set.name, start_score, end_score])
        if adapter_set.name in matching_set_names:
            row_colours[len(table) - 1] = 'green'
    if verbosity > 0:
        print_table(table, print_dest, alignments='LRR', row_colour=row_colours,
                    fixed_col_widths=[35, 8, 8])


def add_full_barcode_adapter_sets(matching_sets):
    """
    This function adds some new 'full' adapter sequences based on what was already found. For
    example, if the ligation adapters and the reverse barcode adapters are found, it assumes we are
    looking at a native barcoding run and so it adds the complete native barcoding adapter
    sequences (with the barcode's upstream and downstream context included).
    """
    matching_set_names = [x.name for x in matching_sets]

    for i in range(1, 97):

        # Native barcode full sequences
        if all(x in matching_set_names
               for x in ['SQK-NSK007', 'Barcode ' + str(i) + ' (reverse)']):
            matching_sets.append(make_full_native_barcode_adapter(i))

        # Rapid barcode full sequences
        if all(x in matching_set_names
               for x in ['Rapid', 'Barcode ' + str(i) + ' (forward)']):
            if 'RBK004_upstream' in matching_set_names:
                matching_sets.append(make_new_full_rapid_barcode_adapter(i))
            elif 'SQK-NSK007' in matching_set_names:
                matching_sets.append(make_old_full_rapid_barcode_adapter(i))

    return matching_sets


def find_adapters_at_read_ends(reads, matching_sets, verbosity, end_size, extra_trim_size,
                               end_threshold, scoring_scheme_vals, print_dest, min_trim_size,
                               threads, check_barcodes, barcode_threshold, barcode_diff,
                               require_two_barcodes, forward_or_reverse_barcodes):
    if verbosity > 0:
        print(bold_underline('Trimming adapters from read ends'),
              file=print_dest)
        name_len = max(max(len(x.start_sequence[0])
                           if x.start_sequence else 0 for x in matching_sets),
                       max(len(x.end_sequence[0])
                           if x.end_sequence else 0 for x in matching_sets))
        for matching_set in matching_sets:
            if matching_set.start_sequence:
                print('  ' + matching_set.start_sequence[0].rjust(name_len) + ': ' +
                      red(matching_set.start_sequence[1]), file=print_dest)
            if matching_set.end_sequence:
                print('  ' + matching_set.end_sequence[0].rjust(name_len) + ': ' +
                      red(matching_set.end_sequence[1]), file=print_dest)
        print('', file=print_dest)

    read_count = len(reads)
    if verbosity == 1:
        output_progress_line(0, read_count, print_dest)

    # Alignments are done natively (in parallel) for a chunk of reads, then applied read by read.
    start_adapters = [x.start_sequence[1] if x.start_sequence else None for x in matching_sets]
    end_adapters = [x.end_sequence[1] if x.end_sequence else None for x in matching_sets]
    for chunk_start in range(0, read_count, NATIVE_CHUNK):
        chunk = reads[chunk_start:chunk_start + NATIVE_CHUNK]
        start_hits = align_regions([read.seq[:end_size] for read in chunk], start_adapters,
                                   scoring_scheme_vals, threads)
        end_hits = align_regions([read.seq[-end_size:] for read in chunk], end_adapters,
                                 scoring_scheme_vals, threads)
        for i, read in enumerate(chunk):
            read.find_start_trim(matching_sets, end_size, extra_trim_size, end_threshold,
                                 scoring_scheme_vals, min_trim_size, check_barcodes,
                                 forward_or_reverse_barcodes, start_hits[i])
            read.find_end_trim(matching_sets, end_size, extra_trim_size, end_threshold,
                               scoring_scheme_vals, min_trim_size, check_barcodes,
                               forward_or_reverse_barcodes, end_hits[i])
            if check_barcodes:
                read.determine_barcode(
                    barcode_threshold, barcode_diff, require_two_barcodes)
            if verbosity == 1:
                output_progress_line(chunk_start+i+1, read_count, print_dest)
            elif verbosity == 2:
                print(read.formatted_start_and_end_seq(end_size, extra_trim_size, check_barcodes),
                      file=print_dest)
            elif verbosity > 2:
                print(read.full_start_end_output(end_size, extra_trim_size, check_barcodes),
                      file=print_dest)

    if verbosity == 1:
        output_progress_line(read_count, read_count,
                             print_dest, end_newline=True)
    if verbosity > 0:
        print('', file=print_dest)


def display_read_end_trimming_summary(reads, verbosity, print_dest):
    if verbosity < 1:
        return
    start_trim_total = sum(x.start_trim_amount for x in reads)
    start_trim_count = sum(1 if x.start_trim_amount else 0 for x in reads)
    end_trim_count = sum(1 if x.end_trim_amount else 0 for x in reads)
    end_trim_total = sum(x.end_trim_amount for x in reads)
    print(int_to_str(start_trim_count).rjust(len(int_to_str(len(reads)))) + ' / ' +
          int_to_str(len(reads)) + ' reads had adapters trimmed from their start (' +
          int_to_str(start_trim_total) + ' bp removed)', file=print_dest)
    print(int_to_str(end_trim_count).rjust(len(int_to_str(len(reads)))) + ' / ' +
          int_to_str(len(reads)) + ' reads had adapters trimmed from their end (' +
          int_to_str(end_trim_total) + ' bp removed)', file=print_dest)
    print('\n', file=print_dest)


def find_adapters_in_read_middles(reads, matching_sets, verbosity, middle_threshold,
                                  extra_trim_good_side, extra_trim_bad_side, scoring_scheme_vals,
                                  print_dest, threads, discard_middle):
    if verbosity > 0:
        verb = 'Discarding' if discard_middle else 'Splitting'
        print(bold_underline(verb + ' reads containing middle adapters'),
              file=print_dest)

    adapters = []
    for matching_set in matching_sets:
        if matching_set.start_sequence:
            adapters.append(matching_set.start_sequence)
        if matching_set.end_sequence:
            if (not matching_set.start_sequence) or \
                    matching_set.end_sequence[1] != matching_set.start_sequence[1]:
                adapters.append(matching_set.end_sequence)

    start_sequence_names = set()
    end_sequence_names = set()
    for matching_set in matching_sets:
        if matching_set.start_sequence:
            start_sequence_names.add(matching_set.start_sequence[0])
        if matching_set.end_sequence:
            end_sequence_names.add(matching_set.end_sequence[0])

    read_count = len(reads)
    if verbosity == 1:
        output_progress_line(0, read_count, print_dest)

    # Alignments are done natively (in parallel) for a chunk of reads, then applied read by read.
    adapter_seqs = [x[1] for x in adapters]
    for chunk_start in range(0, read_count, NATIVE_CHUNK):
        chunk = reads[chunk_start:chunk_start + NATIVE_CHUNK]
        hits = find_middle_hits([read.get_seq_with_start_end_adapters_trimmed() for read in chunk],
                                adapter_seqs, scoring_scheme_vals, middle_threshold, threads)
        for i, read in enumerate(chunk):
            read.find_middle_adapters(adapters, middle_threshold, extra_trim_good_side,
                                      extra_trim_bad_side, scoring_scheme_vals,
                                      start_sequence_names, end_sequence_names, hits[i])
            if verbosity == 1:
                # The original's thread pool loop reported one read ahead.
                done = chunk_start + i + 1
                output_progress_line(done if threads == 1 else done + 1, read_count, print_dest)
            if read.middle_adapter_positions and verbosity > 1:
                print(read.middle_adapter_results(
                    verbosity), file=print_dest, flush=True)

    if verbosity == 1:
        output_progress_line(read_count, read_count,
                             print_dest, end_newline=True)
        print('', flush=True, file=print_dest)


def display_read_middle_trimming_summary(reads, discard_middle, verbosity, print_dest):
    if verbosity < 1:
        return
    middle_trim_count = sum(
        1 if x.middle_adapter_positions else 0 for x in reads)
    verb = 'discarded' if discard_middle else 'split'
    print(int_to_str(middle_trim_count) + ' / ' + int_to_str(len(reads)) + ' reads were ' + verb +
          ' based on middle adapters\n\n', file=print_dest)


def output_reads(reads, out_format, output, read_type, verbosity, discard_middle,
                 min_split_size, print_dest, barcode_dir, input_filename,
                 untrimmed, threads, discard_unassigned):
    if verbosity > 0:
        trimmed_or_untrimmed = 'untrimmed' if untrimmed else 'trimmed'
        if barcode_dir is not None:
            verb = 'Saving '
            destination = 'barcode-specific files'
        elif output is None:
            verb = 'Outputting '
            destination = 'stdout'
        else:
            verb = 'Saving '
            destination = 'file'
        print(bold_underline(verb + trimmed_or_untrimmed + ' reads to ' + destination),
              flush=True, file=print_dest)

    if out_format == 'auto':
        if output is None:
            out_format = read_type.lower()
            if barcode_dir is not None and input_filename.lower().endswith('.gz'):
                out_format += '.gz'
        elif '.fasta.gz' in output.lower():
            out_format = 'fasta.gz'
        elif '.fastq.gz' in output.lower():
            out_format = 'fastq.gz'
        elif '.fasta' in output.lower():
            out_format = 'fasta'
        elif '.fastq' in output.lower():
            out_format = 'fastq'
        else:
            out_format = read_type.lower()

    gzipped_out = False
    gzip_command = 'gzip'
    if out_format.endswith('.gz') and (barcode_dir is not None or output is not None):
        gzipped_out = True
        out_format = out_format[:-3]
        if shutil.which('pigz'):
            if verbosity > 0:
                print('pigz found - using it to compress instead of gzip')
            gzip_command = 'pigz -p ' + str(threads)
        else:
            if verbosity > 0:
                print('pigz not found - using gzip to compress')

    # Output reads to barcode bins.
    if barcode_dir is not None:
        if not os.path.isdir(barcode_dir):
            os.makedirs(barcode_dir)
        barcode_files = {}
        barcode_read_counts, barcode_base_counts = defaultdict(
            int), defaultdict(int)
        for read in reads:
            barcode_name = read.barcode_call
            if discard_unassigned and barcode_name == 'none':
                continue
            if out_format == 'fasta':
                read_str = read.get_fasta(
                    min_split_size, discard_middle, untrimmed)
            else:
                read_str = read.get_fastq(
                    min_split_size, discard_middle, untrimmed)
            if not read_str:
                continue
            if barcode_name not in barcode_files:
                barcode_files[barcode_name] = \
                    open(os.path.join(barcode_dir,
                                      barcode_name + '.' + out_format), 'wt')
            barcode_files[barcode_name].write(read_str)
            barcode_read_counts[barcode_name] += 1
            if untrimmed:
                seq_length = len(read.seq)
            else:
                seq_length = read.seq_length_with_start_end_adapters_trimmed()
            barcode_base_counts[barcode_name] += seq_length
        table = [['Barcode', 'Reads', 'Bases', 'File']]

        for barcode_name in sorted(barcode_files.keys()):
            barcode_files[barcode_name].close()
            bin_filename = os.path.join(
                barcode_dir, barcode_name + '.' + out_format)

            if gzipped_out:
                if not os.path.isfile(bin_filename):
                    continue
                bin_filename_gz = bin_filename + '.gz'
                if os.path.isfile(bin_filename_gz):
                    os.remove(bin_filename_gz)
                try:
                    subprocess.check_output(gzip_command + ' ' + bin_filename,
                                            stderr=subprocess.STDOUT, shell=True)
                except subprocess.CalledProcessError:
                    pass
                bin_filename = bin_filename_gz

            table_row = [barcode_name, int_to_str(barcode_read_counts[barcode_name]),
                         int_to_str(barcode_base_counts[barcode_name]), bin_filename]
            table.append(table_row)

        if verbosity > 0:
            print('')
            print_table(table, print_dest, alignments='LRRL',
                        max_col_width=60, col_separation=2)

    # Output to all reads to stdout.
    elif output is None:
        for read in reads:
            read_str = read.get_fasta(min_split_size, discard_middle) if out_format == 'fasta' \
                else read.get_fastq(min_split_size, discard_middle)
            print(read_str, end='')
        if verbosity > 0:
            print('Done', flush=True, file=print_dest)

    # Output to all reads to file.
    else:
        if gzipped_out:
            out_filename = 'TEMP_' + str(os.getpid()) + '.fastq'
        else:
            out_filename = output
        with open(out_filename, 'wt') as out:
            for read in reads:
                read_str = read.get_fasta(min_split_size, discard_middle) if out_format == 'fasta' \
                    else read.get_fastq(min_split_size, discard_middle)
                out.write(read_str)
        if gzipped_out:
            subprocess.check_output(gzip_command + ' -c ' + out_filename + ' > ' + output,
                                    stderr=subprocess.STDOUT, shell=True)
            os.remove(out_filename)
        if verbosity > 0:
            print('\nSaved result to ' + os.path.abspath(output), file=print_dest)

    if verbosity > 0:
        print('', flush=True, file=print_dest)


def output_progress_line(completed, total, print_dest, end_newline=False, step=10):
    if step > 1 and completed % step != 0 and completed != total:
        return
    progress_str = int_to_str(completed) + ' / ' + int_to_str(total)
    if total > 0:
        percent = 100.0 * completed / total
    else:
        percent = 0.0
    progress_str += ' (' + '%.1f' % percent + '%)'

    end_char = '\n' if end_newline else ''
    print('\r' + progress_str, end=end_char, flush=True, file=print_dest)
