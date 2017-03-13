"""Prost TSV File to Excel Converter.

Usage:
  prost_excel_converter [-hvom] [--in <in_prefix>] [--out <out_prefix>] [TSV_FILES ...]

This converter takes Prost output TSV (tab separated values) files and converts
them to one or more Excel Spreadsheets.  It also formats and conditionally
formats the sheets (for example, adding color to percentage columns).

By default, prost_excel_converter will combine all TSV files into one workbook,
with each TSV file converted into a different sheet within the workbook.

Arguments:
    TSV_FILES   Optional list of TSV input files.

Options:
    -h, --help          Show this screen.
    -v, --version       Show the version.
    -m, --many          Create multiple workbooks, one for each TSV file
    --in <in_prefix>    The TSV input file prefix [default: prost_output]
    --out <out_prefix>  The Excel output file prefix [default: prost_output]

Unimplemented:
    Options for controlling the color (e.g. red), the range (e.g. 0-30%), etc.
    Options for controlling the two color scheme.

"""

# Copyright (C) 2014-2017 Peter Batzel and Jason Sydes
#
# This file is part of Prost!.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the License with this program.
#
# Written by Jason Sydes and Peter Batzel.
# Advised by Thomas Desvignes.

###############
### Imports ###
###############

# Python 3 imports
from __future__ import absolute_import
from __future__ import division

# version
from prost._version import __version__

## Prost imports
import prost.constants

## Other imports

from docopt import docopt
import os
import re
import xlsxwriter


#################
### Constants ###
#################

# Workbook-wide defaults
FORMAT_BOOK_DEFAULTS = {
    'strings_to_numbers': True
}

# Format percentage cells (0.00%)
FMT_PERCENT = {
    'num_format': '0.00%'
}

# Default two-color-conditional formatting for all '%' columns.
FMT_PERCENT_CONDITIONAL_TWO_COLOR = {
    'type': '2_color_scale',
    'min_type': 'percent',
    'max_type': 'percent',
    'min_color': 'white',
    'max_color': 'red',
    'min_value': 0,
    'max_value': 40,
}
FMT_PERCENT_CONDITIONAL_THREE_COLOR = {
    'type': '3_color_scale',
    'min_type': 'num',
    'mid_type': 'num',
    'max_type': 'num',
    'min_color': 'blue',
    'mid_color': 'white',
    'max_color': 'red',
    'min_value': -2,
    'mid_value': 0,
    'max_value': 2,
}


#################
### Functions ###
#################

def to_xlsx_filename(tsv_filename, in_prefix, out_prefix):
    """Convert a tsv filename to an xlsx filename.

    This is only used when in the '--many' mode.  See above.

    Arguments:
        tsv_filename (str): The TSV filename.
        in_prefix (str): The TSV input filename prefix (e.g. 'prost_output').
        out_prefix (str): The TSV output filename prefix (e.g. 'prost_output').

    Returns:
        (str): The output workbook filename.

    """
    name, ext = os.path.splitext(tsv_filename)
    if ext != ".tsv" or name == "":
        raise Exception
    # replace in_prefix with out_prefix in name
    return name.replace(in_prefix, out_prefix) + '.xlsx'

def to_xlsx_sheetname(tsv_filename, in_prefix):
    """Convert a tsv filename to an xlsx sheet name.

    Because MS Excel sheet names have to be short, we also chop out
    "compressed" and "uncompressed".

    Examples:
        This is the naming scheme we've got currently.

        prost_output_candidate_arm_switch.tsv           -> candidate_arm_switch
        prost_output_candidate_mirror-miRNAs.tsv        -> candidate_mirror-miRNAs
        prost_output_compressed_by_annotation.tsv       -> by_annotation
        prost_output_compressed_by_genomic_location.tsv -> by_genomic_location
        prost_output_compressed_by_seed.tsv             -> by_seed
        prost_output_no_genomic_hits.tsv                -> no_genomic_hits
        prost_output_uncompressed_isomiRs.tsv           -> isomiRs

    Arguments:
        tsv_filename (str): The TSV filename.
        in_prefix (str): The TSV input filename prefix (e.g. 'prost_output').

    Returns:
        (str): The shortened sheet name.
    """
    name, ext = os.path.splitext(tsv_filename)
    if ext != ".tsv" or name == "":
        raise Exception
    # Get rid of prefix (e.g. 'prost_output_compressed_by_seed' -> 'compressed_by_seed.tsv')
    name = name.replace("{}_".format(in_prefix), '')
    # Get rid of "compressed" or "uncompressed"
    name = name.replace("uncompressed_", "")
    name = name.replace("compressed_", "")

    return name

def tsv_to_header(tsv_file):
    """Extract the header line from a TSV file.

    Arguments:
        tsv_file (str): The name of the TSV file.

    Returns:
        [str]: An array of strings, where each string represents one cell of
            header.
    """
    line = None
    with open(tsv_file, 'r') as f:
        line = f.readline()
    return line.split("\t")

def header_to_first_and_last_percent_columns(header):
    """Take the header and return 2-tuple representing the range of "percent
    columns" (e.g. %seed_shifted) in the header. For example, if column 18 is
    the first % column, and 29 is the last % column, then this function will
    return (18,29).

    Arguments:
        header ([str]): The header row of the worksheet.

    Returns:
        (int,int): A two tuple representing the range of percent columns in the
            passed header. Note that (None, None) is returned if the header
            does not have any percent columns.
    """
    perc_cols = []
    for cell in header:
        if re.search(r'%', cell):
            perc_cols.append(True)
        else:
            perc_cols.append(False)

    first, last = None, None
    found_first = False
    for i, col in enumerate(perc_cols):
        if not found_first and col is True:
            first = i
            found_first = True
        elif found_first and col is False:
            last = i - 1
            break
        else:
            continue

    return (first, last)

def header_to_first_and_last_log_fold_change_expr_5p_3p_columns(header):
    """Take the header and return 2-tuple representing the range of
    "log_fold_change_expr_5p_3p" columns in the header.  For example, if column
    18 is the first "log_fold_change_expr_5p_3p" column, and 29 is the last
    "log_fold_change_expr_5p_3p" column, then this function will return
    (18,29).

    Arguments:
        header ([str]): The header row of the worksheet.

    Returns:
        (int,int): A two tuple representing the range of
            log_fold_change_expr_5p_3p columns in the passed header. Note that
            (None, None) is returned if the header does not have any
            log_fold_change_expr_5p_3p columns.
    """
    log_fold_change_expr_5p_3p_cols = []
    for cell in header:
        if re.search(r'log_fold_change_expr_5p_3p$', cell):
            log_fold_change_expr_5p_3p_cols.append(True)
        else:
            log_fold_change_expr_5p_3p_cols.append(False)

    first, last = None, None
    found_first = False
    for i, col in enumerate(log_fold_change_expr_5p_3p_cols):
        if not found_first and col is True:
            first = i
            found_first = True
        elif found_first and col is False:
            last = i - 1
            break
        else:
            continue

    return (first, last)

def convert_tsv_to_sheet(book, tsv_file, in_prefix):
    """Convert a TSV file to a worksheet, and add the sheet to the passed
    Workbook object.

    Arguments:
        book (xlsxwriter.workbook.Workbook): The Excel Workbook to which the
            converted TSV file converted to a worksheet will be added.
        tsv_file (str): The TSV filename.
        in_prefix (str): The TSV input filename prefix (e.g. 'prost_output').

    Returns:
        (int, xlsxwriter.worksheet.Worksheet): A tuple containing a) the number
            of rows in the tsv_file, and b) the Excel Worksheet object created
            from the tsv_file.

    """
    header = tsv_to_header(tsv_file)
    sheetname = to_xlsx_sheetname(tsv_file, in_prefix)
    sheet = book.add_worksheet(sheetname)

    num_rows = 0
    with open(tsv_file, 'r') as tsv:
        for row_num, line in enumerate(tsv):
            row = line.rstrip().split('\t')
            sheet.write_row(row_num, 0, row)
            num_rows = row_num

    return (num_rows, sheet)

def perform_formatting(sheet, header, num_rows, fmt_percent,
            fmt_cond_two_color, fmt_cond_three_color):
    """Perform all formatting for a given sheet.

    Arguments:
        sheet (xlsxwriter.worksheet.Worksheet): The Excel Worksheet to be
            formatted.
        header ([str]): The header row of the worksheet.
        num_rows (int): The number of rows in this sheet.
        fmt_percent (xlsxwriter.format.Format): The 0.00% percentage format.
        fmt_cond_two_color (xlsxwriter.format.Format): The two-color conditional
            formatting to apply.
        fmt_cond_three_color (xlsxwriter.format.Format): The three-color conditional
            formatting to apply.
    """
    # Format the percent columns as 0.00%
    sheet_format_percent(sheet, header, fmt_percent)
    # Format the percent columns as conditional two color.
    sheet_format_percent_cond_two_color(sheet, header, num_rows,
            fmt_cond_two_color)
    # Format the percent columns as conditional three color.
    sheet_format_log_fold_change_expr_5p_3p_cond_three_color(sheet, header,
            num_rows, fmt_cond_three_color)
    # Freeze the header row
    sheet.freeze_panes(1, 0)

def sheet_format_percent(sheet, header, fmt):
    """Format an Excel Worksheet such that all columns with a '%' in their
    column header name to be formatted as percentages (i.e. 0.00%).

    Arguments:
        sheet (xlsxwriter.worksheet.Worksheet): The Excel Worksheet to be
            formatted.
        header ([str]): The header row of the worksheet.
        fmt (xlsxwriter.format.Format): The 0.00% percentage format.
    """

    # Determine first and last percent columns
    first, last = header_to_first_and_last_percent_columns(header)
    if first:
        # Convert "percent columns" to type 'percent'
        sheet.set_column(first, last, None, fmt)

def sheet_format_percent_cond_two_color(sheet, header, num_rows, fmt):
    """Format an Excel Worksheet such that all columns with a '%' in their
    column header name are two-color-conditionally formatted.

    Arguments:
        sheet (xlsxwriter.worksheet.Worksheet): The Excel Worksheet to be
            formatted.
        header ([str]): The header row of the worksheet.
        num_rows (int): The number of rows in this sheet.
        fmt (xlsxwriter.format.Format): The conditional formatting to apply.
    """

    # Determine first and last percent columns
    first, last = header_to_first_and_last_percent_columns(header)
    if first:
        # Conditionally format "percent columns"
        sheet.conditional_format(1, first, num_rows, last, fmt)

def sheet_format_log_fold_change_expr_5p_3p_cond_three_color(sheet, header, num_rows, fmt):
    """Format an Excel Worksheet such that all columns with a
    'log_fold_change_expr_5p_3p' in their column header name are
    three-color-conditionally formatted.

    Arguments:
        sheet (xlsxwriter.worksheet.Worksheet): The Excel Worksheet to be
            formatted.
        header ([str]): The header row of the worksheet.
        num_rows (int): The number of rows in this sheet.
        fmt (xlsxwriter.format.Format): The conditional formatting to apply.
    """

    # Determine first and last percent columns
    first, last = header_to_first_and_last_log_fold_change_expr_5p_3p_columns(header)
    if first:
        # Conditionally format "percent columns"
        sheet.conditional_format(1, first, num_rows, last, fmt)

def setup_workbook(filename):
    """Setup the workbook.

    Arguments:
        filename (str): The filename to be given to the workbook.

    Returns:
        (xlsxwriter.workbook.Workbook, xlsxwriter.format.Format): A tuple
            consisting of a) the workbook just created, b) the "percentage
            columns" Format.

    """
    # Create the book
    book = xlsxwriter.Workbook(filename, FORMAT_BOOK_DEFAULTS)

    # Add the percentage columns 0.00% formating
    fmt_percent = book.add_format(FMT_PERCENT)

    return (book, fmt_percent)

def command_line():
    """Entry point when calling this script from the command line.

    This function performs basic parsing of command line options, then hands
    off.

    By default, if no TSV files are passed in, all Prost! output files will be
    processed.

    """
    args = docopt(__doc__, version=__version__)

    # in_prefix is always what the user tells us (out_prefix may not be, see
    # below)
    in_prefix = args['--in']

    # Build the list of tsv_filenames if not passed
    tsv_files = args['TSV_FILES']
    if not tsv_files:
        # Build the full list
        for suffix in prost.constants.OUTPUT_FILE_SUFFIXES:
            if in_prefix:
                tsv_files.append("{}_{}.tsv".format(in_prefix, suffix))
            else:
                tsv_files.append("{}.tsv".format(suffix))

    # Create one or many Workbooks
    if args['--many']:
        # out_prefix is exactly what the user tells us when running --many
        # mode.
        out_prefix = args['--out']

        # Create many workbooks, one for each TSV file.
        filenames = create_many(in_prefix, out_prefix, tsv_files)
        print "Created {}.".format(", ".join(filenames))
    else:
        # Do not allow blank out_prefix (as we would have nothing to name our
        # xlsx file)
        if args['--out']:
            out_prefix = args['--out']
        else:
            out_prefix = prost.constants.DEFAULT_OUTPUT_FILE_PREFIX

        # Create one big workbook
        filename = create_one(in_prefix, out_prefix, tsv_files)
        print "Created {}.".format(filename)

def create_one(in_prefix, out_prefix, tsv_files):
    """Create one workbook containing one worksheet for each TSV file.

    Arguments:
        in_prefix (str): The TSV input filename prefix (e.g. 'prost_output').
        out_prefix (str): The TSV output filename prefix (e.g. 'prost_output').
        tsv_files ([str]): The list of TSV filenames to convert.

    Returns:
        str: The name of the Excel Workbook created.
    """

    book_filename = "{}.xlsx".format(out_prefix)
    book, fmt_percent = setup_workbook(book_filename)
    for tsv_file in tsv_files:
        num_rows, sheet = convert_tsv_to_sheet(book, tsv_file, in_prefix)
        header = tsv_to_header(tsv_file)
        perform_formatting(sheet, header, num_rows, fmt_percent,
                FMT_PERCENT_CONDITIONAL_TWO_COLOR,
                FMT_PERCENT_CONDITIONAL_THREE_COLOR)
    return book_filename

def create_many(in_prefix, out_prefix, tsv_files):
    """Create multiple workbooks, one for each TSV file.

    Arguments:
        in_prefix (str): The TSV input filename prefix (e.g. 'prost_output').
        out_prefix (str): The TSV output filename prefix (e.g. 'prost_output').
        tsv_files ([str]): The list of TSV filenames to convert.

    Returns:
        [str]: The names of the Excel Workbooks created.
    """

    book_filenames = []
    for tsv_file in tsv_files:
        book_filename = to_xlsx_filename(tsv_file, in_prefix, out_prefix)
        book_filenames.append(book_filename)
        book, fmt_percent = setup_workbook(book_filename)
        num_rows, sheet = convert_tsv_to_sheet(book, tsv_file, in_prefix)
        header = tsv_to_header(tsv_file)
        perform_formatting(sheet, header, num_rows, fmt_percent,
                FMT_PERCENT_CONDITIONAL_TWO_COLOR,
                FMT_PERCENT_CONDITIONAL_THREE_COLOR)
    return book_filenames


# vim: softtabstop=4:shiftwidth=4:expandtab
