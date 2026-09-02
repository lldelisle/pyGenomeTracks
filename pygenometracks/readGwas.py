# -*- coding: utf-8 -*-
import collections

from .utilities import InputError, to_string


class ReadGwas(object):
    """
    Reads a GWAS file. The expected fields are:
    chromosome, position, name, and pvalue.

    Example:
    gwas = ReadGwas(open("file.gwas", 'r'))
    for record in gwas:
        print(record.chromosome, record.position, record.pvalue)
    """

    def __init__(self, file_handle, has_header=False):
        """
        :param file_handle: file handle
        """
        self.file_handle = file_handle
        self.line_number = 0

        # Define the fields for GWAS
        self.fields = ['chromosome', 'position', 'name', 'pvalue']
        self.GwasRecord = collections.namedtuple('GwasRecord', self.fields)

        # Skip the header line if present
        if has_header:
            next(self.file_handle)
            self.line_number += 1

    def __iter__(self):
        return self

    def get_no_comment_line(self):
        """
        Skips comment lines starting with '#' or empty lines.
        :return: a valid line
        """
        line = next(self.file_handle)
        line = to_string(line)
        if line.startswith("#") or line.strip() == '':
            line = self.get_no_comment_line()

        self.line_number += 1
        return line

    def __next__(self):
        """
        :return: GwasRecord object
        """
        line = self.get_no_comment_line()
        return self.get_gwas_record(line)

    def get_gwas_record(self, gwas_line):
        """
        Processes each line from a GWAS file and returns a namedtuple object.

        :param gwas_line: a single line from the GWAS file
        :return: GwasRecord object
        """
        line_data = gwas_line.strip()
        line_data = to_string(line_data)
        line_data = line_data.split("\t")

        if len(line_data) < 4:
            raise InputError(f"Line {self.line_number} does not have 4 fields: {gwas_line}."
                             f"We expect at least 4 field, corresponding to: chromosome, position, name, pvalue.")

        try:
            chromosome = line_data[0]
            position = int(line_data[1])
            name = line_data[2]
            pvalue = float(line_data[3])
        except ValueError as e:
            raise InputError(f"Error parsing line {self.line_number}: {gwas_line}\n{e}")

        return self.GwasRecord(chromosome, position, name, pvalue)
