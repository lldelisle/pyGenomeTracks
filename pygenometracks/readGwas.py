# -*- coding: utf-8 -*-
from .readTabular import ReadTabular
from .utilities import InputError


class ReadGwas(ReadTabular):
    """
    Reads a GWAS file. The expected fields are:
    chromosome, position, name, and pvalue.

    Example:
    gwas = ReadGwas("file.gwas")
    for record in gwas:
        print(record.chromosome, record.position, record.pvalue)
    """
    # Define the fields for GWAS
    fields = ['chromosome', 'position', 'name', 'pvalue']

    def __init__(self, file_handle, has_header=False):
        """
        :param file_handle: file handle
        """
        super().__init__(file_handle)

        # Skip the header line if present
        if has_header:
            next(self.file_handle)
            self.line_number += 1
            self.length -= 1

    def get_record(self, gwas_line):
        """
        Processes each line from a GWAS file and returns a namedtuple object.

        :param gwas_line: a single line from the GWAS file
        :return: Record object
        """
        line_data = self.get_line_data(gwas_line)

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

        return self.Record(chromosome, position, name, pvalue)
