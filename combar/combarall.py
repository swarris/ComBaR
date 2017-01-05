from pyPaSWAS import pyPaSWAS
from pyPaSWAS.Core import resource_filename
from pyPaSWAS.Core.Exceptions import InvalidOptionException, ReaderException
from pyPaSWAS.Core.Readers import BioPythonReader
from pyPaSWAS.Core.Scores import BasicScore, CustomScore, DnaRnaScore
from pyPaSWAS.Core.Formatters import SamFormatter, DefaultFormatter
from pyPaSWAS import set_logger, normalize_file_path 

from combar.Core.Programs import ComBaRMapper,  ComBaRIndexer, GenomePlotter 
from combar.Core.Formatters import PlotterFormatter
from combar.Core.HitList import HitListComBaR
from combar import parse_cli


import logging
import os.path


class ComBaR(pyPaSWAS):
    '''
    This class represents the main program. It parses the command line and runs
    one of the programs from the pyPaSWAS suite.
    Settings and arguments are stored in two variables created in parse_cli.

    An effort is made to comply as much as possible with the options and arguments
    as used by NCBI blastall version 2.2.21.
    '''
    def __init__(self):
        pyPaSWAS.__init(self)
        self.config_file = resource_filename(__name__, '/Core/cfg/defaults.cfg')
        self._get_default_logger()
        self.settings = None
        self.arguments = None

    def _get_formatter(self, results):
        '''Sets the format of the output file.'''
        formatter = ''
        if self.output_format == "SAM":
            formatter = SamFormatter(self.logger, results, self.outputfile)
        elif self.output_format == "plot":
            formatter = PlotterFormatter(self.logger, results, self.outputfile)
        else:
            formatter = DefaultFormatter(self.logger, results, self.outputfile)
        return formatter


    def _set_scoring_matrix(self):
        ''' Instantiate the scoring matrix. For more information refer to the
            documentation of the different score types. '''
        matrix_name = self.settings.matrix_name.upper()
        score = None
        if matrix_name == 'DNA-RNA':
            score = DnaRnaScore(self.logger, self.settings)
        elif matrix_name == 'BASIC':
            score = BasicScore(self.logger, self.settings)
        elif matrix_name == 'CUSTOM':
            score = CustomScore(self.logger, self.settings)
        else:
            raise InvalidOptionException(matrix_name + ' is not a valid substitution matrix')
        self.score = score

    def _set_output_format(self):
        '''Determines the output format.
            Currently only TXT for text and SAM for SAM are supported.
        '''
        if self.settings.out_format.upper() == 'SAM':
            self.output_format = 'SAM'
        elif self.settings.out_format.upper() == 'TXT':
            self.output_format = 'TXT'
        elif self.settings.out_format.upper() == 'PLOT':
            self.output_format = 'plot'
        elif self.settings.out_format.upper() == "GRAPH":
            self.output_format = 'GRAPH'
        else:
            raise InvalidOptionException('Invalid output format {0}.'.format(self.settings.out_format))

    def _set_program(self):
        '''Determines what program from the suite should be used and instantiates it'''
        if self.settings.program == 'mapper':
            self.program = ComBaRMapper(self.logger, self.score, self.settings, self.arguments)
            self.logger.warning("Removing limits on length of sequences for ComBaR mapping!")
            self.settings.limit_length = 10**20
        elif self.settings.program == "plotter":
            self.program = GenomePlotter(self.logger, self.score, self.settings, self.arguments)
            self.logger.warning("Removing limits on length of sequences for genome plotter!")
            self.settings.limit_length = 10**20
        else:
            raise InvalidOptionException('Invalid program selected {0}'.format(self.settings.program))

