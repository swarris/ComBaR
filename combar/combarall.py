from pyPaSWAS.pypaswasall import Pypaswas
from pyPaSWAS.Core import resource_filename
from pyPaSWAS.Core.Exceptions import InvalidOptionException, ReaderException
from pyPaSWAS.Core.Readers import BioPythonReader
from pyPaSWAS.Core.Scores import BasicScore, CustomScore, DnaRnaScore
from pyPaSWAS.Core.Formatters import SamFormatter, DefaultFormatter
from pyPaSWAS import set_logger, normalize_file_path 
from pyPaSWAS.Core.HitList import HitList 
from combar.Core.Programs import ComBaRMapper, GenomePlotter,ReadDistance, QGramLinker
from combar.Core.Formatters import PlotterFormatter
from combar import parse_cli


import logging
import os.path


class ComBaR(Pypaswas):
    '''
    This class represents the main program. It parses the command line and runs
    one of the programs from the pyPaSWAS suite.
    Settings and arguments are stored in two variables created in parse_cli.

    An effort is made to comply as much as possible with the options and arguments
    as used by NCBI blastall version 2.2.21.
    '''
    def __init__(self):
        Pypaswas.__init__(self, config = resource_filename(__name__, '/Core/cfg/defaults.cfg'))
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
            self.program = QGramLinker(self.logger, self.score, self.settings, self.arguments)
            self.logger.warning("Removing limits on length of sequences for genome plotter!")
            self.settings.limit_length = 10**20
            self.logger.warning("Forcing output to csv for plotting")
            self.output_format = "plot"
        else:
            raise InvalidOptionException('Invalid program selected {0}'.format(self.settings.program))

    def run(self):
        '''The main program of pyPaSWAS.'''
        # Read command-line arguments
        self.settings, self.arguments = parse_cli(self.config_file)
        self.logger = set_logger(self.settings)
        self.logger.info("Initializing application...")
        self._set_outfile()
        self._set_scoring_matrix()
        self.logger.info('Application initialized.')
        self.logger.info('Setting program...')
        self._set_output_format()
        self._set_program()
        self.logger.info('Program set.')
        
        queriesToProcess = True
        
        query_start = int(self.settings.start_query)
        query_end = int(self.settings.start_query) + int(self.settings.query_step)
        if query_end > int(self.settings.end_query) and int(self.settings.start_query) != int(self.settings.end_query):
            query_end = int(self.settings.end_query)
        
        start_index = int(self.settings.start_target)

        end_index = int(self.settings.start_target) + int(self.settings.sequence_step) 
        if end_index > int(self.settings.end_target) and int(self.settings.start_target) != int(self.settings.end_target):
            end_index = int(self.settings.end_target)
        
        results = HitList(self.logger)
        
        while queriesToProcess:
            self.logger.info('Reading query sequences {} {}...'.format(query_start, query_end))
            try:
                if self.settings.program == "plotter":
                    query_sequences = self._get_target_sequences(self.arguments[0], start=query_start, end=query_end)
                else:
                    query_sequences = self._get_query_sequences(self.arguments[0], start=query_start, end=query_end)
                self.logger.info('Query sequences OK.')
            except ReaderException:
                queriesToProcess = False

            sequencesToProcess = True
            if not self.settings.program == "palindrome":
                start_index = int(self.settings.start_target)
                end_index = int(self.settings.start_target) + int(self.settings.sequence_step) 
                if end_index > int(self.settings.end_target) and int(self.settings.start_target) != int(self.settings.end_target):
                    end_index = int(self.settings.end_target)
            
            while queriesToProcess and sequencesToProcess:
                
                self.logger.info('Reading target sequences {}, {}...'.format(start_index,end_index))
                try:
                    target_sequences = self._get_target_sequences(self.arguments[1], start=start_index, end=end_index)
                    self.logger.info('Target sequences OK.')
                except ReaderException:
                    sequencesToProcess = False
                    
    
                if not sequencesToProcess or not queriesToProcess or len(query_sequences) == 0 or len(target_sequences) == 0:
                    sequencesToProcess = False
                    self.logger.info('Processing done')
                else:
                    self.logger.info('Processing {0}- vs {1}-sequences'.format(len(query_sequences),
                                                                            len(target_sequences)))
                    results.extend(self.program.process(query_sequences, target_sequences, self))
                
                if sequencesToProcess and len(target_sequences) <= end_index:
                    # for palindrome program, skip directly to next
                    sequencesToProcess = False
                    self.logger.info('Processing done')

                start_index = start_index + int(self.settings.sequence_step)
                end_index = end_index + int(self.settings.sequence_step)
                if  self.settings.program == "palindrome" or (int(self.settings.end_target) > 0 and int(self.settings.end_target) < end_index):
                    sequencesToProcess = False

            if int(self.settings.end_query) > 0 and int(self.settings.end_query) < query_end:
                queriesToProcess = False
            query_start = query_start + int(self.settings.query_step)
            query_end = query_end + int(self.settings.query_step)


        nhits = len(results.hits)
        # retrieve and print results!
        self.logger.info('Processing OK ({} hits found).'.format(nhits))
        if nhits > 0:
            self.logger.info('Formatting output...')
            formatter = self._get_formatter(results)
            self.logger.info('Formatting OK.')

            self.logger.info('Writing output...')
            formatter.print_results()
            self.logger.info('Writing OK.')
            self.logger.info('Finished')
        else:
            self.logger.warning('No suitable hits produced, exiting...')
