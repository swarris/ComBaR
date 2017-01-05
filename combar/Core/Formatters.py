'''This module contains the output formatters for combar'''

from pyPaSWAS.Core.Formatters import  DefaultFormatter

        
class PlotterFormatter(DefaultFormatter):
    def __init__(self, logger, hitlist, outputfile):
        DefaultFormatter.__init__(self, logger, hitlist, outputfile)

    def _format_hit(self, hit):
        return '\t'.join([hit.sequence_info.id, hit.target_info.id, str(hit.seq_location[0]),str(hit.target_location[0]), str(hit.target_info.distance)])+"\n"

    def print_results(self, hitlist):
        self.hitlist = hitlist
        '''sets, formats and prints the results to a file.'''
        self.logger.info('plotting results...')
        #format header and hit lines
        for hit in self.hitlist.real_hits.itervalues():
            self.output.write(self._format_hit(hit))

    def open(self):
        self.output = open(self.outputfile, 'w')

    def close(self):
        self.output.close()
        self.logger.debug('finished plotting results')

