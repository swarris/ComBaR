'''This module contains the output formatters for combar'''

from pyPaSWAS.Core.Formatters import  DefaultFormatter

        
class PlotterFormatter(DefaultFormatter):
    def __init__(self, logger, hitlist, outputfile):
        DefaultFormatter.__init__(self, logger, hitlist, outputfile)

    def _format_hit(self, hit):
        seqID = hit.sequence_info.id
        seqLoc =hit.seq_location
        distance = hit.target_info.distance
        
        if seqID[-3:] == "_RC":
            seqID = seqID[:-3]
            seqLoc = (hit.sequence_info.original_length - hit.seq_location[0], hit.sequence_info.original_length - hit.seq_location[1])
            distance = -distance
        targetID = hit.target_info.id
        targetLoc = hit.target_location
        if targetID[-3:] == "_RC":
            targetID = targetID[:-3]
            targetLoc = (hit.target_info.original_length - hit.target_location[0], hit.target_info.original_length - hit.target_location[1])
            distance = -distance
            
        return '\t'.join([seqID, str(seqLoc[0]), str(seqLoc[1]), targetID, str(targetLoc[0]),str(targetLoc[1]), "distance="+str(distance)])+"\n"

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

