''' This module contains the programs from the pyPaSWAS suite '''
from pyPaSWAS.Core.HitList import HitList
from operator import itemgetter,attrgetter
from pyPaSWAS.Core.Programs import Aligner

from pyPaSWAS.Core.SWSeqRecord import SWSeqRecord
from Bio.Seq import Seq
from Hit import Distance
from pyPaSWAS.Core.Hit import Hit


class ComBaRMapper(Aligner):
    
    def __init__(self, logger, score, settings, arguments):
        Aligner.__init__(self, logger, score, settings)
        self.arguments = arguments
        self.qindexerCUDA = False
        self.qindexerOCL = False
        if (self.settings.framework.upper() == 'OPENCL'):
            self.qindexerOCL = True
        else:
            self.qindexerCUDA = True

    def process(self, records_seqs, targets, pypaswas):
        '''This methods creates index files for targets based on the length of the records.
        '''

        
        # step through the targets                                                                                                                                                                           
        self.logger.debug('ComBaR mapping...')

        indexer = None
        while len(records_seqs) > 0:
            prevLength = len(records_seqs[0])
            #create indexer
            
            if self.qindexerOCL:
                from pyPaSWAS.Core.QIndexerOCL import QIndexerOCL
                indexer = QIndexerOCL(self.settings, self.logger, 0.1, records_seqs[0:1], int(self.settings.qgram))
            elif self.qindexerCUDA:
                from pyPaSWAS.Core.QIndexerCUDA import QIndexerCUDA
                indexer = QIndexerCUDA(self.settings, self.logger, 0.1, records_seqs[0:1], int(self.settings.qgram))                
            else:
                from pyPaSWAS.Core.QIndexer import QIndexer
                indexer = QIndexer(self.settings, self.logger, 0.1, records_seqs[0:1], int(self.settings.qgram))
            """
            indexer = QIndexer(self.settings, self.logger, 0.1, records_seqs[0:1], int(self.settings.qgram))
            """
            #while indices to process, process all reads with same length
            # when done, remove these reads
            while indexer.indicesToProcessLeft():
                currentRead = 0
                # only create index with first read of same length 
                if currentRead == 0:
                    indexer.createIndexAndStore(targets, self.arguments[1])
                 
                while currentRead < len(records_seqs) and len(records_seqs[currentRead]) == prevLength:
                    currentBlockOfReads = currentRead
                    
                    while currentRead - currentBlockOfReads < indexer.readsToProcess and currentRead < len(records_seqs) and len(records_seqs[currentRead]) == prevLength:
                        currentRead += 1


                    allLocations = indexer.findIndices(records_seqs[currentBlockOfReads:currentRead])

                    if len(allLocations) > 0 :
                        for read in xrange(currentBlockOfReads, currentRead):
                            
                            firstRead = records_seqs[read]
                            filteredRecordsSeqs = [firstRead]
            
                            locations = allLocations[read-currentBlockOfReads]
                            locs = []
                            if (len(locations) > 0):
                                self.logger.info("Processing seq: " + records_seqs[read].id)
                                for value in locations.itervalues():
                                    locs.extend(value)

                                if len(locs) > float(self.settings.fraction_of_seeds) * indexer.indicesStepSize:
                                    self.logger.warning("Too many locations for fast processing. Sorting and slicing now with fraction_of_seeds = {}.".format(self.settings.fraction_of_seeds))
                                    locs.sort(key = lambda x:x[2]) # sort on distance
                                    locs = locs[:int(float(self.settings.fraction_of_seeds)*len(locs))]

                                splittedTargets = []
        
                                for loc in locs:
                                    swSeqRecord = indexer.getSWSeqRecord(loc, targets)
                                    #if len(swSeqRecord) > 0:
                                    swSeqRecord.distance = loc[2]
                                    swSeqRecord.id = targets[loc[0][1]].id
                                    swSeqRecord.refID = loc[0][1]
                                    splittedTargets.append(swSeqRecord)
        
                                if (len(splittedTargets) > 0 and len(filteredRecordsSeqs) > 0):
                                    splittedTargets.sort(key=lambda seqIO : len(seqIO.seq), reverse=True)
                                    target_index = 0
                                    # process of the seeds:                                                                                                                                                          
                                    while target_index < len(splittedTargets):
                                        last_target_index = self.smith_waterman.set_targets(splittedTargets, target_index, None, filteredRecordsSeqs)
                                        self.logger.debug('At target: {0} of {1}, processing up to {2}'.format(target_index, len(splittedTargets), str(last_target_index)))
                                        results = self.smith_waterman.align_sequences(filteredRecordsSeqs, splittedTargets, target_index)
                                        self.hitlist.extend(results)
                                        target_index = last_target_index
            
            #filter out reads already processed:
            currentRead = 0
            while currentRead < len(records_seqs) and len(records_seqs[currentRead]) == prevLength:
                currentRead += 1
            records_seqs = records_seqs[currentRead:]
        
        if indexer != None and self.qindexerCUDA:
            indexer.pop_context()                             
        self.logger.debug('ComBaR mapping finished.')
        return self.hitlist


class GenomePlotter(Aligner):
    
    def __init__(self, logger, score, settings, arguments):
        Aligner.__init__(self, logger, score, settings)
        self.arguments = arguments
        self.qindexerCUDA = False
        self.qindexerOCL = False
        if (self.settings.framework.upper() == 'OPENCL'):
            self.qindexerOCL = True
        else:
            self.qindexerCUDA = True


    def process(self, records_seqs, targets, pypaswas):
        '''This methods creates indices and determines differences between sequences.
        '''

        
        # step through the targets                                                                                                                                                                           
        # step through the targets                                                                                                                                                                          
        self.logger.debug('Genome plotter...')
        self.logger.info("Clearing output file")
        formatter = pypaswas._get_formatter(self.hitlist)
        formatter.open()
        indexer = None
        plotter = None

        window = int(self.settings.window_length)
        stepSize = int(0.1 * window)
        block = 100
        indexStepSize = 1000

        seq = records_seqs[0]
        dummySeq = [SWSeqRecord(Seq(str(seq.seq[0:window]), seq.seq.alphabet), seq.id, 0, original_length = len(seq.seq), distance = 0, refID = seq.id)]

        if self.qindexerOCL:
            from pyPaSWAS.Core.QIndexerOCL import QIndexerOCL, GenomePlotter
            indexer = QIndexerOCL(self.settings, self.logger, 0.1, dummySeq, int(self.settings.qgram), block, indexStepSize, nAs='A')
            #plotter = GenomePlotter(indexer, dummySeq, block, indexStepSize)                                                                                                                            
        elif self.qindexerCUDA:
            from pyPaSWAS.Core.QIndexerCUDA import QIndexerCUDA, GenomePlotter
            block = 10
            indexStepSize=1000
            indexer = QIndexerCUDA(self.settings, self.logger, 0.1, dummySeq, int(self.settings.qgram), block, indexStepSize, nAs='A')
            #plotter = GenomePlotter(indexer, dummySeq, block, indexStepSize)                                                                                                                            
        else:
            from pyPaSWAS.Core.QIndexer import QIndexer
            indexer = QIndexer(self.settings, self.logger, 0.1, dummySeq, int(self.settings.qgram))

        self.logger.debug("Starting indexer on targets")
        while indexer.indicesToProcessLeft():
            indexer.createIndexAndStore(targets, self.arguments[1])
            plotter = GenomePlotter(indexer, dummySeq, block, indexStepSize, nAs='C')
            self.logger.debug("Starting indexer on sequences")
            while plotter.indicesToProcessLeft():
                plotter.createIndexAndStore(records_seqs, self.arguments[0])
    
                allLocations = plotter.findDistances(indexer)

                for a in xrange(len(allLocations)):
                    locations = allLocations[a]
                    locs = []
                    if (len(locations) > 0):
                        for value in locations.itervalues():
                            locs.extend(value)
                        splittedTargets = []

                        for loc in locs:
                            #self.logger.debug("Loc: {}, len: {}".format(loc, len(records_seqs)))                                                                                                        
                            distance = loc[2]
                            targetID = targets[loc[1][1]].id
                            targetStartIndex = loc[1][0]
                            refID = loc[1][1]
                            seqRefID = loc[0][1]
                            seqID = records_seqs[loc[0][1]].id
                            seqStartIndex = loc[0][0]
                            t = SWSeqRecord(Seq("", seq.seq.alphabet), targetID, targetStartIndex, original_length = 0, distance = distance, refID = targetID)
                            s = SWSeqRecord(Seq("", seq.seq.alphabet), seqID, seqStartIndex, original_length = 0, distance = distance, refID = seqID)
                            self.hitlist.add(Distance(self.logger, s, t, (seqStartIndex,seqStartIndex+window), (targetStartIndex, targetStartIndex+window)))
                            #self.logger.debug("Target {}, read {}, distance {}".format(t, s, distance))                                                                                                 
                formatter.print_results(self.hitlist)

                self.hitlist = HitList(self.logger)

        if indexer != None and self.qindexerCUDA:
            indexer.pop_context()
        
        formatter.close()
        self.logger.debug('Genome plotter finished.')
        return self.hitlist

        
