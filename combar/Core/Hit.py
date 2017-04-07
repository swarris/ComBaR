from pyPaSWAS.Core.Exceptions import CudaException
from Bio.Seq import Seq
from pyPaSWAS.Core.SWSeqRecord import SWSeqRecord
from pyPaSWAS.Core.Hit import Hit

class Distance(Hit):
    def __init__(self, logger, sequence_info, target_info, sequence_location, target_location):
       Hit.__init__(self, logger, sequence_info, target_info, sequence_location, target_location)

    @staticmethod
    def _is_a_valid_location(location, sequence_length):
        ''' Verifies wether or not a location is valid.
        '''
        return isinstance(location, tuple) and len(location) == 2 and location[0] <= location[1] 
    def get_sam_line(self):
        '''Creates and returns an alignment line as described in http://samtools.sourceforge.net/SAM1.pdf
           The following mappings are used for record lines:
                QNAME:  self.get_seq_id()
                FLAG:   only bit 0x2 and 0x10 are set since the other values are either false or unknown
                RNAME:  self.get_target_id()
                POS:    hit.target_location[0] + 1 (+1 since SAM is 1 based!)
                MAPQ:   set to 255 to indicate quality score is not available
                CIGAR:  see _get_sam_cigar
                RNEXT:  set to * (not available)
                PNEXT:  set to 0 (not available)
                TLEN:   set to 0 (not available)
                SEQ:    self.sequence_info.sequence
                QUAL:   set to * (not available)
                AS:     self.score
                AD:f:   computed euclidian distance between sequence and target
                RS:f:   computed relative score
        '''
        rc_seq = False
        rc_target = False

        
        identifier = self.get_seq_id()
        if identifier[-2:] == 'RC':
            identifier = identifier[:-3]
            rc_seq = True
        target_id = self.get_target_id()
        if target_id[-2:] == 'RC':
            target_id = target_id[:-3]
            rc_target = True
            
        #rc = (rc_seq and not rc_target) or (not rc_seq and rc_target)
        self.rc = rc_target or rc_seq 
        
        hit_pos = str(self.target_location[0] + 1)
        sam_sequence = self.sequence_info.seq
        if self.rc:
            self.alignment = self.alignment[::-1]
            self.sequence_match = self.sequence_match[::-1]
            sam_sequence =  Seq(str(self.sequence_info.seq), self.sequence_info.seq.alphabet).reverse_complement()
            hit_pos = str(self.target_info.original_length - self.target_location[1])
        if len(self.sequence_match) == 0 or len(self.target_match) == 0 or len(self.alignment) == 0:
            raise CudaException('sequence_match, target_match and alignment should be set and have lengths > 0.')

        return '\t'.join([identifier, self._get_sam_flag(), target_id,
                          hit_pos, Hit._get_sam_mapq(), self._get_sam_cigar(),
                          '*', str(0), str(0), str(sam_sequence), '*', self._get_sam_alignment_score(),
                          self._get_sam_euclidian_distance(),
                          self._get_sam_relative_score(),
                          self._get_sam_relative_base_score(),
                          self._get_sam_query_coverage(),
                          self._get_sam_query_identity()])

    def _get_sam_euclidian_distance(self):
        '''returns the computed distance as an optional sam field
        '''
        return "AD:f:" + str(self.get_euclidian_distance())

    def get_euclidian_distance(self):
        '''Returns the calculated Euclidian distance between sequence and target'''
        if self.target_info.distance is None:
            return 0.0
        else:
            return self.target_info.distance

    def keys(self):
        '''Returns a list of tuples that will be used as keys for this hit in a hitlist.'''
        if self.get_seq_id is not None and self.get_target_id is not None:
            return [(self.get_seq_id(), self.get_target_id()),
                    (self.get_seq_id(), self.get_target_id()),
                    (self.get_seq_id(), self.get_target_id())]
        else:
            return None


class QGramLink(Distance):
    def __init__(self, logger, sequence_info, target_info, sequence_location, target_location):
       Distance.__init__(self, logger, sequence_info, target_info, sequence_location, target_location)

    def keys(self):
        '''Returns a list of tuples that will be used as keys for this hit in a hitlist.'''
        if self.get_seq_id is not None and self.get_target_id is not None:
            return [(self.get_seq_id(), self.get_target_id(), self.seq_location[0], self.target_location[0]),
                    (self.get_seq_id(), self.get_target_id(), self.seq_location[1], self.target_location[1]),
                    (self.get_seq_id(), self.get_target_id(), self.seq_location[0], self.target_location[0])]
        else:
            return None
