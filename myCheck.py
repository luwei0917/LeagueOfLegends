#!/usr/bin/env python
# Score is:
#   % gene level right or - gene's num  txpts
# + % txpt level right
# + % uniq calls right
# + % definitive rel calls right
import os, re, sys
from os.path  import basename, exists, getsize
from operator import itemgetter
from sys      import stderr
import argparse

class G:
  """Generic container for shared variables."""
  args = None  # set by getParameters()
  gene_transcript_counts = None
  notify_quiet = False
  notify_logFH = None
  def __init__(self):
    __name__ = "G"
    self.notify_logFH = None
pass

## GLOBALS ##
global g
g = G()

def main():  ##########################################################
  args = getParameters()
  notify('#' * 47 + '\n' )
  notify('# Starting evaluation of assignments\n' )
  #'A|B|2:At|Bt|Call|5:Score|Ag|Bg|8:Category|...'
  refD, source_ids = loadValidationFile( args.known )
  loadGFFLoci()  # gets counts of child transcripts for each gene from transcript fasta and gff3
  counts = scoreSolutions( refD, source_ids )
  scores = []
  metrics = [
    ('gene',       'obs_gene_ct', 'max_gene_ct', ' ' * 6 ),
    ('transcript', 'obs_txpt_ct', 'max_txpt_ct', ' ' * 0 ),
    ('unique',     'obs_uniq_ct', 'max_uniq_ct', ' ' * 4 ),
    ('definitive', 'obs_defn_ct', 'definite_ct', ' ' * 0 )
  ]
  notify('#' * 47 + '\n' )
  notify('# Results:\n')
  for A in list( counts.keys() ):
    for B in list( counts[ A ].keys() ):
      for metric, obs_tag, tot_tag, pad in metrics:
        try:
          obs, tot = counts[A][B][ obs_tag ], counts[A][B][ tot_tag ]
        except KeyError as e:
          notify( '# %s metrics were unset for %s to %s\n'%( metric, A,B) )
          continue
#         notify('Obs:%d Tot:%d\n'%(obs, tot) )
        pct = 100.0 * abs( obs ) / tot
        if obs < 0: pct = pct * -1.0
        scores.append( pct )
        notify('# %s to %s %s call accuracy:%s\t% 3.2f%%\n'%(A, B, metric, pad, pct))
  notify('# ' + '-' * 45 + '\n' )
  notify('# Final Score:                       \t% 3.2f\n'%( sum( scores ) ) )
  notify('#' * 47 + '\n' )
# end main()  ######################################################################

def commafy( number ):  ############################################################
  """commafy(number) - add commas to a number"""
  import re
  return re.sub(r'(\d{3})(?=\d)',r'\1,' ,str(number)[::-1])[::-1]
# end commafy()

def getParameters():
  epilog = 'Example:\n' +\
           '  %(prog)s -e 3 A.transcripts.fasta B.transcripts.fasta A.gff3 B.gff3 dev_validation_set.tsv my_solutions.tsv'
  parser = argparse.ArgumentParser(prog='check_min_accuracy.py',
                                   description='Checks a proposed solution file against a check set.\n' +\
                                     'Remember that identifiers should remain exactly as supplied.',
                                   formatter_class=argparse.RawTextHelpFormatter, epilog=epilog )
  parser.add_argument( '-e', '--errors-to-report', metavar='#', default=3, dest='max_error_ct',
                       help='Maximum number of assignment errors to report to STDERR')
  parser.add_argument( 'txpt', metavar="TRANSCRIPT-FASTA", nargs=2,
                       help="Supply space-separated paths to each transcript FASTA file for the genome pair.")
  parser.add_argument( 'gff', metavar="GFF3", nargs=2,
                       help="Supply space-separated paths to each GFF3 file for the genome pair")
  parser.add_argument( 'known', metavar="REFERENCE_FILE",
                       help="Path to the supplied reference validation file used to assess accuracy")
  parser.add_argument( 'test', metavar="SOLUTION_FILE",
                       help='Path to the solution file to assess for accuracy in the format required for the proposal:\n' +\
                            'SourceA|SourceB|A_transcript|[B_transcript]|Call|Score|A_gene|[B_gene]\n'+\
                            'The pipes represent TAB characters\n')

  args = parser.parse_args()
  for f in args.txpt:
    if not exists( f ): sys.exit('ERROR: Failed to find transcript file %s.\n'%( f ) )
    if getsize( f )< 1: sys.exit('ERROR: Transcript file %s is empty.\n'%( f) )
  for f in args.gff:
    if not exists( f ): sys.exit('ERROR: Failed to find GFF3 file %s.\n'%( f ) )
    if getsize( f )< 1: sys.exit('ERROR: GFF3 file %s is empty.\n'%( f) )
  if not exists( args.known ): sys.exit('ERROR: Failed to find the validation file.  You supplied:\n%s\n'%( args.known))
  if getsize( args.known )< 1: sys.exit('ERROR: Supplied validation file appears empty.  You supplied:\n%s\n'%( args.known))
  if not exists( args.test ):  sys.exit('ERROR: Failed to find the solution file.  You supplied:\n%s\n'%( args.test))
  if getsize( args.test ) < 1: sys.exit('ERROR: Supplied solution file appears empty.  You supplied:\n%s\n'%( args.test))
  g.args = args
  return args
# end getParameters()

def loadGFFLoci():  ################################################################
  """ Get counts of transcripts per gene """
  txptNameD =  loadTranscriptSets() # source->set([txptIDs])
  features_to_track = set( ['mRNA', 'miRNA', 'lincRNA'] )
  g.gene_transcript_counts = {}
  gffPat = re.compile( r'Name=((?:transcript|gene)([A-Z])\d+[^;]*);Parent=(gene([A-Z])\d+)' )
  for gff3_file in g.args.gff:
    notify('# Examining feature relationships in %s\n'%( basename( gff3_file ) ) )
    with open( gff3_file, 'r' ) as iFH:
      for line in iFH:
        if line.startswith('#'): continue
        r = line.strip().split( '\t' )
        if len( r ) < 9: continue
        if r[2].strip() in features_to_track:
          m = gffPat.search( r[8] )
          if not m:
            sys.exit('# ERROR: Unhandled GFF3 match type:\n# %s'%(line) )
          t_source = m.group( 2 )
          g_source = m.group( 4 )
          if t_source != g_source:
            sys.exit('# ERROR: Child and parent should never have different sources.\n# %s'%(line) )
          txpt_id, gene_id = m.group( 1 ), m.group( 3 )
          if t_source not in g.gene_transcript_counts: g.gene_transcript_counts[ t_source ] = {}
          if gene_id not in g.gene_transcript_counts[ t_source ]:
            g.gene_transcript_counts[ t_source ][ gene_id ] = 0
          g.gene_transcript_counts[ t_source ][ gene_id ] += 1
          if txpt_id in txptNameD[ t_source ]: txptNameD[ t_source ].remove( txpt_id )
  for g_source in list( g.gene_transcript_counts.keys() ): # GFF known transcripts have been removed
    unplaced_ct = len( txptNameD[ g_source ] )
    if unplaced_ct > 0:
      notify("# Source %s had %s unplaced transcripts, treating them as single genes\n"%(g_source,
                                                                             commafy(unplaced_ct) ) )
    else: notify('# Source %s had no unplaced transcripts\n'%( g_source) )
    for txpt_id in txptNameD[ g_source ]:
      g.gene_transcript_counts[ g_source ][ txpt_id + '_gene' ] = 1
  # now we have enough info to penalize incorrect gene calls
# end loadGFFLoci()  ################################################################

def loadTranscriptSets():  ##########################################################
  """ Get full listing of transcripts which may vary from GFF3 contents """
  namePat = re.compile( r'\>(transcript([A-Z])\d+)' )
  txptNameD = {}  # source->set([txptIDs])
  for txpt_file in g.args.txpt:
    notify('# Getting a list of transcripts from %s\n'%( basename( txpt_file ) ) )
    with open( txpt_file, 'r' ) as iFH:
      for line in iFH:
        m = namePat.search( line )
        if m:
          txpt_id = m.group( 1 )
          source = m.group( 2 )
          if source not in txptNameD: txptNameD[ source ] = set([])
          txptNameD[ source ].add( txpt_id )
  return txptNameD
# end loadTranscriptSets()  #########################################################

def loadValidationFile( known_file ):  ##############################################
  """ this version stores:
    refA->refB->ATxpt = 'c':call, 't':set(Btxpts), 'g':set(Bgenes)
  """
  refD = {}
  source_ids = set( [] )
  uniques = set([ 'absent_gene', 'absent_genome', 'absent_transcript', 'gene_fusion',
                  'unique_transcript' ] )
  dupe_call_check = set([])
  with open( known_file, 'r' ) as iFH:
    for line in iFH:
      if line.startswith( '#' ): continue
      r = line.strip( '\n' ).split( '\t' )
      if len( r ) < 1: continue
      if len( r ) < 3:
        if len( line.strip() ) == 0: continue
        sys.exit( 'Error: validation file row has less than the minimum of 3 values for A_Src|B_Src|A_ID\n"%s"'%line)
      A, B = r[0], r[1]
      if A not in refD:
        refD[ A ] = {}
        source_ids.add( A )
      if B not in refD[ A ]:
        refD[ A ][ B ] = { 'max_gene_ct':0, 'max_txpt_ct':0, 'max_uniq_ct':0, 'definite_ct':0 }
        source_ids.add( B )
      At = r[ 2 ].strip()
      if At not in refD[ A ][ B ]: refD[ A ][ B ][ At ] = { 'c':'','t':set([]),'g':set([]) }
      c_call = None
      if len( r ) < 8:
        notify( line )
        notify( '# %s\n'%( str( r ) ) )
        sys.exit('\n# ERROR: Unexpectedly short line in check file.  8 or more columns expected!\n')
      Bt, call, score = r[3].strip(), r[4].strip().lower(), float( r[5].strip() )
      call_tupe = (A, B, At, Bt, call, r[7].strip() )
      if call_tupe in dupe_call_check:
        notify("# WARNING: Already observed this mapping for %s:%s"%(At, line))
        continue
      else: dupe_call_check.add( call_tupe )
      refD[A][B]['max_txpt_ct'] += 1
      # manage the call for this row
      refD[A][B][At][ 'c' ] = call
      if call != 'no_call': refD[A][B][ 'definite_ct' ] += 1 # track number of possible definite calls
      if call == 'unique_transcript': refD[A][B]['max_uniq_ct'] += 1
      Ag = r[6].strip()  # do nothing with the A gene
      for Bg in r[7].strip().split( ';' ):
        refD[A][B][At]['g'].add( Bg )
        refD[ A ][ B ][ 'max_gene_ct' ] += 1
      if Bt: refD[A][B][At]['t'].add( Bt )
  return refD, source_ids
# end loading the validation file ######################################################

def scoreSolutions( refD, source_ids ):  ###############################################
# refD[ A ][ B ] = { 'max_gene_ct':0, , 'max_uniq_ct':0, 'definite_ct':0 }
#  if At not in refD[ A ][ B ]: refD[ A ][ B ][ At ] = { 'c':'','t':set([]),'g':set([]) }
#
  global g
  # begin loading the solution file
  notify("# Scoring relationships assigned in '%s' by validation set '%s'\n"%( basename( g.args.test ),
                                                                           basename( g.args.known )))
  counts            = {}
  reported_failures = 0
  unknown_gene_ct   = 0
  fails_to_report   = g.args.max_error_ct
  dupe_call_check   = set([])
  DBG = True
  obs_pairs = set([])
  absents = set( ['absent_gene', 'absent_transcript', 'absent_genome', 'gene_fusion'] )
  with open( g.args.test, 'r' ) as iFH:
    for line in iFH:
      if line.startswith( '#' ):     continue
      r = line.strip( '\n' ).split( '\t' )
      if len( r ) < 1 or r[0] == '': continue
      if len( r ) < 8:
        notify( '# ' + line )
        notify( '# %s\n'%( str( r ) ) )
        sys.exit('\n# ERROR: Unexpectedly short line in solution file.  8 or more columns expected!\n')
      A, B = r[0].strip(), r[1].strip()
      At, Bt, call, score = r[2].strip(), r[3].strip(), r[4].strip().lower(), float( r[5].strip() )
      Ag = r[6].strip()
      if A not in refD or B not in refD[ A ]:  # verify known and test file sources match
        msg = 'Error: We expect the source columns to ONLY contain "%s"\n'%('" or "'.join( list( source_ids ) ) )
        msg+= '       You supplied: "%s" and "%s"\n\n'%( A, B )
        sys.exit( msg )

      if ( A , B ) not in obs_pairs:  # init observation trackers
        obs_pairs.add( ( A , B ) )
        refD[A][B][ 'obs_gene_ct' ] = 0
        refD[A][B][ 'obs_txpt_ct' ] = 0
        refD[A][B][ 'obs_uniq_ct' ] = 0
        refD[A][B][ 'obs_defn_ct' ] = 0

      if At not in refD[A][B]:   continue  # only score things in the known set
      FAIL = False
      # check for repeated assertions and ignore them
      call_tupe = (A, B, At, Bt, call, r[7].strip() )
      if call_tupe in dupe_call_check:
        notify( '*' * 30 + '\n' )
        notify( "# ERROR: already observed this mapping for %s.v.%s:\n%s\n"%(A,B,str(call_tupe) ) )
        notify( '*' * 30 + '\n' )
        continue
      else:  dupe_call_check.add( call_tupe )
      CORRECT_CALL = True
      if call != refD[A][B][At]['c']:
        CORRECT_CALL = False
        FAIL = True
      # get B transcript score
      if CORRECT_CALL:
        if Bt:
          if Bt in refD[A][B][At]['t']:
            refD[A][B][ 'obs_txpt_ct' ] += 1
            if call == 'unique_transcript': refD[A][B][ 'obs_uniq_ct' ] += 1
          else: pass # unknown but check set may be incomplete so ignore
        elif call in absents:    refD[A][B][ 'obs_txpt_ct' ] += 1
      elif Bt:                   refD[A][B][ 'obs_txpt_ct' ] -= 1 # got it wrong

      # get B definitive call score
      if CORRECT_CALL and call != 'no_call': refD[A][B][ 'obs_defn_ct' ] += 1
      elif call == 'no_call': pass
      elif CORRECT_CALL == False: refD[A][B][ 'obs_defn_ct' ] -= 1

      # get B gene score
      gScore = 0
      for Bg in r[7].strip().split( ';' ):
        if Bg not in refD[A][B][At]['g']:
          if len( Bg ) < 1: continue
          if unknown_gene_ct < 3:
            notify('# WARNING: Unrecognized gene "%s", possibly for '%(Bg) +\
                                         'transcript missing from GFF. Ignoring\n')
            if unknown_gene_ct == 2: notify('# WARNING: Further unrecognizezd gene examples will not be reported\n')
          unknown_gene_ct += 1
          continue
        if Bg in refD[A][B][At]['g']: gScore += 1
        else: # uh oh, this is wrong penalize by A gene's number of transcripts
          gScore -= g.gene_transcript_counts[ A ][ Ag ]
      refD[A][B]['obs_gene_ct'] += gScore

      if FAIL and reported_failures < int(fails_to_report):
        if reported_failures == 0:
          notify('#'*47 +'\n')
          notify('# Sample failed assignments:\n' )
        eCall = refD[A][B][At]['c']
        eTxpts = '|'.join( refD[A][B][At]['t'] )
        eGenes = '|'.join( refD[A][B][At]['g'] )
        notify("# Expected call '%s' B txpt(s): '%s' B gene(s): '%s'\n"%( eCall, eTxpts, eGenes ) )
        notify("# Your data:%s\n"%( '|'.join(r) ) )
        reported_failures += 1
  return refD
# end scoreSolutions()  ############################################################################

def mixSort(list, key=None):  ######################################################################
  """ Sorts a list in place as humans expect
      This implementation based on Ned Batchelders sort from the web
  """

  numPat = re.compile(r'([0-9]+)')
  if key is None:
    ################################################################################################
    def tryint(s):
      """ Supports mixSort """
      try: return int(s)
      except: return s
    # end tryint()  ################################################################################

    def alphanum_key(s):
      """ Turns a string into a list of string and number chunks.
          "z23a" -> ["z", 23, "a"]
          Supports mixSort
      """
      if isinstance(s, basestring):
        return [ tryint(c) for c in numPat.split(s) ]
      else:
        return s
    # end alphanum_key   ###########################################################################

  else: # supplied a way to retrieve the part of the iterable item to use for sorting
    ''' from  http://stackoverflow.com/questions/6849047 '''
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda item:[ convert( c ) for c in numPat.split( key(item) ) ]
  list.sort(key=alphanum_key)
# end mixSort()  ###################################################################################

####################################################################################################
def notify(msg, suppress=None, logFH=None):
  """determines whether program feedback should be written back to the terminal
  param1: a text message to print if it is ok to do so
  param2: if set to True, notify will not write anything, nor will subsequent calls to notify
  """
  if suppress != None: G.notify_quiet = suppress
  if logFH is not None:
    if G.notify_logFH is None:
      G.notify_logFH = logFH
    elif logFH != G.notify_logFH:
      G.notify_logFH = logFH
  if msg == None: msg = ""
  if G.notify_quiet == False:
    if msg.find("\r") > -1:
      stderr.write( msg )
    else:
      stderr.write( msg )
  if G.notify_logFH != None:
    G.notify_logFH.write( msg )
  return 0
# end notify()  ####################################################################################

if __name__ == '__main__':
  if len(sys.argv) < 2:
    sys.argv.append("--help")
  main()