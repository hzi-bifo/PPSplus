import sys
from Bio import Entrez

#define email for entrez login
db           = "nuccore"
# Entrez.email = "some_email@somedomain.com"
batchSize    = 100
retmax       = 10**9

#load accessions from arguments
if len(sys.argv[1:]) > 1:
  accs = sys.argv[1:]
else: #load accesions from stdin
  accs = [ l.strip() for l in sys.stdin if l.strip() ]
#first get GI for query accesions
sys.stderr.write( "Fetching %s entries from GenBank: %s\n" % (len(accs), ", ".join(accs[:10])))
query  = " ".join(accs)
handle = Entrez.esearch( db=db,term=query,retmax=retmax )
giList = Entrez.read(handle)['IdList']
sys.stderr.write( "Found %s GI: %s\n" % (len(giList), ", ".join(giList[:10])))
#post NCBI query
search_handle     = Entrez.epost(db=db, id=",".join(giList))
search_results    = Entrez.read(search_handle)
webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"]
#fecth all results in batch of batchSize entries at once
for start in range( 0,len(giList),batchSize ):
  sys.stderr.write( " %9i" % (start+1,))
  #fetch entries in batch
  handle = Entrez.efetch(db=db, rettype="gb", retstart=start, retmax=batchSize, webenv=webenv, query_key=query_key)
  #print output to stdout
  sys.stdout.write(handle.read())
