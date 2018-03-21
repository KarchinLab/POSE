'''
Get python modules
'''
import platform, os, re, sys, time, urllib, urllib2
import xml.etree.ElementTree as ET
import cPickle
from random import sample

'''
Get third-party modules
'''

'''
Get POSE modules
'''

def get_user_agent():
    '''
    Makes requests look like they are coming from a person/browser. 
    '''
    
    ClientRevision = '$Revision: 2673 $'
    ClientVersion = '0'
    if len(ClientRevision) > 11:
        ClientVersion = ClientRevision[11:-2]
    UserAgent = 'EBI-Sample-Client/%s (%s; Python %s; %s) %s' % (
        ClientVersion, os.path.basename( __file__ ),
        platform.python_version(), platform.system(),
        'Python-urllib/%s' % urllib2.__version__
    )
    
    return UserAgent

def rest_request(URL):
    '''
    Make requests until you have a response
    '''
    
    try:
        # Set the User-agent.
        UserAgent = get_user_agent()
        Request = urllib2.Request(URL, None, { 'User-Agent' : UserAgent })
        # Make the request (HTTP GET).
        RequestHTTP = urllib2.urlopen(Request)
        Result = RequestHTTP.read()
        RequestHTTP.close()
	# Errors are indicated by HTTP status codes.
    except urllib2.HTTPError, ex:
        # Trap exception and output the document to get error message.
        print >>sys.stderr, ex.read()
        raise

    return Result

def get_sequence(Arguments, Organism="human"):
    '''
    If (e)POSE is enlisted to create the initial sequence pool (i.e., multiple sequence alignment 
    (MSA)), then this function will retrieve the sequence for the target gene. The indentifier for 
    the target gene is passed by the GeneIdentifier argument, and must be either the hugo symbol
    (e.g., PIK3CA) or UniProt ID (e.g., P42336).
    '''

    BaseUniprotURL = "http://www.uniprot.org/uniprot/"

    try:
        URL = BaseUniprotURL + Arguments.GeneIdentifier + ".fasta"
        Request = urllib2.Request(URL)
        Response = urllib2.urlopen(Request)
        FASTA = Response.read()
    except urllib2.HTTPError:
        URL = BaseUniprotURL + "?query=reviewed:yes+gene:" + Arguments.GeneIdentifier + "+organism:" + Organism + "&format=fasta"
        Request = urllib2.Request(URL)
        Response = urllib2.urlopen(Request)
        FASTA = Response.read()

    if not len(FASTA):
        print Arguments.GeneIdentifier, "is not a valid protein/gene identifier. Please set the 'GeneIdentifier' argument to a valid ",
        print "hugo gene symbol or UniProt identifier and rerun"
        exit()
    else:
        ID = FASTA.replace(" ", "$").split()[0].replace("$", " ")
        Sequence = "".join(FASTA.replace(" ", "$").split()[1:])
    
    return ID.split("|")[1], Sequence

def blast_sequence(Arguments, Database="uniprotkb", Program="blastp", SequenceType="protein", \
                   ExpectThreshold=100, MaxTargetSequences=1000):
    '''
    If (e)POSE is enlisted to create the initial sequence pool (i.e., multiple sequence alignment 
    (MSA)), then this function will do a blaspt search to retrieve putative homologs for your 
    target protein sequence. 
    '''

    #Have to provide an email adress for BLAST searches
    if not Arguments.Email or ("@" and "." not in Arguments.Email):
        print "If you're not providing your own MSA, then you need to supply a valid email address so that we can BLASTp",
        print "your sequence. Please provide your own MSA or rerun setting the 'Email' argument with a valid email address"
        print "exiting..."
        exit()


    BaseUniprotURL = "http://www.uniprot.org/uniprot/"
    BaseNCBI_URL = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast'

    #Arguments used to run BLASTp
    BlastArguments = {}

    print "Fetching sequence for", Arguments.GeneIdentifier
    TargetID, TargetSequence = get_sequence(Arguments)
    print "Blasting sequence..."

    BlastArguments['sequence'] = TargetSequence

    #If we want BLASTp to focus on a user-defined domain
    if Arguments.Domain[0]:
        UpStream = int(Arguments.Domain[0]) - 1
        DownStream = int(Arguments.Domain[1])
        BlastArguments['sequence'] = TargetSequence[UpStream:DownStream]

    #Fill in the rest of the blast arguments. For the time being, you have to use MyPOSE to have control over these. 
    BlastArguments['email'] = Arguments.Email
    BlastArguments['database'] = Database
    BlastArguments['program'] = Program
    BlastArguments['stype'] = SequenceType
    BlastArguments['exp'] = ExpectThreshold
    BlastArguments['alignments'] = MaxTargetSequences

    #Submit the blast job...
    RequestURL = BaseNCBI_URL + '/run/'
    RequestData = urllib.urlencode(BlastArguments)
    UserAgent = get_user_agent()
    Request = urllib2.Request(RequestURL, None, { 'User-Agent' : UserAgent })
    RequestHTTP = urllib2.urlopen(Request, RequestData)
    JobID = RequestHTTP.read()
    RequestHTTP.close()

    #...wait until it completes
    Status = 'IN_PROGRESS'
    while Status != 'FINISHED':
        RequestURL = BaseNCBI_URL + '/status/' + JobID
        Status = rest_request(RequestURL)

    #Fetch XML output from job
    RequestURL = BaseNCBI_URL + '/result/' + JobID + '/xml'
    Request = urllib2.Request(RequestURL)
    XML_Response = urllib2.urlopen(Request)
    XML_IDs = ET.fromstring(str(XML_Response.read()))
    
    #Parse the XML output to IDs for each of the blast-retrieved sequences
    Description = {} 
    for e in XML_IDs.find("{http://www.ebi.ac.uk/schema}SequenceSimilaritySearchResult/{http://www.ebi.ac.uk/schema}hits"):
        if e.get("ac") != TargetID: #we'll add the target to the top of the list (below)
            Description[e.get("ac")] = e.get("description")

    UniProtIDs = Description.keys()
    UniProtIDs.insert(0, TargetID)
    Description[TargetID] = "Target"
    
    return UniProtIDs, Description

def make_alignment(Arguments):
    '''
    If (e)POSE is enlisted to create the initial sequence pool (i.e., multiple sequence alignment 
    (MSA)), then this function will align sequences retrieved from a blastp search. 
    '''

    UniProtIDs, Description = blast_sequence(Arguments)
    Target = UniProtIDs[0]
    UniProtIDs = UniProtIDs[1:]
    UniProtIDs = sample(UniProtIDs, min(UniProtIDs, 99))
    UniProtIDs.insert(0, Target)
    BaseAlignURL = 'http://www.uniprot.org/align/' 
    AlignmentArguments = {}
    AlignmentArguments["query"] = " ".join(UniProtIDs)
    print "Sending sequences to UniProt for alignment..."
    Success = False
    TimeKeeper = 0
    while not Success:
        if TimeKeeper == 3:
            UniProtIDs = sample(UniProtIDs, min(UniProtIDs, 99))
            UniProtIDs.insert(0, Target)
            TimeKeeper = 0
            print "Trying different sequences"
            
        try:
            Data = urllib.urlencode(AlignmentArguments)
            Request = urllib2.Request(BaseAlignURL, Data)
            Request.add_header('User-Agent', 'Python %s' % Arguments.Email)
            Response = urllib2.urlopen(Request)
            Success = True
            print "UniProt successfully contacted..."
        except:
            TimeKeeper +=1
            print "Failed to contact UniProt, trying again..."
            
    JobURL = Response.geturl()
    Response.close()
    
    Status = 'NOT_COMPLETED'
    while (Status != 'COMPLETED'):
        RequestURL = JobURL + '.stat'
        Status = rest_request(RequestURL)

    RequestURL = JobURL +'.fasta'
    Request = urllib2.Request(RequestURL)
    Alignment = urllib2.urlopen(Request)

    Genes = [] #We need a list to preserve order (dict won't!)
    MSA = {}
    for Line in Alignment.readlines():
        if ">" in Line:
            Gene = Line.strip(),  Description[Line.replace(">", "").strip()]
            Genes.append(Gene)
            MSA[Gene] = []
        else:
            MSA[Gene].extend(Line.strip())
    
    #return the target gene ID, and MSA
    try:
        return Genes[0], MSA
    except IndexError:
        print "Sorry...UniProt returned an empty list, which is not all that uncommon. Please try again."
        exit()

def write_reduced_alignment(Arguments):
    '''
    Only consider alignment columns (positions) for which the "target" sequence has a residue (i.e., not "-").
    '''

    TargetGene, MSA = make_alignment(Arguments)
    ReducedAlignment = []

    for Position in zip(*MSA.values()):
        #Find position of target gene, and make sure that has a residue
        if Position[MSA.keys().index(TargetGene)] != "-":
            ReducedAlignment.append(Position)

    ReducedAlignment = zip(*ReducedAlignment)

    Sequences = {}
    Genes = MSA.keys()
    for Gene in Genes:
        Sequences[Gene] = "".join(ReducedAlignment[Genes.index(Gene)])

    cPickle.dump([TargetGene, Sequences], open(Arguments.Filename + ".fasta", "wb"), -1)

def pre_process(Arguments):
    '''
    So far just making the alignment...change script name if that remains the case
    '''
    
    if Arguments.PreProcess == "MakeAlignment":
        write_reduced_alignment(Arguments)
        
    return 
