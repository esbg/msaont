from rdflib.namespace import OWL, RDF, RDFS, XSD
from rdflib import Graph, Literal, Namespace, URIRef
from urllib import quote,unquote
from progress.bar import Bar
from biocma import cma, utils
from Bio import SeqIO

def clean_records(block):
    """Removes sequences with duplicate ids. """
    seen = set()
    new_seqs = []
    for i,rec in enumerate(block['sequences']):
        if rec['id'] not in seen:
            seen.add(rec['id'])
            new_seqs.append(rec)

    return {'level':block['level'],
             'one': block['one'],
             'name': block['name'],
             'params': block['params'],
             'query_length': block['query_length'],
             'query_chars': block['query_chars'],
             'sequences': new_seqs}

def manipulate_fasta(ifile, p):
    seqs = SeqIO.parse('temp.fasta','fasta')
    slen = len(seqs.next().seq)
    tempseq = cma.seqrecord2sequence(seqs.next(),0,1)
    tempseq['id'] = 'filler'
    temp = [[] for x in range(slen)]
    seqs = SeqIO.parse('temp.fasta','fasta')
    for seq in seqs:
        for j in range(len(seq.seq)):
            temp[j].append(seq.seq[j])
    caps = []
    tseq = ''
    for x in temp:
        if float(x.count('-'))/len(x) > p:
            caps.append('l')
        else:
            caps.append('u')
            tseq += 'A'
    tempseq['seq'] = tseq
    new_seqs,i = [tempseq],0
    seqs = SeqIO.parse('temp.fasta','fasta')
    for seq in seqs:
        i += 1
        nseq = ''
        cmaseq = cma.seqrecord2sequence(seq,0,i)
        for ind in range(len(cmaseq['seq'])):
            if caps[ind] == 'l': 
                if cmaseq['seq'][ind] != '-':
                    nseq += cmaseq['seq'][ind].lower()
            else:
                nseq += cmaseq['seq'][ind].upper()
        cmaseq['seq'] = nseq
        cmaseq['head_len'] = 0
        new_seqs.append(cmaseq)
    return {'level':0,
             'one': 1,
             'name': ifile,
             'params': '',
             'query_length': caps.count('u'),
             'query_chars': ''.join(['*' for x in range(caps.count('u'))]),
             'sequences': new_seqs}
    
def populate(sequences, name, outfile, namespace):
    MSA = Namespace(namespace)
    graph = Graph()
    graph.bind("rdf", RDF)
    graph.bind("rdfs", RDFS)
    graph.bind("msaont", MSA)
    
    #add aligned fasta support
    #all_aln = cma.read(ifile)
    dedup_aln = clean_records(sequences)
    dedup_eqv = utils.get_equivalent_positions(dedup_aln)

    #MSA instance
    curi = URIRef(MSA[name])	
    graph.add((curi, RDF.type, MSA.msa))
    graph.add((curi, MSA.id, Literal(name)))

    bar = Bar('Scanning sequences', max=len(dedup_aln['sequences']))

    for rec in dedup_aln['sequences'][1:]:
        #assume sequence uri already exists
        tsplit = rec['id'].split('_')
        sequri = MSA['_'.join(tsplit[1:])]
        seq = ''.join((c for c in rec['seq'] if not c.islower()))
        og_seq = rec['seq']
        acc = quote(rec['id'])
        #sequence instance
        sequri = URIRef(MSA[acc])
        graph.add((sequri, RDF.type, MSA.sequence))
        graph.add((sequri, MSA.id, Literal(acc)))
        for i,r in enumerate(seq,start=1):
            ruri = URIRef(MSA[acc+str(i)])
            graph.add((sequri, MSA.has_segment, ruri))
            if og_seq[0] == r:
                og_seq = og_seq[1:]
            else:
                j = 0
                while og_seq[j] != r:
                    j += 1
                insert = og_seq[:j]
                og_seq = og_seq[j+1:]
                rurii = MSA[acc+'_i_'+str(i)]
                graph.add((rurii, RDF.type, MSA.Insertion))
                graph.add((rurii, MSA.hasAlignedPosition, Literal(i)))
                graph.add((rurii, MSA.hasNativeResidue, Literal(insert.upper())))
                graph.add((rurii, MSA.hasLength, Literal(len(insert))))
                if rec['id'] in dedup_eqv[i]:
                    graph.add((rurii, MSA.hasNativePosition, Literal(dedup_eqv[i][rec['id']]+1)))
            if r == '-':
                #add deletion node
                graph.add((ruri, RDF.type, MSA.deletion))
                graph.add((ruri, MSA.deleted_aln_pos, Literal(i)))
            else:
                graph.add((ruri, RDF.type, MSA.aligned_residue))
                graph.add((ruri, MSA.aln_pos, Literal(i)))
                if unquote(acc) in dedup_eqv[i]:
                    graph.add((ruri, MSA.native_pos, Literal(dedup_eqv[i][unquote(acc)])))
                else:
                    print "shouldn't happen" #deletions taken care of above
                graph.add((ruri, MSA.native_residue, Literal(r)))
        bar.next()
    bar.finish()	
    graph.serialize(destination=outfile, format='pretty-xml')
    return

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description = 'Populate the Multiple Sequence Alignment Ontology')
    parser.add_argument('infile', metavar='infile', type=str, help='input aligned sequences')
    parser.add_argument('name', metavar='name', type=str, help='Graph name')
    parser.add_argument('--namespace', metavar='namespace', type=str, default='http://localhost/msaont#', help='Graph namespace')
    parser.add_argument('-p', dest='prop', action='store', default=0.25, help='proportion of inserts allowed in an aligned residue')
    parser.add_argument('-o', dest='outfile', action='store', default='out.rdf', help='outfile')
    args = parser.parse_args()
    ext = args.infile.split('.')[-1]
    if ext != 'cma':
        seqs = manipulate_fasta(args.infile, args.prop)
    else:
        seqs = cma.read(args.infile)
    populate(seqs, args.name, args.outfile, args.namespace)
