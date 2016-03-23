from rdflib.namespace import OWL, RDF, RDFS, XSD
from rdflib import Graph, Literal, Namespace, URIRef
from urllib import quote,unquote
from progress.bar import Bar
from biocma import utils,cma

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

def build_empty_ontology():
	#new namespace definition
        MSA = Namespace("http://localhost/msaont#")

	#instantiate empty graph
	graph = Graph()
	uris = {}

	#namespace binding
	graph.bind("rdf", RDF)
	graph.bind("rdfs", RDFS)
	graph.bind("msaont", MSA)
	
	#MSA class
	msa = URIRef(MSA["msa"])
	uris['msa'] = msa
	graph.add((msa, RDF.type, RDFS.Class))
	graph.add((msa, MSA.id, Literal("Unique identifier for MSA.")))
	
	#Sequence class
	seq = URIRef(MSA["sequence"])
	uris['sequence'] = seq
	graph.add((seq, RDF.type, RDFS.Class))
	graph.add((seq, MSA.id, Literal("Unique identifier for Sequence.")))

	#Segment class
	segment = URIRef(MSA["segment"])
	uris['segment'] = segment
	graph.add((segment, RDF.type, RDFS.Class))

	#subclasses
	ar = URIRef(MSA["aligned_residue"])
	uris['aligned_residue'] = ar
	graph.add((ar, RDFS.subClassOf, segment))
	graph.add((ar, MSA.aln_pos, Literal("Position of residue in profile.")))
	graph.add((ar, MSA.native_pos, Literal("Position of residue in native sequence.")))
	graph.add((ar, MSA.native_residue, Literal("Residue type in native sequence.")))

	insert = URIRef(MSA["insertion"])
	uris['insertion'] = insert
	graph.add((insert, RDFS.subClassOf, segment))
	graph.add((insert, MSA.insert_after_aln_pos, Literal("Insertion occurs after this position in the profile.")))
	graph.add((insert, MSA.native_pos, Literal("Positions (start, end) of insertion in native sequence.")))
	graph.add((insert, MSA.native_residues, Literal("Residues inserted.")))

	delete = URIRef(MSA["deletion"])
	uris['deletion'] = delete
	graph.add((delete, RDFS.subClassOf, segment))
	graph.add((delete, MSA.deleted_aln_pos, Literal("Profile position not represented in sequence.")))

	#has_sequence
	graph.add((MSA.has_sequence, RDF.type, RDF.Property))
	graph.add((MSA.has_sequence, RDFS.domain, msa))
	graph.add((MSA.has_sequence, RDFS.range, seq))
	
	#has_segment
	graph.add((MSA.has_segment, RDF.type, RDF.Property))
	graph.add((MSA.has_segment, RDFS.domain, seq))
	graph.add((MSA.has_segment, RDFS.range, segment))	
	
	return uris, graph

if __name__ == "__main__":
        name = 'pdb_aln'
	MSA = Namespace("http://localhost/msaont#")
	#uris,graph = build_empty_ontology() 
        graph = Graph()
	graph.bind("rdf", RDF)
	graph.bind("rdfs", RDFS)
	graph.bind("msaont", MSA)
	#ifile = 'totnrtxTK_with_consensus.cma'
	#ifile = 'prokino-dedupe.cma'
	#ifile = 'temp.cma'
        #ifile = 'seqs-pressed-cons.cma'
        ifile = '150806_nrtx.unique_is119_is10.cma'
        outfile = 'msaont_nr.rdf'

	all_aln = cma.read(ifile)
	dedup_aln = clean_records(all_aln)
	dedup_eqv = utils.get_equivalent_positions(dedup_aln)
	#MSA instance
	curi = URIRef(MSA[name])	
	graph.add((curi, RDF.type, MSA.msa))
	graph.add((curi, MSA.id, Literal(name)))

	bar = Bar('Scanning sequences', max=len(dedup_aln['sequences']))

        #limited number of records for testing
	for rec in dedup_aln['sequences'][1:6]:
                print rec['id']
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
					#graph.add((ruri, MSA.native_pos, Literal(dedup_eqv[i+1][acc])))
				graph.add((ruri, MSA.native_residue, Literal(r)))
			#need to implement insertions
		bar.next()
						
	bar.finish()	

	quer = '''SELECT 
			?seq ?res ?pos
			WHERE { 
				?seq rdf:type msaont:sequence; 
					msaont:has_segment ?seg . 
				?seg rdf:type msaont:aligned_residue; 
					msaont:aln_pos ?pos; 
					msaont:native_residue ?res . 
				filter(?res = 'R')
				filter(?pos = 59)
			}'''

	quer2 = '''SELECT 
				?pos ?res (Count(*) as ?Count)
				WHERE {
					?seg rdf:type msaont:aligned_residue;
						msaont:aln_pos ?pos;
						msaont:native_residue ?res .
				}
				GROUP BY ?pos ?res'''

        graph.serialize(destination=outfile, format='pretty-xml')
	qres = graph.query(quer)


