"""
Tests of (some) pathway-tools functionality using nosetests.

The EcoCyc database is queried to run these tests. In most cases the
test passes provided the query returns a nonempty list of EcoCyc
frames or a nonempty string. In the few cases where a specific value
is required to pass the test, rather than just the successful return
of something that looks vaguely appropriate, changes in the database
will cause tests to fail even when the code is performing correctly,
but short of installing a purpose-built test PGDB I am not sure how to
avoid this.

"""

import pycyc
from nose.tools import raises, assert_is_instance, assert_raises
from pycyc.frame import Frame

# Specify the values which should be returned for true and false queries
# (which may change in the future)
_NIL = None
_T = True

ec = pycyc.open('ecoli')

sample_reaction = 'trypsyn-rxn' # '1.8.4.13-RXN', 'oropribtrans-rxn'
pathway_info = ('TRPSYN-PWY','RXN0-2381','OROPRIBTRANS-RXN')
(sample_pathway, pathway_reaction, pathway_non_reaction) = pathway_info
sample_compound = 'trp' 
sample_gene = 'EG11024' # should encode at least one enzyme in at least one pwy
                        # should lie in at least one transcription unit 
                        # and have at least one regulator
sample_enzyme = 'TRYPSYN' # should be of types chemical-change, small-molecule,
                          # non-transport
sample_protein = 'TRYPSYN'
sample_transport_reaction='TRANS-RXN-157'
# The sample transport protein should NOT be an enzyme of type non-transport
sample_transport_protein='CPLX-165' 
sample_tf = 'PHOSPHO-CREB'
sample_regulatory_gene = 'EG11218'
sample_non_tf = 'TRYPSYN'
sample_multicomponent_protein = 'TRYPSYN'
sample_frame_label = 'pyrophosphate'
sample_protein_with_nontrivial_container = 'TRYPSYN-APROTEIN'
sample_protein_with_modified_forms = 'CREB-MONOMER'
sample_tu = 'TU00038'
sample_tu_with_mrna_site = 'TU00355'
sample_inhibited_tu = 'TU0-1123'
sample_activated_tu = 'TU0-13935'
multifn_info = ('RIB5PISOMB-CPLX','RXN0-303','RIB5PISOM-RXN')
sample_multienzyme, multienzyme_r1, multienzyme_r2 = multifn_info
sample_label = 'pyrophosphate'
sample_binding_site = 'BS00117'
sample_instance = 'trp'
sample_instance_has_slot = 'MOLECULAR-WEIGHT'
sample_instance_non_slot = 'UTTERLY-FAKE-SLOT-NAME'
sample_class = '|Genes|'

def assert_list(l):
    assert_is_instance(l,list)

def assert_string(s): 
    assert_is_instance(s,str)

def assert_frame_list(l, n=1):
    """ Check that l is a list of valid frames of length >= n. 

    """
    
    assert_list(l)
    assert(len(l) >= n)
    for f in l:
        assert_frame(f)
    
def assert_frame(s):
    """ Check that s is a valid frame.
    
    """
    assert_is_instance(s, Frame)
    assert s._kb == ec
    assert s.test_frame()

# Ensure the exception class is present
@raises(pycyc.PathwayToolsError)
def test_pathway_tools_error():
    raise pycyc.PathwayToolsError

# Test check that a new database is valid
@raises(ValueError)
def test_open_wrong_database():
    pycyc.open('_long_name_that_should_really_not_be_a_valid_pgdb')

### Functions provided directly at module level

def test_all_orgs():
    """ pycyc.all_orgs """
    orgs = pycyc.all_orgs() # list of strings
    for org in orgs:
        assert org.islower()
        assert not org.endswith('base')
        # Do NOT call pycyc.open(org) here as this could hang the
        # server if, eg, one of the databases needs to be updated to a 
        # new Pathway Tools version.

### Functions imported from ptools-fns.html

# create_instance_w_generated_id(self, class_):
#     >>> ecocyc.create_instance_w_generated_id('|Genes|')
# Hard to test, and we'd then need to remove the frame from DB

# def get_name_string(self, item, strip_html=None, direction=None, 
#                     rxn_eqn_as_name=True, name_slot=None):
def test_get_name_string():
    assert_is_instance(ec.get_name_string(sample_reaction),str) 

def test_get_name_string_options():
    assert_is_instance(ec.get_name_string(sample_reaction,strip_html=True,
                                          direction='R2L',
                                          rxn_eqn_as_name=False),str) 

def test_get_name_string_other_slot():
    assert_is_instance(ec.get_name_string(sample_compound,
                                          name_slot='MOLECULAR-WEIGHT'),float) 


# find_pgdbs(self, substring=None):
# Not tested; apparently not available (23 August 2012)

def test_all_rxns():
    l = ec.all_rxns()
    assert_is_instance(l,list)
    assert (len(l) > 0 )

def test_all_rxn_options():
    # Wrap the frame list check and set the wrapper's description 
    # to prevent the long lists from cluttering the verbose nosetest
    # results
    test = lambda l: assert_frame_list(l)
    for option in ['smm','all','protein-small-molecule-reaction',
                   'protein-reaction','trna-reaction','spontaneous',
                   'enzyme','small-molecule','metab-pathways','metab-smm',
                   'metab-all','transport']:
        test.description = 'all_rxns option %s' % option
        
        yield test, ec.all_rxns(type=option)

@raises(pycyc.PathwayToolsError)   
def test_bad_option():
    ec.all_rxns(type='_presumably_bogus_reaction_type_')

# def genes_of_reaction(self, rxn):
def test_genes_of_reaction():
    assert_frame_list(ec.genes_of_reaction(sample_reaction))

# def substrates_of_reaction(self, rxn):
def test_substrates_of_reaction():
    assert_frame_list(ec.substrates_of_reaction(sample_reaction))

# def enzymes_of_reaction(self, rxn):
def test_enzymes_of_reaction():
    assert_frame_list(ec.enzymes_of_reaction(sample_reaction))

# def reaction_reactants_and_products(self, rxn, direction='L2R', pwy=None):
def test_reactants_and_products_default():
    r,p = ec.reaction_reactants_and_products(sample_reaction)
    assert_frame_list(r)
    assert_frame_list(p)

def test_reactants_and_products_by_dir():
    r,p = ec.reaction_reactants_and_products(sample_reaction,
                                             direction='R2L')
    assert_frame_list(r)
    assert_frame_list(p)
    
def test_reactants_and_products_by_pwy():
    r,p = ec.reaction_reactants_and_products(pathway_reaction,direction=None,
                                             pwy=sample_pathway)
    assert_frame_list(r)
    assert_frame_list(p)

@raises(pycyc.PyCycError)
def test_reactants_and_products_by_wrong_pwy():
    r,p = ec.reaction_reactants_and_products(pathway_non_reaction,
                                             direction=None,
                                             pwy=sample_pathway)

@raises(pycyc.PyCycError)
def test_reactants_and_products_bad_arguments():
    r,p = ec.reaction_reactants_and_products(sample_reaction,
                                             direction=None,
                                             pwy=None)

# def reaction_type(self, rxn, type_):
def test_is_reaction_of_type():
    assert (ec.is_reaction_of_type(sample_transport_reaction,
                             'transport') is True)

# def transported_chemicals(self, rxn):
def test_transported_chemicals():
    assert_frame_list(ec.transported_chemicals(sample_transport_reaction))

# def get_predecessors(self, rxn, pwy):
def test_predecessors():
    assert_frame_list(ec.get_predecessors(pathway_reaction, sample_pathway))

# def get_successors(self, rxn, pwy):
def test_predecessors():
    assert_frame_list(ec.get_predecessors(pathway_reaction, sample_pathway))

# def all_pathways(self, base=None, selector='all'):
def test_all_pathways_default():
    assert_frame_list(ec.all_pathways())

def test_all_pathways_options():
    assert_frame_list(ec.all_pathways(base=True, selector='small-molecule'))

# def genes_of_pathway(self, pwy):
def test_genes_of_pathway():
    assert_frame_list(ec.genes_of_pathway(sample_pathway))

# def enzymes_of_pathway(self, pwy):
def test_enzymes_of_pathway():
    assert_frame_list(ec.enzymes_of_pathway(sample_pathway))

# def compounds_of_pathway(self, pwy):
def test_compounds_of_pathway():
    assert_frame_list(ec.compounds_of_pathway(sample_pathway))

# def substrates_of_pathway(self, pwy):
#     """ List reactants and products of pwy.
def test_substrates_of_pathway():
    t = ec.substrates_of_pathway(sample_pathway)
    reactants, proper_reactants, products, proper_products = t
    assert_frame_list(reactants)
    assert_frame_list(products)
    assert_frame_list(proper_reactants)
    assert_frame_list(proper_products)

# def all_transcription_factors(self, allow_modified_forms=True):
def test_all_tfs():
    assert_frame_list(ec.all_transcription_factors())

def test_all_tfs_option():
    assert_frame_list(ec.all_transcription_factors(allow_modified_forms=False))

# def transcription_factor(self, protein):
def test_transcription_factor_true():
    assert(ec.transcription_factor(sample_tf) is _T)

def test_transcription_factor_false():
    assert(ec.transcription_factor(sample_non_tf) is _NIL)

# def all_cofactors(self):
def test_all_cofactors():
    assert_frame_list(ec.all_cofactors())

# def all_modulators(self):
def test_all_modulators():
    # Modulators clearly need not be frames: ecocyc's list of modulators
    # includes such sensible non-frame strings as "macabomycin", 
    # "chlorobiphenyl vancomycin", and "amoxicillin", as well as, 
    # more confusingly, "2".
    modulators = ec.all_modulators()
    assert_is_instance(modulators, list)
    assert len(modulators) > 0 
    for m in modulators:
        assert (isinstance(m,Frame) or isinstance(m,str))

# def monomers_of_protein(self, p):
def test_monomers_of_protein():
    assert_frame_list(ec.monomers_of_protein(sample_multicomponent_protein))

# def components_of_protein(self, p, coefficient=None):
def test_components_of_protein():
    assert_frame_list(ec.components_of_protein(sample_multicomponent_protein))

# def genes_of_protein(self, p):
def test_genes_of_protein():
    assert_frame_list(ec.genes_of_protein(sample_protein))

# def reactions_of_enzyme(self, e):
def test_reactions_of_enzyme():
    assert_frame_list(ec.reactions_of_enzyme(sample_enzyme))

# def enzyme(self, protein, type='any'):
def test_enzyme_default():
    assert (ec.enzyme(sample_enzyme) is True)

def test_enzyme_chemical_change():
    assert_frame_list(ec.enzyme(sample_enzyme, type='chemical-change'))

def test_enzyme_small_molecule():
    assert (ec.enzyme(sample_enzyme, type='small-molecule') is True)

def test_enzyme_nontransport():
    assert (ec.enzyme(sample_transport_protein, type='non-transport') is
            _NIL)

# def transporter(self, protein):
def test_transporter():
    assert (ec.transporter(sample_transport_protein) is _T)
    assert (ec.transporter(sample_enzyme) is _NIL)

# def containers_of(self, protein, exclude_self=None):
def test_containers_of():
    p = sample_protein_with_nontrivial_container
    containers = ec.containers_of(p)
    assert_frame_list(containers,2)
    assert p in containers
    proper_containers = ec.containers_of(p, exclude_self=True)
    assert_frame_list(proper_containers)
    assert p not in proper_containers


# def modified_forms(self, protein, exclude_self=None):
def test_modified_forms():
    p = sample_protein_with_modified_forms
    mod_forms = ec.modified_forms(p)
    assert_frame_list(mod_forms,2)
    assert p in mod_forms
    proper_mod_forms = ec.modified_forms(p, exclude_self=True)
    assert_frame_list(proper_mod_forms)
    assert p not in proper_mod_forms

# def modified_containers(self, protein):
def test_modified_containers():
    p = sample_protein_with_modified_forms
    assert_frame_list(ec.modified_containers(p),2)

# def top_containers(self, protein):
def test_top_containers():
    p = sample_protein_with_nontrivial_container
    assert_frame_list(ec.top_containers(p))

# def reactions_of_protein(self, p):
def test_reactions_of_protein():
    assert_frame_list(ec.reactions_of_protein(sample_enzyme))

# def regulon_of_protein(self, protein):
def test_regulon_of_protein():
    assert_frame_list(ec.regulon_of_protein(sample_tf))

# def transcription_units_of_protein(self, protein):
def test_tus_of_protein():
    assert_frame_list(ec.transcription_units_of_protein(sample_tf))

# def regulator_proteins_of_transcription_unit(self, tu):
def test_regulators_of_tu():
    assert_frame_list(ec.regulator_proteins_of_transcription_unit(sample_tu))

# def full_enzyme_name(self, enzyme):
def test_full_enzyme_name():
    assert_is_instance(ec.full_enzyme_name(sample_enzyme),str)

# def enzyme_activity_name(self, enzyme, reaction=None):
def test_enzyme_activity_name():
    # sample_multienzyme, mulitenzyme_r1, multienzyme_r2 = multifn_info
    name1 = ec.enzyme_activity_name(sample_multienzyme)
    name2 = ec.enzyme_activity_name(multienzyme_r1)
    assert_is_instance(name1, str)
    assert_is_instance(name2, str)
    assert name1
    assert name2
    assert name1 != name2

# def enzymes_of_gene(self, gene):
def test_enzymes_of_gene():
    assert_frame_list(ec.enzymes_of_gene(sample_gene))

# def all_products_of_gene(self, gene):
def test_all_products_of_gene():
    assert_frame_list(ec.all_products_of_gene(sample_gene))

# def reactions_of_gene(self, gene):
def test_reactions_of_gene():
    assert_frame_list(ec.reactions_of_gene(sample_gene))

# def pathways_of_gene(self, gene):
def test_pathways_of_gene():
    assert_frame_list(ec.pathways_of_gene(sample_gene))

# def chromosome_of_gene(self, gene):
def test_chromosome_of_gene():
    assert_frame(ec.chromosome_of_gene(sample_gene))

# def transcription_units_of_gene(self, gene):
def test_transcription_units_of_gene():
    assert_frame_list(ec.transcription_units_of_gene(sample_gene))

# def genes_regulating_gene(self, gene):
def test_genes_regulating_gene():
    assert_frame_list(ec.genes_regulating_gene(sample_gene))

# def genes_regulated_by_gene(self, gene):
def test_genes_regulated_by_gene():
    assert_frame_list(ec.genes_regulated_by_gene(sample_regulatory_gene))

# def terminators_affecting_gene(self, gene):
def test_terminators_affecting_gene():
    assert_frame_list(ec.terminators_affecting_gene(sample_gene))

# def get_gene_sequence(self, gene):
def test_get_gene_sequence():
    sequence = ec.get_gene_sequence(sample_gene)
    assert_is_instance(sequence, str)
    assert sequence
    bases = set('ATGC')
    assert bases.issuperset(set(sequence))

# def nucleotide_to_protein_sequence(self, sequence, code_num=1):
def test_nucleotide_to_protein():
    # The LISP function appears to translate the first three bases to 
    # M and interpret the last three as a silent stop codon regardless of
    # their actual values.
    dna = "ATGGAACGCTACTAA"
    protein = ec.nucleotide_to_protein_sequence(dna)
    assert protein == 'MERY'

# def transcription_unit_promoter(self, tu):
def test_transcription_unit_promoter():
    assert_frame_list(ec.transcription_unit_promoter(sample_tu))

# def transcription_unit_genes(self, tu):
def test_transcription_unit_promoter():
    assert_frame(ec.transcription_unit_promoter(sample_tu))

# def transcription_unit_binding_sites(self, tu):
#     """ List the DNA binding sites within transcription unit tu.
def test_transcription_unit_binding_sites():
    assert_frame_list(ec.transcription_unit_binding_sites(sample_tu))

# def transcription_unit_mrna_binding_sites(self, tu):
def test_transcription_unit_mrna_binding_sites():
    s = sample_tu_with_mrna_site
    assert_frame_list(ec.transcription_unit_mrna_binding_sites(s))

# def transcription_unit_transcription_factors(self, tu):
def test_transcription_unit_transcription_factors():
    assert_frame_list(ec.transcription_unit_transcription_factors(sample_tu))

# def binding_site_transcription_factors(self, bsite):
def test_binding_site_transcription_factors():
    site = sample_binding_site
    assert_frame_list(ec.binding_site_transcription_factors(site))

# def transcription_unit_terminators(self, tu):
def test_transcription_unit_terminators():
    assert_frame_list(ec.transcription_unit_terminators(sample_tu))

# def transcription_unit_activators(self, tu):
def test_transcription_unit_activators():
    assert_frame_list(ec.transcription_unit_activators(sample_tu))

# def transcription_unit_inhibitors(self, tu):
def test_transcription_unit_inhibitors():
    assert_frame_list(ec.transcription_unit_inhibitors(sample_tu))

# def containing_tus(self, site):
def test_containing_tus():
    assert_frame_list(ec.containing_tus(sample_binding_site))

# def all_transported_chemicals(self):
def test_all_transported_chemicals():
    assert_frame_list(ec.all_transported_chemicals())

# def reactions_of_compound(self, cpd):
def test_reactions_of_compound():
    assert_frame_list(ec.reactions_of_compound(sample_compound))

# def direct_activators(self, item):
def test_direct_activators():
    assert_frame_list(ec.direct_activators(sample_activated_tu))

# def direct_inhibitors(self, item):
def test_direct_inhibitors():
    assert_frame_list(ec.direct_inhibitors(sample_inhibited_tu))

def test_reaction_type():
    allowed_values = set(['trna-reaction', 'small-molecule',
                          'transport', 'protein-small-molecule-reaction',
                          'protein-reaction', 'null-reaction',
                          'other'])
    sample_reaction_list = ['trypsyn-rxn', 'TRANS-RXN-112',
                            'ALANINE--TRNA-LIGASE-RXN','TORT-RXN',
                            'UBIQUINOL-OX-RXN','2.7.7.8-RXN',
                            'RXN0-4082']
    for r in sample_reaction_list:
        assert ec.reaction_type(r) in allowed_values

### GFP functions

# def get_slot_values(self, frame, slot):
def test_get_slot_values():
    synonyms = ec.get_slot_values('trp','synonyms')
    assert_is_instance(synonyms,list)
    for word in ['trp', 'W', 'tryptacin', 'trofan', 'tryptophan',              
                 '2-amino-3-indolylpropanic acid']:
        assert word in synonyms

# def get_slot_value(self, frame, slot):
def test_get_slot_value():
    assert ec.get_slot_value('trp','common-name') == 'L-tryptophan'

# def get_class_all_instances(self, class_):
def test_get_class_all_instances():
    assert_frame_list(ec.get_class_all_instances('|Reactions|'))

# def instance_all_instance_of(self, instance, class_):
def test_instance_of():
    assert(ec.instance_all_instance_of(sample_reaction,'|Reactions|') is _T)
    assert(ec.instance_all_instance_of(sample_protein,'|Reactions|') is _NIL)

# def member_slot_value(self, frame, slot, value):
def test_member_slot_value():
    assert(ec.member_slot_value('trp','synonyms','W') == 'W')
    assert(ec.member_slot_value('trp','synonyms', symbol='W') is _NIL)
    assert(ec.member_slot_value('trp','synonyms','leucine') is _NIL)
    assert(ec.member_slot_value('trp','appears-in-left-side-of',
                                symbol='RXN0-287') == 'RXN0-287')
    assert(ec.member_slot_value('trp','appears-in-left-side-of',
                                'RXN0-287') is _NIL)

@raises(ValueError)
def test_member_slot_value_no_value():
    ec.member_slot_value('trp','synonyms')

# def fequal(self, frame1, frame2):
def test_fequal():
    # I don't know how to come up with a case where the result of
    # fequal is not obvious, so these tests are not particularly useful
    assert(ec.fequal(sample_reaction, sample_reaction))
    assert(ec.fequal(sample_reaction, sample_compound) is _NIL)

# def get_frame_labeled(self, label):
def test_get_frame_labeled():
    assert_frame_list(ec.get_frame_labeled(sample_label))

## Additional GFP functions from JavaCyc
def test_coercible_to_frame():
    assert(ec.is_coercible_to_frame('trp') is _T)
    assert(ec.is_coercible_to_frame('not-a-frame') is _NIL)

def test_type_of():
    assert(ec.is_class_all_type_of('|Reactions|',sample_reaction) is _T)
    assert(ec.is_class_all_type_of('|Reactions|',sample_gene) is _NIL)

def test_get_instance_types():
    assert_frame_list(ec.get_instance_direct_types(sample_reaction))
    assert_frame_list(ec.get_instance_all_types(sample_reaction))

def test_get_frame_slots():
    slots = ec.get_frame_slots(sample_compound)
    assert_is_instance(slots, list)
    assert slots


## Additional GFP functions added independently

def test_is_class():
    assert ec.is_class(sample_class)
    assert not ec.is_class(sample_instance)

def test_is_instance():
    assert not ec.is_instance(sample_class)
    assert ec.is_instance(sample_instance)

def test_is_slot():
    assert ec.is_slot(sample_instance, sample_instance_has_slot)
    assert not ec.is_slot(sample_instance, sample_instance_non_slot)

## Python conveniences, object behavior

def test_frame_in():
    sample_frame = 'trp'
    sample_non_frame = '_not_a_frame_'
    assert sample_frame in ec
    assert sample_non_frame not in ec

def test_get_frame():
    sample_frame = 'trp'
    trp = Frame('trp',ec)
    assert ec[sample_frame] == trp

@raises(KeyError)
def test_get_non_frame():
    sample_non_frame = '_not_a_frame_'
    ec[sample_non_frame]

def test_special_attributes():
    special_attributes = {'|Genes|': ec.genes,
                          '|Reactions|': ec.reactions,
                          '|Compounds|': ec.compounds,
                          '|Pathways|': ec.pathways}    
    for k,v in special_attributes.iteritems():
        assert v == ec.get_class_all_instances(k)

@raises(AttributeError)
def test_bad_attribute():
    ec.notanattribute

def test_dir():
    l = dir(ec)
    for a in ['genes', '_org_id', 'is_class', '__doc__']:
        assert a in l
