"""
Test still further pathway-tools functionality using nosetests.

The EcoCyc database is queried to run these tests. In most cases the
successful return of something that looks vaguely appropriate (not
necessarily the correct result) is enough for a test to pass. Changes
in the database will cause tests to fail even when the code is
performing correctly, but short of installing a purpose-built test
PGDB I am not sure how to avoid this.

Docstrings here typically describe the arguments to the function tested,
(for reference in modifying the test) rather than the test itself. For
historical reasons, the arguments may be listed out of order; consult the
interface class method docstrings or the LISP API documentation.

"""

import pycyc
from nose.tools import raises, assert_is_instance, assert_raises
from frame import Frame

# Specify the values which should be returned for true and false queries
# (which may change in the future)
_NIL = None
_T = True

ec = pycyc.open('ecoli')

sample_reaction = 'trypsyn-rxn' # '1.8.4.13-RXN', 'oropribtrans-rxn'
rxn = sample_reaction
# pathway_info = ('TRPSYN-PWY','RXN0-2381','OROPRIBTRANS-RXN')
# (sample_pathway, pathway_reaction, pathway_non_reaction) = pathway_info
pwy = 'TRPSYN-PWY'
# sample_compound = 'trp' 
gene = 'EG11024' # should encode at least one enzyme in at least one pwy
#                         # should lie in at least one transcription unit 
#                         # and have at least one regulator
# sample_enzyme = 'TRYPSYN' # should be of types chemical-change, small-molecule,
#                           # non-transport
protein = 'TRYPSYN'
# sample_transport_reaction='TRANS-RXN-157'
# # The sample transport protein should NOT be an enzyme of type non-transport
# sample_transport_protein='CPLX-165' 
tf = 'PHOSPHO-CREB'
# sample_regulatory_gene = 'EG11218'
# sample_non_tf = 'TRYPSYN'
# sample_multicomponent_protein = 'TRYPSYN'
# sample_frame_label = 'pyrophosphate'
# sample_protein_with_nontrivial_container = 'TRYPSYN-APROTEIN'
# sample_protein_with_modified_forms = 'CREB-MONOMER'
tu = 'TU00038'
# sample_tu_with_mrna_site = 'TU00355'
# sample_inhibited_tu = 'TU0-1123'
# sample_activated_tu = 'TU0-13935'
# multifn_info = ('RIB5PISOMB-CPLX','RXN0-303','RIB5PISOM-RXN')
# sample_multienzyme, multienzyme_r1, multienzyme_r2 = multifn_info
# sample_label = 'pyrophosphate'
binding_site = 'BS00117'
# sample_instance = 'trp'
# sample_instance_has_slot = 'MOLECULAR-WEIGHT'
# sample_instance_non_slot = 'UTTERLY-FAKE-SLOT-NAME'
# sample_class = '|Genes|'
generic_rxn = 'ALCOHOL-DEHYDROG-GENERIC-RXN'
enzrxn = 'TRYPSYN-ENZRXN'
promoter = 'PM721'
reg_frame = 'REG0-11121'


def assert_list(l):
    assert_is_instance(l,list)

def assert_string(s): 
    assert_is_instance(s,str)

def assert_frame_list(l, n=1):
    """ Check that l is a list of valid frames of length >= n. 

    Incompletely implemented! Currently this just checks that L is a
    non-empty list, as coercible_to_frame_p is not available yet.

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

def assert_instance_list(l, *classes):
    """ Check l is a list of frames belonging to at least one of the classes. 

    """
    if not l:
        return
    classes = [ec[class_] for class_ in classes]
    for item in l:
        assert_instance_of_frames(item, classes)

def assert_class_instance(frame, *classes):
    """ Check frame belongs to at least one of the classes. 

    """
    classes = [ec[class_] for class_ in classes]
    assert_instance_of_frames(frame, classes)

def assert_instance_of_frames(frame, classframes):
    """ Verify frame is an instance of at least one of a list of classframes.

    """
    assert any([frame in class_ for class_ in classframes])

############################################################        
def test_all_direct_forms_of_protein():
    """ all_direct_forms_of_protein
  
    protein: An instance of the class |Proteins|.

    """
    result = ec.all_direct_forms_of_protein(protein)
    assert_instance_list(result, "|Proteins|")

def test_rxns_w_isozymes():
    """ rxns_w_isozymes

    rxns: A list of instances of the class |Reactions|. Defaults to
        the result of (all-rxns 'enzyme').

    """
    rxns = ['TRYPSYN-RXN', 'RIBULPEPIM-RXN']
    result = ec.rxns_w_isozymes()
    assert_instance_list(result, "|Reactions|")
    result = ec.rxns_w_isozymes(rxns=rxns)
    assert_instance_list(result, "|Reactions|")

def test_all_substrates():
    """ all_substrates
     Arguments
    ---------
    rxns: A list of reaction frames.

    """
    rxns = ['TRYPSYN-RXN', 'RIBULPEPIM-RXN']
    result = ec.all_substrates(rxns)
    for r in result:
        if isinstance(r,str):
            continue
        else:
            assert_class_instance(r,'|Compounds|') # check

def test_gene_transcription_units():
    """ gene_transcription_units
     Arguments
    ---------
    gene: An instance of class |Genes|.

    """
    result = ec.gene_transcription_units(gene)
    assert_instance_list(result, '|Transcription-Units|')

def test_previous_gene_on_replicon():
    """ previous_gene_on_replicon
     Arguments
    ---------
    g: An instance of class |Genes|.

    """
    result = ec.previous_gene_on_replicon(gene)
    assert_is_instance(result, list)
    assert(result[0] is None or result[0] in ec['|Genes|'])
    # The following condition is never met in ecoli-- no linear replicon?
    if len(result) > 1: 
        assert(result[1] in (None, ':first'))

def test_binding_site_transcription_units():
    """ binding_site_transcription_units
     Arguments
    ---------
    promoter: An instance of class |DNA-Binding-Sites| or |mRNA-
        Binding-Sites|.

    """
    result = ec.binding_site_transcription_units(binding_site)
    assert_instance_list(result, '|Transcription-Units|')

def test_leader_peptide():
    """ leader_peptide
     Arguments
    ---------
    protein: An instance of the class |Proteins|.

    """
    result = ec.leader_peptide(protein)
    assert result in [True, False, None]

def test_chromosome_of_object():
    """ chromosome_of_object
     Arguments
    ---------
    object: An instance of class |All-Genes|, |Transcription-Units|,
        |Promoters|, |Terminators|, |Misc-Features|, or |DNA-
        Binding-Sites|.

    """
    result = ec.chromosome_of_object(gene)
    assert_class_instance(result, "|Genetic-Elements|")

def test_regulators_of_gene():
    """ regulators_of_gene
     Arguments
    ---------
    by_function: If non-None, then return two values: a list of
        activator proteins and a list of inhibitor proteins.
    gene: An instance of the class |Genes|.

    """
    
    # first test with by_function
    result = ec.regulators_of_gene(gene,by_function=True)
    assert_is_instance(result, list)
    assert_is_instance(result[0], list)
    assert_is_instance(result[1], list)
    if result[0]:
        assert_frame_list(result[0])
    if result[1]:
        assert_frame_list(result[1])
    # then without
    result = ec.regulators_of_gene(gene,by_function=False)
    assert_frame_list(result)

def test_tfs_bound_to_compound():
    """ tfs_bound_to_compound
     Arguments
    ---------
    include_inactive: If non-None, then the inactive form of the
        protein is also checked. See the function #'transcription-
        factor? for more information.
    cpd: An instance of class |Compounds|.

    """
    cpd = 'TRP'
    result = ec.tfs_bound_to_compound(cpd)
    assert_instance_list(result, "|Proteins|")
    result = ec.tfs_bound_to_compound(cpd,True)
    assert_instance_list(result, "|Proteins|")

def test_regulon_of_protein():
    """ regulon_of_protein
     Arguments
    ---------
    protein: An instance frame of class |Proteins|.

    """
    result = ec.regulon_of_protein(tf)
    assert_instance_list(result, "|Transcription-Units|")

def test_all_pathways():
    """ all_pathways
     Arguments
    ---------
    base: If this argument evaluates to true, only includes base
        pathways. Otherwise, all pathways, including superpathways,
        will be returned.
    selector: Selects whether all pathways, or just small-molecule
        metabolism base pathways. Can take either a symbol of 'all'
        or 'small-molecule'. Defaults to <code>'all'</code>.

    """
    result = ec.all_pathways()
    assert_instance_list(result, "|Pathways|")
    result = ec.all_pathways(base=True,selector='small-molecule')
    assert_instance_list(result, "|Pathways|")

def test_binding_sites_affecting_gene():
    """ binding_sites_affecting_gene
     Arguments
    ---------
    gene: An instance of the class |Genes|.

    """
    result = ec.binding_sites_affecting_gene(gene)
    assert_instance_list(result, "|DNA-Binding-Sites|")

def test_unmodified_gene_product():
    """ unmodified_gene_product
     Arguments
    ---------
    gene: An instance of class |Genes|.

    """
    result = ec.unmodified_gene_product(gene)
    assert (result == 'RNA' or result in ec['|Polypeptides|'])

def test_gene_clusters():
    """ gene_clusters
     Arguments
    ---------
    genes: A list of instances of class |Genes|.
    max_gap: An integer representing the number of genes any pair
        from genes can be from one another. Default value is 10.

    """
    genes=['EG11010','EG11011','EG11024','EG11025']
    result = ec.gene_clusters(genes)
    assert_is_instance(result, list)
    for l in result:
        assert l[0] in genes
        assert_instance_list(l, "|Genes|")

def test_transcription_unit_all_components():
    """ transcription_unit_all_components
     Arguments
    ---------
    tu: An instance of the class |Transcription-Units|.

    """
    result = ec.transcription_unit_all_components(tu)
    assert_instance_list(result, "|Transcription-Units|",
                         "|mRNA-Binding-Sites|", "|DNA-Binding-Sites|",
                         "|Promoters|", "|Genes|", "|Terminators|")
    

def test_pathway_hole():
    """ pathway_hole
     Arguments
    ---------
    hole_if_any_gene_without_position: If true, then genes without
        specified coordinates for the current organism's genome are
        not counted when determining the status of the reaction.
    rxn: An instance of the class |Reactions|.

    """
    result = ec.pathway_hole(rxn)
    assert result in [True, False, None]
    result = ec.pathway_hole(rxn,hole_if_any_gene_without_position=True)
    assert result in [True, False, None]

def test_substrate_of_generic_rxn():
    """ substrate_of_generic_rxn
     Arguments
    ---------
    cpd: An instance of class |Compounds|.
    rxn: An instance of class |Reactions|.

    """
    rxn = 'ALCOHOL-DEHYDROG-GENERIC-RXN'
    cpd = 'LACTALD'
    result = ec.substrate_of_generic_rxn(cpd, generic_rxn)
    assert (result == True)

def test_cofactors_and_pgroups_of_enzrxn():
    """ cofactors_and_pgroups_of_enzrxn
     Arguments
    ---------
    enzrxn: An instance of the class |Enzymatic-Reactions|.

    """
    result = ec.cofactors_and_pgroups_of_enzrxn(enzrxn)
    assert_is_instance(result, list)
    for c in result:
        assert (isinstance(c,str) or c in ec['|Chemicals|'])

def test_pathways_of_enzrxn():
    """ pathways_of_enzrxn
     Arguments
    ---------
    include_super_pwys: If non-None, then not only will the direct
        pathways in which enzrxn is associated in be returned, but
        also any enclosing super-pathways. If enzrxn is associated
        with a reaction that is directly associated with a super-
        pathway, then the function might return super-pathways even
        if this option is None.
    enzrxn: An instance of the class |Enzymatic-Reactions|.

    """
    result = ec.pathways_of_enzrxn(enzrxn)
    assert_instance_list(result, "|Pathways|")
    result = ec.pathways_of_enzrxn(enzrxn, include_super_pwys=True)
    assert_instance_list(result, "|Pathways|")

def test_transcription_factor_active_forms():
    """ transcription_factor_active_forms
     Arguments
    ---------
    tfs: An instance of the class |Proteins|.

    """
    result = ec.transcription_factor_active_forms(tf)
    assert_instance_list(result, "|Proteins|")

def test_phantom_gene():
    """ phantom_gene
     Arguments
    ---------
    gene: An instance of the class |Genes|.

    """
    result = ec.phantom_gene(gene)
    assert result in [True, False, None]

def test_polypeptide_or_homomultimer():
    """ polypeptide_or_homomultimer
     Arguments
    ---------
    protein: An instance of the class |Proteins|.

    """
    result = ec.polypeptide_or_homomultimer(protein)
    assert result in [True, False, None]

def test_cotranscribed_genes():
    """ cotranscribed_genes
     Arguments
    ---------
    gene: An instance of class |Genes|.

    """
    result = ec.cotranscribed_genes(gene)
    assert_instance_list(result, "|Genes|")

def test_protein_coding_gene():
    """ protein_coding_gene
     Arguments
    ---------
    gene: An instance of the class |Genes|.

    """
    result = ec.protein_coding_gene(gene)
    assert result == True

def test_transcription_unit_first_gene():
    """ transcription_unit_first_gene
     Arguments
    ---------
    tu: An instance of the class |Transcription-Units|.

    """
    result = ec.transcription_unit_first_gene(tu)
    assert_class_instance(result, "|Genes|")

def test_unmodified_form():
    """ unmodified_form
     Arguments
    ---------
    protein: An instance of the class |Proteins|.

    """
    result = ec.unmodified_form(protein)
    assert_class_instance(result, "|Proteins|")

def test_next_gene_on_replicon():
    """ next_gene_on_replicon
     Arguments
    ---------
    g: An instance of class |Genes|.

    """
    result = ec.next_gene_on_replicon(gene)
    assert_is_instance(result, list)
    assert(result[0] in ec['|Genes|'] or result[0] is None)
    if len(result) > 1: # again, doesn't occur in ec
        assert(result[1] in [None, ':last'])

def test_rxn_in_compartment():
    """ rxn_in_compartment
     Arguments
    ---------
    default_ok: If true, then we return true if the reaction has no
        associated compartment information, or one of its associated
        locations is a super-class of one of the members of the
        compartments argument.
    rxn: An instance of the class |Reactions|.
    pwy: If supplied, the search for associated enzymes of the
        argument rxn is limited to the given child of |Pathways|.
    loose: If true, then the compartments 'CCO-CYTOPLASM and 'CCO-
        CYTOSOL are treated as being the same compartment.
    compartments: A list of cellular compartments, as defined in the
        Cellular Components Ontology. See frame 'CCO.

    """
    compartments = ['CCO-CYTOSOL', 'CCO-MIT']
    result = ec.rxn_in_compartment(rxn, compartments)
    # This test may not be exactly correct as a list of compartments
    # could be returned (maybe?) The API document doesn't specify.
    assert (result in compartments or result in [False, None])
    result = ec.rxn_in_compartment(rxn, compartments, pwy=pwy, loose=True,
                                   default_ok=True)
    assert (result in compartments or result in [False, None])

def test_rxn_w_isozymes():
    """ rxn_w_isozymes
     Arguments
    ---------
    rxn: An instance of the class |Reactions|.

    """
    result = ec.rxn_w_isozymes('RIB5PISOM-RXN')
    assert result == True

def test_small_molecule_cplxes_of_prot():
    """ small_molecule_cplxes_of_prot
     Arguments
    ---------
    protein: An instance of the class |Proteins|.

    """
    result = ec.small_molecule_cplxes_of_prot('CPLX0-7759')
    assert_instance_list(result, "|Proteins|")

def test_noncontiguous_pathway():
    """ noncontiguous_pathway
     Arguments
    ---------
    pwy: An instance of the class |Pathways|.

    """
    result = ec.noncontiguous_pathway(pwy)
    assert result in [True, False, None]

def test_autocatalytic_reactions_of_enzyme():
    """ autocatalytic_reactions_of_enzyme
     Arguments
    ---------
    prot: An instance frame of class '|Proteins|.

    """
    result = ec.autocatalytic_reactions_of_enzyme('FERREDOXIN-MONOMER')
    assert_instance_list(result, "|Reactions|")

def test_homomultimeric_containers_of():
    """ homomultimeric_containers_of
     Arguments
    ---------
    protein: An instance of the class |Proteins|.
    exclude_self: If true, then protein will not be included in the
        return value.

    """
    result = ec.homomultimeric_containers_of(protein)
    assert_instance_list(result, "|Proteins|")
    result = ec.homomultimeric_containers_of(protein, exclude_self = True)
    assert_instance_list(result, "|Proteins|")

def test_transcription_factor_ligands():
    """ transcription_factor_ligands
     Arguments
    ---------
    tfs: An instance or a list of instances of the class |Proteins|.
        If tfs is not the active form, then the active form is
        determined automatically.

    mode: One of the following values: 'activator', 'inhibitor', or
        'both'.
    """
    
    tf = 'EG20249-MONOMER' # has ligands
    result = ec.transcription_factor_ligands(tf, mode='both')
    chemicals = ec['|Chemicals|']
    for c in result:
        assert(isinstance(c, str) or c in chemicals)
    tfs = ['PHOSPHO-CREB', 'PHOSPHO-ARCA', 'EG20249-MONOMER']
    result = ec.transcription_factor_ligands(tfs, mode='inhibitor')
    for c in result:
        assert(isinstance(c, str) or c in chemicals)

def test_pseudo_gene():
    """ pseudo_gene
     Arguments
    ---------
    gene: An instance of the class |Genes|.

    """
    result = ec.pseudo_gene(gene)
    assert result in [True, False, None]

def test_protein():
    """ protein
     Arguments
    ---------
    frame: An instance of the class |Proteins|.

    """
    assert ec.protein(protein)
    assert not ec.protein(rxn)

def test_binding_site_to_regulators():
    """ binding_site_to_regulators
     Arguments
    ---------
    bsite: An instance of class |DNA-Binding-Sites|.

    """
    result = ec.binding_site_to_regulators(binding_site)
    assert_instance_list(result, "|Proteins|")

def test_enzrxn_inhibitors():
    """ enzrxn_inhibitors
     Arguments
    ---------
    phys_relevant_only: If true, then only return inhibitors that
        are associated with |Regulation| instances that have the
        'physiologically-relevant? slot set to non-None.
    er: An instance of the class |Enzymatic-Reactions|.

    """
    result = ec.enzrxn_inhibitors(enzrxn)
    assert_instance_list(result, "|Chemicals|")

def test_specific_forms_of_rxn():
    """ specific_forms_of_rxn
     Arguments
    ---------
    rxn: A child of the class |Reactions|.

    """
    result = ec.specific_forms_of_rxn(generic_rxn)
    assert_instance_list(result, "|Reactions|")

def test_reduce_modified_proteins():
    """ reduce_modified_proteins
     Arguments
    ---------
    prots: A list of instances of the class |Proteins|.
    debind: When non-None, the proteins are further simplified by
        obtaining the unbound form of the protein, if it is bound to
        a small molecule.

    """
    prots = ['TRYPSYN', 'PHOSPHO-CREB']
    result = ec.reduce_modified_proteins(prots)
    assert_instance_list(result, "|Proteins|")
    result = ec.reduce_modified_proteins(prots,debind=True)
    assert_instance_list(result, "|Proteins|")

def test_containing_chromosome():
    """ containing_chromosome
     Arguments
    ---------
    site: An instance of class |Transcription-Units|, |mRNA-Binding-
        Sites|, |DNA-Binding-Sites|, |Promoters|, |Genes|, or
        |Terminators|.

    """
    result = ec.containing_chromosome(binding_site)
    assert_class_instance(result, "|Genetic-Elements|")

def test_unmodified_gene_products():
    """ unmodified_gene_products
     Arguments
    ---------
    gene: An instance of class |Genes|.

    """
    result = ec.unmodified_gene_products(gene)
    assert_is_instance(result, list)
    for p in result:
        assert (p == 'RNA' or p in ec['|Polypeptides|'])

def test_genes_in_same_operon():
    """ genes_in_same_operon
     Arguments
    ---------
    gene: An instance of class |Genes|.

    """
    result = ec.genes_in_same_operon(gene)
    assert_instance_list(result, "|Genes|")

def test_all_transporters():
    """ all_transporters
     Arguments
    ---------
    None.

    """
    result = ec.all_transporters()
    assert_instance_list(result, "|Proteins|")

def test_genes_of_proteins():
    """ genes_of_proteins
     Arguments
    ---------
    p: A list of instances of the class |Proteins|.

    """
    prots = ['TRYPSYN', 'PHOSPHO-CREB']
    result = ec.genes_of_proteins(prots)
    assert_instance_list(result, "|Genes|")

def test_regulators_of_operon_transcription():
    """ regulators_of_operon_transcription
     Arguments
    ---------
    by_function: If non-None, then return two values: a list of
        activator proteins and a list of inhibitor proteins.
    operon_list: A list of instances of the class |Transcription-
        Units|.

    """
    operon_list = ['TU0-1123', 'TU0-13935']
    result = ec.regulators_of_operon_transcription(operon_list)
    # In practice, these are not, as advertised, always proteins
    assert_frame_list(result)
    result = ec.regulators_of_operon_transcription(operon_list,
                                                   by_function=True)
    assert_frame_list(result[0])
    assert_frame_list(result[1])

def test_neighboring_genes():
    """ neighboring_genes
     Arguments
    ---------
    g2: An instance of class |Genes|.
    g1: An instance of class |Genes|.
    n: An integer representing the number of genes g1 and g2 can be
        from one another. Default value is 10.

    """
    g2 = 'EG11024'
    g1 = 'EG11027'
    result = ec.neighboring_genes(g2, g1)
    assert result
    result = ec.neighboring_genes(g2, g1, 2)
    assert not result

def test_pathway_components():
    """ pathway_components
     Arguments
    ---------
    rxn_list: The list of reactions to use as the starting list of
        connected component clusters. Defaults to (get-slot-values
        pwy 'reaction-list).
    pred_list: The list of reaction predecessors to iterate from in
        order to cluster the reactions in rxn_list. Defaults to
        (get-slot-values pwy 'predecessors).
    pwy: An instance of the class |Pathways|, which is not a super-
        pathway (i.e., does not have any entries in its 'sub-
        pathways slot).

    """
    # Not testing the options as it's not clear what is going on here
    # anyway. 
    result = ec.pathway_components(pwy)
    assert_is_instance(result, list)
    for l in result[0]:
        assert_instance_list(l,'|Reactions|')
    assert_is_instance(result[1], int)
    assert_is_instance(result[2], int)
    # Test optional arguments?

def test_direct_regulators():
    """ direct_regulators
     Arguments
    ---------
    x: A frame.
    filter_fn: A predicate used to filter the regulation objects
        used to find the regulators.

    """
    x = 'TU0-1123'
    result = ec.direct_regulators(x)
    assert_frame_list(result)
    result = ec.direct_regulators(x, filter_fn = 'coercible-to-frame-p')
    assert_frame_list(result)

def test_all_genetic_regulation_proteins():
    """ all_genetic_regulation_proteins
     Arguments
    ---------
    class: The class |Regulation| or a subclass. It defaults to
        |Regulation-of-Transcription-Initiation|.
    allow_modified_forms: A boolean value. If true, modified and
        unmodified forms of the protein are returned. If false, then
        only unmodified forms of the proteins are returned. The
        default value is True.

    """
    result = ec.all_genetic_regulation_proteins()
    assert_frame_list(result)
    result = ec.all_genetic_regulation_proteins(class_='|Regulation|', 
                                                allow_modified_forms=False)
    assert_instance_list(result, '|Proteins|')

# def test_compartment_of_rxn():
# # Not implemented: LISP function does not exist.
#     """ compartment_of_rxn
#      Arguments
#     ---------
#     default: The default compartment for reactions without any
#         compartment annotations on their substrates. The default
#         value is 'CCO-CYTOSOL.
#     rxn: An instance of the class |Reactions|.

#     """
#     result = ec.compartment_of_rxn(rxn)
#     assert(result in ec['CCO'])

def test_rxn_present():
    """ rxn_present
     Arguments
    ---------
    rxn: An instance of the class |Reactions|.

    """
    result = ec.rxn_present(rxn)
    assert result

def test_promoter_binding_sites():
    """ promoter_binding_sites
     Arguments
    ---------
    promoter: An instance of class |Promoters|.

    """
    result = ec.promoter_binding_sites(promoter)
    assert_instance_list(result, "|DNA-Binding-Sites|")

def test_all_forms_of_protein():
    """ all_forms_of_protein
     Arguments
    ---------
    protein: An instance of the class |Proteins|.

    """
    result = ec.all_forms_of_protein(protein)
    assert_instance_list(result, "|Proteins|")

def test_base_components_of_protein():
    """ base_components_of_protein
     Arguments
    ---------
    p: An instance of the class |Proteins|.
    exclude_small_molecules: If None, then small molecule components
        are also returned. Default value is true.

    """
    result = ec.base_components_of_protein(protein)
    assert_instance_list(result[0], '|Polypeptides|',
                         '|RNAs|', '|Compounds|')
    for c in result[1]:
        assert(isinstance(c, int) or isinstance(c, float))

def test_transcription_unit_activation_frames():
    """ transcription_unit_activation_frames
     Arguments
    ---------
    tu: An instance of the class |Transcription-Units|.

    """
    result = ec.transcription_unit_activation_frames(tu)
    assert_instance_list(result, "|Regulation|")

def test_all_protein_complexes():
    """ all_protein_complexes
     Arguments
    ---------
    filter: The type of protein complexes to return. The argument
        must be one of the following values:
        'hetero': Return all heteromultimers.
        'homo': Return all homomultimers.
        'all' or True: Return all protein complexes.

    """
    result = ec.all_protein_complexes()
    assert_instance_list(result, '|Proteins|') # properly, list of p. complexes

# Not implemented.
# def test_protein_in_compartment():
#     """     Arguments
#     ---------
#     default_ok: If true, then we return true if the reaction has no
#         associated compartment information, or one of its associated
#         locations is a super-class of one of the members of the
#         compartments argument.
#     rxn: An instance of the class |Reactions|.
#     pwy: If supplied, the search for associated enzymes of the
#         argument rxn is limited to the given child of |Pathways|.
#     loose: If true, then the compartments 'CCO-CYTOPLASM and 'CCO-
#         CYTOSOL are treated as being the same compartment.
#     compartments: A list of cellular compartments, as defined in the
#         Cellular Components Ontology. See frame 'CCO.

#     """
#     result = ec.protein_in_compartment(rxn, compartments)
#     assert result in [True, False, None]

def test_reaction_without_sequenced_enzyme():
    """ reaction_without_sequenced_enzyme
     Arguments
    ---------
    complete: If true, the predicate will return true when there is
        any associated gene without a sequence. If None, the
        predicate will return true when all associated genes are
        without a sequence.
    rxn: An instance of the class |Reactions|.

    """
    result = ec.reaction_without_sequenced_enzyme(rxn)
    assert result in [True, False, None]
    result = ec.reaction_without_sequenced_enzyme(rxn, complete=True)
    assert result in [True, False, None]

def test_dna_binding_site():
    """ dna_binding_site
     Arguments
    ---------
    gene: A frame.

    """
    result = ec.dna_binding_site(gene)
    assert result in [True, False, None]

def test_variants_of_pathway():
    """ variants_of_pathway
     Arguments
    ---------
    pwy: An instance of the class |Pathways|.

    """
    result = ec.variants_of_pathway(pwy)
    assert_instance_list(result, "|Pathways|")

def test_all_sigma_factors():
    """ all_sigma_factors
     Arguments
    ---------
    None.

    """
    result = ec.all_sigma_factors()
    assert_instance_list(result, '|Sigma-Factors|')

def test_all_enzymes():
    """ all_enzymes
     Arguments
    ---------
    type: A type as taken from the argument to #'enzyme?. Defaults
        to 'chemical-change'.

    """
    result = ec.all_enzymes()
    assert_instance_list(result, "|Proteins|")
    result = ec.all_enzymes(type='any')
    assert_instance_list(result, "|Proteins|")

def test_all_operons():
    """ all_operons
     Arguments
    ---------
    None.

    """
    result = ec.all_operons()
    for l in result:
        assert_instance_list(l, '|Transcription-Units|')

def test_transcription_unit_inhibition_frames():
    """ transcription_unit_inhibition_frames
     Arguments
    ---------
    tu: An instance of the class |Transcription-Units|.

    """
    result = ec.transcription_unit_inhibition_frames(tu)
    assert_instance_list(result, "|Regulation|")

def test_unmodified_or_unbound_form():
    """ unmodified_or_unbound_form
     Arguments
    ---------
    protein: An instance of the class |Proteins|.

    """
    result = ec.unmodified_or_unbound_form(protein)
    assert_class_instance(result, "|Proteins|")

def test_compartments_of_reaction():
    """ compartments_of_reaction
     Arguments
    ---------
    default_compartment: The default compartment, as determined by
        the function (default_compartment), which currently is set
        to 'CCO-CYTOSOL.
    sides: The slots of the reaction to consider. The default value
        is '(LEFT RIGHT).
    rxn: An instance of the class |Reactions|.

    """
    result = ec.compartments_of_reaction(rxn)
    assert_instance_list(result, 'CCO')
    result = ec.compartments_of_reaction(rxn, sides=['LEFT'])
    assert_instance_list(result, 'CCO')

def test_complex():
    """ complex
     Arguments
    ---------
    frame: An instance of the class |Proteins|.

    """
    result = ec.complex(protein)
    assert result in [True, False, None]

def test_operon_of_gene():
    """ operon_of_gene
     Arguments
    ---------
    gene: An instance of class |Genes|.

    """
    result = ec.operon_of_gene(gene)
    assert_instance_list(result, "|Transcription-Units|")

def test_modified_and_unmodified_forms():
    """ modified_and_unmodified_forms
     Arguments
    ---------
    protein: An instance of the class |Proteins|.

    """
    result = ec.modified_and_unmodified_forms(protein)
    assert_instance_list(result, "|Proteins|")

# def test_rxns_adjacent_in_pwy():
#     """     Arguments
#     ---------
#     rxn1: An instance of the class |Reactions|.
#     rxn2: An instance of the class |Reactions|.
#     pwy: An instance of the class |Pathways|.

#     """
#     pwy = 'TRPSYN-PWY'
#     rxn1 = 'PRTRANS-RXN'
#     rxn2 = 'PRAISOM-RXN'
#     result = ec.rxns_adjacent_in_pwy(rxn1, rxn2, pwy)
#     assert result in [True, False, None]

def test_adjacent_genes():
    """ adjacent_genes
     Arguments
    ---------
    g2: An instance of class |Genes|.
    g1: An instance of class |Genes|.

    """
    g1 = 'EG11024'
    g2 = 'EG11025'
    g3 = 'EG11218'
    assert ec.adjacent_genes(g1,g2)
    assert not ec.adjacent_genes(g2,g3)

def test_regulation_frame_transcription_units():
    """ regulation_frame_transcription_units
     Arguments
    ---------
    reg_frame: An instance of the class |Regulation-of-
        Transcription|.

    """
    result = ec.regulation_frame_transcription_units(reg_frame)
    assert_instance_list(result, "|Transcription-Units|")

def test_transcription_unit_regulation_frames():
    """ transcription_unit_regulation_frames
     Arguments
    ---------
    tu: An instance of the class |Transcription-Units|.

    """
    result = ec.transcription_unit_regulation_frames(tu)
    assert_instance_list(result, "|Regulation|")

def test_dna_binding_sites_of_protein():
    """ dna_binding_sites_of_protein
     Arguments
    ---------
    tf: An instance of the class |Proteins|.
    all_forms: When true, then return the DNA binding sites of
        modified forms and subunits of tf as well.

    """
    result = ec.dna_binding_sites_of_protein(tf)
    assert_instance_list(result, "|DNA-Binding-Sites|")
    result = ec.dna_binding_sites_of_protein(tf, all_forms=True)
    assert_instance_list(result, "|DNA-Binding-Sites|")

def test_gene():
    """ gene
     Arguments
    ---------
    frame: A frame.

    """
    assert ec.gene(gene)
    assert not ec.gene(rxn)

def test_regulator_of_type():
    """ regulator_of_type
     Arguments
    ---------
    protein: An instance frame of class |Proteins|.
    class: A subclass of |Regulation|.

    """
    regulation_subclass = '|Regulation-of-Transcription|'
    assert ec.regulator_of_type(tf, regulation_subclass)
    assert not ec.regulator_of_type('TRYPSYN', regulation_subclass)

def test_binding_site_promoters():
    """ binding_site_promoters
     Arguments
    ---------
    tu: An instance of the class |DNA-Binding-Sites|.

    """
    result = ec.binding_site_promoters(tu)
    assert_instance_list(result, "|Promoters|")

def test_inhibition():
    """ inhibition
     Arguments
    ---------
    reg_frame: An instance of class |Regulation|.

    """
    result = ec.inhibition(reg_frame)
    assert result in [True, False, None]

def test_protein_or_rna_containers_of():
    """ protein_or_rna_containers_of
     Arguments
    ---------
    protein: An instance of the class |Proteins|.
    exclude_self: If true, then protein will not be included in the
        return value.

    """
    result = ec.protein_or_rna_containers_of(protein)
    assert_instance_list(result, "|Proteins|")
    result = ec.protein_or_rna_containers_of(protein, True)
    assert_instance_list(result, "|Proteins|")


regulated_enzrxn = 'PEPCARBOX-ENZRXN'
def test_enzrxn_activators():
    """ enzrxn_activators
     Arguments
    ---------
    phys_relevant_only: If true, then only return activators that
        are associated with |Regulation| instances that have the
        'physiologically-relevant? slot set to non-None.
    er: An instance of the class |Enzymatic-Reactions|.

    """
    # The return should be a list of children of |Chemicals| 
    # but we have no way to test subclasses such as |Alcohols| for
    # descent from |Chemicals| now
    chemicals = ec['|Chemicals|']
    result = ec.enzrxn_activators(regulated_enzrxn)
    assert_frame_list(result)
    for f in result:
        assert f in chemicals
    result = ec.enzrxn_activators(regulated_enzrxn,
                                  phys_relevant_only=True)
    for f in result:
        assert f in chemicals
    assert_frame_list(result)

def test_genes_regulated_by_protein():
    """ genes_regulated_by_protein
     Arguments
    ---------
    protein: An instance of the class |Proteins|.

    """
    result = ec.genes_regulated_by_protein(tf)
    assert_instance_list(result, "|Genes|")

def test_rxn_specific_form_of_rxn():
    """ rxn_specific_form_of_rxn
     Arguments
    ---------
    generic_rxn: A child of the class |Reactions|.
    specific_rxn: A child of the class |Reactions|.

    """
    specific_rxn_true = 'ALCOHOL-DEHYDROG-RXN'
    specific_rxn_false = 'TRYPSYN-RXN'
    assert ec.rxn_specific_form_of_rxn(generic_rxn, specific_rxn_true)
    assert not ec.rxn_specific_form_of_rxn(generic_rxn, specific_rxn_false)

def test_activation():
    """ activation
     Arguments
    ---------
    reg_frame: An instance of class |Regulation|.

    """
    result = ec.activation(reg_frame)
    assert result in [True, False, None]

def test_all_transporters_across():
    """ all_transporters_across
     Arguments
    ---------
    membranes: Either 'all' or a list of instances of the class 'CCO-MEMBRANE.
        Defaults to 'all'.
    method: Either 'location' or 'reaction-compartments'. 'location'
        will check the 'locations slot, while 'reaction-
        compartments' will examine the compartments of reaction
        substrates. Default value is 'location'.

    """
    result = ec.all_transporters_across()
    assert_instance_list(result, "|Proteins|")
    result = ec.all_transporters_across(membranes=['CCO-GOLGI-MEM',
                                                   'CCO-OUTER_MEM'],
                                        method='reaction-compartments')
    assert_instance_list(result, "|Proteins|")

def test_terminator():
    """ terminator
     Arguments
    ---------
    frame: A frame.

    """
    terminator_frame = 'TERM2'
    result = ec.terminator(terminator_frame)
    assert result
    result = ec.terminator(tf)
    assert not result

def test_chromosome_of_operon():
    """ chromosome_of_operon
     Arguments
    ---------
    tu: An instance of the class |Transcription-Units|.

    """
    result = ec.chromosome_of_operon(tu)
    assert_class_instance(result, "|Genetic-Elements|")

def test_rxns_catalyzed_by_complex():
    """ rxns_catalyzed_by_complex
     Arguments
    ---------
    rxns: A list of instances of the class |Reactions|. Defaults to
        the result of (all-rxns 'enzyme').

    """
    rxn_list = ['TRYPSYN-RXN', 'OROPRIBTRANS-RXN']
    result = ec.rxns_catalyzed_by_complex()
    assert_instance_list(result, "|Reactions|")
    result = ec.rxns_catalyzed_by_complex(rxns=rxn_list)
    assert_instance_list(result, "|Reactions|")

def test_transcription_units_of_promoter():
    """ transcription_units_of_promoter
     Arguments
    ---------
    promoter: An instance of class |Promoters|.

    """
    result = ec.transcription_units_of_promoter(promoter)
    assert_instance_list(result, "|Transcription-Units|")

def test_pathway_allows_enzrxn():
    """ pathway_allows_enzrxn
     Arguments
    ---------
    enzrxn: An instance of the class |Enzymatic-Reactions|.
    rxn: An instance of the class |Reactions|.
    single_species: An instance of the class |Organisms| If set,
        then enzrxn has the further stricture that it must be an
        enzymatic reaction present in the organism specified by the
        value passed to single_species.
    pwy: An instance of the class |Pathways|.

    """
    result = ec.pathway_allows_enzrxn(enzrxn, rxn, pwy)
    assert result in [True, False, None]

def test_rna_coding_gene():
    """ rna_coding_gene
     Arguments
    ---------
    gene: An instance of the class |Genes|.

    """
    result = ec.rna_coding_gene(gene)
    assert result in [True, False, None]

def test_connecting_genes():
    """ connecting_genes
     Arguments
    ---------
    g2: An instance of class |Genes|.
    g1: An instance of class |Genes|.

    """
    g1 = 'EG11024'
    g2 = 'EG11025'
    result = ec.connecting_genes(g2, g1)
    if result is None:
        pass
    else:
        assert_frame(result[0])
        assert_is_instance(result[1], int)
        assert_is_instance(result[2], int)

def test_nonspecific_forms_of_rxn():
    """ nonspecific_forms_of_rxn
     Arguments
    ---------
    rxn: An instance of the class |Reactions|.

    """
    result = ec.nonspecific_forms_of_rxn('ALCOHOL-DEHYDROG-RXN')
    assert_instance_list(result, "|Reactions|")

def test_species_of_protein():
    """ species_of_protein
     Arguments
    ---------
    p: A list of instances of the class |Proteins|.

    """
    result = ec.species_of_protein(protein)
    if isinstance(result, str):
        pass
    if isinstance(result, list):
        assert_class_instance(result[0], "|Organisms|")
    else:
        assert_class_instance(result, "|Organisms|")

def test_activated_or_inhibited_by_compound():
    """ activated_or_inhibited_by_compound
     Arguments
    ---------

    cpds: An instance or list of instances of class |Compounds|.
    mode: Represents the type of regulation. Can take on the values
        of "+", "-", or 'None.
    mechanisms: Keywords from the 'mechanism slot of the
        corresponding sub-class of the class |Regulation|. If non-
        None, only regulation objects with mechanisms in this list
        will be explored for regulated objects.
    phys_relevant: If true, then only return inhibitors that are
        associated with |Regulation| instances that have the
        'physiologically-relevant? slot set to non-None.
    slots: A list of enzymatic reaction slots to

    """
    cpds = ['TRP', 'PPI']
    result = ec.activated_or_inhibited_by_compound(cpds)
    assert_instance_list(result, "|Enzymatic-Reactions|")


def test_find_indexed_frame():
    """ find_indexed_frame """
    # Hard to test without a clear description of the behavior, 
    # but this apparently returns a list whose first entry is TRP
    result = ec.find_indexed_frame('trp','|Compounds|')
    assert result[0] == 'TRP'

# def test_pathways_of_organism_in_meta():
#     # No way to test, undocumented.
#     raise NotImplementedError

# def test_enzymes_of_organism_in_meta():
#     # No way to test, undocumented.
#     raise NotImplementedError

def test_lower_taxa_or_species():
    """ lower_taxa_or_species (incomplete test) """
    assert ec.lower_taxa_or_species('ecoli')
    # expected return: (|strain| |subspecies| |varietas|)

def test_get_reaction_list():
    """ get_reaction_list """
    assert_frame_list(ec.get_reaction_list(pwy))

def test_get_class_all_subs():
    """ get_class_all_subs """
    result = ec.get_class_all_subs('|RNAs|')
    assert_frame_list(result)
    for f in result:
        assert ec.is_class(f)

def test_is_class_all_sub_of():
    """ is_class_all_sub_of """
    c1 = '|D-fructofuranuronate|'
    c2 = 'FRUCTURONATE'
    c3 = '|Genes|'
    assert ec.is_class_all_sub_of(c1,c2) == True
    assert not ec.is_class_all_sub_of(c1, c3)
