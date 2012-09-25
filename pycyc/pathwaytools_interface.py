"""
Expose Pathway Tools Internal LISP functions as methods of an interface class.

Most method docstring text is copied with minor modifications from 
http://bioinformatics.ai.sri.com/ptools/ptools-fns.html and
and http://bioinformatics.ai.sri.com/ptools/gfp.html (accessed July 26, 2012.)
and http://bioinformatics.ai.sri.com/ptools/api/ (accessed September 2012.)

Eli Bogart
July-September 2012

"""

from server_interface import BasicPathwayToolsDBInterface
from errors import PyCycError
from frame import Frame

class PathwayToolsInterface(BasicPathwayToolsDBInterface):
    """ Extend the Pathway Tools interface to provide LISP API functions. """

    # Frame access methods

    def __contains__(self, thing):
        return self.is_coercible_to_frame(thing) or False

    def __getitem__(self, item):
        frame = Frame(item, self)
        if frame.test_frame():
            return frame
        else:
            raise KeyError('No frame %s in database %s' % (item, self))

    # Convenience attributes for accessing classes

    special_attributes = {'genes': '|Genes|',
                          'reactions': '|Reactions|',
                          'compounds': '|Compounds|',
                          'pathways': '|Pathways|'}
    
    def __getattr__(self, attr):
        if attr in self.special_attributes:
            return self.get_class_all_instances(self.special_attributes[attr])
        else:
            raise AttributeError('PathwayToolsInterface instance has no '
                                 'attribute %s' % attr)


    def __dir__(self):
        attributes = (dir(self.__class__) + self.__dict__.keys() + 
                      self.special_attributes.keys())
        return attributes

    # Functions automatically translated from ptools-fns.html
    
    # Each method here should look like this:
    # def translated_function_name(self, descriptive_argument_name_1, 
    #                              descriptive_argument_name_2,
    #                              optional_argument_1=None, 
    #                              keyword_argument_1=None,
    #                              etc...):
    #     """ Summary description.
    #
    #     Detailed docstring copied from the ptools-functions page, with names
    #     changed as necessary.
    #
    #     Original LISP function: 'translated-function-name'
    #
    #     """
    #     return self.call('translated-function-name', 
    #                 descriptive_argument_name_1,
    #                 descriptive_argument_name_2, 
    #                 **{'lisp-name-keyword-arugment-1': keyword_argument_1})
    #    
    # 
    # Where appropriate, self.multiple_value_list_call should be used 
    # instead of self.call. 
    def activated_or_inhibited_by_compound(self, cpds, 
                                           mode=None, 
                                           mechanisms=None, 
                                           phys_relevant=False,
                                           slots=[], 
                                           ):
        """ List enzymatic reactions regulated by compound.

        This function is incompletely described in the LISP API: use at
        your own risk!

        Arguments
        ---------
        cpds: An instance or list of instances of class |Compounds|.
        mode: Represents the type of regulation. Can take on the values
            of "+", "-", or None (default)
        mechanisms: Keywords from the 'mechanism' slot of the
            corresponding sub-class of the class |Regulation|. If non-
            None, only regulation objects with mechanisms in this list
            will be explored for regulated objects.
        phys_relevant: If true, then only return inhibitors that are
            associated with |Regulation| instances that have the
            'physiologically-relevant? slot set to non-None.
        slots: A list of enzymatic reaction slots to [API document cuts off
        here.]

        Returns
        -------
        A list of instances of class |Enzymatic-Reactions|.

        LISP equivalent: activated-or-inhibited-by-compound

        """
        kwargs = {}
        if mechanisms != None:
            kwargs['mechanisms'] = mechanisms
        kwargs['mode'] = mode
        kwargs['slots'] = slots
        kwargs['phys-relevant?'] = phys_relevant
        return self.call("activated-or-inhibited-by-compound", cpds, **kwargs)

    def activation(self, reg_frame):
        """ Determine whether a given regulation frame describes activation.

        Arguments
        ---------
        reg_frame: An instance of class |Regulation|.

        Returns
        -------
        A boolean value.

        LISP equivalent: activation-p

        """
        return self.call("activation-p", reg_frame)

    def adjacent_genes(self, g1, g2):
        """ Determine whether two genes are on the same replicon and adjacent.

        Arguments
        ---------
        g1: An instance of class |Genes|.
        g2: An instance of class |Genes|.

        Returns
        -------
        A boolean value.

        LISP equivalent: adjacent-genes?

        """
        return self.call("adjacent-genes?", g2, g1)

    def all_cofactors(self):
        """ List all cofactors used by enzymes in the current organism.

        Returns
        -------
        All cofactors used by enzymes in the current organism.

        LISP equivalent: all-cofactors ()

        """
        return self.call("all-cofactors")

    def all_direct_forms_of_protein(self, protein):
        """ List direct subunits/containers of protein's forms.

        Given a protein, this function will return all of the directly
        related proteins of its modified and unmodified forms, meaning all
        of their direct subunits and all of their direct containers.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.

        Returns
        -------
        A list of instances of the class |Proteins|.

        LISP equivalent: all-direct-forms-of-protein

        """
        return self.call("all-direct-forms-of-protein", protein)

    def all_enzymes(self, type='chemical-change'):
        """ Return all enzymes of a given type.

        Arguments
        ---------
        type: A type of enzyme, as described in the docstring of
             enzyme(). Defaults to 'chemical-change'.

        Returns
        -------
        A list of instances of class |Proteins|.

        LISP equivalent: all-enzymes

        """
        kwargs = {}
        kwargs['type'] = (':' + type if
                          isinstance(type,str) else type)
        return self.call("all-enzymes", **kwargs)

    def all_forms_of_protein(self, protein):
        """ List all related proteins of all forms of protein.

        Given a protein, this function will return all of the related
        proteins of its modified and unmodified forms, meaning all of
        their subunits and all of their containers. Unlike
        all_direct_forms_of_protein, this function is not limited to
        the direct containers only.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.

        Returns
        -------
        A list of instances of the class |Proteins|.

        LISP equivalent: all-forms-of-protein

        """
        return self.call("all-forms-of-protein", protein)

    def all_genetic_regulation_proteins(self, 
                                        class_='default', 
                                        allow_modified_forms=True):
        """ List proteins involved in a class of regulation.

        Enumerates all proteins that are involved in genetic regulation of
        a particular given class. Optionally, just unmodified forms of the
        proteins may be returned.

        Arguments
        ---------
        class_: The class |Regulation| or a subclass. It defaults to
            |Regulation-of-Transcription-Initiation|.
        allow_modified_forms: A boolean value. If true, modified and
            unmodified forms of the protein are returned. If false, then
            only unmodified forms of the proteins are returned. The
            default value is True.

        Returns
        -------
        A list of protein frames that are involved in the specified form
        of regulation.

        LISP equivalent: all-genetic-regulation-proteins

        """
        kwargs = {}
        if class_ != 'default':
            kwargs['class'] = class_
        kwargs['allow-modified-forms?'] = allow_modified_forms
        return self.call("all-genetic-regulation-proteins", **kwargs)


    def all_modulators(self):
        """ List all modulators (activators and inhibitors) of enzymes.

        Returns
        -------
        All modulators (activators and inhibitors) that enzymes in the
        current organism are sensitive to.

        LISP equivalent: all-modulators ()

        """
        return self.call("all-modulators")

    def all_operons(self):
        """ Enumerates all operons

        Enumerates all operons. In this case, an operon is defined as a
        list of overlapping instances of |Transcription-Units|.

        Arguments
        ---------
        None.

        Returns
        -------
        A list of lists of |Transcription-Units|, where all
        |Transcription-Units| in the list belong to the same operon.

        LISP equivalent: all-operons

        """
        return self.call("all-operons")

    def all_pathways(self, selector='all', base=False):
        """ List pathways. 

        Returns a list of pathway instance frames of a specified type.

        Arguments
        ---------
        selector: Selects whether all pathways, or just small-molecule
            metabolism base pathways. Can take values 'all' or
            'small-molecule'. Defaults to 'all'.
        base: If this argument evaluates to true, only includes base
            pathways. Otherwise, all pathways, including
            superpathways, will be returned. Through the Python
            interface, defaults to False.

        Returns
        -------
        A list of instances of class |Pathways|.

        LISP equivalent: all-pathways

        """
        if isinstance(selector,str):
            selector = ':' + selector
        return self.call("all-pathways", selector, base)

    def all_products_of_gene(self, gene):
        """ List gene products of gene (including modified forms, complexes.)

        Arguments
        ---------
        gene: An instance of class genes.

        Returns
        -------
        All gene products of gene (including modified forms and
        complexes), including those that are not enzymes.

        LISP equivalent: all-products-of-gene (gene)

        """
        return self.call("all-products-of-gene", gene)

    def all_protein_complexes(self, filter=True):
        """ Enumerates different types of protein complexes.

        Arguments
        ---------
        filter: The type of protein complexes to return. The argument
            must be one of the following:
            'hetero': Return all heteromultimers.
            'homo': Return all homomultimers.
            'all' or True: Return all protein complexes (the default.)

        Returns
        -------
        A list of protein complex frames.

        LISP equivalent: all-protein-complexes

        """
        kwargs = {}
        kwargs['filter'] = (':' + filter if
                            isinstance(filter,str) else filter)
        return self.call("all-protein-complexes", **kwargs)

    def all_rxns(self, type='metab-smm'):
        """ List all reactions in the database (of a given type.)

        Arguments
        ---------
        type: The type of reaction to return. Defaults to 'metab-smm' .
            The possible values are:
            'all': All reactions.
            'metab-pathways': All reactions found within metabolic
                pathways. Includes reactions that are pathway holes. May
                include a handfull of reactions whose substrates are
                macromolecules, e.g., ACP. Excludes transport reactions.
            'metab-smm': All reactions of small molecule metabolism,
                whether or not they are present in a pathway. Subsumes
                'metab-pathways'.
            'metab-all': All enzyme-catalyzed reactions. Subsumes
                'metab-smm'.
            'enzyme': All enzyme-catalyzed reactions (i.e., instances of
                either EC-Reactions class or Unclassified-Reactions
                class).
            'transport': All transport reactions.
            'small-molecule': All reactions whose substrates are all
                small molecules, as opposed to macromolecules. Excludes
                transport reactions.
            'protein-small-molecule-reaction': One of the substrates of
                the reaction is a macromolecule, and one of the
                substrates of the reaction is a small molecule.
            'protein-reaction': All substrates of the reaction are
                proteins.
            'trna-reaction': One of the substrates of the reaction is a
                tRNA.
            'spontaneous': Spontaneous reactions.
            'non-spontaneous': Non-spontaneous reactions that are likely
                to be enzyme catalyzed. Some reactions will be returned
                for type 'non-spontaneous' that will not be returned by
                'enzyme'.

        Returns
        -------
        A list of reaction frames.

        LISP equivalent: all-rxns

        """
        if isinstance(type,str):
            type = ':' + type
        return self.call("all-rxns", type)

    def all_sigma_factors(self):
        """ Enumerate all RNA polymerase sigma factors.

        Arguments
        ---------
        None.

        Returns
        -------
        A list of all instances of the class |Sigma-Factors|.

        LISP equivalent: all-sigma-factors

        """
        return self.call("all-sigma-factors")

    def all_substrates(self, rxns):
        """ Returns all unique substrates used in a list of reactions.

        Returns all unique substrates used in the reactions specified by
        the argument rxns.

        Arguments
        ---------
        rxns: A list of reaction frames.

        Returns
        -------
        A list of compound frames. There might be strings in the list, as
        the 'left' and 'right' slots of a reaction frame can contain
        strings.

        LISP equivalent: all-substrates

        """
        return self.call("all-substrates", rxns)

    def all_transcription_factors(self, allow_modified_forms=True):
        """ List all transcription factors in the current organism.

        Arguments
        ---------
        allow_modified_forms: if False, we return unmodified forms of the
            proteins only.

        Returns
        -------
        All transcription factors in the current organism (computed as any
        protein that is a substrate in a binding reaction where another
        substrate is a DNA site).

        LISP equivalent: all-transcription-factors (&key (allow-modified-
        forms? t))

        """
        kwargs = {}
        kwargs['allow-modified-forms?'] = allow_modified_forms
        return self.call("all-transcription-factors", **kwargs)

    def all_transported_chemicals(self):
        """ List all chemicals transported by all defined transport reactions.

        Arguments
        ---------
        None.

        Returns
        -------
        All chemicals that are transported by the set of all defined
        transport reactions in the current organism.

        LISP equivalent: all-transported-chemicals ()

        """
        return self.call("all-transported-chemicals")

    def all_transporters(self):
        """ List all transport proteins.

        Arguments
        ---------
        None.

        Returns
        -------
        A list of instances of class |Proteins|.

        LISP equivalent: all-transporters

        """
        return self.call("all-transporters")

    def all_transporters_across(self, membranes='all', method='location'):
        """ Returns a list of transport proteins that transport across o

        Returns a list of transport proteins that transport across one of
        the given membranes.

        Arguments
        ---------   
        membranes: Either 'all' or a list of instances of the
            class 'CCO-MEMBRANE'. Defaults to 'all'.
        method: Either 'location' or 'reaction-compartments'. 'location'
            will check the 'locations slot, while 'reaction-
            compartments' will examine the compartments of reaction
            substrates. Default value is 'location'.


        Returns
        -------
        A list of instances of class |Proteins|.

        LISP equivalent: all-transporters-across

        """
        if isinstance(method,str):
            method = ':' + method
        if isinstance(membranes,str):
            membranes = ':' + membranes
        kwargs = {'membranes': membranes,
                  'method': method}
        return self.call("all-transporters-across", **kwargs)

    def autocatalytic_reactions_of_enzyme(self, prot):
        """ List possibly autocatalytic reactions of protein. 

        Returns a list of reaction frames, where the protein participates
        as a substrate of the reaction, and the reaction has no associated
        Enzymatic Reaction frame. This implies that the protein substrate
        of the reaction might be autocatalyzing the reaction.

        Arguments
        ---------
        prot: An instance frame of class |Proteins|.

        Returns
        -------
        A list of instances of class |Reactions|.

        LISP equivalent: autocatalytic-reactions-of-enzyme

        """
        return self.call("autocatalytic-reactions-of-enzyme", prot)

    def base_components_of_protein(self, p, exclude_small_molecules=True):
        """ List components of protein complex.

        Same as function monomers_of_protein, but also returns
        components of the protein that are RNAs or compounds, not just
        polypeptides.

        Arguments
        ---------
        p: An instance of the class |Proteins|.
        exclude_small_molecules: If None, then small molecule components
            are also returned. Default value is true.

        Returns
        -------
        Two values. The first value is a list of the components, which may
        be instances of the following classes: |Polypeptides|, |RNAs|, and
        |Compounds|. The second value is a list of the corresponding
        coefficients of the components in the first value.

        LISP equivalent: base-components-of-protein

        """
        kwargs = {}
        kwargs['exclude-small-molecules?'] = exclude_small_molecules
        return self.multiple_value_list_call("base-components-of-protein", p,
                                             **kwargs)

    def binding_site_promoters(self, tu):
        """ Returns the promoters of the given DNA binding site.

        Arguments
        ---------
        tu: An instance of the class |DNA-Binding-Sites|.

        Returns
        -------
        A list of instances of class |Promoters|.

        LISP equivalent: binding-site-promoters

        """
        return self.call("binding-site-promoters", tu)

    def binding_site_to_regulators(self, bsite):
        """ List the transcription factors of the given binding site.

        Arguments
        ---------
        bsite: An instance of class |DNA-Binding-Sites|.

        Returns
        -------
        A list of instances of class |Proteins|.

        LISP equivalent: binding-site->regulators

        """
        return self.call("binding-site->regulators", bsite)

    def binding_site_transcription_factors(self, bsite):
        """ List the transcription factors that bind to bsite.

        Arguments
        ---------
        bsite:  An instance of class DNA-Binding-Sites.

        Returns
        -------
        A list of the transcription factors that bind to the DNA binding
        site bsite.

        LISP equivalent: binding-site-transcription-factors (bsite)

        """
        return self.call("binding-site-transcription-factors", bsite)

    def binding_site_transcription_units(self, promoter):
        """ Returns all transcription units of a given binding site.

        Arguments
        ---------
        promoter: An instance of class |DNA-Binding-Sites| or |mRNA-
            Binding-Sites|.

        Returns
        -------
        A list of instances of class |Transcription-Units|.

        LISP equivalent: binding-site-transcription-units

        """
        return self.call("binding-site-transcription-units", promoter)

    def binding_sites_affecting_gene(self, gene):
        """ List binding sites in the same transcription units as the gene.

        Returns all binding sites which are present in the same
        transcription units as the given gene.

        Arguments
        ---------
        gene: An instance of the class |Genes|.

        Returns
        -------
        A list of instances of class |DNA-Binding-Sites|.

        LISP equivalent: binding-sites-affecting-gene

        """
        return self.call("binding-sites-affecting-gene", gene)

    def chromosome_of_gene(self, gene):
        """ Return the chromosome on which gene resides.

        Arguments
        ---------
        gene: An instance of class genes.

        Returns
        -------
        The chromosome on which gene resides.

        LISP equivalent: chromosome-of-gene (gene)

        """
        return self.call("chromosome-of-gene", gene)

    def chromosome_of_object(self, object_):
        """ Given an object, return the replicon where it is located.

        Given an object, the replicon where it is located is returned. If
        there is no associated replicon for the object, None is returned.
        If the object is on more than one replicon, an error is thrown.

        Arguments
        ---------
        object_: An instance of class |All-Genes|, |Transcription-Units|,
            |Promoters|, |Terminators|, |Misc-Features|, or |DNA-
            Binding-Sites|.

        Returns
        -------
        An instance of class |Genetic-Elements|.

        LISP equivalent: chromosome-of-object

        """
        return self.call("chromosome-of-object", object_)

    def chromosome_of_operon(self, tu):
        """ Returns the replicon of the given transcription unit.

        Arguments
        ---------
        tu: An instance of the class |Transcription-Units|.

        Returns
        -------
        An instance of class |Genetic-Elements|.

        LISP equivalent: chromosome-of-operon

        """
        return self.call("chromosome-of-operon", tu)

    def cofactors_and_pgroups_of_enzrxn(self, enzrxn):
        """ List cofactors and prosthetic groups of an enzymatic reaction.

        Arguments
        ---------
        enzrxn: An instance of the class |Enzymatic-Reactions|.

        Returns
        -------
        A list of children of class |Chemicals| or strings, representing
        cofactors and/or prosthetic groups.

        LISP equivalent: cofactors-and-pgroups-of-enzrxn

        """
        return self.call("cofactors-and-pgroups-of-enzrxn", enzrxn)

    def compartments_of_reaction(self, rxn,  
                                 sides=['LEFT', 'RIGHT'],
                                 default_compartment='CCO-CYTOSOL'):
        """ Returns the compartments associated with the given reaction.

        Arguments
        ---------
        rxn: An instance of the class |Reactions|.
        sides: The slots of the reaction to consider. The default value
            is ['LEFT', 'RIGHT'].
        default_compartment: The default compartment. Defaults to 
            'CCO-CYTOSOL'.

        Returns
        -------
        A list of children of the class 'CCO.

        LISP equivalent: compartments-of-reaction

        """
        kwargs = {}
        kwargs['sides'] = sides
        kwargs['default-compartment'] = default_compartment
        return self.call("compartments-of-reaction", rxn, **kwargs)

    def complex(self, frame):
        """ Determine whether a protein frame is a protein complex.

        Arguments
        ---------
        frame: An instance of the class |Proteins|.

        Returns
        -------
        A boolean value.

        LISP equivalent: complex-p

        """
        return self.call("complex-p", frame)

    def components_of_protein(self, p, coefficient=None):
        """ Compute the components of a protein p.

        This function recursively computes all of the components of a
        protein p.

        Arguments
        ---------
        p: An instance of class proteins.
        coefficient: Undocumented.

        Returns 
        ------- 
        components: a list of components (including p)
        coefficients: a list of the coefficients of those components.

        LISP equivalent: components-of-protein (p &optional coefficient)

        """
        return self.call("components-of-protein", p, coefficient)

    def compounds_of_pathway(self, pwy):
        """ List all substrates (reactants and products) of reactions in pwy.

        Arguments
        ---------
        pwy: A pathway frame.

        Returns
        -------
        All substrates (meaning reactants and products) of reactions of
        pwy, with duplicates removed.

        LISP equivalent: compounds-of-pathway (pwy)

        """
        return self.call("compounds-of-pathway", pwy)

    def connecting_genes(self, g1, g2):
        """ If two genes are on a replicon, return it and their positions.

        Given two genes, if they are on the same replicon, return the
        replicon, and their numerical ordering among the genes of the
        replicon, with the smaller number first.

        Arguments
        ---------
        g1: An instance of class |Genes|.
        g2: An instance of class |Genes|.

        Returns
        -------
        Three values: the common replicon, the smaller numerical index,
        and the larger numerical index.

        LISP equivalent: connecting-genes

        """
        return self.multiple_value_list_call("connecting-genes", g2, g1)

    def containers_of(self, protein, exclude_self=None):
        """ List all containers of protein.

        Arguments
        ---------
        protein: An instance of class proteins.
        exclude_self: Boolean indicating whether protein should be
            included in the returned list.

        Returns
        -------
        All containers of protein, including itself (unless exclude_self
        is True).  Follows the Component-Of slot.

        LISP equivalent: containers-of (protein &optional exclude-self?)

        """
        return self.call("containers-of", protein, exclude_self)

    def containing_chromosome(self, site):
        """ Get the regulon corresponding to a site.

        Given a site (whether a DNA binding site, a promoter, a gene, or a
        terminator) along a transcription unit, returns the correspodning
        regulon.

        Arguments
        ---------
        site: An instance of class |Transcription-Units|, |mRNA-Binding-
            Sites|, |DNA-Binding-Sites|, |Promoters|, |Genes|, or
            |Terminators|.

        Returns
        -------
        An instance of class |Genetic-Elements|.

        LISP equivalent: containing-chromosome

        """
        return self.call("containing-chromosome", site)

    def containing_tus(self, site):
        """ List the transcription units containing site.

        Arguments
        ---------
        site: An instance of any of the following classes: Genes,
            Promoters, Terminators, DNA-Binding-sites, mRNA-Binding-sites,
            Transcription-Units.

        Returns
        -------
        A list of transcription units containing site.

        LISP equivalent: containing-tus (site)

        """
        return self.call("containing-tus", site)

    def cotranscribed_genes(self, gene):
        """ List genes cotranscribed with the given gene.

        Returns all co-transcribed genes (i.e., genes which are a part of
        one or more of the same transcription units) of the given gene.

        Arguments
        ---------
        gene: An instance of class |Genes|.

        Returns
        -------
        A list of instances of class |Genes|.

        LISP equivalent: cotranscribed-genes

        """
        return self.call("cotranscribed-genes", gene)

    def create_instance_w_generated_id(self, class_):
        """ Creates a new instance of class.

        Creates a new instance of class whose frame ID is generated
        automatically.

        Arguments
        ---------
        class_: The parent class of the instance to be created.

        Returns
        -------
        The created frame (either its frame name or the frame object).

        Example
        -------
        >>> ecocyc.create_instance_w_generated_id('|Genes|')
        'G0-1'

        LISP equivalent: create-instance-w-generated-id (class)

        """
        return self.call("create-instance-w-generated-id", class_)

    def direct_activators(self, item):
        """ List entities directly activating an item that can be regulated. 

        Arguments
        ---------
        item: any frame that can be regulated, e.g. an instance of
            Enzymatic-Reactions, Transcription-Units, Promoters, Terminators,
            etc.

        Returns
        -------
        The list of entities (e.g. protein, RNA, small molecule) that
        directly activate item.

        LISP equivalent: direct-activators (item)

        """
        return self.call("direct-activators", item)

    def direct_inhibitors(self, item):
        """ List entities directly inhibiting an item that can be regulated. 

        Arguments
        ---------
        item: any frame that can be regulated, e.g. an instance of
            Enzymatic-Reactions, Transcription-Units, Promoters, Terminators,
            etc.

        Returns
        -------
        The list of entities (e.g. protein, RNA, small molecule) that
        directly inhibit item.

        LISP equivalent: direct-inhibitors (item)

        """
        return self.call("direct-inhibitors", item)

    def direct_regulators(self, x, filter_fn=None):
        """ List direct regulators of an object.

        Returns all regulators that are connected to a regulated
        object by a single regulation object.

        Arguments
        ---------
        x: A frame.
        filter_fn: A LISP predicate used to filter the regulation
            objects used to find the regulators (default: no filtering.)

        Returns
        -------
        A list of frames that regulate x.

        LISP equivalent: direct-regulators

        """
        # In practice, filter-fn is a keyword argument
        kwargs = {}
        if filter_fn:
            kwargs['filter-fn'] = filter_fn
        return self.call("direct-regulators", x, **kwargs)

    def dna_binding_site(self, frame):
        """ Test if a frame is an instance of class |DNA-Binding-Sites|.

        Arguments
        ---------
        frame: A frame.

        Returns
        -------
        A boolean value.

        LISP equivalent: dna-binding-site-p

        """
        return self.call("dna-binding-site-p", frame)

    def dna_binding_sites_of_protein(self, tf, all_forms=False):
        """ List all a transcription factor's DNA binding sites.

        Arguments
        ---------
        tf: An instance of the class |Proteins|.
        all_forms: When true, then return the DNA binding sites of
            modified forms and subunits of tf as well (default false).

        Returns
        -------
        A list of instances of the class |DNA-Binding-Sites|.

        LISP equivalent: DNA-binding-sites-of-protein

        """
        kwargs = {}
        kwargs['all-forms?'] = all_forms
        return self.call("DNA-binding-sites-of-protein", tf, **kwargs)

    def enzrxn_activators(self, er, phys_relevant_only=False):
        """ List activators of an enzymatic reaction. 

        Returns the list of activators (generally small molecules) of the
        enzymatic reaction frame.

        Arguments
        ---------
        er: An instance of the class |Enzymatic-Reactions|.
        phys_relevant_only: If true, then only return activators that
            are associated with |Regulation| instances that have the
            'physiologically-relevant?' slot set to non-None.

        Returns
        -------
        A list of children of the class |Chemicals|.

        LISP equivalent: enzrxn-activators

        """
        return self.call("enzrxn-activators", er, phys_relevant_only)

    def enzrxn_inhibitors(self, er, phys_relevant_only=False):
        """ List inhibitors of an enzymatic reaction frame.

        Returns the list of inhibitors (generally small molecules) of the
        enzymatic reaction frame er.

        Arguments
        ---------
        er: An instance of the class |Enzymatic-Reactions|.
        phys_relevant_only: If true, then only return inhibitors that
            are associated with |Regulation| instances that have the
            'physiologically-relevant? slot set to non-None. Through
            the Python interface, defaults to False.

        Returns
        -------
        A list of children of the class |Chemicals|.

        LISP equivalent: enzrxn-inhibitors

        """
        return self.call("enzrxn-inhibitors", er, phys_relevant_only)

    def enzyme(self, protein, type='any'):
        """ Check whether protein is an enzyme (of a particular type.) 

        Arguments
        ---------
        protein: An instance of class proteins.
        type: string, one of four values specifying what is meant by an enzyme:
            'any' Any protein that catalyzes a reaction is considered an
                enzyme.
            'chemical-change' -- If the reactants and products of the
                catalyzed reaction differ (not just cellular location), the
                protein is considered an enzyme.
            'small-molecule' -- If the reactants of the catalyzed
                reaction differ and are all small molecules, the
                protein is considered an enzyme.
            'non-transport' -- Exclude all proteins that only catalyze
                transport reactions.

        Returns 
        ------- 
        If 'type' is 'chemical-change', a list of products and
        reactants of a catalyzed reaction is returned if the protein
        is an enzyme by that criterion, and None is returned if the
        protein is not.

        For other values of 'type', the function returns True if
        protein is an enzyme of the specified type, else None.

        LISP equivalent: enzyme? (protein &optional  (type :any))

        """
        type = ':' + type
        return self.call("enzyme?", protein, type)

    def enzyme_activity_name(self, enzyme, reaction=None):
        """ Return a string describing the name of an enzyme.

        Arguments
        ---------
        enzyme: An instance of class Proteins.
        reaction: An instance of class reactions.

        Returns
        -------
        A string describing the name of an enzyme -- an enzyme-activity
        name.  If reaction is non-None, then we compute the name of the
        enzyme in the context of that reaction, otherwise we return the
        full enzyme name, which describes all of the reactions that it
        catalyzes.  If reaction is None, this function is equivalent to
        full-enzyme-name.

        LISP equivalent: enzyme-activity-name (enzyme &optional reaction)

        """
        return self.call("enzyme-activity-name", enzyme, reaction)

    def enzymes_of_gene(self, gene):
        """ List all enzymes coded for by gene.

        Arguments
        ---------
        gene: An instance of class genes.

        Returns
        -------
        All enzymes coded for by gene (meaning both a monomer product of
        gene, and complexes that contain such monomers), meaning proteins
        that catalyze some reaction.

        LISP equivalent: enzymes-of-gene (gene)

        """
        return self.call("enzymes-of-gene", gene)

    def enzymes_of_organism_in_meta(self, gb_tax_id):
        """ Call the LISP function enzymes-of-organism-in-meta.

        Argument name taken from perlcyc; no other documentation available.

        """
        return self.call('enzymes-of-organism-in-meta', gb_tax_id)

    def enzymes_of_pathway(self, pwy):
        """ List all enzymes that catalyze a reaction in the pathway pwy.

        Arguments
        ---------
        pwy: A pathway frame.

        Returns 
        ------- 
        All enzymes that catalyze a reaction in the pathway pwy.

        LISP equivalent: enzymes-of-pathway (pwy)

        """
        return self.call("enzymes-of-pathway", pwy)

    def enzymes_of_reaction(self, rxn):
        """ Return all enzymes that catalyze rxn.

        Arguments
        ---------
        rxn: A reaction frame.

        Returns
        -------
        All enzymes that catalyze the reaction rxn.

        LISP equivalent: enzymes-of-reaction (rxn)

        """
        return self.call("enzymes-of-reaction", rxn)

    def fequal(self, frame1, frame2):
        """ Test whether frame1 and frame2 are the same frame.

        Arguments
        ---------
        frame1: A frame id or object.
        frame2: A frame id or object.

        Returns
        -------
        True if frame1 and frame2 are the same frame.  Must be used as
        the test function for operations that compare frames.

        LISP equivalent: fequal (frame1 frame2)

        """
        return self.call("fequal", frame1, frame2)

    def find_indexed_frame(self, datum, class_):
        """ Call the LISP function find-indexed-frame.

        No other documentation available. 

        """

        return self.multiple_value_list_call('find-indexed-frame', datum, 
                                             class_)

    def full_enzyme_name(self, enzyme):
        """ Return a string describing the full name of the enzyme.

        Arguments
        ---------
        enzyme: An instance of class Proteins.

        Returns
        -------
        A string describing the full name of enzyme, as the concatenation
        of the common name of the protein followed by the common names of
        its enzymatic reactions.  Note that two enzrxns for the same
        enzyme could have the same common name, so we avoid including the
        same name twice.

        LISP equivalent: full-enzyme-name (enzyme)

        """
        return self.call("full-enzyme-name", enzyme)

    def gene(self, frame):
        """ Determine if the given frame is a gene.

        Arguments
        ---------
        frame: A frame.

        Returns
        -------
        A boolean value.

        LISP equivalent: gene-p

        """
        return self.call("gene-p", frame)

    def gene_clusters(self, genes, max_gap=10):
        """ Group genes by proximity.

        Groups together genes based on whether each gene is a gene
        neighbor with other genes.

        Arguments
        ---------
        genes: A list of instances of class |Genes|.
        max_gap: An integer representing the number of genes any pair
            from genes can be from one another. Default value is 10.

        Returns
        -------
        A list of lists, where the first element of each sub-list is a
        gene from genes, and the rest of the list are all of the gene
        neighbors of the first gene.

        LISP equivalent: gene-clusters

        """
        return self.call("gene-clusters", genes, max_gap)

    def gene_transcription_units(self, gene):
        """ Return all the transcription units which contain the given gene.

        Arguments
        ---------
        gene: An instance of class |Genes|.

        Returns
        -------
        A list of instances of class |Transcription-Units|.

        LISP equivalent: gene-transcription-units

        """
        return self.call("gene-transcription-units", gene)

    def genes_in_same_operon(self, gene):
        """ Given a gene, return all other genes in the same operon.

        Arguments
        ---------
        gene: An instance of class |Genes|.

        Returns
        -------
        A list of instances of class |Genes|.

        LISP equivalent: genes-in-same-operon

        """
        return self.call("genes-in-same-operon", gene)

    def genes_of_pathway(self, pwy):
        """ List all genes that code for enzymes/subunits in pwy.

        Arguments
        ---------
        pwy: A pathway frame.

        Returns
        -------
        All genes that code for enzymes (or subunits of enzymes) that
        catalyze a reaction in the pathway pwy.  The returned list does
        not include genes coding for proteins that are substrates in the
        pathway (if any).

        LISP equivalent: genes-of-pathway (pwy)

        """
        return self.call("genes-of-pathway", pwy)

    def genes_of_protein(self, p):
        """ List all genes that code for protein p and its subunits.

        Arguments
        ---------
        p: An instance of class proteins.

        Returns
        -------
        All genes that code for protein p and all of the subunits of p.

        LISP equivalent: genes-of-protein (p)

        """
        return self.call("genes-of-protein", p)

    def genes_of_proteins(self, p):
        """ List genes encoding proteins and their subunits.

        The same as genes_of_protein, except that it takes a list of
        proteins and returns a set of genes.

        Arguments
        ---------
        p: A list of instances of the class |Proteins|.

        Returns
        -------
        A list of instances of the class |Genes|.

        LISP equivalent: genes-of-proteins

        """
        return self.call("genes-of-proteins", p)

    def genes_of_reaction(self, rxn):
        """ Return all genes encoding enzymes/subunits that catalyze rxn.

        Arguments
        ---------
        rxn: A reaction frame.

        Returns
        -------
        All genes that code for enzymes (or subunits of enzyme complexes)
        that catalyze the reaction rxn.  If multiple enzymes catalyze rxn,
        genes encoding all of the enzymes are returned.

        LISP equivalent: genes-of-reaction (rxn)

        """
        return self.call("genes-of-reaction", rxn)

    def genes_regulated_by_gene(self, gene):
        """ List genes whose transcription or translation this gene regulates.

        Arguments
        ---------
        gene: An instance of class genes.

        Returns
        -------
        A list of genes whose transcription or translation is regulated by
        the product of gene.

        LISP equivalent: genes-regulated-by-gene (gene)

        """
        return self.call("genes-regulated-by-gene", gene)

    def genes_regulated_by_protein(self, protein):
        """ List genes for which protein or modified forms act as a regulator.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.

        Returns
        -------
        A list of instances of the class |Genes|.

        LISP equivalent: genes-regulated-by-protein

        """
        return self.call("genes-regulated-by-protein", protein)

    def genes_regulating_gene(self, gene):
        """ List genes that regulate transcription/translation of this gene.

        Arguments
        ---------
        gene: An instance of class genes.

        Returns
        -------
        A list of genes that regulate transcription or translation of the
        argument gene.

        LISP equivalent: genes-regulating-gene (gene)

        """
        return self.call("genes-regulating-gene", gene)

    def get_class_all_instances(self, class_):
        """ Get all frames that are direct or indirect instances of class.

        Arguments
        ---------
        class_: A class frame.

        Returns
        -------
        A list of all frames that are direct or indirect instances of
        class_.

        Example
        -------
        >>> ecocyc.get_class_all_instances('|Reactions|')
        ['RXN0-6505', 'RXN0-3002', 'RXN0-3141', 'R131-RXN', ... ]

        LISP equivalent: get-class-all-instances (class)

        """
        return self.call("get-class-all-instances", class_)


    def get_class_all_subs(self, class_):
        """ List all subclasses of a class.

        Arguments
        ---------
        class_: a class frame

        Returns 
        -------
        A list of class frames.

        """
        return self.call("get-class-all-subs", class_)

    def get_frame_labeled(self, label):
        """ Get frames whose name or synonym matches label.

        Arguments
        ---------
        label: a string

        Returns
        -------
        A list of frames whose name or synonym matches label.

        Example
        -------
        >>> ecocyc.get_frame_labeled('pyrophosphate')
        ['PPI']

        LISP equivalent: get-frame-labeled (label)

        """
        # Here we want to enclose the argument in double quotes rather than
        # quote it, as a symbol, with ', though the latter approach seems 
        # to be equivalent in cursory testing
        return self.evaluate('(get-frame-labeled "%s")' % label)

    def get_frame_slots(self, frame):
        """ List the slots associated with frame.

        Arguments
        ---------
        frame: A frame.
        
        Returns
        -------
        A list of the slots associated with frame.

        The LISP equivalent offers optional arguments to distinguish
        local from template slots, but these are not yet implemented here.

        """
        return self.call('get-frame-slots',frame)

    def get_gene_sequence(self, gene):
        """ Get the nucleotide sequence of a gene if available.

        Arguments
        ---------
        gene: An instance of class genes.

        Returns
        -------
        A string containing the nucleotide sequence for the gene (if
        available).

        LISP equivalent: get-gene-sequence (gene)

        """
        return self.call("get-gene-sequence", gene)

    def get_instance_all_types(self, instance):
        """ List the all-types of an instance.

        Arguments
        ---------
        instance: A frame.

        Returns
        -------
        A list of the all-types of instance.

        """
        return self.call('get-instance-all-types',instance)

    def get_instance_direct_types(self, instance):
        """ List the direct types of an instance.

        Arguments
        ---------
        instance: A frame.

        Returns
        -------
        A list of the direct types of instance.

        """
        return self.call('get-instance-direct-types',instance)

    def get_name_string(self, item, strip_html=None, direction=None, 
                        rxn_eqn_as_name=True, name_slot=None):
        """ Return a string name for the frame item.

        This is the standard function for computing the biological name of
        a frame.

        Exact operation depends on what type of entity item is.  For
        genes, for example, the fn simply retrieves the value of the
        common-name slot, or else the frame ID, if the common-name slot
        has no value.  For other entities, name generation is more
        complex.

        Arguments
        ---------
        item: The frame whose name is to be computed.
        rxn_eqn_as_name: Influences how names are computed for reaction
            frames.  When True, the reaction equation is returned; else the
            EC number (or frame name if there is no EC number) is returned.
        direction: Influences how names are computed for reaction
            frames.  Possible values: 'r2l' or 'l2r'.
        name_slot: Overrides use of common-name slot as the slot from
            which the name is retrieved.
        strip_html: When True causes this fn to strip HTML tags from the
            resulting names, which by default will sometimes contain HTML
            tags such as italics or subscripts.

        Returns
        -------
        A string name for the frame item.

        LISP equivalent: get-name-string (item &key (rxn-eqn-as-name? t)
        direction name-slot strip-html?)

        """
        if direction:
            direction = ':' + direction
        kwargs = {}
        kwargs['strip-html?'] = strip_html
        kwargs['direction'] = direction
        kwargs['rxn-eqn-as-name?'] = rxn_eqn_as_name
        kwargs['name-slot'] = name_slot
        return self.call("get-name-string", item, **kwargs)

    def get_predecessors(self, rxn, pwy):
        """ List all reactions that are direct predecessors of rxn in pwy.

        Arguments
        ---------
        rxn: A reaction frame.
        pwy: A pathway frame.

        Returns
        -------
        A list of all reactions that are direct predecessors of rxn in
        pwy.

        LISP equivalent: get-predecessors (rxn pwy)

        """
        return self.call("get-predecessors", rxn, pwy)

    def get_reaction_list(self, pwy):
        """ List reactions in pathway. 

        Arguments
        ---------
        pwy: A pathway frame.

        Returns
        -------
        A list of the reactions in the pathway.

        """
        return self.call("get-reaction-list", pwy)

    def get_slot_value(self, frame, slot):
        """ Get the first value of slot of frame.

        Arguments
        ---------
        frame: A frame id or object.
        slot: A slot name (symbol).

        Returns
        -------
        The first value of slot of frame.

        Example
        -------
        >>> ecocyc.get_slot_value('trp','common-name')
        'L-tryptophan'

        LISP equivalent: get-slot-value (frame slot)

        """
        return self.call("get-slot-value", frame, slot)

    def get_slot_values(self, frame, slot):
        """ Get all values of slot of frame.

        Arguments
        ---------
        frame: A frame id or object.
        slot: A slot name (symbol).

        Returns
        -------
        A list of all values of slot of frame.

        Example
        -------
        >>> ecocyc.get_slot_values('trp','synonyms')
        ['trp', 'W', 'tryptacin', 'trofan', 'tryptophan', 
        '2-amino-3-indolylpropanic acid']
        
        LISP equivalent: get-slot-values (frame slot)

        """
        return self.call("get-slot-values", frame, slot)

    def get_successors(self, rxn, pwy):
        """ List all reactions that are direcft successors of rxn in pwy.

        Arguments
        ---------
        rxn: A reaction frame.
        pwy: A pathway frame.

        Returns
        -------
        A list of all reactions that are direct successors of rxn in pwy.

        LISP equivalent: get-successors (rxn pwy)

        """
        return self.call("get-successors", rxn, pwy)

    def homomultimeric_containers_of(self, protein, exclude_self=None):
        """ List all homomultimeric containers of protein.

        This function is the same as the function containers_of, except
        that it only includes containers that are homomultimers.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.
        exclude_self: If true, then protein will not be included in the
            return value.

        Returns
        -------
        A list of instances of the class |Proteins|.

        LISP equivalent: homomultimeric-containers-of

        """
        return self.call("homomultimeric-containers-of", protein, exclude_self)

    def inhibition(self, reg_frame):
        """ Determines if a given regulation frame describes inhibition.

        Arguments
        ---------
        reg_frame: An instance of class |Regulation|.

        Returns
        -------
        A boolean value.

        LISP equivalent: inhibition-p

        """
        return self.call("inhibition-p", reg_frame)

    def instance_all_instance_of(self, instance, class_):
        """ Test whether instance is a direct or indirect child of class.

        Arguments
        ---------
        instance: An instance frame.
        class_: A class frame.

        Returns
        -------
        True if instance is a direct or an indirect child of class, 
        else None.

        Example
        -------
        >>> ecocyc.instance_all_instance_of('PM338','|Promoters|')
        True

        LISP equivalent: instance-all-instance-of-p (instance class)

        """
        return self.call("instance-all-instance-of-p", instance, class_)

    def is_class(self, frame):
        """ Determine whether a frame is a class frame. 

        Arguments
        ---------
        frame: A frame.

        Returns
        -------
        True if the frame is a class frame in the current kb, otherwise None.

        """
        return self.call('class-p',frame)

    def is_class_all_type_of(self, class_, instance):
        """ Test whether instance is an all-instance of class_.

        Arguments
        ---------
        class_: String describing a class, e.g., '|Genes|'
        instance: A frame. 

        Returns
        -------
        True if instance is an all-instance of class_, otherwise None.

        """
        return self.call('class-all-type-of-p', class_, instance)

    def is_class_all_sub_of(self, class_, superclass):
        """ Test whether one class is a subclass of another.

        Arguments
        ---------
        class_: A class frame
        superclass: A class frame

        Returns
        -------
        True if class_ is a subclass of superclass, otherwise None.

        """
        return self.call('class-all-sub-of-p', class_, superclass)

    def is_coercible_to_frame(self, thing):
        """ Test whether a string may be interpreted as a frame ID.

        Arguments
        ---------
        thing: a string 
        
        Returns
        -------
        True if the argument may be used as as a frame ID in this KB, else 
            None.

        """
        # Docstring is my best guess at the behavior.
        return self.call('coercible-to-frame-p',thing)

    def is_instance(self, frame):
        """ Determine whether a frame is an instance frame. 

        Arguments
        ---------
        frame: A frame.

        Returns
        ------- 
        True if the frame is an instance frame in the current kb,
        otherwise None.

        """
        return self.call('instance-p',frame)

    def is_reaction_of_type(self, rxn, type_):
        """ Return True if rxn is of type type_.

        Arguments
        ---------
        rxn: A reaction frame.
        type_: Either 'small-molecule' (meaning all substrates must be
            small molecules rather than macromolecules) or 'transport'
            (meaning the reaction is already a child of class
            Transport-Reactions, or that some substrate of the reaction
            appears as both a reactant and a product of the reaction, but on
            different sides of the reaction).

        Returns
        -------
        True if rxn is of type type_.

        LISP equivalent: reaction-type? (rxn type)

        """
        type_ = ':' + type_
        return self.call("reaction-type?", rxn, type_)
    
    def is_slot(self, frame, slot):
        """ Determine whether a frame has a slot. 

        Arguments
        ---------
        frame: A frame.
        slot: A string labeling a possible slot of frame.

        Returns
        -------
        None if 'frame' has no slot 'slot', otherwise True or, in some
        cases, a nonempty list of the slot label and its values. 

        """
        return self.call('slot-p',frame,slot)

    def leader_peptide(self, protein):
        """ Determine whether the given protein is a leader peptide.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.

        Returns
        -------
        A boolean value.

        LISP equivalent: leader-peptide?

        """
        return self.call("leader-peptide?", protein)

    def lower_taxa_or_species(self, org_frame):
        """ Call the LISP function lower-taxa-or-species-p.

        No further documentation available. Argument name taken from PerlCyc.

        """
        return self.call("lower-taxa-or-species-p", org_frame)

    def member_slot_value(self, frame, slot, value='', symbol=None):
        """ Test whether value is one of the values of slot of frame.

        Arguments
        ---------
        frame: A frame id or object.
        slot: A slot name (symbol).
        value: A slot value such as an integer or string (ignored if symbol
            is specified).
        symbol: A string to be interpreted as a frame or other LISP symbol,
            or None (the default).

        One of 'value' or 'symbol' must be given.

        Returns
        -------
        None if 'value' is not one of the values of slot of frame, otherwise
        'value'.

        Example
        -------
        >>> ecocyc.member_slot_value('trp','common-name','L-tryptophan')
        'L-tryptophan'

        LISP equivalent: member-slot-value-p (frame slot value)

        """
        if not (value or symbol):
            raise ValueError('Either a literal value or a symbol must be '
                             'given as an argument to member_slot_value.')

        if symbol is not None: 
            value = "'" + symbol
        elif isinstance(value, str):
            value = '"%s"' % value

        return self.evaluate("(member-slot-value-p '%s '%s %s)" %
                             (frame, slot, value))

    def modified_and_unmodified_forms(self, protein):
        """ Returns all modified and unmodified forms of the given protein.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.

        Returns
        -------
        A list of instances of the class |Proteins|.

        LISP equivalent: modified-and-unmodified-forms

        """
        return self.call("modified-and-unmodified-forms", protein)
    
    def modified_containers(self, protein):
        """ List all containers of a protein and all its modified forms.

        Arguments
        ---------
        protein: An instance of class proteins.

        Returns
        -------
        Both all containers of a protein (meaning protein complexes of
        which this protein is a part) including itself, and all modified
        forms of a protein.

        LISP equivalent: modified-containers (protein)

        """
        return self.call("modified-containers", protein)

    def modified_forms(self, protein, exclude_self=None):
        """ List the modified forms of protein.

        Arguments
        ---------
        protein: An instance of class proteins.
        exclude_self: Boolean indicating whether protein should be
            included in the returned list.

        Returns
        -------
        A list of the modified forms of protein (including protein, unless
        exclude_self is True).  Follows the Modified-Forms slot.

        LISP equivalent: modified-forms (protein &optional exclude-self?)

        """
        return self.call("modified-forms", protein, exclude_self)

    def monomers_of_protein(self, p):
        """ List the monomers of protein p.

        Arguments
        ---------
        p: An instance of class proteins.

        Returns
        -------
        A list of the monomers (instances of class polypeptides) that are
        subunits of protein p.  If p is itself a monomer, the function
        returns p.

        LISP equivalent: monomers-of-protein (p)

        """
        return self.call("monomers-of-protein", p)

    def neighboring_genes(self, g1, g2, n=10):
        """ Test proximity of two genes.

        Given two genes, this predicate determines if the two genes are
        "neighbors", or within a certain number of genes from one another
        along the replicon.

        Arguments
        ---------
        g1: An instance of class |Genes|.
        g2: An instance of class |Genes|.
        n: An integer representing the number of genes g1 and g2 can be
            from one another. Default value is 10.

        Returns
        -------
        A boolean value.

        LISP equivalent: neighboring-genes-p

        """
        return self.call("neighboring-genes-p", g2, g1, n)

    def next_gene_on_replicon(self, g):
        """ Return the next gene on the replicon.

        Arguments
        ---------
        g: An instance of class |Genes|.

        Returns
        -------
        Returns two values. The first value is the next gene, or None if
        there is not a next gene (i.e., the gene is at the end of a linear
        replicon). The second value is ':last' if the gene is the last gene
        on a linear replicon.

        LISP equivalent: next-gene-on-replicon

        """
        return self.multiple_value_list_call("next-gene-on-replicon", g)

    def noncontiguous_pathway(self, pwy):
        """ Determine if pathway contains multiple connected components.

        See function pathway-components for more explanation.

        Arguments
        ---------
        pwy: An instance of the class |Pathways|.

        Returns
        -------
        A boolean value.

        LISP equivalent: noncontiguous-pathway-p

        """
        return self.call("noncontiguous-pathway-p", pwy)

    def nonspecific_forms_of_rxn(self, rxn):
        """ Return all of the generic forms of the given specific reaction.

        Not every reaction will necessarily have a generic form.

        Arguments
        ---------
        rxn: An instance of the class |Reactions|.

        Returns
        -------
        A list of children of the class |Reactions|.

        LISP equivalent: nonspecific-forms-of-rxn

        """
        return self.call("nonspecific-forms-of-rxn", rxn)

    def nucleotide_to_protein_sequence(self, sequence, code_num=1):
        """ Translate nucleotide sequence to protein sequence. 

        Arguments
        ---------
        sequence: A string representing a nucleotide sequence (such as
            might be returned by get-gene-sequence).
        code_num: The number of the genetic code to use (an integer from 1
            to 15).  Defaults to the standard genetic code.

        Returns
        -------
        A string containing the translated sequence.

        LISP equivalent: nucleotide->protein-sequence (sequence &optional
        (code-num 1))

        """
        query = '(nucleotide->protein-sequence "%s" %d)' % (sequence,
                                                            code_num)
        return self.evaluate(query)

    def operon_of_gene(self, gene):
        """ List transcription units forming the gene's containing operon.

        Arguments
        ---------
        gene: An instance of class |Genes|.

        Returns
        -------
        A list of instances of class |Transcription-Units|.

        LISP equivalent: operon-of-gene

        """
        return self.call("operon-of-gene", gene)

    def pathway_allows_enzrxn(self, pwy, rxn, enzrxn, single_species=None):
        """ Test wheter pathway allows enzrxn to catalyze rxn.

        A predicate which returns a true value if the given pathway
        allows the given enzymatic reaction to catalyze the given
        reaction.  Certain pathways have a list of enzymatic reactions
        that are known not to catalyze certain reactions. See the
        documentation of "slot-unit 'enzyme-use" for more
        information.

        Arguments
        ---------
        pwy: An instance of the class |Pathways|.
        rxn: An instance of the class |Reactions|.
        enzrxn: An instance of the class |Enzymatic-Reactions|.
        single_species: An instance of the class |Organisms|. If set,
            then enzrxn has the further stricture that it must be an
            enzymatic reaction present in the organism specified by the
            value passed to single_species.

        Returns
        -------
        A boolean value.

        LISP equivalent: pathway-allows-enzrxn?

        """
        return self.call("pathway-allows-enzrxn?", enzrxn, rxn, pwy, 
                         single_species)

    def pathway_components(self, pwy, rxn_list=None, 
                           pred_list=None):
        """ List all of the connected components of a pathway.

        A connected component of a pathway is a set of reactions in
        the pathway such that for all reactions R1 in the connected
        component, a predecessor relationship holds between R1 and
        some other reaction R2 in the connected component, and each
        connected component is of maximal size.  Every pathway will
        have from 1 to N connected components, where N is the number
        of reactions in the pathway.  Most pathways have one connected
        component, but not all.

        Arguments
        ---------
        pwy: An instance of the class |Pathways|, which is not a super-
            pathway (i.e., does not have any entries in its 'sub-
            pathways' slot).
        rxn_list: The list of reactions to use as the starting list of
            connected component clusters. If None, defaults to (get-slot-values
            pwy 'reaction-list).
        pred_list: The list of reaction predecessors to iterate from in
            order to cluster the reactions in rxn_list. If None, defaults to
            (get-slot-values pwy 'predecessors). 

        Returns
        -------
        A list of lists and three integers. This is in disagreement
        with the LISP API document, which indicates three values
        should be returned, specifically:

        - the connected components as a list of lists of the form
            [[r1, r2, r3], [r4, r5], [r6, r7, r8]], where each sub-
            list contains all reactions in one connected component
        -  the number of connected components, 
        -  the length of the reaction list.

        LISP equivalent: pathway-components

        """
        args = []
        if rxn_list is not None: 
            args.append(rxn_list)
            if pred_list:
                args.append(pred_list)
        elif pred_list is not None:
            rxn_list = self['pwy']['reaction-list']
            args = [rxn_list, pred_list]
        return self.multiple_value_list_call("pathway-components", pwy, 
                                             *args)

    def pathway_hole(self, rxn, hole_if_any_gene_without_position=False):
        """ Determine if the current reaction is considered a 'pathway hole'.

        True if the current reaction is considered to be a 'pathway
        hole', or without an associated enzyme.

        Arguments
        ---------
        rxn: An instance of the class |Reactions|.
        hole_if_any_gene_without_position: If true, then genes without
            specified coordinates for the current organism's genome are
            not counted when determining the status of the reaction. Defaults
            to False (this may differ from the original LISP version, whose
            default behavior is not documented.)

        Returns
        -------
        A boolean value.

        LISP equivalent: pathway-hole-p

        """
        kwargs = {}
        h = hole_if_any_gene_without_position
        kwargs['hole-if-any-gene-without-position?'] = h
        return self.call("pathway-hole-p", rxn, **kwargs)

    def pathways_of_enzrxn(self, enzrxn, include_super_pwys=False):
        """ List pathways in which the given enzymatic reaction participates.

        Arguments
        ---------
        enzrxn: An instance of the class |Enzymatic-Reactions|.
        include_super_pwys: If true, then not only will the direct
            pathways in which enzrxn is associated in be returned, but
            also any enclosing super-pathways. If enzrxn is associated
            with a reaction that is directly associated with a super-
            pathway, then the function might return super-pathways
            even if this option is false. When called through the
            Python interface, defaults to False.

        Returns
        -------
        A list of instances of class |Pathways|.

        LISP equivalent: pathways-of-enzrxn

        """
        kwargs = {}
        kwargs['include-super-pwys?'] = include_super_pwys
        return self.call("pathways-of-enzrxn", enzrxn, **kwargs)

    def pathways_of_gene(self, gene):
        """ List pathways containing reactions catalyzed by products of gene.

        Arguments
        ---------
        gene: An instance of class genes.

        Returns
        -------
        All pathways containing reactions that are catalyzed by proteins
        that are produts of gene (meaning both monomer products of gene,
        and complexes that contain such monomers).

        LISP equivalent: pathways-of-gene (gene)

        """
        return self.call("pathways-of-gene", gene)

    def phantom_gene(self, gene):
        """ Determine if the given gene is a phantom gene.

        Arguments
        ---------
        gene: An instance of the class |Genes|.

        Returns
        -------
        A boolean value.

        LISP equivalent: phantom-gene-p

        """
        return self.call("phantom-gene-p", gene)

    def polypeptide_or_homomultimer(self, protein):
        """ Determine if the given protein is a polypeptide or homomultimer.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.

        Returns
        -------
        A boolean value.

        LISP equivalent: polypeptide-or-homomultimer-p

        """
        return self.call("polypeptide-or-homomultimer-p", protein)

    def previous_gene_on_replicon(self, g):
        """ Return the previous gene on the replicon.

        Arguments
        ---------
        g: An instance of class |Genes|.

        Returns
        -------
        Returns two values. The first value is the previous gene, or None
        if there is not a previous gene (i.e., the gene is at the
        beginning of a linear replicon). The second value is ':first' if
        the gene is the first gene on a linear replicon.

        LISP equivalent: previous-gene-on-replicon

        """
        return self.multiple_value_list_call("previous-gene-on-replicon", g)

    def promoter_binding_sites(self, promoter):
        """ List the binding sites associated with the given promoter.

        Returns all of the binding sites associated with the given
        promoter, across multiple transcription units.

        Arguments
        ---------
        promoter: An instance of class |Promoters|.

        Returns
        -------
        A list of instances of class |DNA-Binding-Sites|.

        LISP equivalent: promoter-binding-sites

        """
        return self.call("promoter-binding-sites", promoter)

    def protein(self, frame):
        """ Determine whether the given frame is a protein.

        Arguments
        ---------
        frame: An instance of the class |Proteins|.

        Returns
        -------
        A boolean value.

        LISP equivalent: protein-p

        """
        return self.call("protein-p", frame)

    def protein_coding_gene(self, gene):
        """ Determine if a gene encodes a protein (rather than an RNA).

        Arguments
        ---------
        gene: An instance of the class |Genes|.

        Returns
        -------
        A boolean value.

        LISP equivalent: protein-coding-gene?

        """
        return self.call("protein-coding-gene?", gene)

    def protein_or_rna_containers_of(self, protein, exclude_self=False):
        """ Find protein[-RNA] complexes containing a given protein.

        This function is the same as the function containers_of, except
        that it only includes containers that are instances of either
        class |Protein-Complexes|, or class |Protein-RNA-Complexes|.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.
        exclude_self: If true, then protein will not be included in the
            return value.

        Returns
        -------
        A list of instances of the class |Proteins|.

        LISP equivalent: protein-or-rna-containers-of

        """
        return self.call("protein-or-rna-containers-of", protein, exclude_self)

    def pseudo_gene(self, gene):
        """ Determine if the given gene is a pseudo-gene.

        Arguments
        ---------
        gene: An instance of the class |Genes|.

        Returns
        -------
        A boolean value.

        LISP equivalent: pseudo-gene-p

        """
        return self.call("pseudo-gene-p", gene)

    def pwys_of_organism_in_meta(self, gb_tax_id):
        """ Call the LISP function pwys-of-organism-in-meta.

        Argument name taken from perlcyc; no other documentation available.

        """
        return self.call('pwys-of-organism-in-meta', gb_tax_id)

    def reaction_reactants_and_products(self, rxn, direction='L2R', pwy=None):
        """ Return the reactants of rxn and the products of rxn.

        The reactants and products are those determined according to
        either the direction of rxn in pwy, or the direction specified by
        the direction arg (which should be either 'L2R' or 'R2L').  In
        other words, either the pwy or the direction arg should be
        specified, but not both.

        Arguments
        ---------
        rxn: A reaction frame.
        direction: string specifying the sense of the reaction. One of:
            'L2R' (left to right, default) 
            'R2L' (right to left)
            None (only if a pathway frame is specified)
        pwy: A pathway frame or None (default), ignored if direction is given.
            If rxn is not in the specified pathway, Pathway Tools returns no 
            reactants or products, and a PyCycError is raised.


        Returns
        -------
        reactants of rxn : list
        products of rxn : list

        LISP equivalent: reaction-reactants-and-products (rxn &key pwy
        direction)

        """
        kwargs = {}
        if direction:
            kwargs['direction'] = ':' + direction
        elif pwy:
            kwargs['pwy'] = pwy
        else:
            raise PyCycError('Either a reaction direction or a pathway '
                             'containing the reaction must be given.')
        
        r =  self.multiple_value_list_call("reaction-reactants-and-products", 
                                           rxn, **kwargs)
        if r == [None]:
            raise PyCycError('No reactants or products found; '
                             'if a pathway was specified, check that it '
                             'includes the reaction.')
        return r

    def reaction_type(self, rxn):
        """ Returns a keyword describing the type of reaction.

        Arguments
        ---------
        rxn: An instance of the class |Reactions|.

        Returns
        -------
        One of the following:
            'trna-reaction': At least one substrate is a tRNA.
            'small-molecule': All substrates are small molecules, or
                small-molecule classes.
            'transport': A substrate is marked with different
                compartment annotations in the 'left' and 'right' slots.
            'protein-small-molecule-reaction': At least one substrate is
                a protein and at least one is a small molecule.
            'protein-reaction': All substrates are proteins.
            'null-reaction': No substrates or reactants are specified.
            'other': None of the preceding cases apply.

        LISP equivalent: reaction-type

        """
        symbol = self.call("reaction-type", rxn)
        symbol = symbol.lower()
        if symbol.startswith(':'):
            symbol = symbol[1:]
        return symbol

    def reaction_without_sequenced_enzyme(self, rxn, complete=None):
        """ Test if a reaction has genes without associated sequence info.

        Arguments
        ---------
        rxn: An instance of the class |Reactions|.
        complete: If true, the predicate will return true when there
            is any associated gene without a sequence. If None (the
            default), the predicate will return true when all
            associated genes are without a sequence.

        Returns
        -------
        A boolean value.

        LISP equivalent: reaction-without-sequenced-enzyme-p

        """
        # Note this is listed in the API as rxn-without-sequenced-enzyme-p.
        return self.call("reaction-without-sequenced-enzyme-p", rxn, complete)

    def reactions_of_compound(self, cpd):
        """ List the reactions in which cpd occurs as reactant or product.

        Arguments
        ---------
        cpd: An instance of class Chemicals.

        Returns
        -------
        The reactions in which cpd occurs as a reactant or a product.

        LISP equivalent: reactions-of-compound (cpd)

        """
        return self.call("reactions-of-compound", cpd)

    def reactions_of_enzyme(self, e):
        """ List all reactions enzyme e is linked to via enzymatic reactions.

        Arguments
        ---------
        e: An instance of class proteins.

        Returns
        -------
        All reactions that enzyme e is linked to via enzymatic reactions.

        LISP equivalent: reactions-of-enzyme (e)

        """
        return self.call("reactions-of-enzyme", e)

    def reactions_of_gene(self, gene):
        """ List reactions catalyzed by products of gene.

        Arguments
        ---------
        gene: An instance of class genes.

        Returns
        -------
        All reactions catalyzed by proteins that are products of gene
        (meaning both monomer products of gene, and complexes that contain
        such monomers).

        LISP equivalent: reactions-of-gene (gene)

        """
        return self.call("reactions-of-gene", gene)

    def reactions_of_protein(self, p):
        """ List all reactions catalyzed by p or by subunits of p.

        Arguments
        ---------
        protein: An instance of class proteins.

        Returns
        -------
        All reactions catalyzed by p or by subunits of p.

        LISP equivalent: reactions-of-protein (p)

        """
        return self.call("reactions-of-protein", p)

    def reduce_modified_proteins(self, prots, debind=False):
        """ List unmodified forms of proteins without duplicates.

        Given a list of proteins, the function converts all of the
        proteins to their unmodified form, and then removes any duplicates
        from the subsequent list.

        Arguments
        ---------
        prots: A list of instances of the class |Proteins|.
        debind: When non-None, the proteins are further simplified by
            obtaining the unbound form of the protein, if it is bound
            to a small molecule. Through the Python interface,
            defaults to False.

        Returns
        -------
        A list of instances of the class |Proteins|.

        LISP equivalent: reduce-modified-proteins

        """
        kwargs = {}
        kwargs['debind?'] = debind
        return self.call("reduce-modified-proteins", prots, **kwargs)

    def regulation_frame_transcription_units(self, reg_frame):
        """ Find trans. units with promoters/terminators regulated by frame.

        Given a regulation object, return the transcription units when
        one of the regulated entities is a promoter or terminator of
        the transcription unit.

        Arguments
        ---------
        reg_frame: An instance of the class |Regulation-of-
            Transcription|.

        Returns
        -------
        A list of instances of the class |Transcription-Units|.

        LISP equivalent: regulation-frame-transcription-units

        """
        return self.call("regulation-frame-transcription-units", reg_frame)

    def regulator_of_type(self, protein, class_):
        """ Determine if a protein is a regulator of the specified class.

        Arguments
        ---------
        protein: An instance frame of class |Proteins|.
        class_: A subclass of |Regulation|.

        Returns
        -------
        A boolean value.

        LISP equivalent: regulator-of-type?

        """
        return self.call("regulator-of-type?", protein, class_)

    def regulator_proteins_of_transcription_unit(self, tu):
        """ Return all proteins that bind to binding sites within TU.

        Arguments
        ---------
        TU: An instance of class Transcription-Units.

        Returns
        -------
        All proteins that bind to binding sites within TU.

        LISP equivalent: regulator-proteins-of-transcription-unit (tu)

        """
        return self.call("regulator-proteins-of-transcription-unit", tu)

    def regulators_of_gene(self, gene, by_function=None):
        """ List proteins that are regulators of the given gene.

        This information taken from the regulators-of-gene-transcription
        entry in the API file; that function seems not to exist, but this
        does, and has no API file entry, so I suspect a typo.

        Arguments
        ---------
        gene: An instance of the class |Genes|.
        by_function: If true, then return two values: a list of
            activator proteins and a list of inhibitor proteins.

        Returns
        -------
        A list of regulators of the gene. If by_function is non-
        None, then two values are returned. The first value is a list of
        activators, and the second value is a list of inhibitor
        proteins. Note that, contrary to the API document, the return values 
        need not be instances of the class |Proteins|, but may include, for 
        example, class frames. 

        LISP equivalent: regulators-of-gene-transcription

        """
        # Note the regulators-of-gene-transcription entry does not 
        # describe by-function? as a keyword argument, but in practice 
        # it appears to be one.
        fn = "regulators-of-gene"
        if by_function:
            kwargs = {'by-function?': by_function}
            a, i = self.multiple_value_list_call(fn, gene, **kwargs)
            if not a: # () translated to None
                a = []
            if not i:
                i = []
            return [a,i] 
        else:
            return self.call(fn, gene)

    def regulators_of_operon_transcription(self, operon_list, 
                                           by_function=False):
        """ Returns a list of transcription factors of an operon.

        Arguments
        ---------
        operon_list: A list of instances of the class |Transcription-
            Units|.
        by_function: If True, then return two values: a list of
            activator proteins and a list of inhibitor proteins.
            Through the Python interface, defaults to False.

        Returns
        -------
        A list of instances of class |Proteins| (or other
        regulators). If the modified form of the protein is the
        transcription factor, then that is the protein returned.

        LISP equivalent: regulators-of-operon-transcription

        """
        if by_function:
            f = "regulators-of-operon-transcription"
            return self.multiple_value_list_call(f,
                             operon_list, by_function)
        else:
            return self.call("regulators-of-operon-transcription", 
                             operon_list, by_function)

    def regulon_of_protein(self, protein):
        """ Return transcription units regulated by any form of protein.

        Arguments
        ---------
        protein: An instance frame of class |Proteins|.

        Returns
        -------
        A list of instances of the class |Transcription-Units|.

        LISP equivalent: regulon-of-protein

        """
        return self.call("regulon-of-protein", protein)

    def rna_coding_gene(self, gene):
        """ A predicate that determines if the given gene encodes an RNA

        Arguments
        ---------
        gene: An instance of the class |Genes|.

        Returns
        -------
        A boolean value.

        LISP equivalent: rna-coding-gene?

        """
        return self.call("rna-coding-gene?", gene)

    def rxn_in_compartment(self, rxn, compartments, default_ok=False, 
                           pwy=None, loose=False):
        """ Check whether rxn is present in one of a list of compartments.

        Arguments
        ---------
        rxn: An instance of the class |Reactions|.
        compartments: A list of cellular compartments, as defined in the
            Cellular Components Ontology. See frame CCO.
        default_ok: If true, then we return true if the reaction has
            no associated compartment information, or one of its
            associated locations is a super-class of one of the
            members of the compartments argument. 
        pwy: If supplied, the search for associated enzymes of the
            argument rxn is limited to the given child of |Pathways|.
        loose: If true, then the compartments 'CCO-CYTOPLASM' and 'CCO-
            CYTOSOL' are treated as being the same compartment. Through
            the Python interface, defaults to False.

        Returns
        -------
        Experiments suggest this function returns an entry from the
        argument list of compartments or a false value. It could
        return also True or a list of multiple entries from the
        compartments argument.

        LISP equivalent: rxn-in-compartment-p

        """
        kwargs = {}
        kwargs['default-ok?'] = default_ok
        if pwy:
            kwargs['pwy'] = pwy
        kwargs['loose?'] = loose
        return self.call("rxn-in-compartment-p", rxn, compartments, **kwargs)

    def rxn_present(self, rxn):
        """ Determine if there is evidence for occurence of a reaction.

        True if there is evidence for the occurrence of the given
        reaction in the current PGDB.

        Arguments
        ---------
        rxn: An instance of the class |Reactions|.

        Returns
        -------
        A boolean value.

        LISP equivalent: rxn-present?

        """
        return self.call("rxn-present?", rxn)

    def rxn_specific_form_of_rxn(self, specific_rxn, generic_rxn):
        """ Test whether a generic reaction generalizes a specific reaction.

        True if the given generic reaction is a generalized form of
        the given specific reaction.

        Arguments
        ---------
        specific_rxn: A child of the class |Reactions|.
        generic_rxn: A child of the class |Reactions|.

        Returns
        -------
        A boolean value.

        LISP equivalent: rxn-specific-form-of-rxn-p

        """
        return self.call("rxn-specific-form-of-rxn-p", 
                         generic_rxn, specific_rxn)

    def rxn_w_isozymes(self, rxn):
        """ Determine whether a reaction has isozymes.

        True if a given reaction has any associated isozymes (distinct
        proteins or protein classes that catalyze the same reaction).

        Arguments
        ---------
        rxn: An instance of the class |Reactions|.

        Returns
        -------
        A boolean value.

        LISP equivalent: rxn-w-isozymes-p

        """
        return self.call("rxn-w-isozymes-p", rxn)

    def rxns_catalyzed_by_complex(self, rxns='default'):
        """ List reactions catalyzed by a protein complex.

        Enumerates all reactions catalyzed by an enzyme that is a protein
        complex.

        Arguments
        ---------
        rxns: A list of instances of the class |Reactions|. Defaults to
            the result of all_rxns(type='enzyme').

        Returns
        -------
        A list of instances of the class |Reactions| with a protein
        complex as an enzyme.

        LISP equivalent: rxns-catalyzed-by-complex

        """
        kwargs = {}
        if rxns != 'default':
            kwargs['rxns'] = rxns
        return self.call("rxns-catalyzed-by-complex", **kwargs)

    def rxns_w_isozymes(self, rxns='default'):
        """ Enumerate all reactions that have isozymes.

        Enumerate all reactions that have isozymes (distinct proteins or
        protein classes that catalyze the same reaction).

        Arguments
        ---------
        rxns: A list of instances of the class |Reactions|. Defaults to
            the result of all_rxns(type='enzyme').

        Returns
        -------
        A list of instances of the class |Reactions| with isozymes.

        LISP equivalent: rxns-w-isozymes

        """
        kwargs = {}
        if rxns != 'default':
            kwargs['rxns'] = rxns
        return self.call("rxns-w-isozymes", **kwargs)

    def small_molecule_cplxes_of_prot(self, protein):
        """ List forms of protein that are complexes with small molecules.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.

        Returns
        -------
        A list of instances of the class |Proteins|.

        LISP equivalent: small-molecule-cplxes-of-prot

        """
        return self.call("small-molecule-cplxes-of-prot", protein)

    def species_of_protein(self, p):
        """ Get the associated species for the given protein.

        This documentation differs slightly from that presented in the
        LISP API document and reflects empirical testing.

        Arguments
        ---------
        p: Instance or list of instances of the class |Proteins|.

        Returns
        -------
        An instance of the class |Organisms|, list of one or more 
        such instances, or a string.

        LISP equivalent: species-of-protein

        """
        return self.call("species-of-protein", p)

    def specific_forms_of_rxn(self, rxn):
        """ List the specific forms of the given generic reaction.

        Not every reaction will necessarily have a specific form.

        Arguments
        ---------
        rxn: A child of the class |Reactions|.

        Returns
        -------
        A list of instances of the class |Reactions|.

        LISP equivalent: specific-forms-of-rxn

        """
        return self.call("specific-forms-of-rxn", rxn)

    def substrate_of_generic_rxn(self, cpd, rxn):
        """ Determine if a parent of cpd is a substrate of a generic reaction.

        True if a parent of the given compound is a substrate of the
        given generic reaction.

        Arguments
        ---------
        cpd: An instance of class |Compounds|.
        rxn: An instance of class |Reactions|.

        Returns
        -------
        A boolean value.

        LISP equivalent: substrate-of-generic-rxn?

        """
        return self.call("substrate-of-generic-rxn?", cpd, rxn)

    def substrates_of_pathway(self, pwy):
        """ List reactants and products of pwy.

        Arguments
        ---------
        pwy: A pathway frame.

        Returns 
        ------- 
        reactants: list of all reactant compounds (compounds occurring on
            the left side of some reaction in pwy).
        proper_reactants: list of proper reactants (the subset of
            reactants that are not also products).
        products: list of all products.  
        proper_products: The list of proper products.

        This function does take into account the proper direction of each
        reaction within pwy.

        LISP equivalent: substrates-of-pathway (pwy)

        """
        return self.multiple_value_list_call("substrates-of-pathway", pwy)

    def substrates_of_reaction(self, rxn):
        """ Return all substrates of the reaction rxn.

        Arguments
        ---------
        rxn: A reaction frame.

        Returns
        -------
        All substrates of the reaction rxn, meaning the union of the
        reactants and products of rxn.

        LISP equivalent: substrates-of-reaction (rxn)

        """
        return self.call("substrates-of-reaction", rxn)

    def terminator(self, frame):
        """ Test whether frame is an instance of |Terminators|.

        Arguments
        ---------
        frame: A frame.

        Returns
        -------
        A boolean value.

        LISP equivalent: terminatorp

        """
        return self.call("terminatorp", frame)

    def terminators_affecting_gene(self, gene):
        """ List terminators in the same TU as gene, upstream of it.

        Arguments
        ---------
        gene: An instance of class genes.

        Returns
        -------
        A list of terminators that are in the same transcription-unit as
        gene and are upstream of gene.  Anything that regulates these
        terminators (i.e. by attenuation) can be assumed to regulate
        transcription of gene.

        LISP equivalent: terminators-affecting-gene (gene)

        """
        return self.call("terminators-affecting-gene", gene)

    def tfs_bound_to_compound(self, cpd, include_inactive=False):
        """ Get complexes that are transcription factors when bound to cpd.

        Returns a list of protein complexes that, when bound to the given
        compound, act as a transcription factor.

        Arguments
        ---------
        cpd: An instance of class |Compounds|.
        include_inactive: If true, then the inactive form of the
            protein is also checked. See the function
            transcription_factor for more information. Through the
            Python interface, defaults to False.

        Returns
        -------
        A list of instances of class |Proteins|.

        LISP equivalent: tfs-bound-to-compound

        """
        kwargs = {}
        kwargs['include-inactive?'] = include_inactive
        return self.call("tfs-bound-to-compound", cpd, **kwargs)

    def top_containers(self, protein):
        """ List top-level containers of a protein.

        Arguments
        ---------
        protein: An instance of class proteins.

        Returns
        -------
        Those containers of a protein (meaning protein complexes of which
        this protein is a part) that themselves have no containers.

        LISP equivalent: top-containers (protein)

        """
        return self.call("top-containers", protein)

    def transcription_factor(self, protein):
        """ Check whether protein is a transcription factor.

        A protein is a transcription factor if it is a reactant of a
        binding-reaction, and another reactant of the same binding
        reaction is a DNA-Binding-Site.

        Arguments
        ---------
        protein: An instance of class proteins.

        Returns 
        ------- 
        True if protein is a transcription factor in the current organism,
        else None.

        LISP equivalent: transcription-factor? (protein)

        """
        return self.call("transcription-factor?", protein)

    def transcription_factor_active_forms(self, tfs):
        """ Find all active forms of the transcription factor.

        For a given transcription factor, find all active forms (i.e, form
        of the protein that regulates) of the transcription factor.

        Arguments
        ---------
        tfs: An instance of the class |Proteins|.

        Returns
        -------
        A list of instances of the class |Proteins|.

        LISP equivalent: transcription-factor-active-forms

        """
        return self.call("transcription-factor-active-forms", tfs)

    def transcription_factor_ligands(self, tfs, mode='both'):
        """ Get transcription factor ligands for one or more tfs.

        Arguments
        ---------
        tfs: An instance or a list of instances of the class |Proteins|.
            If tfs is not the active form, then the active form is
            determined automatically.

        mode: One of the following values: 'activator', 'inhibitor', or
            'both'. In the Python interface, defaults to 'both'.

        Returns
        -------
        A list of instances or subclasses of the class |Chemicals| or
        strings.

        LISP equivalent: transcription-factor-ligands

        """
        # Note that if tfs is not specified, the LISP function returns
        # what is presumably a list of all ligands of all transcription 
        # factors, but as this is undocumented I don't want to support it 
        # based on a guess at what it does; similarly the default mode
        # is likely 'both' but I won't count on that.
        if isinstance(mode,str):
            mode = ':' + mode
        kwargs = {'tfs': tfs, 'mode': mode}
        return self.call("transcription-factor-ligands", **kwargs)

    def transcription_unit_activation_frames(self, tu):
        """ List regulation frames that activate the transcription unit.

        Arguments
        ---------
        tu: An instance of the class |Transcription-Units|.

        Returns
        -------
        A list of instances of the class |Regulation|.

        LISP equivalent: transcription-unit-activation-frames

        """
        return self.call("transcription-unit-activation-frames", tu)

    def transcription_unit_activators(self, tu):
        """ List regulators activating transcription/translation of tu.

        Arguments
        ---------
        tu: An instance of class Transcription-Units.

        Returns
        -------
         A list of the transcription factors or other regulators
        (e.g. small RNAs) that activate transcription or translation of
        the transcription unit tu.  This does not include regulators of
        attenuation (because they regulate only the portion of a
        transcription unit downstream of the terminator).

        LISP equivalent: transcription-unit-activators (tu)

        """
        return self.call("transcription-unit-activators", tu)

    def transcription_unit_all_components(self, tu):
        """ List binding sites/promoters/genes/terminators of trans. unit.

        Returns all components (binding sites, promoters, genes,
        terminators) of the given transcription unit.

        Arguments
        ---------
        tu: An instance of the class |Transcription-Units|.

        Returns
        -------
        A list of instances of class |Transcription-Units|, |mRNA-Binding-
        Sites|, |DNA-Binding-Sites|, |Promoters|, |Genes|, or
        |Terminators|.

        LISP equivalent: transcription-unit-all-components

        """
        return self.call("transcription-unit-all-components", tu)

    def transcription_unit_binding_sites(self, tu):
        """ List the DNA binding sites within transcription unit tu.

        Arguments
        ---------
        tu: An instance of class Transcription-Units.

        Returns
        -------
        A list of the DNA binding sites within the transcription unit tu.

        LISP equivalent: transcription-unit-binding-sites (tu)

        """
        return self.call("transcription-unit-binding-sites", tu)

    def transcription_unit_first_gene(self, tu):
        """ Returns the first gene of the given transcription unit.

        Arguments
        ---------
        tu: An instance of the class |Transcription-Units|.

        Returns
        -------
        An instance of class |Genes|.

        LISP equivalent: transcription-unit-first-gene

        """
        return self.call("transcription-unit-first-gene", tu)

    def transcription_unit_genes(self, tu):
        """ List the genes within the transcription unit tu.

        Arguments
        ---------
        tu: An instance of class Transcription-Units.

        Returns
        -------
        A list of the genes within the transcription unit tu.

        LISP equivalent: transcription-unit-genes (tu)

        """
        return self.call("transcription-unit-genes", tu)

    def transcription_unit_inhibition_frames(self, tu):
        """ List of regulation frames that inhibit the transcription unit.

        Arguments
        ---------
        tu: An instance of the class |Transcription-Units|.

        Returns
        -------
        A list of instances of the class |Regulation|.

        LISP equivalent: transcription-unit-inhibition-frames

        """
        return self.call("transcription-unit-inhibition-frames", tu)

    def transcription_unit_inhibitors(self, tu):
        """ List regulators inhibiting transcription/translation of tu.

        Arguments
        ---------
        tu: An instance of class Transcription-Units.

        Returns
        -------
        A list of the transcription factors or other regulators
        (e.g. small RNAs) that inhibit transcription or translation of the
        transcription unit tu.  This does not include regulators of
        attenuation (because they regulate only the portion of a
        transcription unit downstream of the terminator).

        LISP equivalent: transcription-unit-inhibitors (tu)

        """
        return self.call("transcription-unit-inhibitors", tu)

    def transcription_unit_mrna_binding_sites(self, tu):
        """ List the mRNA binding sites within transcritpion unit tu.

        Arguments
        ---------
        tu: An instance of class Transcription-Units.

        Returns
        -------
        A list of the mRNA binding sites within the transcription unit tu.

        LISP equivalent: transcription-unit-mrna-binding-sites (tu)

        """
        return self.call("transcription-unit-mrna-binding-sites", tu)

    def transcription_unit_promoter(self, tu):
        """ Get the promoter of the transcription unit TU.

        Arguments
        ---------
        tu: An instance of class Transcription-Units.

        Returns
        -------
        The promoter of the transcription unit tu.

        LISP equivalent: transcription-unit-promoter (tu)

        """
        return self.call("transcription-unit-promoter", tu)

    def transcription_unit_regulation_frames(self, tu):
        """ List of regulation frames that regulate the transcription unit.

        Arguments
        ---------
        tu: An instance of the class |Transcription-Units|.

        Returns
        -------
        A list of instances of the class |Regulation|.

        LISP equivalent: transcription-unit-regulation-frames

        """
        return self.call("transcription-unit-regulation-frames", tu)

    def transcription_unit_terminators(self, tu):
        """ List transcription terminator(s) within transcription unit tu.

        Arguments
        ---------
        tu: An instance of class Transcription-Units.

        Returns
        -------
        A list of the transcription terminator(s) within the transcription
        unit tu.

        LISP equivalent: transcription-unit-terminators (tu)

        """
        return self.call("transcription-unit-terminators", tu)

    def transcription_unit_transcription_factors(self, tu):
        """ List the transcription factors controlling transcription unit tu.

        Arguments
        ---------
        tu: An instance of class Transcription-Units.

        Returns
        -------
        A list of the transcription factors that control the transcription
        unit tu.

        LISP equivalent: transcription-unit-transcription-factors (tu)

        """
        return self.call("transcription-unit-transcription-factors", tu)

    def transcription_units_of_gene(self, gene):
        """ List TUs forming the operon containing gene.

        Arguments
        ---------
        gene: An instance of class genes.

        Returns
        -------
        A list of transcription units that form the operon containing
        gene.

        LISP equivalent: transcription-units-of-gene (gene)

        """
        return self.call("transcription-units-of-gene", gene)

    def transcription_units_of_promoter(self, promoter):
        """ Returns all transcription units of a given promoter.

        Arguments
        ---------
        promoter: An instance of class |Promoters|.

        Returns
        -------
        A list of instances of class |Transcription-Units|.

        LISP equivalent: transcription-units-of-promoter

        """
        return self.call("transcription-units-of-promoter", promoter)

    def transcription_units_of_protein(self, protein):
        """ List transcription units activated or inhibited by protein.

        Does not consider other modified or unmodifed forms of protein as
        does regulon-of-protein.

        Arguments
        ---------
        protein: An instance of class proteins.

        Returns
        -------
        The list of transcription units activated or inhibited by the
        supplied protein or modified protein frame.

        LISP equivalent: transcription-units-of-protein (protein)

        """
        return self.call("transcription-units-of-protein", protein)

    def transported_chemicals(self, rxn):
        """ Return compounds that change compartments in rxn.

        Arguments
        ---------
        rxn: A reaction frame.

        Returns
        -------
        Those compounds in the transport reaction rxn that change
        compartments (meaning a compound is in one compartment as a
        reactant and is in a different compartment as a product).

        LISP equivalent: transported-chemicals (rxn)

        """
        return self.call("transported-chemicals", rxn)

    def transporter(self, protein):
        """ Check whether protein is a transporter.

        Arguments
        ---------
        protein: An instance of class proteins.

        Returns
        -------
        True if protein is a transporter.

        LISP equivalent: transporter? (protein)

        """
        return self.call("transporter?", protein)

    def unmodified_form(self, protein):
        """ Get the unmodified form of the given protein.

        The unmodified form may be the same as the given protein.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.

        Returns
        -------
        An instance of the class |Proteins|.

        LISP equivalent: unmodified-form

        """
        return self.call("unmodified-form", protein)

    def unmodified_gene_product(self, gene):
        """ Get one unmodified product of gene.

        Returns the first element of the list returned by the function
        unmodified_gene_products. This is useful if you are sure that
        there are no alternative splice forms of your gene.

        Arguments
        ---------
        gene: An instance of class |Genes|.

        Returns
        -------
        An instance of either class |Polypeptides| or 'RNA.

        LISP equivalent: unmodified-gene-product

        """
        return self.call("unmodified-gene-product", gene)

    def unmodified_gene_products(self, gene):
        """ List the unmodified gene products of gene.

        Returns all of the unmodified gene products (i.e. alternative
        splice forms) of the given gene.

        Arguments
        ---------
        gene: An instance of class |Genes|.

        Returns
        -------
        A list of instances of either class |Polypeptides| or 'RNA.

        LISP equivalent: unmodified-gene-products

        """
        return self.call("unmodified-gene-products", gene)

    def unmodified_or_unbound_form(self, protein):
        """ Get the unmodified/unbound-to-small-molecule form of protein.

        Returns the unmodified form or unbound (to a small molecule) form
        of the given protein, which might be the same as the given
        protein.

        Arguments
        ---------
        protein: An instance of the class |Proteins|.

        Returns
        -------
        An instance of the class |Proteins|.

        LISP equivalent: unmodified-or-unbound-form

        """
        return self.call("unmodified-or-unbound-form", protein)

    def variants_of_pathway(self, pwy):
        """ Returns all variants of a pathway.

        Arguments
        ---------
        pwy: An instance of the class |Pathways|.

        Returns
        -------
        A list of instance of the class |Pathways|.

        LISP equivalent: variants-of-pathway

        """
        return self.call("variants-of-pathway", pwy)

    ############
    # DATABASE MODIFICATION

    def save_kb(self):
        """ Save this database. 

        LISP equivalent: save-kb

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node31.html#3916
        return self.call("save-kb")

    def revert_kb(self):
        """ Revert this database to previously saved state.
        
        LISP equivalent: revert-kb

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node31.html#3920
        return self.call("revert-kb")
    
    def put_slot_value(self, frame, slot, value):
        """ Set the value of slot of frame to a single element.

        Arguments
        ---------
        frame: a frame in this database
        slot: a string specifying one slot of the frame

        value: the new value (a string, number, True, False, None, or
            list; strings will be interpreted as LISP symbols unless
            they begin and end with double quotes, e.g., '"string
            literal"'.))

        Returns
        -------
        List containing the value.

        LISP equivalent: put-slot-value

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node29.html#3761
        return self.call("put-slot-value", frame, slot, value)

    def put_slot_values(self, frame, slot, values):
        """ Assign a set of values to slot of frame.

        Arguments
        ---------
        frame: a frame in this database
        slot: a string specifying one slot of the frame
        values: the new values (a list of strings, numbers, True,
            False, None, or lists; strings will be interpreted 
            as LISP symbols unless they begin and end with double quotes,
            e.g., '"string literal"'.)

        Returns
        -------
        The new list of slot values (may not equal the 'values' argument,
        if 

        LISP equivalent: put-slot-values

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node29.html#3757
        return self.call("put-slot-values", frame, slot, values)

    def remove_slot_value(self, frame, slot, value):
        """ Remove value from the set of values of slot of frame.

        All occurences, if any, of the value in the slot will be 
        removed. (The 'index' keyword argument is not implemented.)

        Arguments
        ---------
        frame: a frame in this database
        slot: a string specifying one slot of the frame
        value: the value to remove (string, number, True, False, None,
            or list; strings will be interpreted as LISP symbols
            unless they begin and end with double quotes, e.g.,
            '"string literal"'.))

        Returns
        -------
        List containing the value.

        LISP equivalent: remove-slot-value

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node29.html#3773
        return self.call("remove-slot-value", frame, slot, value)

    def replace_slot_value(self, frame, slot, old_value, new_value):
        """ Replace one value of slot of frame with another.

        All occurences, if any, of the value in the set of values of
        the slot of the frame will be replaced. (The 'index' keyword
        argument is not implemented.)

        Arguments
        ---------
        frame: a frame in this database
        slot: a string specifying one slot of the frame
        old_value: the value to remove (string, number, True, False,
            None, or list; strings will be interpreted as LISP symbols
            unless they begin and end with double quotes, e.g.,
            '"string literal"'.))
        new_value: the value to replace it with (same format)

        Returns
        -------
        List containing the value.

        LISP equivalent: replace-slot-value

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node29.html#3769
        return self.call("replace-slot-value", frame, slot, old_value, 
                         new_value)

    def add_slot_value(self, frame, slot, value):
        """ Add a value to the set of values of a slot of a frame.

        Depending on the slot, the value may be added only if it is not
        already present. 

        Arguments
        ---------
        frame: a frame in this database
        slot: a string specifying one slot of the frame
        value: the new value (string, number, True, False, None, or
            list; strings will be interpreted as LISP symbols unless
            they begin and end with double quotes, e.g., '"string
            literal"'.))

        Returns
        -------
        List containing the value.

        LISP equivalent: add-slot-value

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node29.html#3765
        return self.call("add-slot-value", frame, slot, value)

    def create_frame(self, name, types=[], supers=[]):
        """ Create a new frame.

        Arguments
        ---------
        name: the new frame's name (string containing a valid LISP symbol)
        types: a list of class frames of which the new frame should be a direct
            instance
        supers: a list of class frames of which the new frame should
            be a direct subclass (if any, the new frame is a class frame.)

        Exactly one of 'types' and 'supers' should be specified.

        Returns
        -------
        The new frame's ID.

        LISP equivalent: create-frame

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node24.html#3463
        # Note that using supers=[':class'], leading to the LISP 
        # kw arg :direct-supers '(:CLASS), leads to an error in my tests,
        # not the creation of a root class frame.

        kwargs = {}
        if types:
            kwargs['direct-types'] = types
        if supers:
            kwargs['direct-supers'] = supers
        return self.call('create-frame', name, **kwargs)

    def create_class(self, name, supers):
        """ Create a new class.
        
        Arguments
        ---------
        name: the new frame's name (string containing a valid LISP symbol)
        supers: a list of class frames of which the new frame should
            be a direct subclass

        Returns
        -------
        The created class.

        LISP equivalent: create-class

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node26.html#3479
        return self.call('create-class', name, supers)

    def create_instance(self, name, types):
        """ Create a new instance frame.
        
        Arguments
        ---------
        name: the new frame's name (string containing a valid LISP symbol)
        types: one or more class frames of which the new frame should
            be a direct instance

        Returns
        -------
        The created instance.

        LISP equivalent: create-instance

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node26.html#3666
        return self.call('create-instance', name, types)

    def put_instance_types(self, instance, new_types):
        """ Change instance's types.

        Changes instance to be an instance of all the classes in
        new_types.

        Arguments
        ---------
        instance: an instance frame in this database
        new_types: a list of class frames in this database
        
        Returns
        -------
        The new list of types.

        LISP equivalent: put-instance-types

        """
        # http://www.ai.sri.com/~gfp/spec/paper/node28.html#3712
        return self.call('put-instance-types', instance, new_types)

    def delete_frame(self, frame):
        """ Delete a frame from database.

        Arguments
        ---------
        frame: a frame to be deleted.

        Returns
        -------
        None.
        
        LISP equivalent: delete-frame
        
        """
        # http://www.ai.sri.com/~gfp/spec/paper/node24.html#3471
        return self.call('delete-frame', frame)

    #######################
    # SLOT VALUE ANNOTATION

    def get_all_annots(self, frame, slot, value):
        """ Get all annotation labels in use for a value of frame of slot. 

        Arguments
        ---------
        frame: a frame in this database
        slot: string specifying a slot 
        value: a value of the slot

        Returns
        -------
        A list of slot value annotation labels, or None.

        LISP equivalent: get-all-annots

        """
        return self.call('get-all-annots', frame, slot, value)

    def get_value_annots(self, frame, slot, value, label):
        """ List labeled annotations of a value of slot of frame. 

        Arguments
        ---------
        frame: a frame in this database
        slot: string specifying a slot 
        value: a value of the slot
        label: string specifying a label for annotations of the value 

        Returns
        -------
        A list of annotations, or None.

        LISP equivalent: get-value-annots
       
        """
        return self.call('get-value-annots', frame, slot, value, label)


    def get_value_annot(self, frame, slot, value, label):
        """ Get one labeled annotation of a value of slot of frame. 

        Arguments
        ---------
        frame: a frame in this database
        slot: string specifying a slot 
        value: a value of the slot
        label: string specifying a label for annotations of the value 

        Returns
        -------
        The annotation value, or None.

        LISP equivalent: get-value-annot

        """
        return self.call('get-value-annot', frame, slot, value, label)

