PyCyc 0.0.10
September 25, 2012

INTRODUCTION
------------
PyCyc provides a Python interface to Pathway Tools. 

A Pathway Tools database such as EcoCyc or MetaCyc is represented 
by a Python interface object; when a method of the object is called, 
PyCyc submits a query to the Pathway Tools server through the socket
interface to evaluate a corresponding LISP function, then parses the 
result and returns appropriate Python values. 

Frames in the database are represented by a python Frame class which
exposes their slot values as attributes. 

Simple interactions with the database and complex queries may both be
written with a fairly natural Pythonic syntax:

>>> [cpd.common_name for cpd in ecocyc.compounds_of_pathway('PWY0-1487')]
['ATP', 'ADP', 'CreB transcriptional regulator', 'CreC sensory
histidine kinase - phosphorylated', 'CreC', 'CreB-Phosphorylated
DNA-binding transcriptional regulator']

For the complete list of implemented LISP functions and their Python
equivalents, see docs/functions.txt. Optional and keyword arguments are
supported.

INSTALLATION
------------
Run 'python setup.py install'.

USAGE
-----
Before using, ensure pathway tools is running in API mode (e.g., by
starting it with 'pathway-tools -lisp -api').

To interact with a particular organism in Pathway Tools, create an interface
to it with pycyc.open(): 

>>> import pycyc
>>> ecocyc = pycyc.open('ecoli')
>>> metacyc = pycyc.open('meta')

(For a list of available organisms, run pycyc.all_orgs().)

Methods of the interface objects call Pathway Tools functions in the
appropriate databases and translate and return the results:

>>> ecocyc.reactions_of_compound('trp')
[TRYPTOPHAN--TRNA-LIGASE-RXN,
 TRYPTOPHAN-RXN,
 RXN0-287,
 TRANS-RXN-76,
 RXN0-2382,
 RXN0-301,
 TRYPSYN-RXN]
>>> ecocyc.transcription_factor('PHOSPHO-CREB')
True

Note that, somewhat awkwardly, Boolean functions return True or None 
(rather than False):

>>> ecocyc.transcription_factor('TRYPSYN-APROTEIN')
>>> 

For convenience, lists of all compounds, reactions, genes, and pathways
in a database are available as db.compounds, db.reactions, db.genes, 
and db.pathways.

If you know how to write your desired query in LISP, you may evaluate
it directly with, eg, ecocyc.evaluate(lisp_query_string).

For more information on individual functions, see their docstrings.

FRAMES
------
Frames (records) in the databases are represented by Frame objects, 
which expose the frame slots as attributes:

>>> trp = ecocyc['trp']
>>> trp.molecular_weight
204.228

No data is stored locally; a new request to the server is made every time 
an attribute is accessed.

Some LISP slot names may not be easily translated to valid python
attribute names; these may be accessed as frame['lisp-slot-name']. The
complete list of LISP slot names is available from
frame.keys(). Translated python attributes are listed in dir(frame),
and available for IPython tab completion. To avoid ambiguity,
attributes of the Python object may not be set directly; frame slot
values may be changed through frame.put_slot_value(),
frame.put_slot_values(), etc.

(Note that frame.slot returns a list if the slot has more than one
value, otherwise the first value; a slot with two values and a slot
with one value which is a two-element list look the same. Use
frame.slot_values() instead where this is not desirable.)

Frame objects and their label strings may be used interchangeably as
arguments to Pathway Tools functions:

>>> ecocyc.reactions_of_compound(trp)
[TRYPTOPHAN--TRNA-LIGASE-RXN, ...]

Frame objects representing classes within the database support
membership testing:

>>> amino_acids = ecocyc['|Amino-Acids|']
>>> 'trp' in amino_acids
True
>>> glucose = ecocyc['glc']
>>> glucose in amino_acids
False

Classes contain their subclasses:

>>> '|Serines|' in amino_acids
True

Wherever a LISP function returns a LISP symbol not starting with ':', a
corresponding Frame object is returned in Python. Though generally
helpful, this occasionally is an annoyance, if the symbol does not
actually represent a frame in the database, or if it is necessary
to interact with the frame name as a string. (In the latter case,
str(frame) will work.) This feature may be removed or made optional in
later versions if it proves to be more trouble than it's worth; feedback 
is welcome.

KNOWN ISSUES
------------
Numbers with 'd', e.g., 3.141592653589793d0, are parsed as strings, not 
numbers.

Some functions have been implemented based only on an inspection of
the corresponding PerlCyc and/or JavaCyc implementation and are thus
incompletely documented and imperfectly tested; use with
caution. These include:

	 pwys-of-organism-in-meta
	 enzymes-of-organism-in-meta
	 lower-taxa-or-species-p
	 find-indexed-frame
	 get-reaction-list

Interface objects are tested for equality by comparing the org-id of the 
organism they query, which may lead to false negatives ('eco' and 'ecoli'
refer to the same database, for instance, but pycyc.open('eco') and
pycyc.open('ecoli') will yield apparently different interfaces and two
copies of the same frame, one retrieved from each instance, will be unequal.)

TESTING 
------- 
A comprehensive test suite is provided in the subpackage pycyc.test.
However, the tests rely on assumptions about specific properties of
frames in EcoCyc and may fail if that database changes, even though
the PyCyc implementation of the tested functions remains correct!

(To prevent unexpected changes to local databases when running the 
tests, the test script for the functions which write to the database is
not distributed by default.)

Test scripts depend on the nose package.

SOURCES 
------- 
Most docstring text is taken with minor modifications from one of the
following sources:

http://bioinformatics.ai.sri.com/ptools/api/
http://bioinformatics.ai.sri.com/ptools/ptools-fns.html
http://bioinformatics.ai.sri.com/ptools/gfp.html
http://www.ai.sri.com/~gfp/spec/paper/paper.html
http://www.ai.sri.com/pkarp/ocelot/

LICENSE
-------
Distributed under the BSD license; see LICENSE.txt.

AUTHOR
------
Eli Bogart, elb87@cornell.edu

ACKNOWLEDGEMENTS
----------------
Thanks to Chris Myers and Lukas Mueller for helpful comments.
