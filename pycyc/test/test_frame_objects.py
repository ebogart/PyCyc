"""
Test the frame object functionality.

These tests require ecocyc and metacyc. As with the core Pathway Tools
API tests, they will fail if the specific objects used for testing
purposes change in certain ways due to database updates.

"""

import re
import pycyc
from nose.tools import raises, assert_is_instance
from frame import Frame

# This differs slightly from the attribute pattern defined in frame.py
# as here we are validating python attribute names, not finding LISP slot
# names which may be made pythonic
attribute_pattern = re.compile(r'[a-z_][a-z0-9_]*')

ec = pycyc.open('ecoli')
meta = pycyc.open('meta')

# Test frame equality
trp_ecoli_1 = Frame('trp',ec)
trp = trp_ecoli_1 
trp_ecoli_2 = Frame('trp',ec)
trp_meta = Frame('trp',meta)
trp_single_value_slot = 'COMMON-NAME'
trp_many_value_slot = 'SYNONYMS'
trp_untranslatable_slot = 'SCHEMA?'
trp_empty_slot = 'SUPERATOMS'
class_ = Frame('|Genes|',ec)
instance = trp
class_member = 'EG11219' 
class_nonmember = 'not-a-gene'
class_subclass = '|Unclassified-Genes|'
class_nonsubclass = '|Reactions|'

def test_frame_equality():
    assert trp_ecoli_1 == trp_ecoli_2
    assert not trp_ecoli_1 != trp_ecoli_2
    assert trp_ecoli_1 != trp_meta
    assert not trp_ecoli_1 == trp_meta

def test_frame_dir():
    attributes = dir(trp_meta)
    assert all([attribute_pattern.match(s) for s in attributes])
    assert 'synonyms' in attributes
    assert '_kb' in attributes

def test_frame_str_repr():
    assert str(trp_meta) == 'trp'
    assert ('%s' % trp_meta) == 'trp'
    assert repr(trp_meta) == 'trp'

def test_frame_getattr():
    synonyms = trp.synonyms
    assert_is_instance(synonyms,list)
    assert 'tryptophan' in synonyms

@raises(AttributeError)
def test_frame_get_bad_attr():
    trp.utterly_fake_attribute_name

def test_frame_get_element():
    trp[trp_untranslatable_slot] 
    trp[trp_empty_slot]

# To do: automatically check that names of non-slots, non-frames, etc,
# are really not slots/frames/etc

@raises(KeyError)
def test_frame_get_bad_element():
    trp['utterly-fake-slot-name']

@raises(ValueError)
def test_getattr_bad_frame():
    nonframe = Frame('not-a-valid-frame',ec)
    nonframe.synonyms

def test_test_frame():
    frame = Frame('trp',ec)
    nonframe = Frame('not-a-valid-frame',ec)
    assert frame.test_frame()
    assert not nonframe.test_frame()

def test_slot_values():
    assert_is_instance(trp.slot_values(trp_single_value_slot), list)

def test_lists_of_values():
    assert not isinstance(trp[trp_single_value_slot], list)
    assert isinstance(trp[trp_many_value_slot], list)

def test_frame_keys():
    keys = trp.keys()
    for a in [trp_single_value_slot, trp_empty_slot, trp_untranslatable_slot]:
        assert a in keys

@raises(TypeError)
def test_in_instance():
    class_member in instance

def test_in_class():
    assert class_member in class_
    assert class_nonmember not in class_
    assert Frame(class_member, ec) in class_
    assert Frame(class_nonmember, ec) not in class_
    assert class_subclass in class_
    assert class_nonsubclass not in class_
    assert Frame(class_subclass, ec) in class_
    assert Frame(class_nonsubclass, ec) not in class_
