"""
Test LISP parsing with a few examples. 

Testing the database-bound parser depends on the availability of EcoCyc.

"""

import pycyc
from lisp_parser import LispParser, DatabaseSpecificLispParser, LispParserError
from frame import Frame
from nose.tools import raises, assert_is_instance

example = '(LOREM "IPSUM (DOLOR) SIT" 9.8 (AMET () ";?\n# -" 161) PPI)'
result = ['LOREM','IPSUM (DOLOR) SIT',9.8,['AMET', [],';?\n# -',161],'PPI'] 

def test_parser():
    parser = LispParser()
    assert result == parser.parse_string(example)

@raises(LispParserError) # not StopIteration
def test_malformed_string():
    parser = LispParser()
    parser.parse_string('(PPI')

def test_special_values():
    parser = LispParser()    
    assert parser.parse_string('(T NIL)') == [True, None]

ec = pycyc.open('ecoli')

def test_frame_token():
    ec_parser = DatabaseSpecificLispParser(ec)
    trp = ec_parser.parse_single_token('TRP')
    assert_is_instance(trp, Frame)
    assert trp == Frame('TRP', ec)

def test_frame_from_string():
    ec_parser = DatabaseSpecificLispParser(ec)
    trp = ec_parser.parse_string('(TRP)')[0]
    assert_is_instance(trp, Frame)
    assert trp == Frame('TRP', ec)

