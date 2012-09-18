"""
Test interface object behavior not directly related to ptools interaction.

These tests require ecocyc and metacyc. As with the core Pathway Tools
API tests, they will fail if the specific objects used for testing
purposes change in certain ways due to database updates.

"""

import pycyc

ec1 = pycyc.open('ecoli')
ec2 = pycyc.open('ecoli')
meta = pycyc.open('meta')

def test_interface_to_str():
    assert str(ec1) == 'ecoli'

def test_interface_eq():
    assert ec1 == ec2
    assert not ec1 == meta 

def test_interface_neq():
    assert not ec1 != ec2
    assert ec1 != meta 
