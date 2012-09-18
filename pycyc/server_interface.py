"""
Python interface to pathway-tools.

Closely based on JavaCyc. Parser and tokenizer modelled on those in Peter 
Norvig's 'How to Write a LISP Interpreter in Python'.

Eli Bogart
January-July 2012

"""

import os
import socket
from lisp_parser import LispParser, DatabaseSpecificLispParser
from errors import PyCycError, PathwayToolsError
from frame import Frame

##########################################################################
# CONSTANTS

_SOCKET_FILE = '/tmp/ptools-socket'
# List the responses from the server which should raise an exception.
# Further testing needed to confirm the newline is always present.
_ERROR_STRINGS = {':error\n'} 

#########################################################################
# DEBUGGING
#import logging
#logging.basicConfig(filename='socket_log.txt', level=logging.INFO)

##########################################################################
# COMMUNICATION WITH PATHWAY-TOOLS

def get_socket_file():
    """
    Return a file object associated with the ptools Unix socket.

    """
    s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    s.connect(_SOCKET_FILE)
    return s.makefile()

def submit_string(query):
    """
    Send a string as a request to ptools and return the response string.

    """

    f = get_socket_file()
    f.write(query)
#    logging.info('Query: %s' % query)
    f.flush()
    response = f.read()
    f.close()
#    logging.info('Response: %s' % response)
    return response

basic_parser = LispParser()

def evaluate_string(query, parser=basic_parser):
    """
    Evaluate the query in ptools, parsing and returning the response.

    """

    result = submit_string(query)

    if result in _ERROR_STRINGS:
        raise PathwayToolsError('Pathway Tools error evaluating ' + 
                                'LISP query %s.' % query)
    else:
        return parser.parse_string(result)

query = evaluate_string # alias for interactive use

def evaluate_with_organism(organism, query, parser=basic_parser):
    """
    Evaluate a query in the context of a specified organism.
    
    """
    wrapped_query = "(with-organism (:org-id '%s) %s)" % (organism, query)
    try:
        return evaluate_string(wrapped_query, parser)
    except PathwayToolsError:
        raise PathwayToolsError(('In organism %s, Pathway Tools error ' + 
                                 'evaluating query %s.') % (organism, query))



##########################################################################
# PATHWAY TOOLS DB INTERFACE CLASS

class BasicPathwayToolsDBInterface():

    """ Basic tools for querying a specific organism DB within Pathway Tools.

    Child classes implement methods exposing the rest of the Pathway Tools
    internal LISP functions. 

    """

    def __init__(self, organism):
        """ Create an interface to a specific PGDB. """
        try: 
           evaluate_with_organism(organism, '()')
        except PathwayToolsError:
            raise ValueError('Invalid organism id.')
        except socket.error as e:
            raise IOError('Error connecting to Pathway Tools socket file; '
                          'ensure Pathway Tools is running in api mode? '
                          '(socket.error: %s)' % e) 
        self._org_id = organism
        self._parser = DatabaseSpecificLispParser(self)

    def __str__(self):
        return self._org_id

    def __eq__(self, other):
        """ Two interface objects are equal if their org-id strings are equal.

        This is 
        
        """
        if isinstance(other, BasicPathwayToolsDBInterface):
            return self._org_id == other._org_id
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def evaluate(self, query):
        """ Evaluate a query in the context of this PGDB. """
        return evaluate_with_organism(self._org_id, query, self._parser)
        
    def call(self, function_name, *args, **kwargs):
        """
        Call a LISP function with arguments in the context of this PGDB.
        
        Arguments:
        function_name -- function to be called
        *args -- expressions to be used as arguments
        **kwargs -- keyword arguments to LISP functions,

        String arguments will be quoted unless they begin with ':'. Arguments
        which are True will be converted to 'T'; arguments which are None
        or False will be converted to 'NIL'. Non-string arguments, e.g.
        integers, will be converted to strings but not quoted.

        Example:
        >>> call("chromosome-of-gene","g1")
        evaluates the LISP expression "(chromosome-of-gene 'g1)" 
        in the PGDB specified by self._org_id.

        """
        
        call_list = [function_name]
        for arg in args:
            call_list.append(self.translate_argument(arg))
        for keyword, arg in kwargs.iteritems():
            call_list.append(':'  + keyword)
            call_list.append(self.translate_argument(arg))                

        return self.evaluate('(%s)' % " ".join(call_list))

    def translate_argument(self, arg, quote=True):
        """ Translate value of a function argument to appropriate LISP string.

        """
        if arg is True:
            return 'T'
        elif arg is None or arg is False:
            return 'NIL'
        elif isinstance(arg, str):
            if arg.startswith(':'):
                return arg
            else:
                return ("'" if quote else "") + arg
            # String literals in double quotes?
        elif isinstance(arg, Frame):
            return ("'" if quote else "") + str(arg)
        elif isinstance(arg, list):
            # Not entirely clear how quoting should be handled;
            # assume only the outermost list is quoted.
            return (("'" if quote else "") + 
                    "(%s)" % ' '.join(self.translate_argument(i,False)
                                      for i in arg)) 
        else:
            return str(arg)
        

    def multiple_value_list_call(self, function_name, *args, **kwargs):
        """
        Call a LISP function that returns multiple values, in this PGDB.
        
        This calls function_name with args and kwargs appropriately, 
        wrapping the call in a call to multiple-value-list to retrieve
        the function's multiple return values as a list. The somewhat
        cumbersome name of this method is intended to emphasize that it
        is a multiple-value version of the method 'call' while avoiding
        ambiguous resemblance to the LISP macro multiple-value-call, 
        which is not related.

        Arguments:
        function_name -- function to be called
        *args -- expressions to be used as arguments
        **kwargs -- keyword arguments to LISP functions,

        String arguments will be quoted unless they begin with ':'. Arguments
        which are True will be converted to 'T'; arguments which are None
        or False will be converted to 'NIL'. Non-string arguments, e.g.
        integers, will be converted to strings but not quoted.

        Example:
        >>> multiple_value_list_call('reaction-reactants-and-products', \
        'TRYPSYN-RXN') 
        evaluates the LISP expression 
        "(multiple-value-list (reaction-reactants-and-products 'TRYPSYN-RXN))" 
        in the PGDB specified by self._org_id.

        """

        call_list = [function_name]
        for arg in args:
            call_list.append(self.translate_argument(arg))
        for keyword, arg in kwargs.iteritems():
            call_list.append(':'  + keyword)
            call_list.append(self.translate_argument(arg))

        return self.evaluate('(multiple-value-list (%s))' % 
                             " ".join(call_list))

