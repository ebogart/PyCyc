"""
LISP parser modelled on Peter Norvig's lispy/lispy2 (norvig.com/lispy2.html).

This parser is quite simple, as we are not trying to evaluate arbitrary 
LISP expressions, but simply translate (nested) lists of symbols, numbers and
strings to (nested) lists of their Python counterparts.

"""

import re
from frame import Frame

token_pattern = re.compile(r'\s*([()]|"[^"]*"|[^"\s()]*)(.*)', re.DOTALL)

class LispParserError(Exception):
    pass

class LispParser():
    def parse_string(self,expression):
        """
        Tokenize and parse a string as a LISP expression.

        """
        try:
            return self.parse_token_generator(self.token_generator(expression))
        except StopIteration:
            # This should occur only if the expression is malformed
            raise LispParserError('Malformed or unsupported LISP expression.')

    def token_generator(self,expression):
        """ etc """

        while expression:
            token, expression = token_pattern.match(expression).groups()
            if token:
                yield token
            else:
                # Handle pathological cases like expression='"', which 
                # gives token = '' and expression = '"' again
                if expression:                
                    raise LispParserError('Failure to parse LISP string')
                break

    def parse_token_generator(self, token_generator):
        """
        Parse a generator of LISP token strings.

        """
        token = next(token_generator)

        if token == '(':
            elements = [token]
            while elements[-1] != ')':            
                elements.append(self.parse_token_generator(token_generator))
            return elements[1:-1]
        elif token == ')':
            return ')'
        else:
            return self.parse_single_token(token)

    def parse_single_token(self, token):
        """
        Parse a single LISP token string to a number, boolean, None, or string.

        """

        if token.startswith('"'):
            return token.strip('"')

        if token == 'T':
            return True
        elif token == 'NIL':
            return None

        try: 
            return int(token)
        except ValueError:
            pass
        try:
            return float(token)
        except ValueError:
            pass

        return token

class DatabaseSpecificLispParser(LispParser):
    def __init__(self, kb):
        self.kb = kb

    def parse_single_token(self, token):
        """
        Parse a single LISP token as a number, Boolean, string, or Frame.

        """

        if token.startswith('"'):
            return token.strip('"')

        if token.startswith(':'):
            return token

        if token == 'T':
            return True
        elif token == 'NIL':
            return None

        try: 
            return int(token)
        except ValueError:
            pass
        try:
            return float(token)
        except ValueError:
            pass

        return Frame(token, self.kb)
