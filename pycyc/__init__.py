from pathwaytools_interface import PathwayToolsInterface
from errors import PathwayToolsError, PyCycError
from frame import Frame
from server_interface import evaluate_string

# Provide a convenience function for interactively creating
# PathwayToolsInterface objects

def open(database):
    """ Create an interface to the specified Pathway Tools database.

    Arguments:
    database -- string, a Pathway Tools org-id.

    """
    
    return PathwayToolsInterface(database)

def all_orgs():
    """ List Pathway Tools databases available from the server. """
    # Unfortunately, calling all-orgs prints a nicely formatted list
    # to standard output, which we can't capture, and returns (at least
    # approximately) '(#<OCELOT-FILE-KB CASSAVABASE @ #x10220f95b2>
    # #<OCELOT-FILE-KB METABASE @ #x10020727e2> #<OCELOT-FILE-KB ECOBASE @
    # #x1002072802>)' which our parser doesn't handle well; it renders as
    # ['#<OCELOT-FILE-KB',
    #  'CASSAVABASE',
    #  '@',
    #  '#x10220f95b2>',
    #  '#<OCELOT-FILE-KB',
    #  'METABASE',
    #  '@',
    #  '#x10020727e2>',
    #  '#<OCELOT-FILE-KB',
    #  'ECOBASE',
    #  '@',
    #  '#x1002072802>']
    # Rather than improve the parser for this single case, we grab the 
    # second, sixth, etc., entries and remove 'BASE' from the end. So far,
    # this works.

    dbs = evaluate_string('(all-orgs)')
    dbs = [db.lower() for db in dbs[1::4]]
    dbs = [db[:-4] if db.endswith('base') else db for db in dbs]
    return dbs
