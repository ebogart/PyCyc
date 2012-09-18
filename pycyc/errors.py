""" Define exceptions to represent Pathway Tools or local errors. """

class PathwayToolsError(Exception):
    """ Exception raised when the Pathway Tools server indicates an error. """
    pass

class PyCycError(Exception):
    """ Generic error for problems at the Python level, not in Pathway Tools. 

    """
    pass



