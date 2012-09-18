"""
Allow interaction with frames as objects whose attributes are slot values.

""" 
import re
from errors import PathwayToolsError, PyCycError

attribute_pattern = re.compile(r'^[a-zA-Z\-][a-zA-Z0-9\-]*$')

def slot_to_attribute(slot_id):
    """ Translate a LISP slot name into a Python attribute name if possible.

    If the slot name begins with a number or contains any character 
    that is not alphanumeric or '-', return None. Otherwise, 
    translate each hyphen to an underscore, lower-case the string, and 
    return. 

    """
    if attribute_pattern.match(slot_id):
        return slot_id.replace('-','_').lower()
    else:
        return None

def attribute_to_slot(attribute_id):
    """ Translate a python attribute name into a likely LISP slot name. 

    This function reverses the action of slot_to_attribute (note that 
    returning the string to upper case is purely cosmetic, as the LISP
    interpreter is insensitive to case in symbols.)

    """

    return attribute_id.replace('_','-').upper()

class Frame():
    """
    Object representing a frame to allow access to its slots.

    Slots may be accessed as object attributes, if it is possible to 
    determine their LISP symbols from valid Python ids by converting 
    underscores to dashes, e.g.

    >>> eco = pycyc.open('ecoli')
    >>> trp = eco['trp']
    >>> trp.molecular_weight
    204.288

    In some cases (where slot ids contain '?' or ':', for example) 
    slots must be accessed through indexing instead:

    >>> trp['SCHEMA?']
    True
    
    Where a slot has only one value, it is returned; where a slot has
    more than one value, a list is returned. If a list is desired
    regardless of the number of slot values, frame.slot_values(slot)
    may be called (supplying the LISP-formatted slot name):

    >>> trp.slot_values('SCHEMA?')
    ['True']

    The complete list of LISP-formatted slot names may be retrieved by
    frame.keys(). Python-formatted slot names, where available,
    are returned in dir(frame), and should apper in IPython tab completion,
    e.g. 

    To facilitate interactive use and passing frame objects as
    arguments to LISP queries, str(frame) and repr(frame) both return
    the frame's ID string. Note that this violates convention by
    providing an incomplete description of the object in its repr().
    
    Class and instance frames are both represented by python Frame
    objects. One frame may be tested for membership in another;

    >>> frame1 in frame2
    
    is True if frame1 is an instance or subclass of the class frame
    frame2.

    (Specifically, the test raises TypeError if frame2 is not a class
    frame. It returns False if the frames do not belong to the same
    kb; otherwise, it returns
    frame2._kb.instance_all_instance_of(frame1, frame2) if frame1 is
    an instance frame and frame2._kb.is_class_all_sub_of(frame1,
    frame2) if frame1 is a class frame.

    If frame1 is a string, not a Frame object, it is interpreted
    as belonging to the same kb as frame2 automatically; if there is
    no corresponding frame, the test will fail.)

    Slot values are immutable. 

    """

    def __init__(self, name, kb):
        self._name = name
        self._kb = kb
        self._frame_tested = False
        self._known_frame = False
    
    def __str__(self):
        return self._name

    def __repr__(self):
        return self._name

    def __getattr__(self, attribute):
        """ Get slot values from Pathway Tools if possible. """

        slot = attribute_to_slot(attribute)
        try:
            return self._read_slot(slot)
        except KeyError:
            raise AttributeError('Frame %s has no slot %s (translated to '
                                 'LISP format from %s). Try dir(frame) for '
                                 'a list of Python-formatted slot names or '
                                 'frame.keys() for all slot names in LISP '
                                 'format.' % (self, slot, attribute))

    def _read_slot(self, slot):
        try:
            values = self._kb.get_slot_values(self, slot)
        except PathwayToolsError:
            self.test_frame()
            if not self._known_frame:
                raise ValueError('Invalid frame %s '
                                 '(not in pathway tools database %s)' %
                                 (self, self._kb))
            else: 
                raise 

        if not values:
            if not self._kb.is_slot(self, slot):
                raise KeyError('Frame %s in database %s has no slot %s' % 
                               (self, self._kb, slot))
            else:
                return values
        elif len(values) == 1:
            return values[0]
        else:
            return values

    def __dir__(self):
        """ List object attributes and frame slots with Python ids. """

        attributes = dir(self.__class__) + self.__dict__.keys()
        for slot in self.keys():
            # Note that slot will in general be a(n invalid) Frame object
            slot_as_attribute = slot_to_attribute(str(slot))
            if slot_as_attribute:
                attributes.append(slot_as_attribute)
        return attributes

    def keys(self):
        """ List slots of this frame, in LISP format. """
        
        return self._kb.get_frame_slots(self)

    def __getitem__(self, slot):
        """ Get the value(s) of a slot of this frame (in LISP format.) """
        return self._read_slot(slot)
        
    def slot_values(self, slot):
        """ Get values of a LISP-formatted slot of this frame as a list. 

        This function, a wrapper for the LISP function
        get-slot-values, returns a list even if only the slot has only
        one value.

        """
        try:
            values = self._kb.get_slot_values(self, slot)
        except PathwayToolsError:
            self.test_frame()
            if not self._known_frame:
                raise ValueError('Invalid frame %s '
                                 '(not in pathway tools database %s)' %
                                 (self, self._kb))
            else: 
                raise 

        if not values:
            if not self._kb.is_slot(self, slot):
                raise KeyError('Frame %s in database %s has no slot %s' % 
                               (self, self._kb, slot))
        return values

    def __setattr__(self, attribute, value):
        """ Prevent most attributes from being changed. """
        if attribute in ['_name', '_kb', '_frame_tested', '_known_frame']:
            self.__dict__[attribute] = value
        else:
            raise TypeError('Frame attributes may not be changed.')

    # For convenience, wrap the slot editing functions
    def put_slot_value(self, slot, value):
        self._kb.put_slot_value(self, slot, value)

    def put_slot_values(self, slot, values):
        self._kb.put_slot_values(self, slot, values)

    def remove_slot_value(self, slot, value):
        self._kb.remove_slot_value(self, slot, value)
        
    def replace_slot_value(self, slot, old_value, new_value):
        self._kb.replace_slot_value(self, slot, old_value, new_value)

    def add_slot_value(self, slot, value):
        self._kb.add_slot_value(self, slot, value)

    def test_frame(self):
        if self._frame_tested:
            return self._known_frame
        else:
            self._known_frame = self._kb.is_coercible_to_frame(self._name)
            self._frame_tested = True
            return self._known_frame

    def __eq__(self, other):
        if isinstance(other, str):
            return self._name == other
        elif isinstance(other, Frame): 
            if self._kb != other._kb:
                return False
            else:
                return self._kb.fequal(self, other)
        else: 
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, other):
        if not self._kb.is_class(self):
            raise TypeError('Frame %s in kb %s is not a class frame.' %
                            (self, self._kb))
        if isinstance(other, Frame):
            if other._kb != self._kb:
                return False
            if not other.test_frame():
                return False
        elif not self._kb.is_coercible_to_frame(other):
            return False

        if self._kb.is_instance(other):
            return self._kb.instance_all_instance_of(other, self)
        elif self._kb.is_class(other):
            return self._kb.is_class_all_sub_of(other, self)
        else:
            return False
        

