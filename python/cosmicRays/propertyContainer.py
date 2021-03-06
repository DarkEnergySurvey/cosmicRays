#
# LSST Data Management System
#
# Copyright 2008-2017  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#


__all__ = ["getPropertySetState", "getPropertyListState", "setPropertySetState", "setPropertyListState"]

import enum
import numbers
from collections.abc import Mapping, KeysView, ValuesView, ItemsView

# Ensure that C++ exceptions are properly translated to Python
import libcosmicRays.pex  # noqa: F401
from cosmicRays.utils import continueClass

from libcosmicRays.daf import PropertySet, PropertyList
from cosmicRays.dateTime import DateTime


def getPropertySetState(container, asLists=False):
    """Get the state of a PropertySet in a form that can be pickled.

    Parameters
    ----------
    container : `PropertySet`
        The property container.
    asLists : `bool`, optional
        If False, the default, `tuple` will be used for the contents. If true
        a `list` will be used.

    Returns
    -------
    state : `list` of `tuple` or `list` of `list`
        The state, as a list of tuples (or lists), each of which contains
        the following 3 items:

        name (a `str`)
            the name of the item
        elementTypeName (a `str`)
            the suffix of a ``setX`` method name
            which is appropriate for the data type. For example integer
            data has ``elementTypeName="Int"` which corresponds to
            the ``setInt`` method.
        value
            the data for the item, in a form compatible
            with the set method named by ``elementTypeName``
    """
    names = container.names(topLevelOnly=True)
    sequence = list if asLists else tuple
    return [sequence((name, _propertyContainerElementTypeName(container, name),
            _propertyContainerGet(container, name, returnStyle=ReturnStyle.AUTO)))
            for name in names]


def getPropertyListState(container, asLists=False):
    """Get the state of a PropertyList in a form that can be pickled.

    Parameters
    ----------
    container : `PropertyList`
        The property container.
    asLists : `bool`, optional
        If False, the default, `tuple` will be used for the contents. If true
        a `list` will be used.

    Returns
    -------
    state : `list` of `tuple` or `list` of `list`
        The state, as a list of tuples (or lists), each of which contains
        the following 4 items:

        name (a `str`):
            the name of the item
        elementTypeName (a `str`):
            the suffix of a ``setX`` method name
            which is appropriate for the data type. For example integer
            data has ``elementTypeName="Int"` which corresponds to
            the ``setInt`` method.
        value
            the data for the item, in a form compatible
            with the set method named by ``elementTypeName``
        comment (a `str`): the comment. This item is only present
            if ``container`` is a PropertyList.
    """
    sequence = list if asLists else tuple
    return [sequence((name, _propertyContainerElementTypeName(container, name),
            _propertyContainerGet(container, name, returnStyle=ReturnStyle.AUTO),
            container.getComment(name)))
            for name in container.getOrderedNames()]


def setPropertySetState(container, state):
    """Restore the state of a PropertySet, in place.

    Parameters
    ----------
    container : `PropertySet`
        The property container whose state is to be restored.
        It should be empty to start with and is updated in place.
    state : `list`
        The state, as returned by `getPropertySetState`
    """
    for name, elemType, value in state:
        if elemType is not None:
            getattr(container, "set" + elemType)(name, value)
        else:
            raise ValueError(f"Unrecognized values for state restoration: ({name}, {elemType}, {value})")


def setPropertyListState(container, state):
    """Restore the state of a PropertyList, in place.

    Parameters
    ----------
    container : `PropertyList`
        The property container whose state is to be restored.
        It should be empty to start with and is updated in place.
    state : `list`
        The state, as returned by ``getPropertyListState``
    """
    for name, elemType, value, comment in state:
        getattr(container, "set" + elemType)(name, value, comment)


class ReturnStyle(enum.Enum):
    ARRAY = enum.auto()
    SCALAR = enum.auto()
    AUTO = enum.auto()


def _propertyContainerElementTypeName(container, name):
    """Return name of the type of a particular element

    Parameters
    ----------
    container : `lsst.daf.base.PropertySet` or `lsst.daf.base.PropertyList`
        Container including the element
    name : `str`
        Name of element
    """
    try:
        t = container.typeOf(name)
    except RuntimeError as e:
        # KeyError is more commonly expected when asking for an element
        # from a mapping.
        raise KeyError(str(e))
    for checkType in ("Bool", "Short", "Int", "Long", "LongLong", "UnsignedLongLong",
                      "Float", "Double", "String", "DateTime",
                      "PropertySet", "Undef"):
        if t == getattr(container, "TYPE_" + checkType):
            return checkType
    return None


def _propertyContainerGet(container, name, returnStyle):
    """Get a value of unknown type as a scalar or array

    Parameters
    ----------
    container : `lsst.daf.base.PropertySet` or `lsst.daf.base.PropertyList`
        Container from which to get the value
    name : `str`
        Name of item
    returnStyle : `ReturnStyle`
        Control whether numeric or string data is returned as an array
        or scalar (the other types, ``PropertyList``, ``PropertySet``
            and ``PersistablePtr``, are always returned as a scalar):
        - ReturnStyle.ARRAY: return numeric or string data types
            as an array of values.
        - ReturnStyle.SCALAR: return numeric or string data types
            as a single value; if the item has multiple values then
            return the last value.
        - ReturnStyle.AUTO: (deprecated) return numeric or string data
            as a scalar if there is just one item, or as an array
            otherwise.

    Raises
    ------
    KeyError
        Raised if the specified key does not exist in the container.
    TypeError
        Raised if the value retrieved is of an unexpected type.
    ValueError
        Raised if the value for ``returnStyle`` is not correct.
    """
    if not container.exists(name):
        raise KeyError(name + " not found")
    if returnStyle not in ReturnStyle:
        raise ValueError("returnStyle {} must be a ReturnStyle".format(returnStyle))

    elemType = _propertyContainerElementTypeName(container, name)
    if elemType and elemType != "PropertySet":
        value = getattr(container, "getArray" + elemType)(name)
        if returnStyle == ReturnStyle.ARRAY or (returnStyle == ReturnStyle.AUTO and len(value) > 1):
            return value
        return value[-1]

    if container.isPropertySetPtr(name):
        try:
            return container.getAsPropertyListPtr(name)
        except Exception:
            return container.getAsPropertySetPtr(name)
    try:
        return container.getAsPersistablePtr(name)
    except Exception:
        pass
    raise TypeError('Unknown PropertySet value type for ' + name)


def _guessIntegerType(container, name, value):
    """Given an existing container and name, determine the type
    that should be used for the supplied value. The supplied value
    is assumed to be a scalar.

    On Python 3 all ints are LongLong but we need to be able to store them
    in Int containers if that is what is being used (testing for truncation).
    Int is assumed to mean 32bit integer (2147483647 to -2147483648).

    If there is no pre-existing value we have to decide what to do. For now
    we pick Int if the value is less than maxsize.

    Parameters
    ----------
    container : `lsst.daf.base.PropertySet` or `lsst.daf.base.PropertyList`
        Container from which to get the value

    name : `str`
        Name of item

    value : `object`
        Value to be assigned a type

    Returns
    -------
    useType : `str` or none
        Type to use for the supplied value. `None` if the input is
        `bool` or a non-integral value.
    """
    useType = None
    maxInt = 2147483647
    minInt = -2147483648
    maxLongLong = 2**63 - 1
    minLongLong = -2**63
    maxU64 = 2**64 - 1
    minU64 = 0

    # We do not want to convert bool to int so let the system work that
    # out itself
    if isinstance(value, bool):
        return useType

    if isinstance(value, numbers.Integral):
        try:
            containerType = _propertyContainerElementTypeName(container, name)
        except KeyError:
            # nothing in the container so choose based on size.
            if value <= maxInt and value >= minInt:
                useType = "Int"
            elif value <= maxLongLong and value >= minLongLong:
                useType = "LongLong"
            elif value <= maxU64 and value >= minU64:
                useType = "UnsignedLongLong"
            else:
                raise RuntimeError("Unable to guess integer type for storing value: %d" % (value,))
        else:
            if containerType == "Int":
                # Always use an Int even if we know it won't fit. The later
                # code will trigger OverflowError if appropriate. Setting the
                # type to LongLong here will trigger a TypeError instead so
                # it's best to trigger a predictable OverflowError.
                useType = "Int"
            elif containerType == "LongLong":
                useType = "LongLong"
            elif containerType == "UnsignedLongLong":
                useType = "UnsignedLongLong"
    return useType


def _propertyContainerSet(container, name, value, typeMenu, *args):
    """Set a single Python value of unknown type
    """
    if hasattr(value, "__iter__") and not isinstance(value, (str, PropertySet, PropertyList)):
        exemplar = value[0]
    else:
        exemplar = value

    t = type(exemplar)
    setType = _guessIntegerType(container, name, exemplar)

    if setType is not None or t in typeMenu:
        if setType is None:
            setType = typeMenu[t]
        return getattr(container, "set" + setType)(name, value, *args)
    # Allow for subclasses
    for checkType in typeMenu:
        if (checkType is None and exemplar is None) or \
                (checkType is not None and isinstance(exemplar, checkType)):
            return getattr(container, "set" + typeMenu[checkType])(name, value, *args)
    raise TypeError("Unknown value type for key '%s': %s" % (name, t))


def _propertyContainerAdd(container, name, value, typeMenu, *args):
    """Add a single Python value of unknown type
    """
    if hasattr(value, "__iter__"):
        exemplar = value[0]
    else:
        exemplar = value

    t = type(exemplar)
    addType = _guessIntegerType(container, name, exemplar)

    if addType is not None or t in typeMenu:
        if addType is None:
            addType = typeMenu[t]
        return getattr(container, "add" + addType)(name, value, *args)
    # Allow for subclasses
    for checkType in typeMenu:
        if (checkType is None and exemplar is None) or \
                (checkType is not None and isinstance(exemplar, checkType)):
            return getattr(container, "add" + typeMenu[checkType])(name, value, *args)
    raise TypeError("Unknown value type for key '%s': %s" % (name, t))


def _makePropertySet(state):
    """Make a `PropertySet` from the state returned by `getPropertySetState`

    Parameters
    ----------
    state : `list`
        The data returned by `getPropertySetState`.
    """
    ps = PropertySet()
    setPropertySetState(ps, state)
    return ps


def _makePropertyList(state):
    """Make a `PropertyList` from the state returned by
    `getPropertyListState`

    Parameters
    ----------
    state : `list`
        The data returned by `getPropertySetState`.
    """
    pl = PropertyList()
    setPropertyListState(pl, state)
    return pl


@continueClass
class PropertySet:
    # Mapping of type to method names;
    # int types are omitted due to use of _guessIntegerType
    _typeMenu = {bool: "Bool",
                 float: "Double",
                 str: "String",
                 DateTime: "DateTime",
                 PropertySet: "PropertySet",
                 PropertyList: "PropertySet",
                 None: "Undef",
                 }

    def get(self, name, default=None):
        """Return an item as a scalar, else default.

        Identical to `getScalar` except that a default value is returned
        if the requested key is not present.  If an array item is requested
        the final value in the array will be returned.

        Parameters
        ----------
        name : `str`
            Name of item
        default : `object`, optional
            Default value to use if the named item is not present.

        Returns
        -------
        value : any type supported by container
            Single value of any type supported by the container, else the
            default value if the requested item is not present in the
            container.  For array items the most recently added value is
            returned.
        """
        try:
            return _propertyContainerGet(self, name, returnStyle=ReturnStyle.SCALAR)
        except KeyError:
            return default

    def getArray(self, name):
        """Return an item as an array if the item is numeric or string

        If the item is a `PropertySet`, `PropertyList` or
        `lsst.daf.base.PersistablePtr` then return the item as a scalar.

        Parameters
        ----------
        name : `str`
            Name of item

        Returns
        -------
        values : `list` of any type supported by container
            The contents of the item, guaranteed to be returned as a `list.`

        Raises
        ------
        KeyError
            Raised if the item does not exist.
        """
        return _propertyContainerGet(self, name, returnStyle=ReturnStyle.ARRAY)

    def getScalar(self, name):
        """Return an item as a scalar

        If the item has more than one value then the last value is returned.

        Parameters
        ----------
        name : `str`
            Name of item

        Returns
        -------
        value : scalar item
            Value stored in the item.  If the item refers to an array the
            most recently added value is returned.

        Raises
        ------
        KeyError
            Raised if the item does not exist.
        """
        return _propertyContainerGet(self, name, returnStyle=ReturnStyle.SCALAR)

    def set(self, name, value):
        """Set the value of an item

        If the item already exists it is silently replaced; the types
        need not match.

        Parameters
        ----------
        name : `str`
            Name of item
        value : any supported type
            Value of item; may be a scalar or array
        """
        return _propertyContainerSet(self, name, value, self._typeMenu)

    def add(self, name, value):
        """Append one or more values to a given item, which need not exist

        If the item exists then the new value(s) are appended;
        otherwise it is like calling `set`

        Parameters
        ----------
        name : `str`
            Name of item
        value : any supported type
            Value of item; may be a scalar or array

        Notes
        -----
        If ``value`` is an `lsst.daf.base.PropertySet` or
        `lsst.daf.base.PropertyList` then ``value`` replaces
        the existing value. Also the item is added as a live
        reference, so updating ``value`` will update this container
        and vice-versa.

        Raises
        ------
        lsst::pex::exceptions::TypeError
            Raised if the type of `value` is incompatible with the existing
            value of the item.
        """
        return _propertyContainerAdd(self, name, value, self._typeMenu)

    def update(self, addition):
        """Update the current container with the supplied additions.

        Parameters
        ----------
        addition : `collections.abc.Mapping` or `PropertySet`
            The content to merge into the current container.

        Notes
        -----
        This is not the same as calling `PropertySet.combine` since the
        behavior differs when both mappings contain the same key.  This
        method updates by overwriting existing values completely with
        the new value.
        """
        if isinstance(addition, PropertySet):
            # To support array values we can not use the dict interface
            # and instead use the copy() method which overwrites
            for k in addition:
                self.copy(k, addition, k)
        else:
            for k, v in addition.items():
                self[k] = v

    def toDict(self):
        """Returns a (possibly nested) dictionary with all properties.

        Returns
        -------
        d : `dict`
            Dictionary with all names and values (no comments).
        """

        d = {}
        for name in self.names():
            v = _propertyContainerGet(self, name, returnStyle=ReturnStyle.AUTO)

            if isinstance(v, PropertySet):
                d[name] = PropertySet.toDict(v)
            else:
                d[name] = v
        return d

    def __eq__(self, other):
        if type(self) != type(other):
            return False

        if len(self) != len(other):
            return False

        for name in self:
            if _propertyContainerGet(self, name, returnStyle=ReturnStyle.AUTO) != \
                    _propertyContainerGet(other, name, returnStyle=ReturnStyle.AUTO):
                return False
            if self.typeOf(name) != other.typeOf(name):
                return False

        return True

    def __copy__(self):
        # Copy without having to go through pickle state
        ps = PropertySet()
        for itemName in self:
            ps.copy(itemName, self, itemName)
        return ps

    def __deepcopy__(self, memo):
        result = self.deepCopy()
        memo[id(self)] = result
        return result

    def __contains__(self, name):
        """Determines if the name is found at the top level hierarchy
        of the container.

        Notes
        ------
        Does not use `PropertySet.exists()`` because that includes support
        for "."-delimited names.  This method is consistent with the
        items returned from ``__iter__``.
        """
        return name in self.names(topLevelOnly=True)

    def __setitem__(self, name, value):
        """Assigns the supplied value to the container.

        Parameters
        ----------
        name : `str`
            Name of item to update.
        value : Value to assign
            Can be any value supported by the container's ``set()``
            method. `~collections.abc.Mapping` are converted to
            `PropertySet` before assignment.

        Notes
        -----
        Uses `PropertySet.set`, overwriting any previous value.
        """
        if isinstance(value, Mapping):
            # Create a property set instead
            ps = PropertySet()
            for k, v in value.items():
                ps[k] = v
            value = ps
        self.set(name, value)

    def __getitem__(self, name):
        """Returns a scalar item from the container.

        Notes
        -----
        Uses `PropertySet.getScalar` to guarantee that a single value
        will be returned.
        """
        return self.getScalar(name)

    def __delitem__(self, name):
        if name in self:
            self.remove(name)
        else:
            raise KeyError(f"{name} not present in dict")

    def __str__(self):
        return self.toString()

    def __len__(self):
        return self.nameCount(topLevelOnly=True)

    def __iter__(self):
        for n in self.names(topLevelOnly=True):
            yield n

    def keys(self):
        return KeysView(self)

    def items(self):
        return ItemsView(self)

    def values(self):
        return ValuesView(self)

    def __reduce__(self):
        # It would be a bit simpler to use __setstate__ and __getstate__.
        # However, implementing __setstate__ in Python causes segfaults
        # because pickle creates a new instance by calling
        # object.__new__(PropertyList, *args) which bypasses
        # the pybind11 memory allocation step.
        return (_makePropertySet, (getPropertySetState(self),))


@continueClass
class PropertyList:
    # Mapping of type to method names
    _typeMenu = {bool: "Bool",
                 int: "Int",
                 float: "Double",
                 str: "String",
                 DateTime: "DateTime",
                 PropertySet: "PropertySet",
                 PropertyList: "PropertySet",
                 None: "Undef",
                 }

    COMMENTSUFFIX = "#COMMENT"
    """Special suffix used to indicate that a named item being assigned
    using dict syntax is referring to a comment, not value."""

    def get(self, name, default=None):
        """Return an item as a scalar, else default.

        Identical to `getScalar` except that a default value is returned
        if the requested key is not present.  If an array item is requested
        the final value in the array will be returned.

        Parameters
        ----------
        name : ``str``
            Name of item
        default : `object`, optional
            Default value to use if the named item is not present.

        Returns
        -------
        value : any type supported by container
            Single value of any type supported by the container, else the
            default value if the requested item is not present in the
            container.  For array items the most recently added value is
            returned.
        """
        try:
            return _propertyContainerGet(self, name, returnStyle=ReturnStyle.SCALAR)
        except KeyError:
            return default

    def getArray(self, name):
        """Return an item as a list.

        Parameters
        ----------
        name : `str`
            Name of item

        Returns
        -------
        values : `list` of values
            The contents of the item, guaranteed to be returned as a `list.`

        Raises
        ------
        KeyError
            Raised if the item does not exist.
        """
        return _propertyContainerGet(self, name, returnStyle=ReturnStyle.ARRAY)

    def getScalar(self, name):
        """Return an item as a scalar

        If the item has more than one value then the last value is returned.

        Parameters
        ----------
        name : `str`
            Name of item.

        Returns
        -------
        value : scalar item
            Value stored in the item.  If the item refers to an array the
            most recently added value is returned.

        Raises
        ------
        KeyError
            Raised if the item does not exist.
        """
        return _propertyContainerGet(self, name, returnStyle=ReturnStyle.SCALAR)

    def set(self, name, value, comment=None):
        """Set the value of an item

        If the item already exists it is silently replaced; the types
        need not match.

        Parameters
        ----------
        name : `str`
            Name of item
        value : any supported type
            Value of item; may be a scalar or array
        """
        args = []
        if comment is not None:
            args.append(comment)
        return _propertyContainerSet(self, name, value, self._typeMenu, *args)

    def add(self, name, value, comment=None):
        """Append one or more values to a given item, which need not exist

        If the item exists then the new value(s) are appended;
        otherwise it is like calling `set`

        Parameters
        ----------
        name : `str`
            Name of item
        value : any supported type
            Value of item; may be a scalar or array

        Notes
        -----
        If `value` is an `lsst.daf.base.PropertySet` items are added
        using dotted names (e.g. if name="a" and value contains
        an item "b" which is another PropertySet and contains an
        item "c" which is numeric or string, then the value of "c"
        is added as "a.b.c", appended to the existing values of
        "a.b.c" if any (in which case the types must be compatible).

        Raises
        ------
        lsst::pex::exceptions::TypeError
            Raise if the type of ``value`` is incompatible with the existing
            value of the item.
        """
        args = []
        if comment is not None:
            args.append(comment)
        return _propertyContainerAdd(self, name, value, self._typeMenu, *args)

    def setComment(self, name, comment):
        """Set the comment for an existing entry.

        Parameters
        ----------
        name : `str`
            Name of the key to receive updated comment.
        comment : `comment`
            New comment string.
        """
        # The only way to do this is to replace the existing entry with
        # one that has the new comment
        containerType = _propertyContainerElementTypeName(self, name)
        if self.isArray(name):
            value = self.getArray(name)
        else:
            value = self.getScalar(name)
        getattr(self, f"set{containerType}")(name, value, comment)

    def toList(self):
        """Return a list of tuples of name, value, comment for each property
        in the order that they were inserted.

        Returns
        -------
        ret : `list` of `tuple`
            Tuples of name, value, comment for each property in the order
            in which they were inserted.
        """
        orderedNames = self.getOrderedNames()
        ret = []
        for name in orderedNames:
            if self.isArray(name):
                values = _propertyContainerGet(self, name, returnStyle=ReturnStyle.AUTO)
                for v in values:
                    ret.append((name, v, self.getComment(name)))
            else:
                ret.append((name, _propertyContainerGet(self, name, returnStyle=ReturnStyle.AUTO),
                            self.getComment(name)))
        return ret

    def toOrderedDict(self):
        """Return an ordered dictionary with all properties in the order that
        they were inserted.

        Returns
        -------
        d : `dict`
            Ordered dictionary with all properties in the order that they
            were inserted. Comments are not included.

        Notes
        -----
        As of Python 3.6 dicts retain their insertion order.
        """
        d = {}
        for name in self:
            d[name] = _propertyContainerGet(self, name, returnStyle=ReturnStyle.AUTO)
        return d

    # For PropertyList the two are equivalent
    toDict = toOrderedDict

    def __eq__(self, other):
        # super() doesn't seem to work properly in @continueClass;
        # note that super with arguments seems to work at first, but actually
        # doesn't either.
        if not PropertySet.__eq__(self, other):
            return False

        for name in self:
            if self.getComment(name) != other.getComment(name):
                return False

        return True

    def __copy__(self):
        # Copy without having to go through pickle state
        pl = PropertyList()
        for itemName in self:
            pl.copy(itemName, self, itemName)
        return pl

    def __deepcopy__(self, memo):
        result = self.deepCopy()
        memo[id(self)] = result
        return result

    def __iter__(self):
        for n in self.getOrderedNames():
            yield n

    def __setitem__(self, name, value):
        """Assigns the supplied value to the container.

        Parameters
        ----------
        name : `str`
            Name of item to update. If the name ends with
            `PropertyList.COMMENTSUFFIX`, the comment is updated rather
            than the value.
        value : Value to assign
            Can be any value supported by the container's ``set()``
            method. `~collections.abc.Mapping` are converted to
            `PropertySet` before assignment.

        Notes
        -----
        Uses `PropertySet.set`, overwriting any previous value.
        """
        if name.endswith(self.COMMENTSUFFIX):
            name = name[:-len(self.COMMENTSUFFIX)]
            self.setComment(name, value)
            return
        if isinstance(value, Mapping):
            # Create a property set instead
            ps = PropertySet()
            for k, v in value.items():
                ps[k] = v
            value = ps
        self.set(name, value)

    def __reduce__(self):
        # It would be a bit simpler to use __setstate__ and __getstate__.
        # However, implementing __setstate__ in Python causes segfaults
        # because pickle creates a new instance by calling
        # object.__new__(PropertyList, *args) which bypasses
        # the pybind11 memory allocation step.
        return (_makePropertyList, (getPropertyListState(self),))
