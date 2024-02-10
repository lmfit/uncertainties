import sys
import itertools

from math import sqrt, log, isnan, isinf  # Optimization: no attribute look-up
# The following restricts the local function getargspec() to the common
# features of inspect.getargspec() and inspect.getfullargspec():
if sys.version_info < (3,):  # !! Could be removed when moving to Python 3 only
    from inspect import getargspec
else:
    from inspect import getfullargspec as getargspec

from .. util import (set_doc, covariance_matrix, partial_derivative,
NumericalDerivatives, IndexableIter, basestring,
)
from . base import AffineScalarFuncBase, LinearCombination, NotUpcast
from .. compat import FLOAT_LIKE_TYPES

########################################

# Definition of boolean operators, that assume that self and
# y_with_uncert are AffineScalarFunc.

# The fact that uncertainties must be small is used, here: the
# comparison functions are supposed to be constant for most values of
# the random variables.

# Even though uncertainties are supposed to be small, comparisons
# between 3+/-0.1 and 3.0 are handled correctly (even though x == 3.0 is
# not a constant function in the 3+/-0.1 interval).  The comparison
# between x and x is handled too, when x has an uncertainty.  In fact,
# as explained in the main documentation, it is possible to give a
# useful meaning to the comparison operators, in these cases.

def eq_on_aff_funcs(self, y_with_uncert):
    """
    __eq__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """
    difference = self - y_with_uncert
    # Only an exact zero difference means that self and y are
    # equal numerically:
    return not(difference._nominal_value or difference.std_dev)

def ne_on_aff_funcs(self, y_with_uncert):
    """
    __ne__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """

    return not eq_on_aff_funcs(self, y_with_uncert)

def gt_on_aff_funcs(self, y_with_uncert):
    """
    __gt__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """
    return self._nominal_value > y_with_uncert._nominal_value

def ge_on_aff_funcs(self, y_with_uncert):
    """
    __ge__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """

    return (gt_on_aff_funcs(self, y_with_uncert)
            or eq_on_aff_funcs(self, y_with_uncert))

def lt_on_aff_funcs(self, y_with_uncert):
    """
    __lt__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """
    return self._nominal_value < y_with_uncert._nominal_value

def le_on_aff_funcs(self, y_with_uncert):
    """
    __le__ operator, assuming that both self and y_with_uncert are
    AffineScalarFunc objects.
    """

    return (lt_on_aff_funcs(self, y_with_uncert)
            or eq_on_aff_funcs(self, y_with_uncert))

def req_on_aff_funcs(self, y_with_uncert):
    return eq_on_aff_funcs(y_with_uncert, self)

def rne_on_aff_funcs(self, y_with_uncert):
    return ne_on_aff_funcs(y_with_uncert, self)

def rgt_on_aff_funcs(self, y_with_uncert):
    return gt_on_aff_funcs(y_with_uncert, self)

def rge_on_aff_funcs(self, y_with_uncert):
    return ge_on_aff_funcs(y_with_uncert, self)

def rlt_on_aff_funcs(self, y_with_uncert):
    return lt_on_aff_funcs(y_with_uncert, self)

def rle_on_aff_funcs(self, y_with_uncert):
    return le_on_aff_funcs(y_with_uncert, self)


class AffineScalarFuncOps(AffineScalarFuncBase):

    ############################################################
    # Operators: operators applied to AffineScalarFunc and/or
    # float-like objects only are supported.  This is why methods
    # from float are used for implementing these operators.

    # Operators with no reflection:
    ########################################

    # __nonzero__() is supposed to return a boolean value (it is used
    # by bool()).  It is for instance used for converting the result
    # of comparison operators to a boolean, in sorted().  If we want
    # to be able to sort AffineScalarFunc objects, __nonzero__ cannot
    # return a AffineScalarFunc object.  Since boolean results (such
    # as the result of bool()) don't have a very meaningful
    # uncertainty unless it is zero, this behavior is fine.
    def __bool__(self):
        """
        Equivalent to self != 0.
        """
        #! This might not be relevant for AffineScalarFunc objects
        # that contain values in a linear space which does not convert
        # the float 0 into the null vector (see the __eq__ function:
        # __nonzero__ works fine if subtracting the 0 float from a
        # vector of the linear space works as if 0 were the null
        # vector of that space):
        return self != 0.  # Uses the AffineScalarFunc.__ne__ function
    
    def force_aff_func_args(self, func):
        """
        Takes an operator op(x, y) and wraps it.

        The constructed operator returns func(x, to_affine_scalar(y)) if y
        can be upcast with to_affine_scalar(); otherwise, it returns
        NotImplemented.

        Thus, func() is only called on two AffineScalarFunc objects, if
        its first argument is an AffineScalarFunc.
        """
        def op_on_upcast_args(x, y):
            """
            Return %s(self, to_affine_scalar(y)) if y can be upcast
            through to_affine_scalar.  Otherwise returns NotImplemented.
            """ % func.__name__

            try:
                y_with_uncert = self._to_affine_scalar(y)
            except NotUpcast:
                # This module does not know how to handle the comparison:
                # (example: y is a NumPy array, in which case the NumPy
                # array will decide that func() should be applied
                # element-wise between x and all the elements of y):
                return NotImplemented
            else:
                return func(x, y_with_uncert)

        return op_on_upcast_args
    ## Logical operators: warning: the resulting value cannot always
    ## be differentiated.

    # The boolean operations are not differentiable everywhere, but
    # almost...

    # (1) I can rely on the assumption that the user only has "small"
    # errors on variables, as this is used in the calculation of the
    # standard deviation (which performs linear approximations):

    # (2) However, this assumption is not relevant for some
    # operations, and does not have to hold, in some cases.  This
    # comes from the fact that logical operations (e.g. __eq__(x,y))
    # are not differentiable for many usual cases.  For instance, it
    # is desirable to have x == x for x = n+/-e, whatever the size of e.
    # Furthermore, n+/-e != n+/-e', if e != e', whatever the size of e or
    # e'.

    # (3) The result of logical operators does not have to be a
    # function with derivatives, as these derivatives are either 0 or
    # don't exist (i.e., the user should probably not rely on
    # derivatives for his code).

    # !! In Python 2.7+, it may be possible to use functools.total_ordering.

    # __eq__ is used in "if data in [None, ()]", for instance.  It is
    # therefore important to be able to handle this case too, which is
    # taken care of when force_aff_func_args(eq_on_aff_funcs)
    # returns NotImplemented.
    ########################################
    @classmethod
    def _add_comparative_ops(cls):
        cls.__eq__ = cls.force_aff_func_args(cls, eq_on_aff_funcs)

        cls.__ne__ = cls.force_aff_func_args(cls, ne_on_aff_funcs)
        cls.__gt__ = cls.force_aff_func_args(cls, gt_on_aff_funcs)

        # __ge__ is not the opposite of __lt__ because these operators do
        # not always yield a boolean (for instance, 0 <= numpy.arange(10)
        # yields an array).
        cls.__ge__ = cls.force_aff_func_args(cls, ge_on_aff_funcs)

        cls.__lt__ = cls.force_aff_func_args(cls, lt_on_aff_funcs)
        cls.__le__ = cls.force_aff_func_args(cls, le_on_aff_funcs)

        cls.__req__ = cls.force_aff_func_args(cls, req_on_aff_funcs)
        cls.__rne__ = cls.force_aff_func_args(cls, rne_on_aff_funcs)
        cls.__rgt__ = cls.force_aff_func_args(cls, rgt_on_aff_funcs)
        cls.__rge__ = cls.force_aff_func_args(cls, rge_on_aff_funcs)
        cls.__rlt__ = cls.force_aff_func_args(cls, rlt_on_aff_funcs)
        cls.__rle__ = cls.force_aff_func_args(cls, rle_on_aff_funcs)
        
    @classmethod
    def wrap(cls, f, derivatives_args=[], derivatives_kwargs={}):
        """
        Wraps a function f into a function that also accepts numbers with
        uncertainties (UFloat objects); the wrapped function returns the
        value of f with the correct uncertainty and correlations. The
        wrapped function is intended to be used as a drop-in replacement
        for the original function: they can be called in the exact same
        way, the only difference being that numbers with uncertainties can
        be given to the wrapped function where f accepts float arguments.

        Doing so may be necessary when function f cannot be expressed
        analytically (with uncertainties-compatible operators and
        functions like +, *, umath.sin(), etc.).

        f must return a float-like (i.e. a float, an int, etc., not a
        list, etc.), unless when called with no number with
        uncertainty. This is because the wrapped function generally
        returns numbers with uncertainties: they represent a probability
        distribution over the real numbers.

        If the wrapped function is called with no argument that has an
        uncertainty, the value of f is returned.

        Parameters: the derivatives_* parameters can be used for defining
        some of the partial derivatives of f. All the (non-None)
        derivatives must have the same signature as f.

        derivatives_args --

            Iterable that, when iterated over, returns either derivatives
            (functions) or None. derivatives_args can in particular be a
            simple sequence (list or tuple) that gives the derivatives of
            the first positional parameters of f.

            Each function must be the partial derivative of f with respect
            to the corresponding positional parameters.  These functions
            take the same arguments as f.

            The positional parameters of a function are usually
            positional-or-keyword parameters like in the call func(a,
            b=None). However, they also include var-positional parameters
            given through the func(a, b, *args) *args syntax. In the last
            example, derivatives_args can be an iterable that returns the
            derivative with respect to a, b and then to each optional
            argument in args.

            A value of None (instead of a function) obtained when
            iterating over derivatives_args is automatically replaced by
            the relevant numerical derivative. This derivative is not used
            if the corresponding argument is not a number with
            uncertainty. A None value can therefore be used for non-scalar
            arguments of f (like string arguments).

            If the derivatives_args iterable yields fewer derivatives than
            needed, wrap() automatically sets the remaining unspecified
            derivatives to None (i.e. to the automatic numerical
            calculation of derivatives).

            An indefinite number of derivatives can be specified by having
            derivatives_args be an infinite iterator; this can for
            instance be used for specifying the derivatives of functions
            with a undefined number of argument (like sum(), whose partial
            derivatives all return 1).

        derivatives_kwargs --

            Dictionary that maps keyword parameters to their derivatives,
            or None (as in derivatives_args).

            Keyword parameters are defined as being those of kwargs when f
            has a signature of the form f(..., **kwargs). In Python 3,
            these keyword parameters also include keyword-only parameters.

            Non-mapped keyword parameters are replaced automatically by
            None: the wrapped function will use, if necessary, numerical
            differentiation for these parameters (as with
            derivatives_args).

            Note that this dictionary only maps keyword *parameters* from
            the *signature* of f. The way the wrapped function is called
            is immaterial: for example, if f has signature f(a, b=None),
            then derivatives_kwargs should be the empty dictionary, even
            if the wrapped f can be called a wrapped_f(a=123, b=42).

        Example (for illustration purposes only, as
        uncertainties.umath.sin() runs faster than the examples that
        follow): wrap(math.sin) is a sine function that can be applied to
        numbers with uncertainties.  Its derivative will be calculated
        numerically.  wrap(math.sin, [None]) would have produced the same
        result.  wrap(math.sin, [math.cos]) is the same function, but with
        an analytically defined derivative.

        Numerically calculated derivatives are meaningless when the
        function is not differentiable (e.g., math.hypot(x, y) in (x, y) =
        (0, 0), and sqrt(x) in x = 0). The corresponding uncertainties are
        either meaningless (case of hypot) or raise an exception when
        calculated (case of sqrt). In such cases, it is recommended (but
        not mandatory) to supply instead a derivative function that
        returns NaN where the function is not differentiable. This
        function can still numerically calculate the derivative where
        defined, for instance by using the
        uncertainties.core.partial_derivative() function.

        The correctness of the supplied analytical derivatives an be
        tested by setting them to None instead and comparing the
        analytical and the numerical differentiation results.

        Note on efficiency: the wrapped function assumes that f cannot
        accept numbers with uncertainties as arguments. If f actually does
        handle some arguments even when they have an uncertainty, the
        wrapped function ignores this fact, which might lead to a
        performance hit: wrapping a function that actually accepts numbers
        with uncertainty is likely to make it slower.
        """
        
        derivatives_args_index = IndexableIter(
            # Automatic addition of numerical derivatives in case the
            # supplied derivatives_args is shorter than the number of
            # arguments in *args:
            itertools.chain(derivatives_args, itertools.repeat(None)))


        # Derivatives for keyword arguments (includes var-keyword
        # parameters **kwargs, but also var-or-keyword parameters, and
        # keyword-only parameters (Python 3):

        derivatives_all_kwargs = {}

        for (name, derivative) in derivatives_kwargs.items():

            # Optimization: None keyword-argument derivatives are converted
            # right away to derivatives (instead of doing this every time a
            # None derivative is encountered when calculating derivatives):

            if derivative is None:
                derivatives_all_kwargs[name] = partial_derivative(f, name)
            else:
                derivatives_all_kwargs[name] = derivative

        # When the wrapped function is called with keyword arguments that
        # map to positional-or-keyword parameters, their derivative is
        # looked for in derivatives_all_kwargs.  We define these
        # additional derivatives:

        try:
            argspec = getargspec(f)
        except TypeError:
            # Some functions do not provide meta-data about their
            # arguments (see PEP 362). One cannot use keyword arguments
            # for positional-or-keyword parameters with them: nothing has
            # to be done:
            pass
        else:
            # With Python 3, there is no need to handle keyword-only
            # arguments (and therefore to use inspect.getfullargspec())
            # because they are already handled by derivatives_kwargs.

            for (index, name) in enumerate(argspec.args):

                # The following test handles the case of
                # positional-or-keyword parameter for which automatic
                # numerical differentiation is used: when the wrapped
                # function is called with a keyword argument for this
                # parameter, the numerical derivative must be calculated
                # with respect to the parameter name. In the other case,
                # where the wrapped function is called with a positional
                # argument, the derivative with respect to its index must
                # be used:

                derivative = derivatives_args_index[index]

                if derivative is None:
                    derivatives_all_kwargs[name] = partial_derivative(f, name)
                else:
                    derivatives_all_kwargs[name] = derivative

        # Optimization: None derivatives for the positional arguments are
        # converted to the corresponding numerical differentiation
        # function (instead of doing this over and over later every time a
        # None derivative is found):

        none_converter = lambda index: partial_derivative(f, index)

        for (index, derivative) in enumerate(
            derivatives_args_index.returned_elements):
            if derivative is None:
                derivatives_args_index.returned_elements[index] = (
                    none_converter(index))

        # Future None values are also automatically converted:
        derivatives_args_index.none_converter = none_converter


        ## Wrapped function:

        #! Setting the doc string after "def f_with...()" does not
        # seem to work.  We define it explicitly:
        @set_doc("""\
        Version of %s(...) that returns an affine approximation
        (AffineScalarFunc object), if its result depends on variables
        (Variable objects).  Otherwise, returns a simple constant (when
        applied to constant arguments).

        Warning: arguments of the function that are not AffineScalarFunc
        objects must not depend on uncertainties.Variable objects in any
        way.  Otherwise, the dependence of the result in
        uncertainties.Variable objects will be incorrect.

        Original documentation:
        %s""" % (f.__name__, f.__doc__))
        def f_with_affine_output(*args, **kwargs):
            
            ########################################
            # The involved random variables must first be gathered, so
            # that they can be independently updated.

            # The arguments that contain an uncertainty (AffineScalarFunc
            # objects) are gathered, as positions or names; they will be
            # replaced by their nominal value in order to calculate
            # the necessary derivatives of f.

            pos_w_uncert = [index for (index, value) in enumerate(args)
                            if isinstance(value, cls)]
            names_w_uncert = [key for (key, value) in kwargs.items()
                            if isinstance(value, cls)]

            ########################################
            # Value of f() at the nominal value of the arguments with
            # uncertainty:

            # The usual behavior of f() is kept, if no number with
            # uncertainty is provided:
            if (not pos_w_uncert) and (not names_w_uncert):
                return f(*args, **kwargs)

            ### Nominal values of the (scalar) arguments:

            # !! Possible optimization: If pos_w_uncert is empty, there
            # is actually no need to create a mutable version of args and
            # one could do args_values = args.  However, the wrapped
            # function is typically called with numbers with uncertainties
            # as positional arguments (i.e., pos_w_uncert is not emtpy),
            # so this "optimization" is not implemented here.

            ## Positional arguments:
            args_values = list(args)  # Now mutable: modified below
            # Arguments with an uncertainty are converted to their nominal
            # value:
            for index in pos_w_uncert:
                args_values[index] = args[index].nominal_value

            ## Keyword arguments:

            # For efficiency reasons, kwargs is not copied. Instead, its
            # values with uncertainty are modified:

            # The original values with uncertainties are needed: they are
            # saved in the following dictionary (which only contains
            # values with uncertainty):

            kwargs_uncert_values = {}

            for name in names_w_uncert:
                value_with_uncert = kwargs[name]
                # Saving for future use:
                kwargs_uncert_values[name] = value_with_uncert
                # The original dictionary is modified (for efficiency reasons):
                kwargs[name] = value_with_uncert.nominal_value

            f_nominal_value = f(*args_values, **kwargs)
        
            # If the value is not a float, then this code cannot provide
            # the result, as it returns a UFloat, which represents a
            # random real variable. This happens for instance when
            # ufloat()*numpy.array() is calculated: the
            # AffineScalarFunc.__mul__ operator, obtained through wrap(),
            # returns a NumPy array, not a float:
            if not isinstance(f_nominal_value, FLOAT_LIKE_TYPES):
                print(args, kwargs)
                print("467", f_nominal_value, type(f_nominal_value))
                return NotImplemented

            ########################################

            # Calculation of the linear part of the function value,
            # defined by (coefficient, argument) pairs, where 'argument'
            # is an AffineScalarFunc (for all AffineScalarFunc found as
            # argument of f):
            linear_part = []

            for pos in pos_w_uncert:
                linear_part.append((
                    # Coefficient:
                    derivatives_args_index[pos](*args_values, **kwargs),
                    # Linear part of the AffineScalarFunc expression:
                    args[pos]._linear_part))

            for name in names_w_uncert:

                # Optimization: caching of the automatic numerical
                # derivatives for keyword arguments that are
                # discovered. This gives a speedup when the original
                # function is called repeatedly with the same keyword
                # arguments:
                derivative = derivatives_all_kwargs.setdefault(
                    name,
                    # Derivative never needed before:
                    partial_derivative(f, name))

                linear_part.append((
                    # Coefficient:
                    derivative(*args_values, **kwargs),
                    # Linear part of the AffineScalarFunc expression:
                    kwargs_uncert_values[name]._linear_part))

            # The function now returns the necessary linear approximation
            # to the function:
            return cls(
                f_nominal_value, LinearCombination(linear_part))

        f_with_affine_output = set_doc("""\
        Version of %s(...) that returns an affine approximation
        (AffineScalarFunc object), if its result depends on variables
        (Variable objects).  Otherwise, returns a simple constant (when
        applied to constant arguments).

        Warning: arguments of the function that are not AffineScalarFunc
        objects must not depend on uncertainties.Variable objects in any
        way.  Otherwise, the dependence of the result in
        uncertainties.Variable objects will be incorrect.

        Original documentation:
        %s""" % (f.__name__, f.__doc__))(f_with_affine_output)
        
        # It is easier to work with f_with_affine_output, which represents
        # a wrapped version of 'f', when it bears the same name as 'f':
        # ! __name__ is read-only, in Python 2.3:
        f_with_affine_output.name = f.__name__

        return f_with_affine_output
    
    @classmethod
    def _add_arithmetic_ops(cls):
        """
        Adds many operators (__add__, etc.) to the AffineScalarFunc class.
        """

        ########################################

        #! Derivatives are set to return floats.  For one thing,
        # uncertainties generally involve floats, as they are based on
        # small variations of the parameters.  It is also better to
        # protect the user from unexpected integer result that behave
        # badly with the division.

        ## Operators that return a numerical value:

        def _simple_add_deriv(x):
            if x >= 0:
                return 1.
            else:
                return -1.
        
        cls.__abs__ = cls.wrap(float.__abs__, [_simple_add_deriv])
        cls.__neg__ = cls.wrap(float.__neg__, [lambda x: -1.])
        cls.__pos__ = cls.wrap(float.__pos__, [lambda x: 1.])
        cls.__trunc__ = cls.wrap(float.__trunc__, [lambda x: 0.])

        ########################################
        # Final definition of the operators for AffineScalarFunc objects:

        cls.__add__ = cls.wrap(float.__add__, [lambda x, y: 1., lambda x, y: 1.])
        cls.__radd__ = cls.wrap(float.__radd__, [lambda x, y: 1., lambda x, y: 1.])
        cls.__sub__ = cls.wrap(float.__sub__, [lambda x, y: 1., lambda x, y: -1.])
        cls.__rsub__ = cls.wrap(float.__rsub__, [lambda y, x: -1., lambda y, x: 1.])
        cls.__mul__ = cls.wrap(float.__mul__, [lambda x, y: y, lambda x, y: x])
        cls.__rmul__ = cls.wrap(float.__rmul__, [lambda y, x: x, lambda y, x: y])
        # cls.__div__ = cls.wrap(float.__div__, ['1/y', '-x/y**2'])
        # cls.__rdiv__ = cls.wrap(float.__rdiv__, ['-1/x', 'y/x**2'])
        cls.__truediv__ = cls.wrap(float.__truediv__, [lambda x, y: 1/y, lambda x, y: -x/y**2])
        cls.__rtruediv__ = cls.wrap(float.__rtruediv__, [lambda y, x: -x/y**2, lambda y, x: 1/y])
        cls.__floordiv__ = cls.wrap(float.__floordiv__, [lambda x, y: 0., lambda x, y: 0.])
        cls.__rfloordiv__ = cls.wrap(float.__rfloordiv__, [lambda y, x: 0., lambda y, x: 0.])
        cls.__mod__ = cls.wrap(float.__mod__, [lambda x, y: 1., lambda x, y: partial_derivative(float.__mod__, 1)(x, y)])
        cls.__rmod__ = cls.wrap(float.__rmod__, [lambda y, x: partial_derivative(float.__mod__, 1)(x, y), lambda y, x: 1.])
                
        def pow_deriv_0(x, y):
            if y == 0:
                return 0.
            elif x != 0 or y % 1 == 0:
                return y*x**(y-1)
            else:
                return float('nan')

        def pow_deriv_1(x, y):
            if x == 0 and y > 0:
                return 0.
            else:
                return log(x)*x**y
        
        pow_deriv_0 = nan_if_exception(pow_deriv_0)
        pow_deriv_1 = nan_if_exception(pow_deriv_1)

        # This module does not handle uncertainties on complex numbers:
        # complex results for the nominal value of some operations cannot
        # be calculated with an uncertainty:
        def no_complex_result(func):
            '''
            Return a function that does like func, but that raises a
            ValueError if the result is complex.
            '''
            def no_complex_func(*args, **kwargs):
                '''
                Like %s, but raises a ValueError exception if the result
                is complex.
                ''' % func.__name__

                value = func(*args, **kwargs)
                if isinstance(value, complex):
                    raise ValueError('The uncertainties module does not handle'
                                    ' complex results')
                else:
                    return value

            return no_complex_func

        cls.__pow__ = cls.wrap(no_complex_result(float.__pow__), [pow_deriv_0, pow_deriv_1])
        cls.__rpow__ = cls.wrap(no_complex_result(float.__rpow__), [lambda y, x: pow_deriv_1(x, y),
                                   lambda y, x: pow_deriv_0(x, y)])


    
###############################################################################

# Some operators can have undefined derivatives but still give
# meaningful values when some of their arguments have a zero
# uncertainty. Such operators return NaN when their derivative is
# not finite. This way, if the uncertainty of the associated
# variable is not 0, a NaN uncertainty is produced, which
# indicates an error; if the uncertainty is 0, then the total
# uncertainty can be returned as 0.

# Exception catching is used so as to not slow down regular
# operation too much:

def nan_if_exception(f):
    '''
    Wrapper around f(x, y) that let f return NaN when f raises one of
    a few numerical exceptions.
    '''

    def wrapped_f(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except (ValueError, ZeroDivisionError, OverflowError):
            return float('nan')

    return wrapped_f

