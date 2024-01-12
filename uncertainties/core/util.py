import math 
import sys

if sys.version_info < (3,):
     from past.builtins import basestring
else:
     # Avoid importing from past in Python 3 since it utilizes the builtin
     # 'imp' module, which is deprecated as of Python 3.4, see
     # https://docs.python.org/3/library/imp.html. The 2to3 tool replaces
     # basestring with str, so that's what we effectively do here as well:
     basestring = str


def set_doc(doc_string):
    """
    Decorator function that sets the docstring to the given text.

    It is useful for functions whose docstring is calculated
    (including string substitutions).
    """
    def set_doc_string(func):
        func.__doc__ = doc_string
        return func
    return set_doc_string


try:
    import numpy
except ImportError:
    pass
else:

    # Entering variables as a block of correlated values.  Only available
    # if NumPy is installed.

    #! It would be possible to dispense with NumPy, but a routine should be
    # written for obtaining the eigenvectors of a symmetric matrix.  See
    # for instance Numerical Recipes: (1) reduction to tri-diagonal
    # [Givens or Householder]; (2) QR / QL decomposition.

    def correlated_values(nom_values, covariance_mat, tags=None):
        """
        Return numbers with uncertainties (AffineScalarFunc objects)
        that correctly reproduce the given covariance matrix, and have
        the given (float) values as their nominal value.

        The correlated_values_norm() function returns the same result,
        but takes a correlation matrix instead of a covariance matrix.

        The list of values and the covariance matrix must have the
        same length, and the matrix must be a square (symmetric) one.

        The numbers with uncertainties returned depend on newly
        created, independent variables (Variable objects).

        nom_values -- sequence with the nominal (real) values of the
        numbers with uncertainties to be returned.

        covariance_mat -- full covariance matrix of the returned numbers with
        uncertainties. For example, the first element of this matrix is the
        variance of the first number with uncertainty. This matrix must be a
        NumPy array-like (list of lists, NumPy array, etc.). 

        tags -- if 'tags' is not None, it must list the tag of each new
        independent variable.
        """

        # !!! It would in principle be possible to handle 0 variance
        # variables by first selecting the sub-matrix that does not contain
        # such variables (with the help of numpy.ix_()), and creating 
        # them separately.
        
        std_devs = numpy.sqrt(numpy.diag(covariance_mat))

        # For numerical stability reasons, we go through the correlation
        # matrix, because it is insensitive to any change of scale in the
        # quantities returned. However, care must be taken with 0 variance
        # variables: calculating the correlation matrix cannot be simply done
        # by dividing by standard deviations. We thus use specific
        # normalization values, with no null value:
        norm_vector = std_devs.copy()
        norm_vector[norm_vector==0] = 1

        return correlated_values_norm(
            # !! The following zip() is a bit suboptimal: correlated_values()
            # separates back the nominal values and the standard deviations:
            list(zip(nom_values, std_devs)),
            covariance_mat/norm_vector/norm_vector[:,numpy.newaxis],
            tags)

    def correlated_values_norm(values_with_std_dev, correlation_mat,
                               tags=None):
        '''
        Return correlated values like correlated_values(), but takes
        instead as input:

        - nominal (float) values along with their standard deviation, and
        - a correlation matrix (i.e. a normalized covariance matrix).

        values_with_std_dev -- sequence of (nominal value, standard
        deviation) pairs. The returned, correlated values have these
        nominal values and standard deviations.

        correlation_mat -- correlation matrix between the given values, except
        that any value with a 0 standard deviation must have its correlations
        set to 0, with a diagonal element set to an arbitrary value (something
        close to 0-1 is recommended, for a better numerical precision).  When
        no value has a 0 variance, this is the covariance matrix normalized by
        standard deviations, and thus a symmetric matrix with ones on its
        diagonal.  This matrix must be an NumPy array-like (list of lists,
        NumPy array, etc.).

        tags -- like for correlated_values().
        '''

        # If no tags were given, we prepare tags for the newly created
        # variables:
        if tags is None:
            tags = (None,) * len(values_with_std_dev)

        (nominal_values, std_devs) = numpy.transpose(values_with_std_dev)

        # We diagonalize the correlation matrix instead of the
        # covariance matrix, because this is generally more stable
        # numerically. In fact, the covariance matrix can have
        # coefficients with arbitrary values, through changes of units
        # of its input variables. This creates numerical instabilities.
        #
        # The covariance matrix is diagonalized in order to define
        # the independent variables that model the given values:
        (variances, transform) = numpy.linalg.eigh(correlation_mat)

        # Numerical errors might make some variances negative: we set
        # them to zero:
        variances[variances < 0] = 0.

        # Creation of new, independent variables:

        # We use the fact that the eigenvectors in 'transform' are
        # special: 'transform' is unitary: its inverse is its transpose:

        variables = tuple(
            # The variables represent "pure" uncertainties:
            Variable(0, sqrt(variance), tag)
            for (variance, tag) in zip(variances, tags))

        # The coordinates of each new uncertainty as a function of the
        # new variables must include the variable scale (standard deviation):
        transform *= std_devs[:, numpy.newaxis] 
        
        # Representation of the initial correlated values:
        values_funcs = tuple(
            AffineScalarFunc(
                value,
                LinearCombination(dict(zip(variables, coords))))
            for (coords, value) in zip(transform, nominal_values))

        return values_funcs
        
    def correlation_matrix(nums_with_uncert):
        '''
        Return the correlation matrix of the given sequence of
        numbers with uncertainties, as a NumPy array of floats.
        '''

        cov_mat = numpy.array(covariance_matrix(nums_with_uncert))

        std_devs = numpy.sqrt(cov_mat.diagonal())

        return cov_mat/std_devs/std_devs[numpy.newaxis].T


def covariance_matrix(nums_with_uncert):
    """
    Return a matrix that contains the covariances between the given
    sequence of numbers with uncertainties (AffineScalarFunc objects).
    The resulting matrix implicitly depends on their ordering in
    'nums_with_uncert'.

    The covariances are floats (never int objects).

    The returned covariance matrix is the exact linear approximation
    result, if the nominal values of the numbers with uncertainties
    and of their variables are their mean.  Otherwise, the returned
    covariance matrix should be close to its linear approximation
    value.

    The returned matrix is a list of lists.
    """
    # See PSI.411 in EOL's notes.

    covariance_matrix = []
    for (i1, expr1) in enumerate(nums_with_uncert, 1):
        derivatives1 = expr1.derivatives  # Optimization
        vars1 = set(derivatives1)  # !! Python 2.7+: viewkeys() would work
        coefs_expr1 = []

        for expr2 in nums_with_uncert[:i1]:
            derivatives2 = expr2.derivatives  # Optimization
            coefs_expr1.append(sum(
                ((derivatives1[var]*derivatives2[var]*var._std_dev**2)
                # var is a variable common to both numbers with
                # uncertainties:
                for var in vars1.intersection(derivatives2)),
                # The result is always a float (sum() with no terms
                # returns an integer):
                0.))

        covariance_matrix.append(coefs_expr1)

    # We symmetrize the matrix:
    for (i, covariance_coefs) in enumerate(covariance_matrix):
        covariance_coefs.extend([covariance_matrix[j][i]
                                 for j in range(i+1, len(covariance_matrix))])

    return covariance_matrix


# Step constant for numerical derivatives in
# partial_derivative(). Value chosen to as to get better numerical
# results:
STEP_SIZE = math.sqrt(sys.float_info.epsilon)

# !! It would be possible to split the partial derivative calculation
# into two functions: one for positional arguments (case of integer
# arg_ref) and one for keyword arguments (case of string
# arg_ref). However, this would either duplicate the code for the
# numerical differentiation, or require a call, which is probably more
# expensive in time than the tests done here.
def partial_derivative(f, arg_ref):
    """
    Return a function that numerically calculates the partial
    derivative of function f with respect to its argument arg_ref.

    arg_ref -- describes which variable to use for the
    differentiation. If f is called with f(*args, **kwargs) arguments,
    an integer represents the index of an argument in args, and a
    string represents the name of an argument in kwargs.
    """

    # Which set of function parameter contains the variable to be
    # changed? the positional or the optional keyword arguments?
    change_kwargs = isinstance(arg_ref, basestring)

    def partial_derivative_of_f(*args, **kwargs):
        """
        Partial derivative, calculated with the (-epsilon, +epsilon)
        method, which is more precise than the (0, +epsilon) method.
        """

        # args_with_var contains the arguments (either args or kwargs)
        # that contain the variable that must be shifted, as a mutable
        # object (because the variable contents will be modified):

        # The values in args need to be modified, for the
        # differentiation: it is converted to a list:
        if change_kwargs:
            args_with_var = kwargs
        else:
            args_with_var = list(args)

        # The step is relative to the parameter being varied, so that
        # shifting it does not suffer from finite precision limitations:
        step = STEP_SIZE*abs(args_with_var[arg_ref])
        if not step:
            # Arbitrary, but "small" with respect to 1:
            step = STEP_SIZE

        args_with_var[arg_ref] += step

        if change_kwargs:
            shifted_f_plus = f(*args, **args_with_var)
        else:
            shifted_f_plus = f(*args_with_var, **kwargs)

        args_with_var[arg_ref] -= 2*step  # Optimization: only 1 list copy

        if change_kwargs:
            shifted_f_minus = f(*args, **args_with_var)
        else:
            shifted_f_minus = f(*args_with_var, **kwargs)

        return (shifted_f_plus - shifted_f_minus)/2/step

    return partial_derivative_of_f



class NumericalDerivatives(object):
    """
    Convenient access to the partial derivatives of a function,
    calculated numerically.
    """
    # This is not a list because the number of arguments of the
    # function is not known in advance, in general.

    def __init__(self, function):
        """
        'function' is the function whose derivatives can be computed.
        """
        self._function = function

    def __getitem__(self, n):
        """
        Return the n-th numerical derivative of the function.
        """
        return partial_derivative(self._function, n)

class IndexableIter(object):
    '''
    Iterable whose values can also be accessed through indexing.

    The input iterable values are cached.

    Some attributes:

    iterable -- iterable used for returning the elements one by one.

    returned_elements -- list with the elements directly accessible.
    through indexing. Additional elements are obtained from self.iterable.

    none_converter -- function that takes an index and returns the
    value to be returned when None is obtained form the iterable
    (instead of None).
    '''

    def __init__(self, iterable, none_converter=lambda index: None):
        '''
        iterable -- iterable whose values will be returned.

        none_converter -- function applied to None returned
        values. The value that replaces None is none_converter(index),
        where index is the index of the element.
        '''
        self.iterable = iterable
        self.returned_elements = []
        self.none_converter = none_converter

    def __getitem__(self, index):

        returned_elements = self.returned_elements

        try:

            return returned_elements[index]

        except IndexError:  # Element not yet cached

            for pos in range(len(returned_elements), index+1):

                value = next(self.iterable)

                if value is None:
                    value = self.none_converter(pos)

                returned_elements.append(value)

            return returned_elements[index]

    def __str__(self):
        return '<%s: [%s...]>' % (
            self.__class__.__name__,
            ', '.join(map(str, self.returned_elements)))

