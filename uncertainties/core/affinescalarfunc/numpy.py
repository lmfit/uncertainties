
import numpy as np
# ufuncs are listed at https://numpy.org/doc/stable/reference/ufuncs.html


ufunc_derivatives = {
    'add': [lambda x, y: 1., lambda x, y: 1.],
    'subtract': [lambda x, y: 1., lambda x, y: -1.],
    'multiply': [lambda x, y: y, lambda x, y: x],
    'divide': [lambda x, y: 1./y, lambda x, y: -x/y**2],
    'true_divide': [lambda x, y: 1./y, lambda x, y: -x/y**2],
    'floor_divide': [lambda x, y: 0., lambda x, y: 0.],
    }

def is_upcast_type(t):
    # This can be used to allow downstream modules to overide operations; see pint
    # TODO add upcast_type list or dict to a public interface
    return False

def implements(numpy_func_string, func_type):
    """Register an __array_function__/__array_ufunc__ implementation for UArray
    objects.
    """
    print(numpy_func_string, func_type)

    def decorator(func):
        if func_type == "function":
            HANDLED_FUNCTIONS[numpy_func_string] = func
        elif func_type == "ufunc":
            HANDLED_UFUNCS[numpy_func_string] = func
        else:
            raise ValueError(f"Invalid func_type {func_type}")
        return func

    return decorator

HANDLED_FUNCTIONS = {}
HANDLED_UFUNCS = {}


def apply_func_elementwise(func, inputs, kwargs, result_dtype="object"):
    if isinstance(inputs[0], np.ndarray):
        result = np.empty_like(inputs[0], dtype=result_dtype)
        for index, x in np.ndenumerate(inputs[0]):
            inputs_ = [x if i == 0 else inputs[i] for i in range(len(inputs))]
            result[index] = func(*inputs_, **kwargs)
    elif isinstance(inputs[1], np.ndarray):
        result = np.empty_like(inputs[1], dtype=result_dtype)
        for index, x in np.ndenumerate(inputs[1]):
            inputs_ = [x if i == 1 else inputs[i] for i in range(len(inputs))]
            result[index] = func(*inputs_, **kwargs)
    else:
        result = func(*inputs, **kwargs)
    return result

def numpy_wrap(func_type, func, args, kwargs, types):
    """Return the result from a NumPy function/ufunc as wrapped by uncertainties."""

    if func_type == "function":
        handled = HANDLED_FUNCTIONS
        # Need to handle functions in submodules
        name = ".".join(func.__module__.split(".")[1:] + [func.__name__])
    elif func_type == "ufunc":
        handled = HANDLED_UFUNCS
        # ufuncs do not have func.__module__
        name = func.__name__
    else:
        raise ValueError(f"Invalid func_type {func_type}")

    if name not in handled or any(is_upcast_type(t) for t in types):
        print("NotImplemented L54")
        return NotImplemented
    return handled[name](*args, **kwargs)

class UFloatNumpy(object):
    # NumPy function/ufunc support
    __array_priority__ = 17

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__":
            # Only handle ufuncs as callables
            return NotImplemented
        
        # Replicate types from __array_function__
        types = {
            type(arg)
            for arg in list(inputs) + list(kwargs.values())
            if hasattr(arg, "__array_ufunc__")
        }

        return numpy_wrap("ufunc", ufunc, inputs, kwargs, types)

    def __array_function__(self, func, types, args, kwargs):
        return numpy_wrap("function", func, args, kwargs, types)
  
    # original code for _add_numpy_ufuncs. may be helpful for writing a generic wraps
    # this can be deleted:
    # @classmethod
    # def _add_numpy_ufuncs(cls):
    #     def implement_ufunc(func_str, derivatives):
    #         func = getattr(np, func_str)
    #         @implements(func_str, "ufunc")
    #         def implementation(*inputs, **kwargs):
    #             if isinstance(inputs[0], np.ndarray):
    #                 result = np.empty_like(inputs[0], dtype=object)
    #                 for index, x in np.ndenumerate(inputs[0]):
    #                     inputs_ = (x if i == 0 else inputs[i] for i in range(len(inputs)))
    #                     result[index] = cls.wrap(func, derivatives)(*inputs_, **kwargs)
    #             elif isinstance(inputs[1], np.ndarray):
    #                 result = np.empty_like(inputs[1], dtype=object)
    #                 for index, x in np.ndenumerate(inputs[1]):
    #                     inputs_ = (x if i == 1 else inputs[i] for i in range(len(inputs)))
    #                     result[index] = cls.wrap(func, derivatives)(*inputs_, **kwargs)
    #             else:
    #                 result = cls.wrap(func, derivatives)(*inputs, **kwargs)
    #             return result
            
    #         return implementation

    #     for func_str, derivatives in ufunc_derivatives.items():
    #         implement_ufunc(func_str, derivatives)

    
    @classmethod
    def _add_numpy_ufuncs(cls):
        def implement_ufunc(func_str, derivatives):
            func = getattr(np, func_str)
            @implements(func_str, "ufunc")
            def implementation(*inputs, **kwargs):
                return apply_func_elementwise(
                    cls.wrap(func, derivatives), inputs, kwargs)
            return implementation

        for func_str, derivatives in ufunc_derivatives.items():
            implement_ufunc(func_str, derivatives)

    @classmethod
    def _add_numpy_comparative_ufuncs(cls):
        def implement_ufunc(func_str, dunder):
            func = getattr(np, func_str)
            @implements(func_str, "ufunc")
            def implementation(*inputs, **kwargs):
                if isinstance(inputs[0], cls):
                    func = getattr(inputs[0], dunder)
                    inputs_ = [inputs[1]]
                elif isinstance(inputs[1], cls):
                    func = getattr(inputs[1],  '__r' + dunder[2:])
                    inputs_ = [inputs[0]]
                
                return apply_func_elementwise(func, inputs_, kwargs, result_dtype=bool)
            
            return implementation

        ufunc_comparatives = {
            'equal': '__eq__',
            'not_equal': '__ne__',
            'less': '__lt__',
            'less_equal': '__le__',
            'greater': '__gt__',
            'greater_equal': '__ge__',
        }
        for func_str, dunder in ufunc_comparatives.items():
            implement_ufunc(func_str, dunder)
