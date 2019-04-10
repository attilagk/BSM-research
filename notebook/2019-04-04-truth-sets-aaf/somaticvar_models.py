
def get_yset(Y, xset, x_max, return_params):
    xsubset = [y for y in xset if y < x_max]
    n = len(xsubset)
    xsum = sum(xsubset)
    a = Y / (n - xsum / x_max)
    b = - a / x_max
    if return_params:
        return({'a': a, 'b': b})
    def M_function(x):
        y = a + b * x
        return(y)
    yset = [M_function(z) for z in xsubset]
    return(yset)
