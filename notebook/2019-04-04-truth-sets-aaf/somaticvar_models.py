
def get_yset(Y, xset, x_max, only_params):
    xsubset = [y for y in xset if y < x_max]
    n = len(xsubset)
    xsum = sum(xsubset)
    a = Y / (n - xsum / x_max)
    b = - a / x_max
    if only_params:
        return({'a': a, 'b': b})
    def M_function(x):
        y = max(a + b * x, 0)
        return(y)
    yset = [M_function(z) for z in xset]
    return(yset)
