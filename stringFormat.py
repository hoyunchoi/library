import numpy as np

def latexString(t_string):
    return "$" + t_string + "$"

def latexSci(t_value):
    exponent = int(np.log10(t_value))
    residual = int(np.power(10, np.log10(t_value)-exponent))
    if (residual == 1):
        return "$10^{" + str(exponent) + "}$"
    else:
        return "$" + str(residual) +  "\\times 10^{" + str(exponent) + "}$"

def latexFloat(t_value, digit=None):
    if (digit):
        string = "{:.{}f}".format(t_value, digit)
    else:
        string = str(t_value)
    return "$" + string + "$"

def list_latexSci(list_value):
    result = np.array([])
    for t_value in list_value:
        result = np.append(result, latexSci(t_value))
    return result

def list_latexFloat(list_value, digit=None):
    result = np.array([])
    for t_value in list_value:
        result = np.append(result, latexFloat(t_value, digit))
    return result

def latexBoldString(t_string):
    return "$\\boldsymbol{" + t_string + "}$"

def latexBoldSci(t_value):
    exponent = int(np.log10(t_value))
    residual = int(np.power(10, np.log10(t_value)-exponent))
    if (residual == 1):
        return "$\\boldsymbol{10^{" + str(exponent) + "}}$"
    else:
        return "$\\boldsymbol{" + str(residual) +  "\\times 10^{" + str(exponent) + "}}$"

def latexBoldFloat(t_value, digit=None):
    if (digit):
        string = "{:.{}f}".format(t_value, digit)
    else:
        string = str(t_value)
    return "$\\boldsymbol{" + string + "}$"

def list_latexBoldSci(list_value):
    result = np.array([])
    for t_value in list_value:
        result = np.append(result, latexSci(t_value))
    return result

def list_latexBoldFloat(list_value, digit=None):
    result = np.array([])
    for t_value in list_value:
        result = np.append(result, latexFloat(t_value, digit))
    return result

if __name__=="__main__":
    print("This is a module latexString.py")