import numpy as np

def latexString(t_string):
    return "$" + t_string + "$"

def latexSci(t_value):
    if (t_value==0):
        return ""
    exponent = int(np.floor(np.log10(t_value)))
    residual = round(np.power(10.0, np.log10(t_value)-exponent), 10)
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

def list_latexString(list_string):
    result = np.array([])
    for string in list_string:
        result = np.append(result, latexString(string))
    return result

def latexBoldString(t_string):
    return "$\\boldsymbol{" + t_string + "}$"

def latexBoldSci(t_value):
    exponent = int(np.floor(np.log10(t_value)))
    residual = round(np.power(10.0, np.log10(t_value)-exponent), 10)
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
        result = np.append(result, latexBoldSci(t_value))
    return result

def list_latexBoldFloat(list_value, digit=None):
    result = np.array([])
    for t_value in list_value:
        result = np.append(result, latexBoldFloat(t_value, digit))
    return result

def list_latexBoldString(list_string):
    result = np.array([])
    for string in list_string:
        result = np.append(result, latexBoldString(string))
    return result

if __name__=="__main__":
    print("This is a module latexString.py")