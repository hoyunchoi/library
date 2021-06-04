import math


def _sci(value, digit):
    if value == 0:
        return ""
    exponent = math.floor(math.log10(value))
    residual = round(math.pow(10.0, math.log10(value) - exponent), 10)
    if (residual == 1):
        return "10^{" + str(exponent) + "}"
    elif (digit):
        return "{:.{}f}".format(residual, digit) + "\\times 10^{" + str(exponent) + "}"
    else:
        return "{}\\times 10^".format(residual) + "{" + str(exponent) + "}"


def _float(value, digit):
    if (digit):
        return "{:.{}f}".format(value, digit)
    else:
        return str(value)


def latex_string(string):
    return "$" + string + "$"


def latex_bold_string(string):
    return "$\\boldsymbol{" + string + "}$"


def latex_sci(value, digit=None):
    return latex_string(_sci(value, digit))


def latex_bold_sci(value, digit=None):
    return latex_bold_string(_sci(value, digit))


def latex_float(value, digit=None):
    return latex_string(_float(value, digit))


def latex_bold_float(value, digit=None):
    return latex_bold_string(_float(value, digit))


def list_latex_string(string_list):
    result = []
    for string in string_list:
        result.append(latex_string(string))
    return result


def list_latex_bold_string(string_list):
    result = []
    for string in string_list:
        result.append(latex_bold_string(string))
    return result


def list_latex_sci(value_list, digit=None):
    result = []
    for value in value_list:
        result.append(latex_sci(value, digit))
    return result


def list_latex_bold_sci(value_list, digit=None):
    result = []
    for value in value_list:
        result.append(latex_bold_sci(value, digit))
    return result


def list_latex_float(value_list, digit=None):
    result = []
    for value in value_list:
        result.append(latex_float(value, digit))
    return result


def list_latex_bold_float(value_list, digit=None):
    result = []
    for value in value_list:
        result.append(latex_bold_float(value, digit))
    return result


if __name__ == "__main__":
    print("This is a module latex_string.py")
