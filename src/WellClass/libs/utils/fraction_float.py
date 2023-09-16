
# handle type hints problem for python version < 3.10
from typing import Union

def fraction_float(frac_str: str) -> Union[float, int]:
    ''' Evaluates numbers as e.g. 12 1/4 
    '''
    frac_str = frac_str.split()    # "12"  "1/4"
    integer = frac_str[0]          # "12"
    rational = eval(integer)       # 12
    if len(frac_str)>1:
        fraction = frac_str[1]     #1/4
        rational += eval(fraction) #12 + 0.25
    return rational


def float_to_fraction_inches(d_in: float) -> str:
        """Converts float number to integer fraction string
        """
        int_d = int(d_in)                                          #int(12.25) = 12
        remainder = d_in - int_d                                   #0.25
        if remainder > 0:
                fraction = remainder.as_integer_ratio()            #0.25 -> 1/4 ->  (1, 4)
                fraction_str = f'{fraction[0]}/{fraction[1]}'      #"1/4"
                d_fmt = f'{int_d} {fraction_str}"'                 #"12 1/4"  
                return d_fmt
        else:
                d_fmt = f'{int_d}"'                                #"12"
                return d_fmt
