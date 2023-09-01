
import re
import io

import pandas as pd

from .fraction_float import fraction_float

def csv_parser(csv_file):
    with open(csv_file, encoding='utf-8-sig') as f:
        lines = f.read()
    lines = re.sub('\,+\n', '\n', lines)

    well = dict()
    for section in lines.split('\n\n'):
        sec_lines = section.strip().split('\n')
       
        if len(sec_lines)>1:
            table_title = sec_lines[0]
            table_str   = '\n'.join(sec_lines[1:])

            if table_title.lower() == 'well_header':
                table = pd.read_csv(io.StringIO(table_str), names=['attribute', 'value'], header=None, index_col=0, usecols=[0,1]).squeeze()
                for attribute, value in table.items():
                    try:
                        table[attribute] = float(value)
                    except:
                        table[attribute] = value

            elif table_title.lower() == 'reservoir_pressure':
                table = pd.read_csv(io.StringIO(table_str), index_col=False).squeeze()
                # print(table)

            elif table_title.lower() == 'barrier_permability':
               table = pd.read_csv(io.StringIO(table_str), index_col=False).squeeze()

            elif table_title.lower() in ['main_barrier', 'co2_datum']:        
                table = pd.read_csv(io.StringIO(table_str), header= None, index_col=0).squeeze()

            else:
                table = pd.read_csv(io.StringIO(table_str))
                if table_title in ['drilling', 'casing_cement']:
                    table['diameter_in'] = table['diameter_in'].map(fraction_float)
            
            try:
                well[table_title] = table.to_dict()
                # well[table_title] = table
            except:
                well[table_title] = table

    return well

