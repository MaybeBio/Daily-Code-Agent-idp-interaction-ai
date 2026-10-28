#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modules used to:

    1. Parse input parameters
    2. Parse param files
    3. Write out list of simulation parameters to text
"""

def write_input_params(file_output,input_params):
    """
    write_input_params writes all the input information from
    *input_params* into *file_output* in the same syntax
    """

    with open(file_output,'a') as f:
        for key in input_params.keys():
            f.write( ''.join(key)+','+str(input_params[key])+'\n');
    f.close()


def input_parse(filename,params_flag=False):
    """
    Parses input files (filename is path to input parameters or param_list file)

    params_flag toggles between input_params vs param_list
    """

    input_parameters  ={};
    with open(filename, 'r') as f:
        count = 0;

        for line in f:
            line=line.strip();
            if line:
                if line.find('#') == -1:
                    if not params_flag:
                        var_name,var_value = line.split(',');
                        try:
                            input_parameters[var_name]=float(var_value);
                        except:
                            input_parameters[var_name] = var_value;
                    else:
                        if count==0:
                            var_name = line.strip('\n');
                            input_parameters[var_name] = [];
                            count+=1;
                        else:
                            input_parameters[var_name].append(float(line.strip('\n')))
    return input_parameters;
