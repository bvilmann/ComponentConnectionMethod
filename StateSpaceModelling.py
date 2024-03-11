import pandas as pd
import numpy as np
from numpy import sin, cos

class StateSpaceVariableStore:
    def __init__(self):
        self.data = {}

        # Add numpy functions or any other required functions to the namespace
        self.data.update({k: getattr(np, k) for k in dir(np) if not k.startswith('_')})

    def set(self, variable, value):
        self.data[variable] = value

    def get_binding(self):
        return self.data

class StateSpaceModelling:
    def __init__(self,filepath,var_sheet,ss_sheet):

        # Load data
        vars = pd.read_excel(filepath,sheet_name=var_sheet,header=0)
        ss = pd.read_excel(filepath,sheet_name=ss_sheet,header=0,index_col=0)

        # Store variables
        self.store = store = StateSpaceVariableStore()
        for i, row in vars.iterrows():
            # store.set(row['group'], row['variable'], row['value'])
            store.set(row['variable'], row['value'])

        # Create symbolic state space representation
        # self.ss_sym = ss.replace('\.', '_', regex=True)
        self.ss_sym = ss

        return

    def evaluate_ssm_cell(self,cell):
        try:
            return eval(cell, self.store.get_binding()) if cell else None
        except Exception as e:
            if (isinstance(cell, float) or isinstance(cell, int)) and not np.isnan(cell):
                try:
                    return float(cell)
                except Exception as e_:
                    print(f'"{cell}" was not recognized as a number, i.e. evaluated as 0!')
                    return 0
            elif (isinstance(cell, float) or isinstance(cell, int)) and np.isnan(cell):
                return 0

            raise ValueError(f'"{cell}" ({type(cell)}) was not recognized as a variable, i.e. evaluated as 0!\nCheck if you miss to initialize and load the variables into the CCM call.')
            return 0

    def eval_num_ss(self,format = 'numpy'):
        # input validation
        assert format in ['numpy','dataframe','df','pandas'], "Expected 'format' to be in ['numpy','dataframe','df','pandas']"

        # Evaluate each cell of state space model
        self.ss_eval = self.ss_sym.applymap(self.evaluate_ssm_cell)

        # Convert to numpy
        self.ss_num = self.ss_eval.to_numpy()

        # Return depending on requested format
        if format == 'numpy':
            return self.ss_num
        else:
            return self.ss_eval

# =================== SCRIPT ===================
# file = r'statespacemodels\ssm_test.xlsx'

# ssm = StateSpaceModelling(file)

# print(ssm.form_ss())


# =================== SCRIPT ===================

