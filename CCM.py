# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 08:50:37 2023

@author: Benjamin Vilmann
"""

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
from scipy.linalg import block_diag
from numpy.linalg import inv,eig
from numpy import sin, cos

import control
import matplotlib.cm as cm

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




#%%
class StateSpaceSystem:
    def __init__(self,file_path:str,ss,vals:dict={},verbose=False,comp_id=1):
        self.ssm = ssm = StateSpaceModelling(file_path,f'var',f'SS{ss}')
        self.L = df = pd.read_excel(file_path, sheet_name=f'L_map')

        self.name = f'{comp_id}'

        self.u = u = df[(df.group == 'comp') & (df.vector == 'input') & (df.comp == ss)].sort_values(['comp','idx'], inplace=False)
        self.y = y = df[(df.group == 'comp') & (df.vector == 'output') & (df.comp == ss)].sort_values(['comp','idx'], inplace=False)
        self.x = x = df[(df.group == 'comp') & (df.vector == 'state') & (df.comp == ss)].sort_values(['comp','idx'], inplace=False)
        self.a = a = df[(df.group == 'sys') & (df.vector == 'input') & (df.comp == ss)].sort_values(['comp','idx'], inplace=False)
        self.b = b = df[(df.group == 'sys') & (df.vector == 'output') & (df.comp == ss)].sort_values(['comp','idx'], inplace=False)

        # Change values based on input
        for k, v in vals.items():
            ssm.store.data[k] = v
        
        # Evaluate the numerical state space model
        ssm.eval_num_ss()

        # Identifying A, B, C, D matrices based on column and row names
        self.idx_row = idx_row = len([c for c in ssm.ss_eval.index if c[0].lower() == 'd' or c[0].lower() == 'x'])
        self.idx_col = idx_col = len([c for c in ssm.ss_eval.columns if c[0].lower() == 'x'])
        self.A = ssm.ss_num[:idx_row,:idx_col]
        self.B = ssm.ss_num[:idx_row,idx_col:]
        self.C = ssm.ss_num[idx_row:,:idx_col]
        self.D = ssm.ss_num[idx_row:,idx_col:]

        # Verbose print statement
        if verbose:
            print(f'A{ss}',self.A.shape,idx_row,idx_col,'\n',self.A)
            print(f'B{ss}',self.B.shape,idx_row,idx_col,'\n',self.B)
            print(f'C{ss}',self.C.shape,idx_row,idx_col,'\n',self.C)
            print(f'D{ss}',self.D.shape,idx_row,idx_col,'\n',self.D)
        
        return



#%%


class ComponentConnectionMethod:
    def __init__(self,file_path,params:dict={},comp_id:str='CCM',system_input:list=None, system_output:list=None):
        # -------- Store data --------
        # Count sheets starting with "SS"
        ss_sheets_count = len([1 for sheet in pd.ExcelFile(file_path).sheet_names if sheet.startswith('SS')])

        # Read file
        subsystems = [StateSpaceSystem(file_path, i, vals=params, comp_id=comp_id) for i in range(0,ss_sheets_count)]

        print(subsystems)

        # -------- Store data --------
        self.name = self.suffix = comp_id
        self.subsystems = subsystems
        I = lambda x: np.diag(np.ones(x))
        self.L = L_map = subsystems[0].L
        self.L_map = subsystems[0].L
        self.a = subsystems[0].a
        self.b = subsystems[0].b
        if system_input is None:
            self.system_input = list(L_map[(L_map['group'] == 'sys') & (L_map['vector'] == 'input')].name.values)
        else:
            self.system_input = system_input
        if system_output is None:
            self.system_output = list(L_map[(L_map['group'] == 'sys') & (L_map['vector'] == 'output')].name.values)
        else:
            self.system_output = system_output

        # -------- Composite system state model (CSSM) --------
        # Prepare dictionary for subsystem matrices
        for M in ['A','B','C','D']:
            setattr(self,M,{})
        
        # Get state space models of subsystems
        for i, ss in enumerate(subsystems):     
            self.A[i+1] = ss.A
            self.B[i+1] = ss.B
            self.C[i+1] = ss.C
            self.D[i+1] = ss.D

        self.u = u = pd.concat([ss.u for ss in subsystems], axis=0)
        self.x = x = pd.concat([ss.x for ss in subsystems], axis=0)
        self.y = y = pd.concat([ss.y for ss in subsystems], axis=0)

        # Create block diagonal matrix
        for M in ['A','B','C','D']:
            matrices = [matrix for key, matrix in getattr(self,M).items()]
            # print(M,'\n',getattr(self,M))
            getattr(self,M)[0] = block_diag(*matrices)
        
        # -------- L-map --------
        if system_input is not None and system_output is not None:
            L_map = L_map[((L_map['name'].isin(system_input)) & (L_map['group']=='sys') & (L_map['vector']=='input')) \
                          |((L_map['name'].isin(system_output)) & (L_map['group']=='sys') & (L_map['vector']=='output')) \
                          |(L_map['group']=='comp')]
        # Get L-map
        self.L = L = self.get_interconnection_matrices(L_map)
        for k,v in L.items():
            setattr(self,f'L{k}',v)
        
        # -------- Composite system state model (CSSM) --------
        # Ease notation
        A = self.A; B = self.B; C = self.C; D = self.D
        nx = len(A[0])
        ny = len(D[0] @ L[1])

        # Calcualte the composite system state model (CSSM)
        self.F =  F =  A[0] + B[0] @ L[1] @ inv(I(ny)-D[0] @ L[1]) @ C[0]
        self.G =  G =  B[0] @ L[1] @ inv(I(ny)-D[0] @ L[1]) @ D[0] @ L[2]+B[0] @ L[2]
        self.H = H = L[3] @ inv(I(ny)-D[0] @ L[1]) @ C[0]
        self.J = J = L[3] @ inv(I(ny)-D[0] @ L[1]) @ D[0] @ L[2]+L[4]

        # Store data in dictionary
        self.sys = {'F':F,
                    'G':G,
                    'H':H,
                    'J':J}

        return

    def __repr__(self):
        data = {
            'F': pd.DataFrame(self.F,index=self.x.name.values,columns=self.x.name.values),
            'G': pd.DataFrame(self.G,index=self.x.name.values,columns=self.a.name.values),
            'H': pd.DataFrame(self.H,index=self.b.name.values,columns=self.x.name.values),
            'J': pd.DataFrame(self.J,index=self.b.name.values,columns=self.a.name.values)
        }
        print_str = ''
        for k, v in data.items():
            print_str += f'# ============== [{k}] ============== #\n{v}\n\n'
        return print_str

    def show(self,save:bool = False,fontsize = 20,ccsm_tick_fontsize= 16,cssm_tick_fontsize=16,boldfacecolor='grey',text_alpha = 0.25):
        # ====================== CCSM ======================
        fig, ax = plt.subplots(1,1,dpi=150)
        M = np.vstack([np.hstack([self.A[0],self.B[0]]),
                      np.hstack([self.C[0],self.D[0]])])

        ax.imshow(np.where(M==0,np.nan,0.75), cmap='gray', vmin=0, vmax=2)

        xlabels = ['$' + n + '$' for n in self.x.latex_name] + ['$' + n + '$' for n in self.u.latex_name]
        ylabels = ['$' + n + '$' for n in self.x.latex_name] + ['$' + n + '$' for n in self.y.latex_name]

        ax.set_xticks([i for i in range(len(xlabels))])
        ax.set_yticks([i for i in range(len(ylabels))])
        ax.set_xticklabels(xlabels,fontsize=ccsm_tick_fontsize)
        ax.set_yticklabels(ylabels,fontsize=ccsm_tick_fontsize)

        # Minor ticks
        ax.set_xticks(np.arange(-.5, len(xlabels), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(ylabels), 1), minor=True)

        ax.axvline(self.x.shape[0]-.5,color='magenta',zorder=6)
        ax.axhline(self.x.shape[0]-.5,color='magenta',zorder=6)


        x_, y_, u_ = 0,0,0
        for ss in self.subsystems[:-1]:
            x_ += len(ss.x)
            u_ += len(ss.u)
            y_ += len(ss.y)

            # x lines
            ax.axhline(x_ - .5, color='cyan',ls = '--')
            ax.axvline(x_ - .5, color='cyan',ls = '--')

            # u line
            ax.axvline(self.x.shape[0] + u_ - .5, color='cyan',ls = '--')

            # y line
            ax.axhline(self.x.shape[0] + y_ - .5, color='cyan',ls = '--')

        ax.text(self.x.shape[0]/2-.5,self.x.shape[0]/2-.5, "A",ha='center',va='center', fontweight="bold",color=boldfacecolor,fontsize=fontsize,alpha=text_alpha,zorder=1)
        ax.text(self.x.shape[0]+self.u.shape[0]/2-.5,self.x.shape[0]/2-.5, "B",ha='center',va='center', fontweight="bold",color=boldfacecolor,fontsize=fontsize,alpha=text_alpha,zorder=1)
        ax.text(self.x.shape[0]/2-.5,self.x.shape[0]+self.y.shape[0]/2-.5, "C",ha='center',va='center', fontweight="bold",color=boldfacecolor,fontsize=fontsize,alpha=text_alpha,zorder=1)
        ax.text(self.x.shape[0]+self.u.shape[0]/2-.5,self.x.shape[0]+self.y.shape[0]/2-.5, "D",ha='center',va='center', fontweight="bold",color=boldfacecolor,fontsize=fontsize,alpha=text_alpha,zorder=1)
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()

        # Gridlines based on minor ticks
        ax.grid(which='minor', color='lightgrey', linestyle=':', linewidth=0.5)
        fig.tight_layout()

        if not save:
            plt.show()
        else:
            plt.savefig('CCSM_matrix.pdf',bbox_inches='tight')
        plt.close()

        # ====================== CSSM SYSTEM ======================
        fig, ax = plt.subplots(1,1,dpi=150)
        M = np.vstack([np.hstack([self.F,self.G]),
                      np.hstack([self.H,self.J])])

        ax.imshow(np.where(M==0,np.nan,0.75), cmap='gray', vmin=0, vmax=2)

        u = ['$' + n[1].latex_name + '$' for n in self.L_map.iterrows() if n[1].group =='sys' and n[1]['name'] in self.system_input]
        y = ['$' + n[1].latex_name + '$' for n in self.L_map.iterrows() if n[1].group =='sys' and n[1]['name'] in self.system_output]
        xlabels = ['$' + n + '$' for n in self.x.latex_name] + u
        ylabels = ['$' + n + '$' for n in self.x.latex_name] + y

        ax.set_xticks([i for i in range(len(xlabels))])
        ax.set_yticks([i for i in range(len(ylabels))])
        ax.set_xticklabels(xlabels,fontsize=cssm_tick_fontsize)
        ax.set_yticklabels(ylabels,fontsize=cssm_tick_fontsize)

        # Minor ticks
        ax.set_xticks(np.arange(-.5, len(xlabels), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(ylabels), 1), minor=True)

        ax.axvline(self.x.shape[0]-.5,color='magenta')
        ax.axhline(self.x.shape[0]-.5,color='magenta')

        ax.text(self.x.shape[0]/2-.5,self.x.shape[0]/2-.5, "F",ha='center',va='center', fontweight="bold",color=boldfacecolor,fontsize=fontsize,alpha=text_alpha,zorder=1)
        ax.text(self.x.shape[0]+len(u)/2-.5,self.x.shape[0]/2-.5, "G",ha='center',va='center', fontweight="bold",color=boldfacecolor,fontsize=fontsize,alpha=text_alpha,zorder=1)
        ax.text(self.x.shape[0]/2-.5,self.x.shape[0]+len(y)/2-.5, "H",ha='center',va='center', fontweight="bold",color=boldfacecolor,fontsize=fontsize,alpha=text_alpha,zorder=1)
        ax.text(self.x.shape[0]+len(u)/2-.5,self.x.shape[0]+len(y)/2-.5, "J",ha='center',va='center', fontweight="bold",color=boldfacecolor,fontsize=fontsize,alpha=text_alpha,zorder=1)
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()

        # Gridlines based on minor ticks
        ax.grid(which='minor', color='lightgrey', linestyle=':', linewidth=0.5)

        fig.tight_layout()

        if not save:
            plt.show()
        else:
            plt.savefig('CSSM_matrix.pdf',bbox_inches='tight')
        plt.close()

        return



    def get_interconnection_matrices(self,L_map):
        
        # Read the specific sheet 'L_map'
        df = L_map

        # df now contains the data from the 'L_map' sheet
        u_c = df[(df.group == 'comp') & (df.vector == 'input')].sort_values(['comp','idx'], inplace=False)
        y_c = df[(df.group == 'comp') & (df.vector == 'output')].sort_values(['comp','idx'], inplace=False)
        u_s = df[(df.group == 'sys') & (df.vector == 'input')].sort_values(['comp','idx'], inplace=False)
        y_s = df[(df.group == 'sys') & (df.vector == 'output')].sort_values(['comp','idx'], inplace=False)

        # Create a directed graph
        G = nx.DiGraph()

        # Add any additional vertices if needed
        G.add_nodes_from([ f"{n['group'][0]}{n['vector'][0]}_{n['name']}" for i, n in u_c.iterrows()])
        G.add_nodes_from([ f"{n['group'][0]}{n['vector'][0]}_{n['name']}" for i, n in y_c.iterrows()])
        G.add_nodes_from([ f"{n['group'][0]}{n['vector'][0]}_{n['name']}" for i, n in u_s.iterrows()])
        G.add_nodes_from([ f"{n['group'][0]}{n['vector'][0]}_{n['name']}" for i, n in y_s.iterrows()])

        L = {}
        edges = []
        edge_clrs = []
        # Set the seed for reproducibility
        for i, pair in enumerate([(y_c,u_c),(u_s,u_c),(y_c,y_s),(u_s,y_s)]):
            from_,to_ = pair   
            from_ = from_.reset_index()
            to_ = to_.reset_index()
            L[int(f'{i+1}')] = np.zeros((len(to_),len(from_)))
            # print(L[int(f'{i+1}')].shape)
            for j, node_j in to_.iterrows():
                # print(j,node_j['name'])
                for k, node_k in from_.iterrows():
                    # print(k,node_k['name'])
                    # print(j,k)
                    if node_j['name'] == node_k['name']:
                        L[int(f'{i+1}')][j,k] = 1

        return L
