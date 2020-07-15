#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 15:14:04 2018

@author: VS
"""    

class RamaGrid():
    def __init__(self,phis,psis,verbose=True):
        '''
        Args:
            phis: list of phi values
            psis: list of psi values
        '''
        self.verbose = verbose
        # Round values to 2 decimal places
        self.phi = []
        self.psi = []
        for phi,psi in zip(phis,psis):
            self.phi.append(round(phi,2))
            self.psi.append(round(psi,2))
    
    def make_grid(self,bin_size=1):
        '''
        Make a 2d-grid and populate it by phi and psi values
        Args:
            bin_size: `int` value (bin_size=1 makes grid of 360*360)
        '''
        
        bin_size = float(bin_size)
        self.bin_size = bin_size
        
        dim = int(360/(bin_size))
        
        ## initialize 2d grid
        self.grid = [[0]*dim for i in range(dim)]
        if self.verbose:
            print('-'*100)
            print('Initialised a 2d grid of dimension - {}*{}'.format(dim,dim))
        ##
        
        self.total_bins = dim*dim
        self.population = len(self.phi)
        
        if self.verbose:
            print('Placing {} (phi,psi) pairs in the grid'.format(self.population))
        # time.sleep(1)
        count = 0
        ##
        
        ## for loop - populate grid
        for phi,psi in zip(self.phi,self.psi):
            ## transform (-180,180) to (0,360)
            phi = phi+180
            psi = psi+180
            # print(phi,psi)
            ##
            ## place in appropriate grid
            if phi==360:
                i=360//bin_size - 1
            else:
                i = phi//(bin_size)
            if psi==360:
                j=360//bin_size - 1
            else:
                j = psi//(bin_size)
            i=int(i)
            j=int(j)
            self.grid[j][i]+=1
            # print(i,j)
            ##
            ## the progress bar again
            count+=1
            ##
        ## end of for loop
        
        ## count non_empty bins
        self.non_empty_bins = self.total_bins
        for i in range(dim):
            for j in range(dim):
                # if zero
                if not self.grid[i][j]:
                    self.non_empty_bins-=1
        ## end of for loop
        
        # time.sleep(1)
        
        if self.verbose:
            print('\n{} of {} bins occupied i.e., {}% of total grid'.format(self.non_empty_bins,self.total_bins,
                    100*self.non_empty_bins/float(self.total_bins)))
            print('-'*100)

        return self.grid
        
    def _is_grid_made(self):
        ## check if grid has been made
        try:
            grid = self.grid
        except NameError:
            print('Grid has not been made yet! Will make now!')
            while True:
                try:
                    bin_size = float(input('Enter bin_size:'))
                    break
                except ValueError:
                    print('Invalid value for bin_size!')
                    print('Try again!')
                    continue
            grid = self.make_grid(bin_size)
        return grid
        
    def plot_grid(self):
        import numpy as np
        import plotly
        import plotly.graph_objs as go
        
        grid = self._is_grid_made()
        grid = np.array(grid)
        axisrange = np.arange(-180,180,self.bin_size)
        trace = go.Heatmap(z=(grid),x=axisrange,y=axisrange)
        data = [trace]
        plotly.offline.plot(data)

    def write_grid_tocsv(self,out_name):
        f = open(out_name,'w')
        for i in range(len(self.grid)):
            for j in range(len(self.grid)):
                f.write('{},'.format(self.grid[i][j]))
            f.write('\n')
        f.close()
        print('Wrote {}'.format(out_name))
    
    def count(self,thresholds=[0,1,2]):
        import numpy as np
        thresholds.sort()
        grid = self._is_grid_made()
        grid = np.array(grid)
        
        flat = grid.flatten()
        flat.sort()
        
        results = {}
        i = 0
        for th in thresholds:
            results[th]=0
            for f in flat[i:]:
                i+=1
                if f>th:
                    results[th]+=len(flat[i:])
                    break
                
        return results
  

def main(args):
    '''
    Things this function will do:
        1. Load coordinates and compute dihedrals
        2. Make atom index dic
        3. Make RamaGrids
        4. Compute flex index
    '''

    # PART 1
    
    import sys
    try:
        import mdtraj as md
        import pandas as pd
        import numpy as np
    except ImportError:
        print('!'*100)
        print('Essential libraries: mdtraj or pandas or numpy is missing !!')
        print('Install these libraries and try again.')
        print('Bye!')
        print('!'*100)
        sys.exit(0)
        
    print('-'*100)
    print('Loading trajectory...')
    try:
        t = md.load(args['pdb'])
    except OSError:
        print('!'*100)
        print('Input pdb file: {0} not found!'.format(args['pdb']))
        print('Bye!')
        print('!'*100)
        
    top = t.topology
    
    print('Calculating dihedrals...')
    
    temp_phi = md.compute_phi(t)
    phi_ind, phi_vals = temp_phi
    phi_res_list = []
    
    temp_psi = md.compute_psi(t)
    psi_ind, psi_vals = temp_psi
    psi_res_list = []
    
    atom_index_dic = {}
    
    # PART 2
    
    for a in top.atoms:
        items = str(a).split('-')
        res = items[0][:3]
        resnum = int(items[0][3:])
        # atom = items[-1]
        atom_index = a.index
    
        atom_index_dic[atom_index] = (res,resnum)
    
    # PART 3
    
    for inds in phi_ind:
        temp = []
        for i in inds:
            try:
                temp.append(atom_index_dic[i])
            except KeyError:
                temp.append((None,None))  
                continue
        phi_res_list.append(temp[-1])
        
    for inds in psi_ind:
        temp = []
        for i in inds:
            try:
                temp.append(atom_index_dic[i])
            except KeyError:
                temp.append((None,None))  
                continue
        psi_res_list.append(temp[0])
    
    # PART 4
    
    frame = 0
    phi_res_dic = {}
    psi_res_dic = {}
    
    for phi_list, psi_list in zip(phi_vals,psi_vals):
    #     dihedrals_by_frame[frame] = []
        
        ind = 0
        for phi in phi_list:
            phi = np.rad2deg(phi)
            # print(phi)
            phi_res = phi_res_list[ind]
            if phi_res in phi_res_dic:
                phi_res_dic[phi_res].append(phi)
            else:
                phi_res_dic[phi_res] = [phi]
            ind+=1
            
        ind = 0
        for psi in psi_list:
            psi = np.rad2deg(psi)
            # print(psi)
            psi_res = psi_res_list[ind]
            if psi_res in psi_res_dic:
                psi_res_dic[psi_res].append(psi)
            else:
                psi_res_dic[psi_res] = [psi]
            ind+=1
    
        frame+=1
    
    # print(psi_res_dic.keys())
    
    residues = []
    for res in phi_res_dic:
        if not res in residues:
            residues.append(res)
    for res in psi_res_dic:
        if not res in residues:
            residues.append(res)
    
    phi_psi_dic = {}
    for residue in residues:
        try:
            phis = phi_res_dic[residue]
        except KeyError:
            phis = None
        try:
            psis = psi_res_dic[residue]
        except KeyError:
            psis = None
        phi_psi_dic[residue] = {}
        if phis == None:
            phis = [180]*len(psis)
        elif psis == None:
            psis = [180]*len(phis)
        phi_psi_dic[residue]['phi'] = phis
        phi_psi_dic[residue]['psi'] = psis
    
    thresholds = args['thresholds']
    
    dic = {'resnum':[]}
    for th in thresholds:
        dic.update({'f_index_th={}'.format(th):[],\
                    'pop_th={}'.format(th):[]})
        
    f_is = {}
    pops = {}
    
    i = 0
    for res in residues:
        resnum = res[-1]
        dic['resnum'].append(resnum)
        dihs = phi_psi_dic[res]
        rg = RamaGrid(dihs['phi'],dihs['psi'],verbose=args['verbose'])
        rg.make_grid(bin_size=int(args['bin_size']))
        temp = rg.count(thresholds=thresholds)
        temp2 = {}
        temp3 = {}
        for key,val in temp.items():
            temp2[key] = val/float(rg.total_bins)
            temp3[key] = val/float(rg.population)
        f_is[resnum] = temp2
        pops[resnum] = temp3
        i+=1
        # rg.plot_grid()
        # rg.write_grid_tocsv('grid_resnum_{0}.csv'.format(resnum))
        # if i==2:
        #     break
    
    # print(f_is)
    for th in thresholds:
        for res in residues:
            resnum = res[-1]
            dic['f_index_th={}'.format(th)].append(f_is[resnum][th])
            dic['pop_th={}'.format(th)].append(pops[resnum][th])
    
    df = pd.DataFrame.from_dict(dic)
    df.to_csv(args['out_name'])
    print('-'*100)
    print('Results of flexibility index written to {}'\
          .format(args['out_name']))
    print('-'*100)

    
if __name__=='__main__':
    import sys
        
    flags = input('Enter flags file:')
    # flags = 'flexindex_pdb_sample.flags'
    print('-'*100)
    try:
        f = open(flags)
        print('Reading flags..')
    except IOError:
        print('Could not open flags file!')
        print('Bye!')
        sys.exit(0)
            
    args = {}
    # Read flags file line by line
    for line in f:
        if line.startswith('#'):
            #ignore commented lines
            pass
        elif line.startswith('-'):
            #find the quotes
            quote_is = line.index("'")
            qst = line.strip()[0:quote_is].strip()
            ans = line.strip()[quote_is+1:].strip()[:-1]
            args[qst[1:]] = ans
        else:
            #ignore blank lines
            pass
        
    if 'pdb' in args:
        pdb = args['pdb']
    if 'res_range' in args:
        a,b = map(int,args['res_range'].split('-'))
        res_range = range(a,b+1)
        args['res_range'] = res_range
    if 'thresholds' in args:
        thresholds = list(map(int,args['thresholds'].split(',')))
        args['thresholds'] = thresholds
    if 'verbose' in args:
        if args['verbose'].lower()=='true':
            args['verbose'] = True
        elif args['verbose'].lower()=='false':
            args['verbose'] = False
    
    print('-'*100)
    print('Arguments given:')
    for arg,val in args.items():
        print('{} : {}'.format(arg,val))
    print('-'*100)
    
    print('Calculating flexibility index now.. This might take a while..')
    main(args)

