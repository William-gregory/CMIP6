import numpy as np
from scipy import stats
import itertools
import operator

class Network:
    def __init__(self,data,dimX=0,dimY=0,dimT=0,nodes={},corrs=[],tau=0,gridcells=[],unavail=[],anomaly={},links={},strength={},strengthmap=[]):
        """
        The input 'data' are expected to be de-trended (zero-mean)
        and in the format x,y,t if an area grid, or lat,lon,t for
        a lat-lon grid.
        """
        self.data = data
        self.dimX,self.dimY,self.dimT = self.data.shape
        self.nodes = nodes
        self.corrs = corrs
        self.tau = tau
        self.gridcells = gridcells
        self.unavail = unavail
        self.anomaly = anomaly
        self.links = links
        self.strength = strength
        self.strengthmap = strengthmap
    
    def get_threshold(self, significance=0.01):
        """
        Compute pairwise correlations between all grid cells.
        The average of all correlations which are positive and
        below a specified significance level will determine the
        threshold which is used to cluster cells to form network
        nodes in the function get_nodes().
        """
        ID = np.where(np.abs(np.nanmax(self.data,2))>0)
        R = np.corrcoef(self.data[ID])
        np.fill_diagonal(R,np.nan)
        self.corrs = np.copy(R)
        self.gridcells = ID[0]*self.dimY + ID[1]
        
        df = self.dimT - 2
        R = R[R>=0]
        T = R*np.sqrt(df/(1 - R**2))
        P = stats.t.sf(T,df)
        R = R[P<significance]

        self.tau = np.mean(R)
    
    def get_nodes(self, latlon=False):
        """
        cluster grid cells together to from nodes of the
        complex network. Clustering is based on a greedy
        algorithm, and the threshold for clustering two 
        grid cells together is defined by self.tau
        """
        ids = np.where(np.isnan(self.data[:,:,:]))
        i_nan = ids[0][0] ; j_nan = ids[1][0]
        
        def area_neighbours(Area, i_nan, j_nan):
            rows = np.array(Area)[:,0]
            cols = np.array(Area)[:,1]
            rows_m = rows-1
            cols_m = cols-1
            rows_p = rows+1
            cols_p = cols+1
            
            p1 = np.array([rows_m,cols]).ravel().reshape(len(rows),2,order='F')
            p2 = np.array([rows_p,cols]).ravel().reshape(len(rows),2,order='F')
            p3 = np.array([rows,cols_m]).ravel().reshape(len(rows),2,order='F')
            p4 = np.array([rows,cols_p]).ravel().reshape(len(rows),2,order='F')
            cond1 = p1[:,0]<0
            cond2 = p2[:,0]>self.dimX-1
            cond3 = p3[:,1]<0
            cond4 = p4[:,1]>self.dimY-1
            if latlon:
                p3[:,1][cond3] = self.dimY-1
                p4[:,1][cond4] = 0
            else:
                p3[:,0][cond3] = i_nan
                p3[:,1][cond3] = j_nan
                p4[:,0][cond4] = i_nan
                p4[:,1][cond4] = j_nan
            p1[:,0][cond1] = i_nan
            p1[:,1][cond1] = j_nan
            p2[:,0][cond2] = i_nan
            p2[:,1][cond2] = j_nan
            p = np.concatenate((p1,p2,p3,p4)).tolist()
            return [i for i in p if i not in self.unavail]

        def area_max_correlation(Area, neighbours):
            Rmean = [] ; X = []
            for cell in neighbours:
                R = []
                new_cell = cell[0]*self.dimY + cell[1]
                if new_cell in self.gridcells:
                    X.append(cell)
                    IDm = np.where(self.gridcells==new_cell)
                    Rmean.append(np.nanmean(self.corrs[cells_in_k,IDm]))
            try:
                Rmax = np.nanmax(Rmean)
            except ValueError:
                Rmax = np.nan
            return np.array(X), Rmean, Rmax
        
        def diag_indices(a, k):
            rows, cols = np.diag_indices_from(a)
            if k < 0:
                return rows[-k:], cols[:k]
            elif k > 0:
                return rows[:-k], cols[k:]
            else:
                return rows, cols

        #S T E P   1   (C R E A T E   N O D E S)

        self.nodes = {}
        self.unavail = []
        if latlon:
            neighbour_corrs1 = self.corrs.diagonal(offset=1)
            neighbour_corrs2 = self.corrs.diagonal(offset=self.dimY-1)
            subset = np.arange(0,len(neighbour_corrs2),self.dimY)
            neighbour_corrs2 = neighbour_corrs2[subset]
            neighbour_corrs = np.concatenate((neighbour_corrs1,neighbour_corrs2))

            cellIDs1 = diag_indices(self.corrs,1)
            cellIDs2 = diag_indices(self.corrs,self.dimY-1)

            cellIDs = (np.concatenate((cellIDs1[0],cellIDs2[0][subset])),\
                       np.concatenate((cellIDs1[1],cellIDs2[1][subset])))
        else:
            neighbour_corrs = self.corrs.diagonal(offset=1)
            cellIDs = diag_indices(self.corrs,1)
            
        cellIDs = (self.gridcells[cellIDs[0]],self.gridcells[cellIDs[1]])
        k = 0
        neighbour_corrs,cellIDs1,cellIDs2 = list(zip(*sorted(zip(neighbour_corrs,cellIDs[0],cellIDs[1]),reverse=True)))
        cell_IDs = (cellIDs1,cellIDs2)
        np.random.seed(2)
        for it in range(len(neighbour_corrs)):
            cells_in_k = []
            i = cell_IDs[0][it]
            j = cell_IDs[1][it]
            r = neighbour_corrs[it]
            
            row_i = int(np.floor(i/self.dimY)) ; col_i = int(i % self.dimY)
            row_j = int(np.floor(j/self.dimY)) ; col_j = int(j % self.dimY)
            
            if ([row_i,col_i] not in self.unavail) & ([row_j,col_j] not in self.unavail):
                if r>self.tau:
                    self.nodes.setdefault(k, []).append([row_i,col_i])
                    self.nodes.setdefault(k, []).append([row_j,col_j])
                    self.unavail.append([row_i,col_i])
                    self.unavail.append([row_j,col_j])
                    cells_in_k.extend(np.where(self.gridcells==i)[0])
                    cells_in_k.extend(np.where(self.gridcells==j)[0])

                    while True: #expand
                        neighbours = area_neighbours(self.nodes[k], i_nan, j_nan)
                        X, Rmean, Rmax = area_max_correlation(Area=self.nodes[k], neighbours=neighbours)
                        if Rmax > self.tau:
                            m = X[Rmean==Rmax].tolist()
                            if len(m)>1:
                                m = m[np.random.randint(low=0,high=len(m))]
                            else:
                                m = m[0]
                            self.nodes.setdefault(k, []).append(m)
                            self.unavail.append(m)
                            cells_in_k.extend(np.where(self.gridcells==m[0]*self.dimY+m[1])[0])
                        else:
                            break
                    if len(self.nodes[k]) <= 2:
                        del self.nodes[k]
                    k += 1
                else:
                    break
        
        #S T E P   2   (M E R G E   N O D E S)
        
        self.unavail = []
        while True:
            Rs = {}
            unavail_neighbours = {}
            num_cells = dict([(area,len(self.nodes[area])) if self.nodes[area] not in self.unavail else (area,np.inf) for area in self.nodes.keys()])
            maxID = min(num_cells.items(), key=operator.itemgetter(1))[0]
            if num_cells[maxID] > 175: #arbitrary choice?
                break
            else:
                cells_in_k = [np.where(self.gridcells==cell[0]*self.dimY+cell[1])[0] for cell in self.nodes[maxID]]
                neighbours = area_neighbours(self.nodes[maxID], i_nan, j_nan)
                for cell in neighbours:
                    gcell = cell[0]*self.dimY + cell[1]
                    Rmean = []
                    cond1 = gcell in self.gridcells
                    cond2 = cell not in self.nodes[maxID]
                    cond3 = cell not in [k for k, g in itertools.groupby(sorted(itertools.chain(*unavail_neighbours.values())))]
                    cond4 = len([area for area, cells in self.nodes.items() if cell in cells]) > 0
                    if (cond1) & (cond2) & (cond3) & (cond4):
                        nID = [area for area, cells in self.nodes.items() if cell in cells][0]
                        unavail_neighbours[nID] = self.nodes[nID]
                        X, Rmean, Rmax = area_max_correlation(Area=self.nodes[nID]+self.nodes[maxID], neighbours=self.nodes[nID]+self.nodes[maxID])
                        if nID not in Rs: 
                            Rs[nID] = np.nanmean(Rmean)
                try:
                    Rs_maxID = max(Rs.items(), key=operator.itemgetter(1))[0]
                    if Rs[Rs_maxID] > self.tau:
                        for cell in self.nodes.pop(Rs_maxID, None):
                            self.nodes.setdefault(maxID, []).append([cell[0],cell[1]])
                    else:
                        self.unavail.append(self.nodes[maxID])
                except ValueError:
                    self.unavail.append(self.nodes[maxID])
        
    def get_links(self, area=None, lat=None):
        """
        compute the anomaly time series associated with
        every node of the network, and subsequently compute
        weighted links (based on covariance) between all of
        these nodes. The strength of each node (also known as
        the weighted degree), is defined as the sum of the
        absolute value of each nodes links. Here the network
        is fully connected, so every node connects to every other
        node
        """
        self.anomaly = {}
        self.links = {}
        self.strength = {}
        self.strengthmap = np.zeros((self.dimX,self.dimY))*np.nan
        if lat is not None:
            scale = np.sqrt(np.cos(np.radians(lat)))
        elif area is not None:
            scale = np.sqrt(area)
        else:
            scale = np.ones((self.dimX,self.dimY))
            
        for A in self.nodes:
            temp_array = np.zeros(self.data.shape)*np.nan
            for cell in self.nodes[A]:
                temp_array[cell[0],cell[1],:] = np.multiply(self.data[cell[0],cell[1],:],scale[cell[0],cell[1]])
            self.anomaly[A] = np.nansum(temp_array, axis=(0,1))
            
        for A in self.anomaly:
            sdA = np.std(self.anomaly[A])
            for A2 in self.anomaly:
                sdA2 = np.std(self.anomaly[A2])
                if A2 != A:
                    self.links.setdefault(A, []).append(stats.pearsonr(self.anomaly[A],self.anomaly[A2])[0]*(sdA*sdA2))
                elif A2 == A:
                    self.links.setdefault(A, []).append(0)
            
        for A in self.links:
            absolute_links = []  
            for link in self.links[A]:
                absolute_links.append(abs(link))
            self.strength[A] = np.nansum(absolute_links)
            for cell in self.nodes[A]:
                self.strengthmap[cell[0],cell[1]] = self.strength[A]