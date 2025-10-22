"""
Cumulative2D
=============
Class created mainly to analyze lipid membranes in different ways

Classes
-------

.. autoclass:: Cumulative2D
    :members:
    :undoc-members:
    :show-inheritance:


"""


import MDAnalysis as mda
import numpy as np
import pandas as pd
from twodanalysis import MembProp
from twodanalysis.analysis import OrderParameters


class Cumulative2D(MembProp):
    def __init__(self,
                 universe,
                lipid_list = None,
                verbose = False,
                edges = None,
                nbins = None,
                periodic = False,
                working_lip = None,
                connection_chains = None,
                forcefield = "charmm",
                ):
        """ Class to compute the Voronoi tessalation for a given universe. It uses the Voronoi class from scipy.spatial

        Parameters
        ----------
        universe : MDAnalysis.core.universe.Universe
            Universe to be used for the analysis
        lipid_list : list, optional
            list of lipids to be considered for Voronoi tesselation, by default None
        verbose : bool, optional
            Descriptive tool, by default False
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
            If None, it uses the default value of 180
        periodic : bool, optional
            If True, it uses the periodicity of the system to extend the data, by default False
        edges : list, optional
            A list with the limits of the grid [xmin,xmax,ymin,ymax], by default None
            If None, it uses the min and max values of the atoms in the first frame.
        connection_chains : dict, optional
            dictionary with the connection chains for the lipids, by default None
            It is a dictionary with the lipid name as key and a list of lists/tuples as value. The lists/tuples contains the
            atoms that connect the headgroup and the lipid tails. For more details refer to When problems arise section.
        working_lip : dict, optional
            dictionary with the working lipids, by default None
            It is a dictionary with the lipid name as key and a dictionary as value. The dictionary contains the atoms that
            are used as the head reference (Commonly P for phospholipids). The minimum you should pass to this dictionary is the key
            "head". However, you can also pass other keys such as "last_c" to define the last and first carbon of the tail or charge.
            For example, for a phospholipid you can pass the following dictionary:
            working_lip = {"DOPC": {"head": "P",
                                    "charge": 0,}}
            The only key needed is "head". Further entries are inferef with MembProp class.


        """

        super().__init__(universe,
                         lipid_list=lipid_list,
                         connection_chains=connection_chains,
                         working_lip=working_lip,
                         verbose=verbose, forcefield=forcefield)

        if edges is None:
            positions = self.all_head.positions[:,:2]
            print(self.all_head)
            v_min = np.min(positions)
            v_max = np.max(positions)
            self.edges = [v_min, v_max, v_min, v_max]
        else:
            self.edges = edges


        self.start = 0
        self.final = -1
        self.step = 1

        self.nbins = nbins if nbins is not None else 50
        self.periodic = periodic

        self.guess_chain_lenght()



    ############## Order parameters related code ####################
    def order_matrix(self,
                        lipid = "DOPC",
                        layer = "top",
                        nbins = None,
                        n_chain = None,
                        edges = None,
                        start = None,
                        final = None,
                        step = 1,
                        method = "numpy",
                        ):
        """Generates 2D SCD for specific lipids

        Parameters
        ----------
        lipid : str, optional
            Lipid to compute SCD, by default "DOPC"
        layer : str, optional
            Layer to project in the 2D plane, by default "top"
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
        n_chain : list(int), optional
            list with the number of carbons in each chain i.e., [16,18], by default None
        edges : list(int), optional
            Edges for the 2D grid in the shape [xmin,xmax,ymin,ymax], by default None
        start : int, optional
            Frame to start, by default None
        final : int, optional
            Final frame to consider, by default None
        step : int, optional
            Frames to skip, by defaulonnection_chains=connection_chains,t 1
        method : str, optional
            Method to do the averages, by default "numpy"

        Returns
        -------
        ndarray(n_grid,n_grid), list
            matrix containing the 2D SCD and list containing edges in
            the following disposition [v_min,v_max,v_min,_vmax]
            (Can be used to plot directly with extent)
        """


        nbins = self.nbins if nbins is None else nbins
        start = self.start if start is None else start
        final = self.final if final is None else final
        all_head = self.all_head

        n_chain = self.chain_info[lipid] if n_chain is None else n_chain
        chain_structure = self.extract_chain_info(lipid)
        #print(chain_structure)

        try:
            n_chain1 = n_chain[0]
            n_chain2 = n_chain[1]

        except:
            n_chain1 = n_chain
            n_chain2 = 0



        matrix = [] # this will store a matrix of the shape (2+n_chain,
        for ts in self.u.trajectory[start:final:step]:
            z = all_head.positions[:,2]
            z_mean = z.mean() # get middle of the membrane
            #Pick atoms in the layer
            if layer == "both":
                layer_at = self.memb.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}))")


            else:
                layer_at = self.memb.select_atoms(f"byres ((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {self.map_layers[layer]} {z_mean})")

            lipid_resids = layer_at.residues.resids

            only_p = layer_at.select_atoms(f"name {self.working_lip[lipid]['head']}")
            positions = only_p.positions[:,:2]
            if n_chain1 !=0:
                layer_at1 = layer_at
                if self.forcefield == "amber":
                    chain1_resids = lipid_resids - 1
                    layer_at1 = self.u.select_atoms(f"resid {' '.join(map(str, chain1_resids))}")
                angles_sn1 = OrderParameters.individual_order_sn1(layer_at1, n_chain1, atoms_inv=chain_structure[1][:n_chain1])
                angles_sn1 = angles_sn1.T
                positions = np.concatenate([positions, angles_sn1], axis = 1)
            if n_chain2 != 0:
                layer_at1 = layer_at
                if self.forcefield == "amber":
                    chain1_resids = lipid_resids + 1
                    layer_at1 = self.u.select_atoms(f"resid {' '.join(map(str, chain1_resids))}")

                angles_sn2 = OrderParameters.individual_order_sn2(layer_at1, n_chain2, atoms_inv=chain_structure[0][:n_chain2])
                angles_sn2 = angles_sn2.T
                positions = np.concatenate([positions, angles_sn2], axis = 1)
            if self.periodic:
                head_pos = positions[:,:3]
                carbons = [positions[:,i+3] for i in range(positions.shape[1] - 3)]
                dimensions = ts.dimensions[:2]

                head_pos, carbons = self.extend_data(head_pos, dimensions, self.periodicity, others = carbons)
                positions = np.concatenate([head_pos, np.array(carbons).T], axis = 1)

            matrix.append(positions) # Expect dim (n_lipids, 2+n_chain1+n_chain2)


        #matrix = np.array(matrix) # Expect dim (frames, n_lipids, 2+n_chain1+n_chain2)
        matrix = np.concatenate(matrix, axis = 0) # Expect dim (n_lipids*frames, 2+n_chain1+n_chain2)



        edges = self.edges if edges is None else edges
        print(n_chain)
        if method == "numpy":
            H, edges = self.numpyhistogram2D(matrix[:,:2],
                                             matrix[:,2:],
                                             n_chain,
                                             bins = nbins,
                                             edges = edges)
        else:
            H, edges = self.histogram2D(matrix[:,:2],
                                        matrix[:,2:],
                                        n_chain,
                                        bins = nbins,
                                        edges = edges)

        H = np.rot90(H)
        H[H==0] = np.nan

        return H, edges



    # Computes teh histogram of the average order parameters in each bin
    def histogram2D(self,
                    sample1,
                    weights,
                    n_chain,
                    bins = None,
                    edges = None):
        """ Computes the 2D histogram of 2D data with various values taking an average of them

        Parameters
        ----------
        sample1 : np.array(n,2)
            2D data information
        weights : np.array(n,m)
            m values can be attached to the data (usually lenght of the tails)
        n_chain : int or list
            Number of carbons in each chain, e.g, 16 or [16,16] for both chains
        bins : int, optional
            Number of bins to split the space, by default 10
        edges : list(float)
            Edges for the 2D grid in the shape [xmin,xmax,ymin,ymax]

        Returns
        -------
        np.array, list
            matrix containining the averaged 2D histogram, edges corresponding to te matrix
        """
        if edges is None:
            limits = [[self.edges[0],self.edges[1]], [self.edges[2], self.edges[2]]]
        else:
            limits = [[edges[0],edges[1]], [edges[2], edges[3]]]
        nbin = np.empty(2,np.intp)
        edges = 2*[None]
        for i in range(2):
            edges[i] = np.linspace(limits[i][0], limits[i][1], bins +1)
            nbin[i] = len(edges[i]) + 1
        Ncount = (tuple(np.searchsorted(edges[i], sample1[:,i], side = "right") for i in range(2)))
        for i in range(2):
            on_edge = (sample1[:,i] == edges[i][-1])
            Ncount[i][on_edge] -= 1
        xy = np.ravel_multi_index(Ncount, nbin)
        xy_test = xy.reshape(-1,1)
        xy_test = np.concatenate((xy_test, weights), axis = 1)
        hist = self.count_order(xy_test, nbin.prod(), n_chain)
        hist = hist.reshape(nbin)
        hist = hist.astype(float, casting = "safe")
        core = 2*(slice(1,-1),)
        hist = hist[core]

        return hist, limits

    # Computes the histogram of the average order parameters in each bin
    def numpyhistogram2D(self,
                         sample1,
                         weights,
                         n_chain,
                         bins = 10,
                         edges = None):
        """ Computes the 2D histogram of 2D data with various values taking an average of them

        Parameters
        ----------
        sample1 : np.array(n,2)
            2D data information
        weights : np.array(n,m)
            m values can be attached to the data (usually length of the tails)
        n_chain : int or list
            Number of carbons in each chain, e.g, 16 or [16,16] for both chains
        bins : int, optional
            Number of bins to split the space, by default 10
        edges : list(float)
            Edges for the grid [xmin,xmax,ymin,ymax]
        Returns
        -------
        np.array, np.array
            matrix containing the averaged 2D histogram, edges corresponding to te matrix
        """
        if edges is None:
            edges = [[self.edges[0], self.edges[1]],
                     [self.edges[2], self.edges[3]]]
        else:
            edges = [[edges[0], edges[1]],
                     [edges[2], edges[3]]]

        #print(sample1.shape, weights.shape)
        if isinstance(n_chain, int):
            n_chain = [n_chain]

        hist, xedges,yedges = np.histogram2d(sample1[:,0],
                                             sample1[:,1],
                                             bins = bins,
                                             range = edges)
        hist[hist == 0] = 1

        count = 0
        mat_chain = []

        n_feat = 0
        for chain in n_chain:
            n_feat += chain

        #print("here", n_feat, weights.shape[1]

        weights = 1.5*weights-0.5

        for chain in n_chain:
            if chain == 0:
                continue
            matrix = np.zeros(hist.shape)

            for i in range(chain):
                temp, xedges, yedges = np.histogram2d(sample1[:,0],
                                                      sample1[:,1],
                                                      weights = weights[:,count * n_chain[0]+ i],
                                                      bins = bins,
                                                      range = edges)
                matrix += np.abs(temp)
            count += 1

            matrix = matrix/(chain*hist)
            matrix[matrix == 0] = np.nan
            mat_chain.append(matrix)

        matrix = np.mean(np.array(mat_chain), axis = 0)



        return matrix, [xedges[0],xedges[-1], yedges[0], yedges[-1]]





    @staticmethod
    def count_order(data, min_lenght, n_chain):
        """ Function used to count and average the order parameter in each grid square

        Parameters
        ----------
        data : ndarray(n,2 + len(n_chain[0]) or 2 + len(n_chain[0]) +len(n_chain[1]))
            Array that contains the order parameters data to be averaged.
        min_lenght : int
            size of the data
        n_chain : int or list
            Sets the number of carbons in the fatty acida chain

        Returns
        -------
        _np.array
            Array containing the mean SCD for each grid square
        """
        columns = ["index"]
        carbons_sn2 = False
        try:
            carbons_sn1 = [f"sn1-{i+2}" for i in range(n_chain[0]) ]
            carbons_sn2 = [f"sn2-{i+2}" for i in range(n_chain[1])]
            columns = columns + carbons_sn1 + carbons_sn2
        except:
            carbons_sn1 = [f"sn1-{i+2}" for i in range(n_chain) ]
            columns = columns + carbons_sn1


        df = pd.DataFrame(data, columns = columns)
        result = []

        for i in range(min_lenght):
            temp = df[df["index"] == i]
            if len(temp) > 0 :
                sn1 = temp[carbons_sn1]
                sn1 = sn1.mean()
                sn1 = 1.5 * sn1 - 0.5
                sn1 = sn1.abs()
                sn1 = sn1.mean()
                final = sn1
                if carbons_sn2:
                    sn2 = temp[carbons_sn2]
                    sn2 = sn2.mean()
                    sn2 = 1.5 * sn2 - 0.5
                    sn2 = sn2.abs()
                    sn2 = sn2.mean()

                    final = (final + sn2)*0.5
                result.append([i,final])
            else:
                result.append([i,0])
        result = np.array(result)
        result = result[:,1]
        return result


    def all_lip_order(self, layer, nbins = None,
                        edges = None,
                        start = None,
                        final = None,
                        step = None,
                        chain = "both"):
        """all_lip_order Find the 2D order parameters for all lipids



        Parameters
        ----------
        layer : (str)
            Layer, can be top, bot, both
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
        edges : list(float)
            Edges for the grid in the shape [xmin,xmax,ymin,ymax]
        start : (int, optional), optional
            start frame, by default None
        final : (int, optional), optional
            final frame, by default None
        step : (int, optional), optional
            step, by default 1
        chain : str, optional
            If "sn1" computes order parameters for sn1 tail. If "sn2" computes order parameters for
            sn2 tail. If "both" computes order parameters for both tails as an average of both.

        Returns
        -------
        ndarray(n,n), ndarray(n+1)
            matrix containing the 2d order, edges of the matrix
        """
        if nbins is None:
            nbins = self.nbins

        lipid_list = list(self.lipid_list)
        if "CHL1" in lipid_list:
            lipid_list.remove("CHL1")
        if "CHL" in lipid_list:
            lipid_list.remove("CHL")

        lipids = self.chain_info

        matrices = []
        for key in lipid_list:
            n_chain = lipids[key].copy()
            print(n_chain)
            if chain == "sn2":
                n_chain[0] = 0
            elif chain == "sn1":
                n_chain[1] = 0

            H, edges = self.order_matrix(key, layer, nbins, n_chain,edges,
                        start = start,
                        final = final,
                        step = step)

            matrices.append(H)
        matrices = np.array(matrices)
        matrices = np.nanmean(matrices, axis = 0)


        return matrices, edges


    ############## End of order parameters related code ############################3



    def surface(self,
                start = None,
                final = None,
                step = None,
                lipid_list = "DSPC",
                layer = 'top',
                function = None,
                splay = False,):
        """Code to loop over the trajectory and print [x,y,z(referenced to zmean), charge] in a file.


        Parameters
        ----------
        start : int, optional
            Start Frame, by default None
        final : int, optional
            Final frame, by default None
        step : int, optional
            Frames to skip, by default None
        lipid_list : str or list, optional
             Lipid to work, by default "DSPC"
        layer : str, optional
            Layer to work, by default 'top'
        function : function, optional
            Function that returns two arrays, the first with atomgroup.residues.resids, and the second with any
            property value.
            Defaults to None
        splay : bool, optional
            Include or not splay angle computations, by default False

        Returns
        -------
        pd.DataFrame
            Dataframe containing x,y,z,id,splay,function for each lipid
        """



        start = self.start if start is None else start
        final = self.final if final is None else final
        step = self.step if step is None else step


        sign = self.map_layers[layer]


        ##### Select all the P atoms to find the middle of the membrane
        all_p = self.all_head
        if isinstance(lipid_list, list):
            names = f" {self.build_name(self.working_lip[lipid_list[0]]['head'])}"
            selection_string = f"(resname {lipid_list[0]} and {names}) "
            for lipid in lipid_list[1:]:
                selection_string += f" or (resname {lipid} and {self.build_name(self.working_lip[lipid]['head'])}) "

        else:
            lipid_list = [lipid_list]
            selection_string = f"(resname {lipid_list} and {self.build_name(self.working_lip[lipid]['head'])}) "

        #print(self.working_lip, selection_string)
        columns = ["x", "y", "z", "id"]




        if splay:
            self.guess_last_cs()
            carbons1 = [self.working_lip[lipid]["last_c"][0] for lipid in lipid_list]
            carbons2 = [self.working_lip[lipid]["last_c"][1] for lipid in lipid_list]
            heads = [self.working_lip[lipid]["head"] for lipid in lipid_list]
            heads = self.build_resname_atom(lipid_list, heads)
            carbons1 = self.build_resname_atom(lipid_list, carbons1)
            carbons2 = self.build_resname_atom(lipid_list, carbons2)
            columns.append("splay")

        if function:
            columns.append("property")

        if self.verbose:
            print("######### Running surface function ########## ")
            print(f"We will compute the surface for a membrane with these lipids {lipid_list}")
            print(f"Currently working on: {lipid_list}")
            print(f"Layer: {layer}")



        #### Loop over trajectory to find the lipids in the requested membrane
        pos_data = []
        for ts in self.u.trajectory[start:final:step]:
            positions = all_p.positions[:,2]
            mean_z = positions.mean()

            # Selects the lipid head and of the working lipid
            if layer == "both":
                final_selection_string = selection_string
            else:
                #print(selection_string)
                final_selection_string = f"({selection_string}) and prop z {sign} {str(mean_z)}"

            # Select the atoms in the head
            if splay:
                final_selection_byres = f"byres {final_selection_string}"
                lipid_ats = self.memb.select_atoms(final_selection_byres)
                head_p = lipid_ats.select_atoms(heads)
                c1 = lipid_ats.select_atoms(carbons1)
                c2 = lipid_ats.select_atoms(carbons2)
                v1 = c1.positions - head_p.positions
                v2 = c2.positions - head_p.positions


                # Compute the cos of splay angle, must bet lenmght nlipids
                costheta = (np.sum(v1 * v2, axis=1)/
                            (np.linalg.norm(v1, axis = 1)* np.linalg.norm(v2, axis = 1)))
                costheta = np.arccos(costheta)
                costheta = np.rad2deg(costheta[:,np.newaxis])
                atoms = lipid_ats.select_atoms(selection_string)
            else:
                atoms = self.memb.select_atoms(final_selection_string)

            ### Get resids
            atom_resid = atoms.residues.resids





            ### Get positions # Maybe have to check what happens with masses 0
            atom_pos = atoms.center_of_mass(compound="residues")



            atom_pos = np.concatenate((atom_pos, atom_resid[:,np.newaxis]), axis = 1)
            if splay:
                atom_pos = np.concatenate((atom_pos, costheta), axis = 1)
            atom_pos[:,2] = np.abs(atom_pos[:,2]-mean_z)



            if function is not None:
                final_selection_byres = f"byres {final_selection_string}"
                lipid_ats = self.memb.select_atoms(final_selection_byres)
                ids, values = function(lipid_ats)


                #print(ids.shape, values.shape, atom_resid.shape)
                mapped_array = np.full_like(atom_resid, np.nan, dtype=float)

                id_to_value = dict(zip(ids, values))

                for i, id_ in enumerate(atom_resid):
                    mapped_array[i] = id_to_value.get(id_, np.nan)
                mapped_array = mapped_array[:,np.newaxis]
                atom_pos = np.concatenate((atom_pos, mapped_array), axis=1)


            if self.periodic:
                dimensions = self.u.trajectory.ts.dimensions[:2]
                head_pos = atom_pos[:,:3]
                others = [atom_pos[:,i + 3] for i in range(len(columns) - 3)]
                head_pos, others = self.extend_data(head_pos, dimensions, self.periodicity, others = others)
                others = np.array(others).T
                atom_pos = np.concatenate([head_pos, others], axis = 1)



            pos_data.append(atom_pos)



        pos_data = np.concatenate(pos_data, axis = 0)
        df_data = pd.DataFrame(pos_data, columns = columns)
        df_data["id"] = df_data["id"].astype(int)


        return df_data



    @staticmethod
    def build_resname_atom(resnames, atomsnames):
        resnames = list(resnames)
        string = f"( (resname {resnames[0]}  and name {atomsnames[0]} ) "
        for i, resname in enumerate(resnames[1:]):
            string += f" or (resname {resname}  and name {atomsnames[ i + 1]}) "
        string +=  " ) "
        return string




    def height_matrix(self,
                      lipid_list,
                      layer,
                      edges = None,
                      start = None,
                      final = None,
                      step = None,
                      nbins = None):
        """ Code to divide the space in a 2D grid and compute the height referenced to zmean

        Parameters
        ----------
        lipid_list : list(str)
            Lipids to include in the height analysis
        layer : str
            Working layer for heigt. It can be bot/top
        edges : list
            Edges for the grid in the shape [xmin,xmax,ymin,ymax]
        start : int, optional
            Frame to start analysis, by default None
        final : int, optional
            Final frame for the analysis, by default None
        step : int, optional
            Steps to skip, by default None
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]


        Returns
        -------
        ndarray(nbins,nbins), edges
            Retun a matrix with the height information
        """
        nbins = self.nbins if nbins is None else nbins

        # Assign values if any is None
        if edges is None:
            edges = [[self.edges[0], self.edges[1]],
                     [self.edges[2], self.edges[3]]]
        else:
            edges = [[edges[0], edges[1]],
                     [edges[2], edges[3]]]


        print(f"Computing matrix for {layer} in frames {start}-{final}")

        #### Obtain data
        data = self.surface(lipid_list= lipid_list,
                            layer = layer,
                            start = start,
                            final = final,
                            step = step)

        ### Obtain weighted 2D histogram
        H_height, x_edges, y_edges = np.histogram2d(x = data["x"],
                                                    y = data["y"],
                                                    weights = data["z"],
                                                    bins = nbins,
                                                    range = edges)
        ### Obtain count 2D histogram
        H_count, x_edges, y_edges = np.histogram2d(x = data["x"],
                                                   y = data["y"],
                                                   bins = nbins,
                                                   range = edges)

        # Replace 0 counts with 1 to avoind division by zero
        H_count[H_count == 0] = 1.

        # Obtain the average for each grid square
        H_avg = H_height/H_count

        # Replace values with zero with nan
        H_avg[H_avg == 0] =  np.nan

        # Rotate to recover the original shape of the x,y plane
        H_avg = np.rot90(H_avg)

        return H_avg, [x_edges[0],x_edges[-1],y_edges[0], y_edges[-1]]


    def splay_matrix(self,
                     lipid_list,
                     layer,
                     edges = None,
                     start = None,
                     final = None,
                     step = None,
                     nbins = None):
        """ Code to divide the space in a 2D grid and compute the splay angle

        Parameters
        ----------
        lipid_list : list(str)
            Lipids to include in the analysis
        layer : str
            Working layer. Can be top/bot
        edges : list
            Edges for the grid in the shape [xmin,xmax,ymin,ymax]
        start : int, optional
            Frame to start analysis, by default None
        final : int, optional
            Final frame for the analysis, by default None
        step : int, optional
            Steps to skip, by default None
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]


        Returns
        -------
        ndarray(nbins,nbins), edges
            Retun a matrix with the height information
        """

        nbins = self.nbins if nbins is None else nbins

        if edges is None:
            edges = [[self.edges[0], self.edges[1]],
                     [self.edges[2], self.edges[3]]]
        else:
            edges = [[edges[0], edges[1]],
                     [edges[2], edges[3]]]




        print(f"Computing matrix for {layer} in frames {start}-{final}")

        data = self.surface(lipid_list=lipid_list,
                            layer = layer,
                            start = start,
                            final = final,
                            step = step,
                            splay=True)
        H_splay, x_edges, y_edges = np.histogram2d(x = data["x"],
                                                   y = data["y"],
                                                   weights = data["splay"],
                                                   bins = nbins,
                                                   range = edges)
        H_count, x_edges, y_edges = np.histogram2d(x = data["x"],
                                                   y = data["y"],
                                                   bins = nbins,
                                                   range = edges)
        H_count[H_count == 0] = 1.
        H_avg = H_splay/H_count
        H_avg[H_avg == 0] =  np.nan
        H_avg = np.rot90(H_avg)
        return H_avg, [x_edges[0],x_edges[-1],y_edges[0], y_edges[-1]]



    def thickness(self,
                  nbins = None,
                  edges = None,
                  lipid_list = None,
                  start = None,
                  final=None,
                  step = None):
        """Find the thickness mapped in a 2d grid

        Parameters
        ----------
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
        edges : list
            Edges for the grid in the shape [xmin,xmax,ymin,ymax]
        lipid_list : list(str)
            Lipids to be considered in the thickness computation
        start : int, optional
            Start frame, by default None
        final : int, optional
            Final frame, by default None
        step : int, optional
            Step frame, by default None

        Returns
        -------
        np.array(nbins,nbins), edges
            Matrix with the thickness, edeges for the matrix
        """

        nbins = self.nbins if nbins is None else nbins
        if lipid_list is None:
            lipid_list = list(self.lipid_list)
            if "CHL1" in lipid_list:
                lipid_list.remove("CHL1")

        matrix_up, edges = self.height_matrix(lipid_list,
                        "top",
                        edges = edges,
                        start = start,
                        final = final,
                        step = step,
                        nbins = nbins)
        matrix_bot, edges = self.height_matrix(lipid_list,
                        "bot",
                        edges = edges,
                        start = start,
                        final = final,
                        step = step,
                        nbins = nbins)

        mat_thickness = matrix_bot + matrix_up

        mat_thickness[mat_thickness == 0] = np.nan

        return mat_thickness, edges

    def project_property(self,
                        function,
                      layer = "top",
                      distribution = False,
                      edges = None,
                      start = None,
                      final = None,
                      step = None,
                      nbins = None,
                      ):
        """ Code to divide the space in a 2D grid and project the value of a property

        Parameters
        ----------
        function : function
            Function that takes mda.Atomgroup and returns two arrays: (1) array containing ids
            (2) array containing values for those ids
        layer : str
            Working layer can be bot/top, defaults to "top"
        distribution : bool
            If true compute the histogram instead of the averages. This is useful when someone wants
            to visualize the lipid lateral distribution over time.
        edges : list(float)
            Edges for the grid in the shape [xmin,xmax,ymin,ymax]
        start : int, optional
            Frame to start analysis, by default None
        final : int, optional
            Final frame for the analysis, by default None
        step : int, optional
            Steps to skip, by default None
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]


        Returns
        -------
        ndarray(nbins,nbins), edges
            Retun a matrix with the height information
        """

        print(f"Computing matrix for {layer} in frames {start}-{final}")

        data = self.surface(function = function, lipid_list=self.lipid_list,
                            layer = layer,
                            start = start,
                            final = final,
                            step = step)



        nbins = self.nbins if nbins is None else nbins
        edges = [self.edges[:2], self.edges[2:]] if edges is None else [edges[:2], edges[2:]]


        H_property, x_edges, y_edges = np.histogram2d(x = data["x"],
                                                    y = data["y"],
                                                    weights = data["property"],
                                                    bins = nbins,
                                                    range = edges)
        H_out = H_property
        if not distribution:
            H_count, x_edges, y_edges = np.histogram2d(x = data["x"],
                                                   y = data["y"],
                                                   bins = nbins,
                                                   range = edges)

            H_count[H_count == 0] = 1.

            H_out = H_property/H_count

            H_out[H_out == 0] =  np.nan


        H_out = np.rot90(H_out)

        return H_out, [x_edges[0],x_edges[-1],y_edges[0], y_edges[-1]]
