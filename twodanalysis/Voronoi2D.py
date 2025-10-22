"""
Voronoi2D
=============
Class created mainly to analyze lipid membranes in different ways

Classes
-------

.. autoclass:: Voronoi2D
    :members:
    :undoc-members:
    :show-inheritance:


"""


import MDAnalysis as mda
import numpy as np
import pandas as pd
from twodanalysis import MembProp
from scipy.spatial import Voronoi
from scipy.spatial import ConvexHull
from twodanalysis.analysis import OrderParameters

class Voronoi2D(MembProp):

    def __init__(self,
                 universe,
                lipid_list = None,
                verbose = False,
                nbins = None,
                edges = None,
                connection_chains = None,
                working_lip = None,
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
                         verbose=verbose,
                         connection_chains=connection_chains,
                         working_lip=working_lip,
                         forcefield=forcefield)

        self.forcefield = forcefield

        if edges is None:
            positions = self.all_head.positions[:,:2]
            self.v_min = np.min(positions)
            self.v_max = np.max(positions)
            self.edges = [self.v_min, self.v_max, self.v_min, self.v_max]
        else:
            self.edges = edges
            self.v_min = edges[0]
            self.v_max = edges[1]

        self.start = 0
        self.final = -1
        self.step = 1
        self.nbins = 180
        if nbins is not None:
            self.nbins = nbins


        self.guess_chain_lenght()




    ####### Code related to APL using voronoi tesselation and delunay triangulations
    # This part needs scypy to process data
    def voronoi_properties(self,
                    layer = 'top',
                        working_lip = None,
                        lipid_list = None,
                        splay = False,
                        function = None,
                        ):
        """Computes the APL for membranes with different lipids
        Parameters
        ----------
        layer : str, optional
            layer to compute the apl. It can be top/bot, by default 'top'
        working_lip : dict, optional
            dictionary mapping lipid and atoms to work for APL, by default None
        lipid_list : list, optional
            list of lipids to be considered for APL, by default None
        splay : bool,
            If True, it computes the splay angle for the lipids.
        function: function:
            Function that takes mda. Atomgroup and returns two arrays: (1) array with ids
            (2) array with value of properties
        Returns
        -------
        dict
            dictionary with vertices, areas, apl (if true), function values (if provided)
        """


        sign = self.map_layers[layer]

        if lipid_list is None:
            lipid_list = list(self.lipid_list)

        if working_lip is None:
            working_lip = self.working_lip
        #print(self.lipid_list)


        all_p = self.all_head
        positions = all_p.positions


        mean_z = positions[:,2].mean()

        selection_string = f"(((resname {lipid_list[0]} and name {working_lip[lipid_list[0]]['head']}) and prop z {sign} {mean_z}))"
        for lipid in lipid_list[1:]:
            selection_string += f" or (((resname {lipid} and name {working_lip[lipid]['head']}) and prop z {sign} {mean_z}))"


        heads = self.memb.select_atoms(selection_string)
        heads_pos = heads.positions[:,:2]
        height_pos =np.abs(heads.positions[:,2] - mean_z)
        resnames_pos = heads.resnames
        orig_len = len(heads_pos)

        ## Extent data
        dimensions = self.u.trajectory.ts.dimensions[:2]
        columns_others = ["resnames", "heights"]
        others = [resnames_pos, height_pos]


        if splay:
            carbons1 = [self.working_lip[lipid]["last_c"][0] for lipid in lipid_list]
            carbons2 = [self.working_lip[lipid]["last_c"][1] for lipid in lipid_list]
            heads = [self.working_lip[lipid]["head"] for lipid in lipid_list]
            heads = self.build_resname_atom(lipid_list, heads)
            carbons1 = self.build_resname_atom(lipid_list, carbons1)
            carbons2 = self.build_resname_atom(lipid_list, carbons2)
            selection_string_byres = f"byres {selection_string}"
            lipid_ats = self.memb.select_atoms(selection_string_byres)


            head_p = lipid_ats.select_atoms(heads)
            c1 = lipid_ats.select_atoms(carbons1)
            c2 = lipid_ats.select_atoms(carbons2)
            v1 = c1.positions - head_p.positions
            v2 = c2.positions - head_p.positions
            #print(heads,carbons2,carbons1, head_p.n_atoms, c1.n_atoms)
            costheta = np.sum(v1 * v2, axis=1)/(np.linalg.norm(v1, axis = 1)* np.linalg.norm(v2, axis = 1))# Compute the cos of splay angle, must bet lenmght nlipids
            costheta = np.arccos(costheta)
            costheta = np.rad2deg(costheta)
            others.append(costheta)
            columns_others.append("splay")

        if function is not None:
            selection_string_byres = f"byres {selection_string}"
            lipid_ats = self.memb.select_atoms(selection_string_byres)
            resids = lipid_ats.residues.resids

            ids, values = function(lipid_ats)

            second_dim = 1
            if isinstance(values, np.ndarray) and np.atleast_2d(values).shape[0] != 1:
                second_dim = values.shape[1]


            mapped_array = np.full((len(resids), second_dim), np.nan, dtype=float)

            id_to_value = dict(zip(ids, values))

            for i, id_ in enumerate(resids):
                mapped_array[i,:] = id_to_value.get(id_, np.nan)
            others.append(mapped_array)
            columns_others.append("function")


        heads_pos, others = self.extend_data(heads_pos, dimensions, self.periodicity, others = others)
        resnames_pos = others[0]
        height_pos = others[1]


        voronoi_dict = {"vertices":[],
                        "points":heads_pos,
                        "areas":[],
                        "orig_len":orig_len
                         }

        for column, other in zip(columns_others, others):
            voronoi_dict[column] = other


        voronoi = Voronoi(heads_pos)
        vertices = voronoi.vertices

        resnames = list(set(voronoi_dict["resnames"]))
        result_dict = {}
        for lipid in resnames:
            result_dict[lipid] = list()

        #for i, region in enumerate(voronoi.point_region[:orig_len]):
        update_points = list()
        for i, region in enumerate(voronoi.point_region):
            #print(i, region, len(voronoi.point_region), len(voronoi_dict["resnames"]))
            if -1 in voronoi.regions[region]:
                #print("here")
                voronoi_dict["areas"].append(np.nan)
                continue

            vertex = vertices[voronoi.regions[region]]
            hull = ConvexHull(vertex)
            area = hull.volume
            voronoi_dict["areas"].append(area)
            #update_points.append(voronoi_dict["points"][i])
            if i < orig_len:
                result_dict[voronoi_dict["resnames"][i]].append(area)

            voronoi_dict["vertices"].append(vertex)
        #voronoi_dict["points"] = np.array(update_points)
        #print(len(voronoi_dict["areas"]))

        for lipid in resnames:
            result_dict[lipid] = np.mean(np.array(result_dict[lipid]))

        voronoi_dict["apl"] = result_dict


        return voronoi_dict


    def map_voronoi(self,
                    voronoi_points,
                    voronoi_property,
                    nbins = None,
                    edges = None):
        """ Function to map voronoi diagram to a 2D plane

        Parameters
        ----------
        voronoi_points : ndarray(:,2)
            [x,y] positions of the points to be considered in the voronoi plot
        voronoi_property : ndarray(:)
            Areas corresponding to the points
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
            If None, it uses the default value of 180
        edges : list(float)
            A list with the lipids of the grid [xmin,xmax,ymin,ymax], defaults to None

        Returns
        -------
        ndarray, edges
            numpy array (nbins,nbins), adn edges of this array
        """

        edges = self.edges if edges is None else edges
        nbins = self.nbins if nbins is None else nbins
        xmin, xmax, ymin, ymax = edges

        if isinstance(nbins, list):
            nbins1  = nbins[0]
            nbins2 =  nbins[1]

        elif isinstance(nbins, int):
            nbins1  = nbins
            nbins2 =  nbins
        xcoords = np.linspace(xmin, xmax, nbins1)
        ycoords = np.linspace(ymin, ymax, nbins2)


        xx, yy = np.meshgrid(xcoords, ycoords)
        grid_points = np.vstack([xx.ravel(), yy.ravel()]).T
        points = voronoi_points


        distances = np.linalg.norm(grid_points[:,None, :]- points[None,:,:], axis = 2)


        closest_seed_indices = np.argmin(distances, axis=1).astype(int)


        voronoi_property = np.atleast_2d(voronoi_property)
        if voronoi_property.shape[0] == 1:
            voronoi_property = voronoi_property.T
            grid = voronoi_property[closest_seed_indices, :].reshape(nbins2, nbins1)
        else:
            grid = voronoi_property[closest_seed_indices, :].reshape(nbins2, nbins1,voronoi_property.shape[1])

        return grid, edges



    def voronoi_apl(self,layer = "top",
                    start = 0,
                    final = -1,
                    step = 1,
                    lipid_list = None,
                    nbins = None,
                    edges = None):
        """Function to compute and map the grid APL for several frames, map them to a 2D grid and average them

        Parameters
        ----------
        layer : str, optional
            working lipid layer, by default "top"
        start : int, optional
            Frame to start, by default None. If None uses all the trajectory
        final : int, optional
            final frame, by default None. If None uses all the trajectory
        step : int, optional
            Frames to skip, by default None. If None uses all the trajectory
        lipid_list : list, optional
            lipids involved in the computation, by default None
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
            If None, it uses the default value of 180
        edges : list, optional
            A list with the limits of the grid [xmin,xmax,ymin,ymax]

        Returns
        -------
        ndarray
            Array with the averaged 2D APL, edges
        """
        lipid_list = list(self.lipid_list) if lipid_list is None else lipid_list



        nbins = self.nbins if nbins is None else nbins

        edges = self.edges if edges is None else edges

        xmin, xmax, ymin, ymax = edges

        grid_size = (xmax-xmin)/nbins
        ng = int(nbins*0.05)
        self.guess_last_cs()
        no_present = [lipid for lipid in list(self.lipid_list) if lipid not in lipid_list]
        matrices = []
        for _ in self.u.trajectory[start:final:step]:

            voronoi_dict = self.voronoi_properties(layer = layer)

            areas = np.array(voronoi_dict["areas"])
            #print(voronoi_dict["resnames"])
            #print(type(areas),
            #       voronoi_dict["points"].shape,
            #       areas.shape,
            #       voronoi_dict["resnames"].shape,
            #       voronoi_dict["splay"].shape)
            for lipid in no_present:
                areas[voronoi_dict["resnames"] == lipid] = np.nan

            points = voronoi_dict["points"]


            selection_x = (points[:,0] > xmin - ng*grid_size) & (points[:,0] < xmax + ng*grid_size)
            selection_y = (points[:,1] > ymin - ng*grid_size) & (points[:,1] < ymax + ng*grid_size)
            selection = selection_x & selection_y

            points = points[selection]
            areas = areas[selection]


            matrix,_ = self.map_voronoi(points, areas, nbins, [xmin, xmax, ymin, ymax])
            matrices.append(matrix)
            #print(matrix.shape)




        final_mat = np.nanmean(np.array(matrices), axis = 0)
        final_mat = np.flipud(final_mat)

        return final_mat, edges

    def windows_apl(self,
                    layer = "top",
                    start = 0,
                    final = -1,
                    step = 1,
                    w_size = 10,
                    nbins = 180,):
        """Function to compute APL and map it to a 2D grid in windows of time defined by the user.

        Parameters
        ----------
        layer : str, optional
            Wroking layer, by default "top"
        start : int, optional
            Start Frame, by default 0
        final : int, optional
            Final frame, by default -1
        step : int, optional
            Frames to skip, by default 1
        w_size : int, optional
            windows size (number of frames), by default 10
        lipid_list : list, optional
            lipids to be included in the analysis, by default None
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
            If None, it uses the default value of 180
        edges : list, optional
            list with the edges of the grid [xmin,xmax,ymin,ymax], by default None

        Returns
        -------
        list(ndarrays(nbins,bins)), list
            list with the windows averages between the start time and the final time, and edges of these matrices
        """


        if final == -1:
            final = len(self.u.trajectory)
        n_windows = (final-start)//w_size
        matrices = []
        for i in range(n_windows):
            matrix = self.voronoi_apl(layer = layer,
                                    start = i*w_size + start,
                                    final =(i+1)*w_size + start ,
                                    step = step,
                                    lipid_list = None,
                                    nbins = nbins)
            matrices.append(matrix)
        return matrices


    @classmethod
    def create_sel_string(self,lipid_list, sign, mean_z):
        selection_string = f"(((resname {lipid_list[0]} and name {self.working_lip[lipid_list[0]]['head']}) and prop z {sign} {mean_z}))"
        for lipid in lipid_list[1:]:
            selection_string += f" or (((resname {lipid} and name {self.working_lip[lipid]['head']}) and prop z {sign} {mean_z}))"
        return selection_string




    def voronoi_thickness(self,
                          start = None,
                          final = None,
                          step = None,
                          lipid_list = None,
                          nbins = None,
                          edges = None):
        """Function to compute and map the grid APL for several frames, map them to a 2D grid and average them


        Parameters
        ----------
        start : int, optional
            Frame to start, by default None. If None uses all the trajectory
        final : int, optional
            final frame, by default None. If None uses all the trajectory
        step : int, optional
            Frames to skip, by default None. If None uses all the trajectory
        lipid_list : list, optional
            lipids involved in the computation, by default None
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
            If None, it uses the default value of 180
        edges : list, optional
            A list with the limits of the grid [xmin,xmax,ymin,ymax]

        Returns
        -------
        ndarray
            Array with the averaged 2D APL, edges
        """

        lipid_list = list(self.lipid_list) if lipid_list is None else lipid_list
        start = self.start if start is None else start
        final = self.final if final is None else final
        step = self.step if step is None else step
        nbins = self.nbins if nbins is None else nbins
        edges = self.edges if edges is None else edges




        no_present = [lipid for lipid in list(self.lipid_list) if lipid not in lipid_list]
        matrices = []
        for _ in self.u.trajectory[start:final:step]:





            voronoi_dict = self.voronoi_properties(layer = "top",
                                            lipid_list=self.lipid_list)
            heights = voronoi_dict["heights"]
            for lipid in no_present:
                heights[voronoi_dict["resnames"] == lipid] = np.nan




            matrix_top,_ = self.map_voronoi(voronoi_dict["points"],
                                         heights,
                                         nbins,
                                         edges,

                                         )

            voronoi_dict = self.voronoi_properties(layer = "bot",
                                            lipid_list=self.lipid_list)

            heights = voronoi_dict["heights"]
            for lipid in no_present:
                heights[voronoi_dict["resnames"] == lipid] = np.nan
            matrix_bot,_ = self.map_voronoi(voronoi_dict["points"],
                                         voronoi_dict["heights"],
                                         nbins,
                                         edges,
                                         )

            matrix_thickness = matrix_top + matrix_bot


            matrices.append(matrix_thickness)

        final_mat = np.nanmean(np.array(matrices), axis = 0)
        final_mat = np.flipud(final_mat)
        return final_mat, edges


    def voronoi_height(self,
                       layer = "top",
                       start = None,
                       final = None,
                       step = None,
                       lipid_list = None,
                       nbins = None,
                       edges = None):
        """Function to compute and map the grid height for several frames, map them to a 2D grid and average them

        Parameters
        ----------
        layer : str, optional
            working lipid layer, by default "top"
        start : int, optional
            Frame to start, by default None. If None uses all the trajectory
        final : int, optional
            final frame, by default None. If None uses all the trajectory
        step : int, optional
            Frames to skip, by default None. If None uses all the trajectory
        lipid_list : list, optional
            lipids involved in the computation, by default None
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
            If None, it uses the default value of 180
        edges : list, optional
            A list with the limits of the grid [xmin,xmax,ymin,ymax]

        Returns
        -------
        ndarray
            Array with the averaged 2D height, edges
        """
        lipid_list = list(self.lipid_list) if lipid_list is None else lipid_list
        start = self.start if start is None else start
        final = self.final if final is None else final
        step = self.step if step is None else step
        nbins = self.nbins if nbins is None else nbins
        edges = self.edges if edges is None else edges


        no_present = [lipid for lipid in list(self.lipid_list) if lipid not in lipid_list]
        matrices = []
        for _ in self.u.trajectory[start:final:step]:

            voronoi_dict = self.voronoi_properties(layer = layer)

            heights = voronoi_dict["heights"]
            for lipid in no_present:
                heights[voronoi_dict["resnames"] == lipid] = np.nan
            matrix_height,_ = self.map_voronoi(voronoi_dict["points"],
                                         heights,
                                         nbins,
                                         edges,
                                         )
            matrices.append(matrix_height)
        final_mat = np.nanmean(np.array(matrices), axis = 0)
        final_mat = np.flipud(final_mat)
        return final_mat, edges


    def voronoi_splay(self,
                      layer = "top",
                      start = None,
                      final = None,
                      step = None,
                      lipid_list = None,
                      nbins = None,
                      edges = None):
        """Function to compute and map the grid splay angle for several frames, map them to a 2D grid and average them

        Parameters
        ----------
        layer : str, optional
            working lipid layer, by default "top"
        start : int, optional
            Frame to start, by default None. If None uses all the trajectory
        final : int, optional
            final frame, by default None. If None uses all the trajectory
        step : int, optional
            Frames to skip, by default None. If None uses all the trajectory
        lipid_list : list, optional
            lipids involved in the computation, by default None
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
            If None, it uses the default value of 180
        edges : list, optional
            A list with the limits of the grid [xmin,xmax,ymin,ymax]

        Returns
        -------
        ndarray
            Array with the averaged 2D splay angle, edges
        """
        if lipid_list is None:
            lipid_list = list(self.lipid_list)

        nbins = self.nbins if nbins is None else nbins
        edges = self.edges if edges is None else edges


        self.guess_last_cs()

        no_present = [lipid for lipid in list(self.lipid_list) if lipid not in lipid_list]
        matrices = []



        for _ in self.u.trajectory[start:final:step]:

            voronoi_dict = self.voronoi_properties(layer = layer, splay=True)


            splay_vect = voronoi_dict["splay"]
            for lipid in no_present:
                splay_vect[voronoi_dict["resnames"] == lipid] = np.nan


            matrix_height,_ = self.map_voronoi(voronoi_dict["points"],
                                         splay_vect,
                                         nbins,
                                         edges,

                                         )




            matrices.append(matrix_height)

        final_mat = np.nanmean(np.array(matrices), axis = 0)
        final_mat = np.flipud(final_mat)
        return final_mat, edges


    def voronoi_order(self,
                     lipid,
                     layer,
                     n_chain = None,
                     start = None,
                     final = None,
                     step = None,
                     nbins = None,
                     edges = None,
                     ):
        """Function to compute the order parameter for a given lipid in a given layer. It uses the voronoi tessalation
        to compute the order parameter and then maps it to a 2D grid. The order parameter is computed using the angles
        between the two chains of the lipid. The order parameter is defined as S = 1/2 * (3 * cos^2(theta) - 1)
        where theta is the angle between the two chains of the lipid.
        The order parameter is computed for each chain and then averaged to obtain the final order parameter.

        Parameters
        ----------
        lipid : str
            lipid to compute the order parameter
        layer : str
            layer to compute the order parameter
        n_chain : list(int), optional
            length of the chain involved in the computation, e.g., [16,18], by default None
        start : int, optional
            Frame to start, by default None. If None uses all the trajectory
        final : int, optional
            final frame, by default None. If None uses all the trajectory
        step : int, optional
            Frames to skip, by default None. If None uses all the trajectory
        lipid_list : list, optional
            lipids involved in the computation, by default None
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
            If None, it uses the default value of 180
        edges : list, optional
            A list with the limits of the grid [xmin,xmax,ymin,ymax], by default None


        Returns
        -------
        ndarray
            Array with the 2D order parameter, edges
        """
        start = self.start if start is None else start
        final = self.final if final is None else final
        step = self.step if step is None else step
        nbins = self.nbins if nbins is None else nbins
        edges = self.edges if edges is None else edges

        n_chain1 = self.chain_info[lipid][0] if n_chain is None else n_chain[0]
        n_chain2 = self.chain_info[lipid][1] if n_chain is None else n_chain[1]

        chain_structure = self.extract_chain_info(lipid)
        print("Chain strcucture: ", chain_structure)
        def get_ids(lipids):
            layer_lip = lipids.select_atoms(f"resname {lipid}")
            lip_resids = layer_lip.residues.resids

            angles = []
            if n_chain1 !=0:
                layer_lip1 = layer_lip
                if self.forcefield == "amber":
                    layer_lip1 = self.u.select_atoms(f"resid {' '.join(map(str, lip_resids - 1))}")
                angles_sn1 = OrderParameters.individual_order_sn1(layer_lip1, n_chain1, atoms_inv=chain_structure[1])
                angles.append(angles_sn1.T)
            if n_chain2 !=0:
                layer_lip1 = layer_lip
                #print("Here::::::", layer_lip1)
                if self.forcefield == "amber":
                    layer_lip1 = self.u.select_atoms(f"resid {' '.join(map(str, lip_resids + 1))}")
                #print("Here::::::", layer_lip1)
                angles_sn2 = OrderParameters.individual_order_sn2(layer_lip1, n_chain2, atoms_inv=chain_structure[0])
                angles.append(angles_sn2.T)
            #print(angles)
            angles = np.concatenate(angles, axis=1)

            layer_ids = layer_lip.residues.resids




            return [layer_ids,angles]

        matrices = []
        for _ in self.u.trajectory[start:final:step]:
            voronoi_dict = self.voronoi_properties(layer = layer, function = get_ids)

            order_matrix, _ = self.map_voronoi(voronoi_points=voronoi_dict["points"],
                                            voronoi_property=voronoi_dict["function"],
                                            nbins=nbins,
                                            edges=edges,
                                            )


            matrices.append(order_matrix)
        order = np.nanmean(np.array(matrices), axis=0)
        chains_order = []
        old_chain = 0
        for chain in n_chain:
            print(chain, old_chain, lipid, n_chain1, n_chain2)
            if chain != 0:
                chains_order.append(np.nanmean(order[:,:,old_chain:old_chain + chain], axis=2))
            old_chain += chain
        final_order = np.nanmean(np.array(chains_order), axis = 0)
        final_order = np.flipud(np.abs(1.5*final_order-0.5))
        return final_order, edges

    def voronoi_all_lip_order(self,
                     lipid_list,
                     layer,
                     chain = "both",
                     start = None,
                     final = None,
                     step = None,
                     nbins = None,
                     edges = None,
                     ):
        """Function to compute the order parameter for all specified lipids (average them) in a given layer.
        It uses the voronoi tessalation to compute the order parameter and then maps it to a 2D grid. The order
        parameter is computed using the angles between the two chains of the lipid. The order parameter is
        defined as S = 1/2 * (3 * cos^2(theta) - 1) where theta is the angle between the two chains of the lipid.

        Parameters
        ----------
        lipid_list : list, optional
            lipids involved in the computation, by default None
        layer : str
            layer to compute the order parameter
        chain : str
            chain to compute the order parameter, by default "both". It can be sn1/sn2/both
        start : int, optional
            Frame to start, by default None. If None uses all the trajectory
        final : int, optional
            final frame, by default None. If None uses all the trajectory
        step : int, optional
            Frames to skip, by default None. If None uses all the trajectory
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
            If None, it uses the default value of 180
        edges : list, optional
            A list with the limits of the grid [xmin,xmax,ymin,ymax], by default None
        Returns
        -------
        ndarray
            Array with the 2D order parameter, edges
        """
        order_matrices = []


        for lipid in lipid_list:

            n_chain = self.chain_info[lipid].copy()



            if chain == "sn2":
                n_chain[0] = 0
            elif chain == "sn1":
                n_chain[1] = 0

            order, edges = self.voronoi_order(lipid, layer=layer,n_chain=n_chain, start = start, final=final, step=step,
                                              nbins=nbins,
                                              edges=edges)
            order_matrices.append(order)
        order_matrices = np.nanmean(np.array(order_matrices), axis = 0)
        return order_matrices, edges


    def project_property(self, function,
                      layer = "top",
                      start = None,
                      final = None,
                      step = None,
                      nbins = None,
                      edges = None):
        """Function to compute and map the grid APL for several frames, map them to a 2D grid and average them

        Parameters
        ----------
        function : function
            Function that takes mda. Atomgroup and returns two arrays: (1) array with ids
            (2) array with value of properties
        layer : str, optional
            working lipid layer, by default "top"
        start : int, optional
            Frame to start, by default None. If None uses all the trajectory
        final : int, optional
            final frame, by default None. If None uses all the trajectory
        step : int, optional
            Frames to skip, by default None. If None uses all the trajectory
        lipid_list : list, optional
            lipids involved in the computation, by default None
        nbins : int or list, optional
            number of bins for the grid, by default None
            If int, it is the same for both dimensions. If list, it is a list with the number of bins for each dimension
            e.g., [nbins1, nbins2]
            If None, it uses the default value of 180
        edges : list, optional
            A list with the limits of the grid [xmin,xmax,ymin,ymax]

        Returns
        -------
        ndarray
            Array with the averaged 2D properties (coming from fuction), edges
        """


        start = self.start if start is None else start
        final = self.final if final is None else final
        step = self.step if step is None else step
        nbins = self.nbins if nbins is None else nbins
        edges = self.edges if edges is None else edges
        matrices = []
        for _ in self.u.trajectory[start:final:step]:
            voronoi_dict = self.voronoi_properties(layer = layer, function = function)
            property_vect = voronoi_dict["function"]
            matrix,_ = self.map_voronoi(voronoi_dict["points"],
                                         property_vect,
                                         nbins,
                                         edges,
                                         )
            matrices.append(matrix)

        final_mat = np.nanmean(np.array(matrices), axis = 0)
        final_mat = np.flipud(final_mat)
        return final_mat, edges



    @staticmethod
    def build_resname_atom(resnames, atomsnames):
        resnames = list(resnames)
        string = f"( (resname {resnames[0]}  and name {atomsnames[0]} ) "

        for i, resname in enumerate(resnames[1:]):
            string += f" or (resname {resname}  and name {atomsnames[ i + 1]}) "

        string +=  " ) "
        return string









