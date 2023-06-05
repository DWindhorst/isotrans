""""

"""
# -*- coding: utf-8 -*-

import iso_flux_connections, iso_Storages


class Project(object):

    def __init__(self):
        """
        Constructor of iso_project.

        Projects are the owner of cells and cells are the owner of storages and storages are the owner of fluxconnections.

        """
        #self.__cells = []  # list holding all cells associated with this project
        #self.__connection_matrix = [1, 5]
        self.storage_instance = iso_Storages.Storage

    def new_storage(self):

        """
        Layer storing isotopes, states
         """

        return self.storage_instance

""""
    def get_flux_nodes(self):
        """
        Returns all flux nodes (storages, soil layers, atmospheres,...) associated with this project.
        """
        return LSI_storages.flux_node.Dict_of_all_flux_nodes.values()

    def get_iso_storages(self):
        """
        Returns all isotope storages (storages,pond, soil layers, ,...) associated with this project.
        """
        return LSI_storages.iso_storage.Dict_of_all_iso_storages.values()

    def get_flux_connections(self):
        """
        Returns all flux connections associated with this project.
        """
        return LSI_flux_connections.flux_connection.Dict_of_all_flux_connections.values()

    def get_cells(self):
        """
        Returns all cells associated with this project
        """
        return self.__cells

    def new_cell(self, area, atmosphere, x, y, z, ):
        """
        Creates a new cell and registers it in the project

        @param x:
        @type x:
        """
        i_point = LSI_storages.Point(x, y, z)
        i_cell = LSI_cell.iso_cell(location, atmosphere, area)
        self.__cells.append(i_cell)

    def run(self, Isotopologue, delta_time):
        """
        Runs SLI for one time step. All state variables (Temperatures) and flux parameters (q_l, q_v) need to be updated before running SLI.
        After one run cycle the states of the sytem will be updated

        B_t+1 * C_t+1 + B_in(x) * C_t+1(x) - B_out * C_t+1 = B * C_t --> x = source

        @param delta_time:
        @type delta_time:
        """

        iso_storages = self.get_iso_storages()
        # matrix_A = numpy.zeros(len(iso_storages),len(iso_storages)) #future Storage according to movements [row storages, columns concentrations]
        matrix_A = lil_matrix((len(iso_storages), len(iso_storages)))  # create a n x n matrix
        matrix_C = numpy.zeros(len(iso_storages))  # current storage

        for i_iso_storage in iso_storages:

            index_iso_storage = iso_storages.index(i_iso_storage)

            matrix_A[
                index_iso_storage, index_iso_storage] = i_iso_storage.liquid_H2O_volume + i_iso_storage.gas_H2O_volume

            matrix_C[index_iso_storage] = i_iso_storage.get_storage_i(
                Isotopologue)  # current storage of isotopes in i_iso_storage

            connections_to_iso_storages = i_iso_storage.connections_to_iso_storages

            for i_connections_to_iso_storage in connections_to_iso_storages:
                if isinstance(i_connections_to_iso_storage, LSI_flux_connections.boundary_connection):
                    # flow from boundary_connection just leave or enter the system at a known concentration and don't move within the system.
                    matrix_C[index_iso_storage] += i_connections_to_iso_storage.calc_flux_i(Isotopologue)
                else:
                    i_calc_flux_liquid = i_connections_to_iso_storage.calc_flux_liquid(Isotopologue)

                    index_righ_node = iso_storages.index(i_connections_to_iso_storage.right_node)
                    index_left_node = iso_storages.index(i_connections_to_iso_storage.left_node)

                    if i_calc_flux_liquid >= 0.0:  # flow from left node to right node with the concentration of the left node
                        matrix_A[index_righ_node, index_left_node] += i_calc_flux_liquid  # flux into right node
                        matrix_A[index_left_node, index_left_node] += i_calc_flux_liquid * -1  # flux from left node

                    else:  # flow from right node to left node with the concentration of the right node
                        matrix_A[index_left_node, index_right_node] += i_calc_flux_liquid * -1  # flux into left node
                        matrix_A[index_right_node, index_right_node] += i_calc_flux_liquid  # flux from right node

        """
        to solve sparse matrix with a small diagonal band use
        http://docs.scipy.org/doc/scipy/reference/sparse.html


        http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.solve.html

        """
        matrix_A = matrix_A.tocsr()
        c_i_layers_t1 = spsolve(A, B)

        # update c_i

        print
        c_i_layers_t1

    def print_hello(self):
        print('hello!')
        return

    def print_hello2(self):
        print('hello!')
        return

    """