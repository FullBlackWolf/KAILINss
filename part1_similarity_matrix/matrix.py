import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import time

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as ssp

from cospar.tmap import _tmap_core as tmap_core
from cospar.tmap import _utils as tmap_util

from cospar import help_functions as hf
from cospar import logging as logg
from cospar import settings
from cospar import tool as tl
from types import MappingProxyType
from typing import (
    Union,
    Optional,
    Any,
    Mapping,
    Callable,
    NamedTuple,
    Generator,
    Tuple,
    Literal,
)
import warnings

import numpy as np
import scipy
from anndata import AnnData
from scipy.sparse import issparse, coo_matrix, csr_matrix
from sklearn.utils import check_random_state

def Hematopoiesis_distance_matrix(
    adata,
    round_of_smooth=20,
    neighbor_N=20,
    beta=0.1,
    truncation_threshold=0.001,
    save_subset=True,
    use_existing_KNN_graph=False,
    compute_new_Smatrix=True,
    use_full_Smatrix = True,
):
    """
    Generate similarity matrix (Smatrix) through graph diffusion

    It generates the similarity matrix via iterative graph diffusion.
    Similarity matrix from each round of diffusion will be saved, after truncation
    to promote sparsity and save space. If save_subset is activated, only save
    Smatrix for smooth rounds at the multiples of 5 (like 5,10,15,...). If a Smatrix is pre-computed,
    it will be loaded directly if compute_new_Smatrix=Flase.

    Parameters
    ----------
    adata: :class:`~anndata.AnnData` object
    file_name: str
        Filename to load pre-computed similarity matrix or save the newly
        computed similarity matrix.
    round_of_smooth: `int`, optional (default: 10)
        The rounds of graph diffusion.
    neighbor_N: `int`, optional (default: 20)
        Neighber number for constructing the KNN graph, using the UMAP method.
    beta: `float`, option (default: 0.1)
        Probability to stay at the origin in a unit diffusion step, in the range [0,1]
    truncation_threshold: `float`, optional (default: 0.001)
        At each iteration, truncate the similarity matrix using
        truncation_threshold. This promotes the sparsity of the matrix,
        thus the speed of computation. We set the truncation threshold to be small,
        to guarantee accracy.
    save_subset: `bool`, optional (default: True)
        If true, save only Smatrix at smooth round [5,10,15,...];
        Otherwise, save Smatrix at each round.
    use_existing_KNN_graph: `bool`, optional (default: False)
        If true and adata.obsp['connectivities'], use the existing knn graph to build
        the similarity matrix, regardless of neighbor_N.
    compute_new_Smatrix: `bool`, optional (default: False)
        If true, compute a new Smatrix, even if there is pre-computed Smatrix with the
        same parameterization.

    Returns
    -------
        similarity_matrix: `sp.spmatrix`
    """




            
    if (1):  # compute now
        from cospar import logging as logg
        logg.hint(f"Compute similarity matrix: computing new; beta={beta}")

        # add a step to compute PCA in case this is not computed

        if (not use_existing_KNN_graph) or ("connectivities" not in adata.obsp.keys()):
            # here, we assume that adata already has pre-computed PCA
            n_neighbors=neighbor_N
            sc.pp.neighbors(adata, n_neighbors=neighbor_N)###################################################################
        
        else:
            logg.hint(
                "Use existing KNN graph at adata.obsp['connectivities'] for generating the smooth matrix"
            )
        adjacency_matrix = adata.obsp["connectivities"]

        ############## The new method
        adjacency_matrix = (adjacency_matrix + adjacency_matrix.T) / 2
        ##############

        adjacency_matrix = hf.sparse_rowwise_multiply(
            adjacency_matrix, 1 / adjacency_matrix.sum(1).A.squeeze()
        )
        nrow = adata.shape[0]
        similarity_matrix = ssp.lil_matrix((nrow, nrow))
        similarity_matrix.setdiag(np.ones(nrow))
        transpose_A = adjacency_matrix.T

        if round_of_smooth == 0:
            SM = 0
            similarity_matrix = ssp.csr_matrix(similarity_matrix)
            ##ssp.save_npz(file_name + f"_SM{SM}.npz", similarity_matrix)

        for iRound in range(round_of_smooth):
            SM = iRound + 1

            logg.info("Smooth round:", SM)
            t = time.time()
            similarity_matrix = (
                beta * similarity_matrix + (1 - beta) * transpose_A * similarity_matrix
            )
            # similarity_matrix =beta*similarity_matrix+(1-beta)*similarity_matrix*adjacency_matrix
            # similarity_matrix_array.append(similarity_matrix)

            logg.hint("Time elapsed:", time.time() - t)

            t = time.time()
            sparsity_frac = (similarity_matrix > 0).sum() / (
                similarity_matrix.shape[0] * similarity_matrix.shape[1]
            )
            if sparsity_frac >= 0.1:
                # similarity_matrix_truncate=similarity_matrix
                # similarity_matrix_truncate_array.append(similarity_matrix_truncate)

                logg.hint(f"Orignal sparsity={sparsity_frac}, Thresholding")
                similarity_matrix = hf.matrix_row_or_column_thresholding(
                    similarity_matrix, truncation_threshold
                )
                sparsity_frac_2 = (similarity_matrix > 0).sum() / (
                    similarity_matrix.shape[0] * similarity_matrix.shape[1]
                )
                # similarity_matrix_truncate_array.append(similarity_matrix_truncate)

                logg.hint(f"Final sparsity={sparsity_frac_2}")

                logg.info(
                    f"similarity matrix truncated (Smooth round={SM}): ",
                    time.time() - t,
                )

            # logg.info("Save the matrix")
            # file_name=f'data/20200221_truncated_similarity_matrix_SM{round_of_smooth}_kNN{neighbor_N}_Truncate{str(truncation_threshold)[2:]}.npz'
            similarity_matrix = ssp.csr_matrix(similarity_matrix)

            ############## The new method
            # similarity_matrix=similarity_matrix.T.copy()
            ##############

            if save_subset:
                if SM % 5 == 0:  # save when SM=5,10,15,20,...

                    logg.hint("Save the matrix at every 5 rounds")
                    ##ssp.save_npz(file_name + f"_SM{SM}.npz", similarity_matrix)
            else:  # save all

                logg.hint("Save the matrix at every round")
                ##ssp.save_npz(file_name + f"_SM{SM}.npz", similarity_matrix)
    distance_matrix = similarity_matrix
    return distance_matrix



def Reprogramming_distance_matrix(
    adata,
    round_of_smooth=15,
    neighbor_N=20,
    beta=0.1,
    truncation_threshold=0.001,
    save_subset=True,
    use_existing_KNN_graph=False,
    compute_new_Smatrix=True,
    use_full_Smatrix = True,
):
    """
    Generate similarity matrix (Smatrix) through graph diffusion

    It generates the similarity matrix via iterative graph diffusion.
    Similarity matrix from each round of diffusion will be saved, after truncation
    to promote sparsity and save space. If save_subset is activated, only save
    Smatrix for smooth rounds at the multiples of 5 (like 5,10,15,...). If a Smatrix is pre-computed,
    it will be loaded directly if compute_new_Smatrix=Flase.

    Parameters
    ----------
    adata: :class:`~anndata.AnnData` object
    file_name: str
        Filename to load pre-computed similarity matrix or save the newly
        computed similarity matrix.
    round_of_smooth: `int`, optional (default: 10)
        The rounds of graph diffusion.
    neighbor_N: `int`, optional (default: 20)
        Neighber number for constructing the KNN graph, using the UMAP method.
    beta: `float`, option (default: 0.1)
        Probability to stay at the origin in a unit diffusion step, in the range [0,1]
    truncation_threshold: `float`, optional (default: 0.001)
        At each iteration, truncate the similarity matrix using
        truncation_threshold. This promotes the sparsity of the matrix,
        thus the speed of computation. We set the truncation threshold to be small,
        to guarantee accracy.
    save_subset: `bool`, optional (default: True)
        If true, save only Smatrix at smooth round [5,10,15,...];
        Otherwise, save Smatrix at each round.
    use_existing_KNN_graph: `bool`, optional (default: False)
        If true and adata.obsp['connectivities'], use the existing knn graph to build
        the similarity matrix, regardless of neighbor_N.
    compute_new_Smatrix: `bool`, optional (default: False)
        If true, compute a new Smatrix, even if there is pre-computed Smatrix with the
        same parameterization.

    Returns
    -------
        similarity_matrix: `sp.spmatrix`
    """




            
    if (1):  # compute now
        from cospar import logging as logg
        logg.hint(f"Compute similarity matrix: computing new; beta={beta}")

        # add a step to compute PCA in case this is not computed

        if (not use_existing_KNN_graph) or ("connectivities" not in adata.obsp.keys()):
            # here, we assume that adata already has pre-computed PCA
            n_neighbors=neighbor_N
            sc.pp.neighbors(adata, n_neighbors=neighbor_N)###################################################################
        
        else:
            logg.hint(
                "Use existing KNN graph at adata.obsp['connectivities'] for generating the smooth matrix"
            )
        adjacency_matrix = adata.obsp["connectivities"]

        ############## The new method
        adjacency_matrix = (adjacency_matrix + adjacency_matrix.T) / 2
        ##############

        adjacency_matrix = hf.sparse_rowwise_multiply(
            adjacency_matrix, 1 / adjacency_matrix.sum(1).A.squeeze()
        )
        nrow = adata.shape[0]
        similarity_matrix = ssp.lil_matrix((nrow, nrow))
        similarity_matrix.setdiag(np.ones(nrow))
        transpose_A = adjacency_matrix.T

        if round_of_smooth == 0:
            SM = 0
            similarity_matrix = ssp.csr_matrix(similarity_matrix)
            ##ssp.save_npz(file_name + f"_SM{SM}.npz", similarity_matrix)

        for iRound in range(round_of_smooth):
            SM = iRound + 1

            logg.info("Smooth round:", SM)
            t = time.time()
            similarity_matrix = (
                beta * similarity_matrix + (1 - beta) * transpose_A * similarity_matrix
            )
            # similarity_matrix =beta*similarity_matrix+(1-beta)*similarity_matrix*adjacency_matrix
            # similarity_matrix_array.append(similarity_matrix)

            logg.hint("Time elapsed:", time.time() - t)

            t = time.time()
            sparsity_frac = (similarity_matrix > 0).sum() / (
                similarity_matrix.shape[0] * similarity_matrix.shape[1]
            )
            if sparsity_frac >= 0.1:
                # similarity_matrix_truncate=similarity_matrix
                # similarity_matrix_truncate_array.append(similarity_matrix_truncate)

                logg.hint(f"Orignal sparsity={sparsity_frac}, Thresholding")
                similarity_matrix = hf.matrix_row_or_column_thresholding(
                    similarity_matrix, truncation_threshold
                )
                sparsity_frac_2 = (similarity_matrix > 0).sum() / (
                    similarity_matrix.shape[0] * similarity_matrix.shape[1]
                )
                # similarity_matrix_truncate_array.append(similarity_matrix_truncate)

                logg.hint(f"Final sparsity={sparsity_frac_2}")

                logg.info(
                    f"similarity matrix truncated (Smooth round={SM}): ",
                    time.time() - t,
                )

            # logg.info("Save the matrix")
            # file_name=f'data/20200221_truncated_similarity_matrix_SM{round_of_smooth}_kNN{neighbor_N}_Truncate{str(truncation_threshold)[2:]}.npz'
            similarity_matrix = ssp.csr_matrix(similarity_matrix)

            ############## The new method
            # similarity_matrix=similarity_matrix.T.copy()
            ##############

            if save_subset:
                if SM % 5 == 0:  # save when SM=5,10,15,20,...

                    logg.hint("Save the matrix at every 5 rounds")
                    ##ssp.save_npz(file_name + f"_SM{SM}.npz", similarity_matrix)
            else:  # save all

                logg.hint("Save the matrix at every round")
                ##ssp.save_npz(file_name + f"_SM{SM}.npz", similarity_matrix)
    distance_matrix = similarity_matrix
    return distance_matrix


