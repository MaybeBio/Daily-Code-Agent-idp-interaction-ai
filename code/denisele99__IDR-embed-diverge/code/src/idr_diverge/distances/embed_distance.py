from sklearn.metrics import pairwise_distances, pairwise
import pandas as pd
import numpy as np

#Based on implementation from:
#Littmann et al. (2021) "Embeddings from deep learning transfer GO annotations beyond homology"
#https://github.com/Rostlab/goPredSim/
#Modified by Denise Le
#modifications include addition of add_to_distance_matrix function


class EmbedDistanceMatrix:
    '''Modified from Littmann et al. 2021 goPredSim'''
    def __init__(self, embedding_db: pd.DataFrame):
        """
        Initialize the embedding distance matrix object with a lookup database.

        Args:
            embedding_db (pd.DataFrame): A DataFrame where the first column contains IDs 
                                         and the remaining columns contain embedding vectors.
        """
        # Validate input DataFrame
        if embedding_db.shape[1] < 2:
            raise ValueError("embedding_db must have at least one feature column.")
        if embedding_db.iloc[:, 0].duplicated().any():
            raise ValueError("Duplicate IDs found in the embedding database.") #Can comment out
            embedding_db=embedding_db.drop_duplicates(subset=embedding_db.columns[0])
            
        
        # Store embeddings in a dictionary for easy lookup
        self.embedding_db = dict(zip(embedding_db.iloc[:, 0], embedding_db.iloc[:, 1:].to_numpy()))
        
        # Prepare data for distance matrix computation
        self.ids, self.raw_data = zip(*self.embedding_db.items())
        self.raw_data = np.array(self.raw_data)  # Convert to NumPy array
        
        # Compute the initial distance matrix using all embeddings TODO is this necessary?Delete?
        self.distance_matrix, self.distance_ids = self.run_embedding_lookup_distance(
            queries=self.embedding_db, metric='cosine'
        )

    def run_embedding_lookup_distance(self, queries, metric='cosine'):
        """
        Calculate pairwise distances between query embeddings and the lookup database.

        Args:
            queries (dict or np.ndarray): A dictionary of query IDs to embeddings, 
                                          or a NumPy array of query embeddings.
            metric (str): The distance metric to use (e.g., 'cosine', 'euclidean').

        Returns:
            tuple:
                - distances (np.ndarray): A 2D array of distances between queries and the database.
                - query_ids (list): A list of query IDs.
        """
        # Validate distance metric
        if metric not in pairwise.PAIRWISE_DISTANCE_FUNCTIONS.keys():
            raise ValueError(f"{metric} is not a valid distance metric. Use valid metrics from sklearn.")
        
        # Handle queries as dictionary or NumPy array
        if isinstance(queries, dict):
            query_ids, raw_query_data = zip(*queries.items())
        else:
            raw_query_data = queries
            query_ids = range(len(queries))
        
        raw_query_data = np.array(raw_query_data).squeeze()

        # Ensure 2D array for single query
        if len(raw_query_data.shape) == 1:
            raw_query_data = raw_query_data.reshape(1, -1)

        # Compute pairwise distances
        distances = pairwise_distances(raw_query_data, self.raw_data, metric=metric)
        return distances, query_ids

    def add_to_distance_matrix(self, query_embeddings: pd.DataFrame):
        """
        Extend the existing distance matrix by adding query embeddings.

        Args:
            query_embeddings (pd.DataFrame): A DataFrame where the first column contains query IDs 
                                             and the remaining columns are embedding vectors.

        Returns:
            tuple:
                - final_dist_matrix (np.ndarray): Updated distance matrix.
                - updated_ids (list): List of all IDs (existing and new).
        """
        # Validate input DataFrame
        if query_embeddings.shape[1] < 2:
            raise ValueError("query_embeddings must have at least one feature column.")
        if query_embeddings.iloc[:, 0].duplicated().any():
            #raise ValueError("Duplicate IDs found in query_embeddings.") #Can comment out
            print("Duplicate IDs found in query_embeddings.")
            print(list(query_embeddings[query_embeddings.duplicated(subset=query_embeddings.columns[0])][0]))
            query_embeddings=query_embeddings.drop_duplicates(subset=query_embeddings.columns[0])

        # Extract query IDs and vectors
        query_ids = query_embeddings.iloc[:, 0]
        query_vectors = query_embeddings.iloc[:, 1:].to_numpy()
        
        # Number of query embeddings
        n_query = len(query_vectors)
        
        # Calculate distances between query embeddings and existing embeddings
        dist_to_base = pairwise_distances(query_vectors, np.concatenate((self.raw_data,query_vectors),axis=0), metric='cosine')
        
        # Update the distance matrix
        # 1. Add columns for distances from the base embeddings to the queries
        extended_columns = np.concatenate((self.distance_matrix, dist_to_base.T[:-n_query]), axis=1)
        
        # 2. Add rows for distances from the queries to all embeddings (base + query)
        final_dist_matrix = np.concatenate((extended_columns, dist_to_base), axis=0)
        
        # Update IDs
        updated_ids = list(self.ids) + list(query_ids)
        return final_dist_matrix, updated_ids
