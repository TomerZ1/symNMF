import sys
import math
import sklearn as sk
import sklearn.metrics as skm
import symnmf
import symnmf_ext
import kmeans
import numpy as np

def read_cmd():
    """
    Reads command line arguments and returns the values of k and file_name.
    
    Returns:
        tuple: A tuple containing:
            - k (int): The number of clusters to create.
            - file_name (str): The name of the file to read data from.
    
    Raises:
        SystemExit: If the number of command-line arguments is incorrect.
    """
    if len(sys.argv) != 3:
        print("An Error Has Occurred")
        sys.exit(1)
    k = int(sys.argv[1])
    file_name = sys.argv[2]
    return k, file_name

def euclidean_distance(x, y):
    """
    Computes the Euclidean distance between two points.
    
    Parameters:
        x (list): The first point.
        y (list): The second point.

    Returns:
        float: The Euclidean distance between the two points.
    """
    return np.linalg.norm(np.array(x) - np.array(y))

def process_nmf_result(H):
    """
    Assigns each row in the NMF result matrix to a cluster based on the highest value.
    
    Parameters:
        H (list): The NMF result matrix.

    Returns:
        list: A list of cluster assignments.
    """
    result = [np.argmax(row) for row in H]
    return result

def assign_kmeans_clusters(points, centroids):
    """
    Assigns each point to the nearest KMeans cluster based on Euclidean distance.

    Parameters:
        points (list): A list of points.
        centroids (list): A list of centroids.

    Returns:
        list: A list of cluster assignments
    """
    return [np.argmin([euclidean_distance(point, centroid) for centroid in centroids]) for point in points]

def main():
    # Read cmd arguments and create matrix of points
    k, filename = read_cmd()
    data = np.loadtxt(filename, delimiter=",", dtype=np.float64)
    data = data.tolist()

    # Create normalized matrix, initialize H, and run SymNMF
    norm_matrix = symnmf_ext.norm(data)
    initial_h = symnmf.initialize_matrix_h(k, norm_matrix).tolist()
    symnmf_clustering = process_nmf_result(symnmf_ext.symnmf(norm_matrix, initial_h, k))

    original_argv = sys.argv  # Backup original arguments
    sys.argv = ["kmeans.py", str(k), filename]  # Set new args for kmeans.py
    
    # Run KMeans and assign points to clusters
    kmeans_clustering = kmeans.main()
    kmeans_clustering = assign_kmeans_clusters(data, kmeans_clustering)
    
    sys.argv = original_argv # Restore original sys.argv

    # Compute silhouette scores
    kmeans_sil = skm.silhouette_score(data, kmeans_clustering)
    nmf_sil = skm.silhouette_score(data, symnmf_clustering)

    # Print silhouette scores
    print(f"nmf: {format(nmf_sil, '.4f')}")
    print(f"kmeans: {format(kmeans_sil, '.4f')}")

if __name__ == "__main__":
    main()