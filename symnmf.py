import numpy as np
import sys 
import symnmf_ext
import math 

def main():
    np.random.seed(1234)

    # Read command line arguments
    k, goal, file_name = read_cmd_args() 
    
    # Load data from file and convert to matrix
    data = np.loadtxt(file_name, delimiter=",", dtype=np.float64)
    data = data.tolist()

    # Run the appropriate function based on the goal
    if goal == "sym":
        result = symnmf_ext.sym(data)  # Compute similarity matrix
    elif goal == "ddg":
        result = symnmf_ext.ddg(data)  # Compute diagonal degree matrix
    elif goal == "norm":
        result = symnmf_ext.norm(data)  # Compute normalized similarity matrix
    elif goal == "symnmf":
        norm = symnmf_ext.norm(data)
        initial_h = (initialize_matrix_h(k, norm)).tolist()
        result = symnmf_ext.symnmf(norm, initial_h, k)  # Run full SymNMF algorithm
    else:
        print("An Error Has Occurred")
        sys.exit(1)

    # Print the output matrix, formatted to 4 decimal places
    np.savetxt(sys.stdout, result, fmt="%.4f", delimiter=",")

def read_cmd_args():
    """
    Reads command line arguments and returns the values of k, goal, and file_name.
    
    Returns:
        tuple: A tuple containing:
            - k (int): The number of clusters to create.
            - goal (str): The goal of the program.
            - file_name (str): The name of the file to read data from.
    
    Raises:
        SystemExit: If the number of command-line arguments is incorrect.
    """
    if len(sys.argv) != 4:
        print("An Error Has Occurred")
        sys.exit(1)
    k = int(sys.argv[1])
    goal = sys.argv[2]
    file_name = sys.argv[3]
    return k, goal, file_name

def initialize_matrix_h(k, result):
    """
    Initialize a random matrix H with values between 0 and 2*sqrt(avg/k).
    
    Creates a random matrix with the same number of rows as the result matrix and k columns.
    The range of random values is determined based on the average value of the result matrix.
    
    Parameters:
        k (int): The number of clusters to create (number of columns in the output matrix).
        result (np.ndarray): The result matrix used to determine dimensions and value range.
    
    Returns:
        np.ndarray: A random matrix of shape (rows, k) with values uniformly distributed 
                   between 0 and 2*sqrt(avg/k), where avg is the mean of the result matrix.
    """
    # Calculate the average of the result (norm) matrix and set the range of the random matrix
    avg = np.mean(result)
    matrix_range = 2 * math.sqrt(avg / k)

    # Create a random matrix with the same number of rows as the result matrix
    rows = len(result)
    matrix_h = np.random.uniform(0, matrix_range, (rows, k))
    return matrix_h

if (__name__ == "__main__"):
    main()