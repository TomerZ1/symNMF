import sys
import numpy as np

def main():
    if (len (sys.argv) > 4 or len(sys.argv) <= 2): #added =
        print("An Error Has Occurred")
        sys.exit(1)
    if not is_valid_int(sys.argv[1]):
            print("Invalid number of clusters!") 
            sys.exit(1)
    if (len (sys.argv) == 4):
        if (not is_valid_int(sys.argv[2])):
            print("Invalid maximum iteration!")
            sys.exit(1)   
        iter = int(float(sys.argv[2])) 
    else:
        iter = 200    
    k = int(float(sys.argv[1]))
    file_path = sys.argv[-1]
    matrix = input_to_list(file_path)
    if (not (k> 1 and k < len(matrix))):
        print("Invalid number of clusters!")
        sys.exit(1)
    if (not (iter > 1 and iter < 1000)):
        print("Invalid maximum iteration!")
        sys.exit(1)
    i = 0
    flag = True
    centroids = matrix[:k]

    while ((i < iter) and (flag)):
        cluster_ass = cluster_assign(centroids, matrix)
        flag = update_centroids(centroids, cluster_ass, matrix)
        i += 1
    
    #for i in range(len(centroids)):
     #   for j in range(len(centroids[i])-1):
      #     print("%.4f" % centroids[i][j],end=',')
       # print("%.4f" % centroids[i][j+1])
    #np.savetxt(sys.stdout, centroids, fmt="%.4f", delimiter=", ")
    return centroids
            
def is_valid_int(value):
    for i in range(len(value)):
        if (not value[i].isdigit()) and (value[i] != "."):
            return False
        elif (value[i] == "."):
            if(i == (len(value) -1)):
                return False
            i+=1
            while (i < len(value)):
                if (not (value[i] == '0')):
                    return False
                i+=1

            break
    return True

def input_to_list(path):
    matrix = []
    with open(path, 'r') as file:
        for line in file:
            matrix.append([float(x) for x in line.strip().split(',')])
    return matrix

def distance(point1, point2):
    return sum([(x-y)**2 for x,y in zip(point1, point2)])**0.5

def cluster_assign(centroids, points):
    minDis = float('inf')
    cluster_assignments = [0 for i in range(len(points))]
    for j in range(len(points)):
        for k in range(len(centroids)):
            curr_dis = distance(points[j], centroids[k])
            if (curr_dis < minDis):
                minDis = curr_dis
                cluster_assignments[j] = k
        minDis = float('inf')
    return cluster_assignments

def update_centroids(centroids, cluster_assignments, points):
    flag = False
    for i in range(len(centroids)):
        cluster_points = [points[j] for j in range(len(points)) if cluster_assignments[j] == i]
        if (len(cluster_points) == 0):
            continue
        temp = [0 for k in range(len(centroids[0]))]
        for j in range(len(cluster_points)):
            for k in range(len(centroids[0])):
                temp[k] += cluster_points[j][k]
        for k in range(len(centroids[0])):
            temp[k] = temp[k]/len(cluster_points)
        if (distance(temp, centroids[i]) >= 0.001):
            flag = True
        centroids[i] = temp
    return flag



if (__name__ == "__main__"):
    main()