import sys
import operator
import numpy as np
import pandas as pd
import mykmeanssp as kmeans


DEFAULT_MAX_ITER = 300

def get_arguments():

    # function must get 4 or 5 arguments (max_iter is optional)
    if len(sys.argv) > 6 or len(sys.argv) < 5:
        raise Exception

    # number of arguments is valid
    else:

        try:

            # first argument is always k
            k = int(sys.argv[1])

        # error in converting k to an int
        except:
            raise Exception

        # max iter is not given
        if len(sys.argv) == 5:
            max_iter = DEFAULT_MAX_ITER  # default value

            try:
                epsilon = float(sys.argv[2])

            # error in converting epsilon to a float
            except:
                raise Exception

            input_file1 = sys.argv[3]
            input_file2 = sys.argv[4]

        # max iter is given
        else:

            try:
                max_iter = int(sys.argv[2])
                epsilon = float(sys.argv[3])

            # error in converting max iter to an int or epsilon to a float
            except:
                raise Exception

            input_file1 = sys.argv[4]
            input_file2 = sys.argv[5]

        if (k < 0 or max_iter < 0 or epsilon < 0):
            raise Exception

    return k, max_iter, epsilon, input_file1, input_file2

def calculate_distance(vector1, vector2):
    dist = 0
    for i in range(len(vector1)):
        sub = vector1[i] - vector2[i]
        dist += sub ** 2
    return dist

def update_distances(distances, data, chosen_indexes, num_of_data_points):

    for i in range(num_of_data_points):
        min_dist = float('inf')
        for j in chosen_indexes:
            temp_dist = calculate_distance(data[i].tolist(), data[j].tolist())
            if temp_dist < min_dist:
                min_dist = temp_dist
        distances[i] = min_dist

def update_probabilities(probabilities, distances):

    sum_of_min_distances = np.sum(distances)
    for i in range(len(distances)):
        probabilities[i] = distances[i] / sum_of_min_distances

def choose_centroids(data, k, num_of_data_points):

    chosen_indexes = list()
    distances = np.empty(num_of_data_points)
    probabilities = np.empty(num_of_data_points)

    np.random.seed(0)
    chosen_indexes.append(np.random.choice(num_of_data_points))
    update_distances(distances, data, chosen_indexes, num_of_data_points)
    update_probabilities(probabilities, distances)

    i = 1
    while i < k:
        chosen_indexes.append(np.random.choice(num_of_data_points, p=probabilities))
        update_distances(distances, data, chosen_indexes, num_of_data_points)
        update_probabilities(probabilities, distances)
        i += 1
    return chosen_indexes


if __name__ == "__main__":

    # step 1: get arguments from command line and validate them
    try:
        k, max_iter, epsilon, input_file1, input_file2 = get_arguments()

    # arguments are missing or invalid
    except:
        print("Invalid Input!\n")
        sys.exit()

    # step 2: create the intersection of the input files
    df1 = pd.DataFrame(np.loadtxt(input_file1, delimiter=','))
    df2 = pd.DataFrame(np.loadtxt(input_file2, delimiter=','))
    joined_df = pd.merge(df1, df2, how='inner', on=0)
    joined_df = joined_df.rename(columns={0: 'index'})
    joined_df = joined_df.sort_values(by=['index'])
    joined_df = joined_df.drop(columns=['index'])
    joined_df = joined_df.to_numpy()
    dimensions = np.shape(joined_df)
    num_of_data_points = dimensions[0]
    data_point_size = dimensions[1]

    if k > num_of_data_points:
        print("Invalid Input!\n")
        sys.exit()

    # step 3: kmeans++ algorithm
    chosen_centroids_indexes = choose_centroids(joined_df, k, num_of_data_points)
    initial_centroids = [joined_df[idx].tolist() for idx in chosen_centroids_indexes]

    # step 4: pass data from Python to C and run kmeans algorithm
    # arguments passed to C: data, centroids, k, max iter, num of data points, size of data point, epsilon.
    final_centroids = kmeans.fit(joined_df.tolist(), initial_centroids, k, max_iter, num_of_data_points, data_point_size, epsilon)

    # print results: first indexes, then centroids
    print(",".join(str(x) for x in chosen_centroids_indexes))

    result = np.ndarray(shape=(k, data_point_size))
    for i in range(k):
        for j in range(data_point_size):
            result[i][j] = final_centroids[(i * data_point_size) + j]

    res = ""
    for i in result:
        for j in i:
            res += str(np.round(j, 4)) + ','
        print(res[:-1])
        res = ""
