import sys
import numpy as np
import mySpkmeansModule as mspkm

from importlib import reload


def get_arguments():

    # args should be: file name, k (optional), goal, input file
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Invalid Input!")
        sys.exit()

    # num of args is valid, validate args values
    else:

        try:

            if len(sys.argv) == 4:

                # check that k is an int
                k = int(sys.argv[1])

                # check k is not a negative number
                if k < 0:
                    raise Exception

                # check that goal is valid
                if sys.argv[2] not in ["spk", "wam", "ddg", "lnorm", "jacobi"]:
                    raise Exception

                # check that the file name ends with .txt or .csv
                if not (sys.argv[3].endswith(".txt") or sys.argv[3].endswith(".csv")):
                    raise Exception

                if sys.argv[2] == "spk" and k == 1:
                    print("An Error Has Occurred", end="")
                    sys.exit()

                return k, sys.argv[2], sys.argv[3]

            else:  # num of args is 3, k was not passed

                # check that goal is valid
                if sys.argv[1] not in ["spk", "wam", "ddg", "lnorm", "jacobi"]:
                    raise Exception

                # check that the file name ends with .txt or .csv
                if not (sys.argv[2].endswith(".txt") or sys.argv[2].endswith(".csv")):
                    raise Exception

                return 0, sys.argv[1], sys.argv[2]

        except Exception:
            print("Invalid Input!")
            sys.exit()


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
        if sum_of_min_distances == 0:
            probabilities[i] = 0
        else:
            probabilities[i] = distances[i] / sum_of_min_distances


def choose_centroids(data, k, num_of_data_points):

    chosen_indexes = list()
    distances = np.empty(num_of_data_points)
    probabilities = np.empty(num_of_data_points)

    # np.random.seed(0)
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

def format_decimal(num):
    digits_after_point = num.split('.')[1]
    num_of_digit = len(digits_after_point)
    return num + ((4 - num_of_digit) * "0")

if __name__ == '__main__':

    # get CMD arguments
    k, goal, file_path = get_arguments()

    # load data from input file
    input_matrix = np.loadtxt(file_path, delimiter=',')

    # input matrix metadata
    num_of_data_points = np.shape(input_matrix)[0]  # N
    data_point_size = np.shape(input_matrix)[1]  # number of entries in each data point

    # k must be less than or equal to the number of data points
    if num_of_data_points < k:
        print("Invalid Input!")
        sys.exit()

    matrix_T = mspkm.calculate_mat_T(input_matrix.flatten().tolist(), k, num_of_data_points, data_point_size, goal)

    if goal == 'spk':
        np.random.seed(0)
        np.set_printoptions(formatter={'float_kind': '{:.4f}'.format})

        if (k == 0):
            if len(matrix_T) == 0:
                print("An Error Has Occurred", end="")
                sys.exit()
            k = int(matrix_T.pop())
        else:
            matrix_T.pop()

        matrix_T = np.array(matrix_T).reshape(-1, k)
        num_of_data_points = np.shape(matrix_T)[0]
        data_point_size = np.shape(matrix_T)[1]
        reload(mspkm)

        # step 3: kmeans++ algorithm
        chosen_centroids_indexes = choose_centroids(matrix_T, k, num_of_data_points)
        initial_centroids = [matrix_T[idx].tolist() for idx in chosen_centroids_indexes]
        final_centroids = mspkm.fit(matrix_T.tolist(), initial_centroids, k, 300, num_of_data_points, data_point_size, 0)

        # print results: first indexes, then centroids
        print(",".join(str(x) for x in chosen_centroids_indexes))

        result = np.ndarray(shape=(k, data_point_size))
        for i in range(k):
            for j in range(data_point_size):
                result[i][j] = final_centroids[(i * data_point_size) + j]

        res = ""
        for i in result:
            for j in i:
                res += format_decimal(str(np.round(j, 4))) + ','
            print(res[:-1])
            res = ""
