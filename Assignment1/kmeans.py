#!/usr/bin/env python3

import sys
import math


# constants
EPSILON = 0.001


class DataPoint:

    def __init__(self, points: list):
        self.points = points

    def get_closest_centroid(self, centroids: list):
        min_distance = calculate_distance(self.points, centroids[0].center)
        closest_centroid = centroids[0]
        for i in range (1, len(centroids)):
            distance = calculate_distance(self.points, centroids[i].center)
            if distance < min_distance:
                min_distance = distance
                closest_centroid = centroids[i]
        return closest_centroid

    def __repr__(self):
        return str(self.points)


class Centroid:

    def __init__(self, datapoint: DataPoint):
        self.center = datapoint.points
        self.related_datapoints = list()

    def update_centroid(self):
        num_of_related_datapoints = len(self.related_datapoints)
        sums = [0 for i in range(len(self.center))]
        for dp in self.related_datapoints:
            for idx in range(len(sums)):
                sums[idx] += dp.points[idx]
        new_center = [x/num_of_related_datapoints for x in sums]
        changed = calculate_distance(self.center, new_center) >= EPSILON
        self.center = new_center
        self.related_datapoints = list()
        return changed

    def add_related_datapoint(self, datapoint: DataPoint):
        self.related_datapoints.append(datapoint)

    def __repr__(self):
        return f'Center: {str(self.center)}'


def calculate_distance(vector1, vector2):
    dist = 0
    for i in range(len(vector1)):
        sub = vector1[i] - vector2[i]
        dist += sub ** 2
    return dist


if __name__ == "__main__":

    # function must get 4 or 5 arguments
    if len(sys.argv) > 5 or len(sys.argv) < 4:
        print("Invalid Input!")
        sys.exit()

    else:

        K = int(sys.argv[1])  # first argument is always K

        if len(sys.argv) == 4:
            max_iter = 200  # default value
            input_filename = sys.argv[2]
            output_filename = sys.argv[3]
        else:
            max_iter = int(sys.argv[2])
            input_filename = sys.argv[3]
            output_filename = sys.argv[4]

        if (K < 0 or max_iter < 0):
            print("An Error Has Occurred")
            sys.exit()

        # open file and read datapoints
        try:
            file = open(input_filename, "r")
        except:
            print("An Error Has Occurred")
            sys.exit()

        lines = file.readlines()
        file.close()
        datapoints = [DataPoint([float(x) for x in line.replace("\n", "").split(',')]) for line in lines if line != '']

        centroids = [Centroid(datapoints[i]) for i in range(K)]

        stable = False
        iteration_number = 0

        while not stable or iteration_number < max_iter:

            # update vectors clusters according to minimum distance from centroids
            for datapoint in datapoints:
                closest_centroid = datapoint.get_closest_centroid(centroids)
                closest_centroid.add_related_datapoint(datapoint)

            # update centroids
            centroids_are_stable = [centroid.update_centroid() for centroid in centroids]
            stable = True not in centroids_are_stable

            iteration_number += 1

        # create output file
        output_file = open(output_filename, "w")
        for centroid in centroids:
            for i in range(len(centroid.center)):
                temp_num = math.floor(centroid.center[i] * 10000)
                centroid.center[i] = format(float(temp_num / 10000), '.4f')

            output_file.write(f'{",".join(centroid.center)}\n')

        output_file.close()
