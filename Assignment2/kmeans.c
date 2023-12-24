#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>

/* structures definition */
typedef struct data_point {
    double* coordinates; /* (x1,x2,..,xn) */
    int index;           /* the index of closest centroid */
} data_point;

typedef struct centroid {
    double* curr_centroid;
    double* new_centroid;
    int num_of_related_points;
    double epsilon;
} centroid;

/* functions signatures*/
static PyObject* fit(PyObject *self, PyObject *args);
static double* listToPointer(PyObject* vector, int data_point_size);
double calculate_distance(double* point1, double* point2, int data_point_size);
data_point update_min_distance(data_point point, int data_point_size, int k);
centroid update_centroid(centroid c, int data_point_size);
int get_data_point_memory_size(int data_point_size);
int get_centroid_memory_size(int data_point_size);

/* variables declartions */
centroid *centroids;
data_point *points;
int i, j, r, iteration_number, has_changed, data_point_memory_size, centroid_memory_size, idx;
double distance, min_distance, diff, dist;
double *point;

static double* listToPointer(PyObject* vector, int data_point_size) {

    Py_ssize_t l;

    point = (double *)calloc(data_point_size, sizeof(double));
    if (point == NULL) {
        return NULL;
    }
    for (l = 0; l < data_point_size; l++) {
        PyObject *coordinate = PyList_GetItem(vector, l);
        point[l] = PyFloat_AsDouble(coordinate);
    }

    return point;
}

/* calculate the distance between points */
double calculate_distance(double* point1, double* point2, int data_point_size) {
    dist = 0;
    diff = 0;
    for (r = 0; r < data_point_size; r++) {
		diff = (point1[r] - point2[r]);
        dist += diff * diff;
	}
	return dist;
}

/* 
   Calculates the distance between a data point and each of the centroids.
   Updates the centroid (num of releted points and new center) and the data point (closest centroid)
*/
data_point update_min_distance(data_point datapoint, int data_point_size, int k) {
	
    idx = 0;
    min_distance = calculate_distance(datapoint.coordinates, centroids[0].curr_centroid, data_point_size);
    distance = min_distance;
    
    /* find the closest centroid index */ 
    for (j = 1; j < k; j++) {
        distance = calculate_distance(datapoint.coordinates, centroids[j].curr_centroid, data_point_size);
        if (distance < min_distance) {
			min_distance = distance;
            idx = j;
		}
	}
    datapoint.index = idx;
    
    /* add point coordinates to new centroid coordinates */
    for (j = 0; j < data_point_size; j++) {
        centroids[idx].new_centroid[j] += datapoint.coordinates[j];
    }
    
    /* increase the num of related points */
    centroids[idx].num_of_related_points += 1;
    return datapoint;
}

/*
    divide each centroid coordinates by the num of releted points,
    calculate EPSILON for every centroid (difference between new center to old center),
    copy the values from new_centroid to curr_centroid,
    set new_centroid values to 0,
    set the centroid num_of_related_points to 0.
*/
centroid update_centroid(centroid c, int data_point_size) {
    
    /* update each centroid / num of releted points */
    for (j = 0; j < data_point_size; j++) {
        c.new_centroid[j] = c.new_centroid[j] / c.num_of_related_points;
    }

    /* update epsilon for every centroid */ 
    c.epsilon = calculate_distance(c.new_centroid, c.curr_centroid, data_point_size);
    
    /* set new centroid values (curr_centroid = new_centroid) */
    for (j = 0; j < data_point_size; j++) {
        c.curr_centroid[j] = c.new_centroid[j];
    }

    c.num_of_related_points = 0;

    for (j = 0; j < data_point_size; j++) {
        c.new_centroid[j] = 0;
    }
    
    return c;
}

int get_data_point_memory_size(int data_point_size) {
    int size = data_point_size * sizeof(double);    /* one point size */
    size += sizeof(int);                            /* index size */
    return size;
}

int get_centroid_memory_size(int data_point_size) {
    int size = data_point_size * sizeof(double);    /* one point size */
    size = size * 2;                                /* we have two points */
    size += sizeof(double);                         /* epsilon size */
    size += sizeof(int);                            /* num_of_related_points size */
    return size;
}

static PyObject* fit(PyObject *self, PyObject *args) {
    PyObject *py_data;
    PyObject *py_centroids;
    int k;
    int max_iter;
    int num_of_data_points;
    int data_point_size;
    double EPSILON;
    if (!PyArg_ParseTuple(args, "OOiiiid", &py_data, &py_centroids, &k, &max_iter, &num_of_data_points, &data_point_size, &EPSILON)) {
        return NULL;
    }

    /* calculate memory for each struct */
    data_point_memory_size = get_data_point_memory_size(data_point_size);
    centroid_memory_size = get_centroid_memory_size(data_point_size);

    /* Alocate memory for k centroids */
    centroids = (centroid *)calloc(k, centroid_memory_size);
    if (centroids == NULL) {
        return NULL;
    }

    /* Alocate memory for num_of_data_points data points */
    points = (data_point *)calloc(num_of_data_points, data_point_memory_size);
    if (points == NULL) {
        return NULL;
    }

    /* initiate all centroid structs */
    for (i = 0; i < k; i++) {
        double *cur = (double *)calloc(data_point_size, sizeof(double));
        if (cur == NULL) {
            return NULL;
        }

        double *zeros = (double *)calloc(data_point_size, sizeof(double));
        if (zeros == NULL) {
            return NULL;
        }

        cur = listToPointer(PyList_GetItem(py_centroids, i), data_point_size);
        if (cur == NULL) {
            return NULL;
        }

        centroids[i].curr_centroid = cur;
        centroids[i].new_centroid = zeros;
        centroids[i].num_of_related_points = 0;
        centroids[i].epsilon = 2 * EPSILON;
    }

    /* crate an array of the data points */
    for (i = 0; i < num_of_data_points; i++) {

        double *coor = (double *)calloc(data_point_size, sizeof(double));
        if (coor == NULL) {
            return NULL;
        }
        coor = listToPointer(PyList_GetItem(py_data, i), data_point_size);
        if (coor == NULL) {
            return NULL;
        }

        points[i].coordinates = coor;
        points[i].index = 0;
    }

    /* main loop: kmeans algorithm */
    iteration_number = 0;
    has_changed = 1; /* initiated to True */
    while (iteration_number < max_iter && has_changed == 1) {

        /* update closest centroid index in data point */
        for (i = 0; i < num_of_data_points; i++) {
            points[i] = update_min_distance(points[i], data_point_size, k);
        }

        /* add data points to closest centroid and process centroid's data  */
        for (i = 0; i < k; i++) {
            centroids[i] = update_centroid(centroids[i], data_point_size);
        }

        /* check if all centroids' epsilons are smaller than EPSILON */
        has_changed = 0;
        for (i = 0; i < k; i++) {
            if (centroids[i].epsilon >= EPSILON) {
                has_changed = 1;
            }
        }
        iteration_number++;
    }

    /* create results object */
    PyObject *result = PyList_New(0);
    if (!result) {
        return NULL;
    }
    for (i = 0; i < k; i++) {
        for (j = 0; j < data_point_size; j++) {
            PyObject *num = PyFloat_FromDouble((double)centroids[i].curr_centroid[j]);
            if (!num) {
                Py_DECREF(result);
                return NULL;
            }
            PyList_Append(result, num);
        }
    }

    /* free data points memory */
    for (i = 0; i < num_of_data_points; i++) {
        free(points[i].coordinates);
    }

    /* free centroids memory */
    for (i = 0; i < k; i++) {
        free(centroids[i].curr_centroid);
        free(centroids[i].new_centroid);
    }

    free(points);
    free(centroids);
    return result;
}

static PyMethodDef mykmeansspMethods[] = {
        {"fit",
            (PyCFunction) fit,
                METH_VARARGS,
                    PyDoc_STR("run kmeans algorithm with C to get centroids")},
        {NULL,  NULL, 0, NULL}
};


static struct PyModuleDef _moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        mykmeansspMethods
};


PyMODINIT_FUNC
PyInit_mykmeanssp(void) {
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
