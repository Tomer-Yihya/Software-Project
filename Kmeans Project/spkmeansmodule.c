#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include "spkmeans.h"


static PyObject* calculate_mat_T(PyObject *self, PyObject *args) {
    PyObject *python_matrix;
    int k;
    int num_of_dp;
    int dp_size;
    char *goal;
    if (!PyArg_ParseTuple(args, "Oiiis", &python_matrix, &k, &num_of_dp, &dp_size, &goal)) {
        return NULL;
    }

    /* create matrix T */
    dp_matrix.data = (double**) calloc(num_of_dp, sizeof(double*));
    assert(dp_matrix.data != NULL);

    for (i = 0; i < num_of_dp; i++) {
        dp_matrix.data[i] = (double*) calloc(dp_size, sizeof(double));
        assert(dp_matrix.data[i] != NULL);
    }

    dp_matrix.rows = num_of_dp;
    dp_matrix.columns = dp_size;

    /* convert py object to matrix struct */
    for (i = 0; i < num_of_dp; i++) {
        for (j = 0; j < dp_size ; j++) {
            PyObject *num = PyList_GetItem(python_matrix, (dp_size * i) + j);
            dp_matrix.data[i][j] = PyFloat_AsDouble(num);
        }
    }
    PyObject *lst = PyList_New(0);
    if (!lst) {
        return NULL;
    }

    if (strcmp(goal, "wam") == 0) {
        handle_wam(dp_matrix);
    } else if (strcmp(goal, "ddg") == 0) {
        handle_ddg(dp_matrix);
    } else if (strcmp(goal, "lnorm") == 0) {
        handle_lnorm(dp_matrix);
    } else if (strcmp(goal, "jacobi") == 0) {
        handle_jacobi(dp_matrix);
    } else if (strcmp(goal, "spk") == 0) {
        create_diagonal_matrix(&w, dp_matrix.rows, 0);
        get_weighted_adjacency_matrix(dp_matrix, &w);
        create_diagonal_matrix(&ddg, w.columns, 0);
        get_diagonal_degree_matrix(w, &ddg);
        sqrt_diag_matrix(&ddg);
        lnorm = get_lnorm_matrix(w, ddg);
        jacobi = get_jacobi(lnorm);
        eigen_values = get_eigen_values(jacobi);

        /* get k */
        sorted_eigenvalues = (double*) calloc(jacobi.rows, sizeof(double));
        assert(sorted_eigenvalues != NULL);
        for (i = 0; i < jacobi.rows; i++) {
            sorted_eigenvalues[i] = eigen_values[i];
        }
        bubble_sort(sorted_eigenvalues, jacobi.rows);

        if (k == 0) {
            k = calculate_k(sorted_eigenvalues, jacobi.rows);
        }

        if (k == 0 || k == 1) {
            return lst;
        }

        mat_T = get_T_mat(eigen_values, v, k);
        for (i = 0; i < mat_T.rows; i++) {
            for (j = 0; j < mat_T.columns; j++) {
                PyObject *num = PyFloat_FromDouble((double)mat_T.data[i][j]);
                if (!num) {
                    Py_DECREF(lst);
                    return NULL;
                }
                PyList_Append(lst, num);
            }
        }

        PyObject *num = PyFloat_FromDouble((double) k);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_Append(lst, num);
    }
    free(eigen_values);
    free_matrix_memory(&w);
    free_matrix_memory(&ddg);
    free(sorted_eigenvalues);
    free_matrix_memory(&lnorm);
    free_matrix_memory(&jacobi);
    free_matrix_memory(&dp_matrix);
    free_matrix_memory(&mat_T);
    return lst;
}

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

/* Python C-API */

static PyMethodDef Module_Methods[] = {
        {"fit", (PyCFunction) fit, METH_VARARGS, "run kmeans++ algorithm"},
        {"calculate_mat_T", (PyCFunction) calculate_mat_T, METH_VARARGS, "perform the request action by goal"},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef Moduledef = {
        PyModuleDef_HEAD_INIT,
        "mySpkmeansModule",
        NULL,
        -1,
        Module_Methods
};


PyMODINIT_FUNC PyInit_mySpkmeansModule(void) {
    PyObject *m;
    m = PyModule_Create(&Moduledef);
    if(!m) {
        return NULL;
    }
    return m;
}
