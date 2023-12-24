#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "spkmeans.h"

/* constants */

const double EPSILON = 0.00001;
const int MAX_ITERATION_NUMBER = 100;

/* functions implementations */

/* ----------------------------- handle functions ----------------------------- */
/* handle wam goal */
void handle_wam(matrix input_mat) {

    /* calculate affinity matrix and print it */
    create_diagonal_matrix(&w, input_mat.rows, 0);
    get_weighted_adjacency_matrix(input_mat, &w);
    print_matrix(w);
    free_matrix_memory(&w);
}

/* handle ddg goal */
void handle_ddg(matrix input_mat) {

    /* calculate diagonal degree matrix and print it */
    create_diagonal_matrix(&w, input_mat.rows, 0);
    get_weighted_adjacency_matrix(input_mat, &w);
    if (is_diagonal(w) == 1) {
        print_matrix(w);
    } else {
        create_diagonal_matrix(&ddg, w.rows, 0);
        get_diagonal_degree_matrix(w, &ddg);
        print_matrix(ddg);
        free_matrix_memory(&ddg);
    }
    free_matrix_memory(&w);
}

/* handle lnorm goal */
void handle_lnorm(matrix input_mat) {

    /* calculate ddg matrix and use it to calc lnorm, then print lnorm */
    create_diagonal_matrix(&w, input_mat.rows, 0);
    get_weighted_adjacency_matrix(input_mat, &w);
    if (is_diagonal(w) == 1) {
        print_matrix(w);
    } else {
        create_diagonal_matrix(&ddg, w.rows, 0);
        get_diagonal_degree_matrix(w, &ddg);
        sqrt_diag_matrix(&ddg);
        lnorm = get_lnorm_matrix(w, ddg);
        print_matrix(lnorm);
        free_matrix_memory(&ddg);
        free_matrix_memory(&lnorm);
    }
    free_matrix_memory(&w);
}

/* handle jacobi goal */
void handle_jacobi(matrix input_mat) {
    if (input_mat.rows != input_mat.columns) {
        printf("Invalid Input!\n");
    }
    else if (is_symmetric(input_mat) == 0) {
        printf("Invalid Input!\n");
    } else {
        jacobi = get_jacobi(input_mat);

        /* print eigenvalues */
        for (i = 0; i < jacobi.rows-1; i++) {
            if (jacobi.data[i][i] > -0.00005 && jacobi.data[i][i] <= 0) {
                printf("0.0000,");
            } else {
                printf("%.4f,", jacobi.data[i][i]);
            }
        }
        if (jacobi.data[jacobi.rows-1][jacobi.rows-1] > -0.00005 && jacobi.data[jacobi.rows-1][jacobi.rows-1] <= 0) {
            printf("0.0000\n");
        } else {
            printf("%.4f\n", jacobi.data[jacobi.rows-1][jacobi.rows-1]);
        }

        /* print eigenvectors (matrix v = P1xP2XP3...) */
        print_matrix(v);
        free_matrix_memory(&v);
        free_matrix_memory(&jacobi);
    }
}

/* ----------------------------- calculation functions ----------------------------- */
/* calculate affinity matrix */
void get_weighted_adjacency_matrix(matrix input_mat, matrix *w) {
    for (i = 0; i < input_mat.rows; i++) {
        for (j = i+1; j < input_mat.rows; j++) {
            weight = calculate_weight(input_mat.data[i], input_mat.data[j], input_mat.columns);
            w->data[i][j] = weight;
            w->data[j][i] = weight;
        }
    }
    w->rows = input_mat.rows;
    w->columns = input_mat.rows;
}

/* sum the squares of all off-diagonal elements of a matrix */
double calculate_matrix_off(matrix mat) {
    squares_sum = 0;
    for (i = 0; i < mat.rows; i++) {
        for (j = 0; j < mat.columns; j++) {
            if (i != j) {
                squares_sum += mat.data[i][j] * mat.data[i][j];
            }
        }
    }
    return squares_sum;
}

void sqrt_diag_matrix(matrix *mat) {
    for (i = 0; i < mat->rows; i++) {
        if (mat->data[i][i] != 0) {
            mat->data[i][i] = (1/(sqrt(mat->data[i][i])));
        }
    }
}

/* receive two matrices: A, A' and calculates off^2(A)-off^2(A') */
/* returns 1 (true) if the result is less than or equal to epsilon, or 0 (false) if it is greater than epsilon */
int does_converge(matrix a, matrix a_tag, double epsilon) {
    off_a = calculate_matrix_off(a);
    off_a_tag = calculate_matrix_off(a_tag);
    if (fabs(off_a - off_a_tag) <= epsilon) {
        return 1;
    }
    return 0;
}

/* receive a matrix and a row and returns the sum of the row's elements of this matrix */
double sum_matrix_row(matrix mat, int row) {
    sum = 0;
    for (j = 0; j < mat.columns; j++) {
        sum += mat.data[row][j];
    }
    return sum;
}

/* receive a matrix and a row and returns the sum of the row's square elements */
double sum_square_row(matrix mat, int row) {
    sum = 0;
    for (index_1 = 0; index_1 < mat.columns; index_1++) {
        sum += mat.data[row][index_1] * mat.data[row][index_1];
    }
    return sqrt(sum);
}

/* calculate D^(-1/2), the diagonal degree matrix */
void get_diagonal_degree_matrix(matrix w, matrix *ddg) {
    for (i = 0; i < w.rows; i++) {
        temp_num = sum_matrix_row(w, i);
        ddg->data[i][i] = temp_num;
    }
}

/* calculate lnorm  = I - D^-1/2 * w * D^-1/2 */
matrix get_lnorm_matrix(matrix w, matrix ddg) {
    create_diagonal_matrix(&unity, w.rows, 1);
    create_diagonal_matrix(&multiplication_result, w.rows, 0);
    create_diagonal_matrix(&subtraction_result, w.rows, 0);
    multiply_matrices(&ddg, &w, &multiplication_result);
    lnorm = multiplication_result;
    multiply_matrices(&lnorm, &ddg, &multiplication_result);
    lnorm = multiplication_result;
    subtract_matrices(&unity, &lnorm, &subtraction_result);
    lnorm = subtraction_result;

    free_matrix_memory(&unity);
    free_matrix_memory(&multiplication_result);

    return lnorm;
}

/* save pivot indexes- the off-diagonal element with the highest absolute value */
void get_pivot(matrix *a) {
    pivot_value = 0;
    for (i = 0; i < a->rows; i++) {
        for (j = i+1; j < a->columns; j++) {
            if (i != j) {
                if (fabs(a->data[i][j]) > pivot_value) {
                    pivot_i = i;
                    pivot_j = j;
                    pivot_value = fabs(a->data[i][j]);
                }
            }
        }
    }
}

/* calculate the rotation matrix P */
void get_p_matrix(matrix *a, matrix *p) {
    for (i = 0; i < p->rows; i++) {
        for (j = 0; j < p->columns; j++) {
            p->data[i][j] = i == j ? 1 : 0;
        }
    }
    get_pivot(a);
    theta = (a->data[pivot_j][pivot_j] - a->data[pivot_i][pivot_i]) / (2 * a->data[pivot_i][pivot_j]);
    sign_theta = theta < 0 ? -1 : 1;
    t = sign_theta / (fabs(theta) + pow((theta * theta + 1), 0.5));
    c = 1 / pow((t * t + 1), 0.5);
    s = c * t;
    p->data[pivot_i][pivot_i] = c;
    p->data[pivot_j][pivot_j] = c;
    p->data[pivot_j][pivot_i] = -s;
    p->data[pivot_i][pivot_j] = s;
}

/* run jacobi procedure */
matrix get_jacobi(matrix a) {

    create_diagonal_matrix(&v, a.rows, 1);
    create_diagonal_matrix(&a_tag, a.rows, 1);
    if (is_diagonal(a) == 1) {
        copy_matrix_values(&a, &a_tag);
        jacobi = a_tag;
    } else {
        create_diagonal_matrix(&p_t, a.rows, 1);
        create_diagonal_matrix(&p, a.rows, 1);
        create_diagonal_matrix(&multiplication_result, a.rows, 0);
        iteration_number = 0;
        while (iteration_number < MAX_ITERATION_NUMBER) {
            get_p_matrix(&a, &p);
            multiply_matrices(&v, &p, &multiplication_result);
            copy_matrix_values(&multiplication_result, &v);
            transpose_matrix(&p, &p_t);
            multiply_matrices(&p_t, &a, &multiplication_result);
            multiply_matrices(&multiplication_result, &p, &a_tag);
            if (does_converge(a, a_tag, EPSILON)) {
                break;
            }
            copy_matrix_values(&a_tag, &a);
            iteration_number++;
        }
        jacobi = a_tag;
        free_matrix_memory(&p);
        free_matrix_memory(&p_t);
        free_matrix_memory(&multiplication_result);
    }
    return jacobi;
}

/* return a double array of eigenvalues */
double * get_eigen_values(matrix mat) {
    eigen_values = (double*)calloc(mat.rows, sizeof(double));
    assert(eigen_values != NULL);
    for (i = 0; i < mat.rows; i++) {
        eigen_values[i] = mat.data[i][i];
    }
    return eigen_values;
}

/* use heuristic to calculate k when k = 0 */
int calculate_k(double *eigen_values, int size) {
    k_diff = (double*) calloc(size/2, sizeof(double));
    assert(k_diff != NULL);
    for (i = 0; i < (size/2); i++) {
        k_diff[i] = fabs(eigen_values[i+1] - eigen_values[i]);
    }
    max_k_diff = get_max_k_diff(k_diff, size/2);
    free(k_diff);
    return max_k_diff + 1;
}

/* helper function for calculate k */
int get_max_k_diff(double *k_diff, int size) {
    k_index = 0;
    for (i = 0; i < size; i++)
        if (k_diff[i] > k_diff[k_index]) {
            k_index = i;
        }
    return k_index;
}

/* implementation of bubble sort algorithm: sorting both arrays according to main array */
void two_arrays_bubble_sort(double *main_arr, int size, double *secondary_arr) {
    for (i = 0; i < size-1; i++) {
        for (j = 0; j < size-i-1; j++) {
            if (main_arr[j] > main_arr[j+1]) {
                swap(&main_arr[j], &main_arr[j + 1]);
                swap(&secondary_arr[j], &secondary_arr[j + 1]);
            }
        }
    }
}

/* implementation of bubble sort algorithm */
void bubble_sort(double *main_arr, int size) {
    for (i = 0; i < size-1; i++) {
        for (j = 0; j < size-i-1; j++) {
            if (main_arr[j] > main_arr[j+1]) {
                swap(&main_arr[j], &main_arr[j + 1]);
            }
        }
    }
}

/* bubble sort helper function */
void swap(double *num1, double *num2) {
    temp_num = *num1;
    *num1 = *num2;
    *num2 = temp_num;
}

/* create matrix U and normalize it to get matrix T */
matrix get_T_mat(double *eigenvalues, matrix v, int k) {

    /* create indexes array for sorting */
    indexes_arr = (double*)calloc(v.rows, sizeof(double));
    assert(indexes_arr != NULL);
    int_indexes_arr = (int*)calloc(v.rows, sizeof(int));
    assert(int_indexes_arr != NULL);
    for (i = 0; i < v.rows; i++) {
        indexes_arr[i] = i;
    }

    /* memory for matrices U, T  */
    mat_U.data = (double**)calloc(v.rows, sizeof(double*));
    assert(mat_U.data != NULL);
    mat_T.data = (double**)calloc(v.rows, sizeof(double*));
    assert(mat_T.data != NULL);
    for (i = 0; i < v.rows; i++) {
        mat_U.data[i] = (double*)calloc(k, sizeof(double));
        assert(mat_U.data[i] != NULL);
        mat_T.data[i] = (double*)calloc(k, sizeof(double));
        assert(mat_T.data[i] != NULL);
    }

    /* initialize matrices U, T parameters */
    mat_U.rows = v.rows;
    mat_U.columns = k;
    mat_T.rows = v.rows;
    mat_T.columns = k;

    /* sort eigenvalues and eigenvectors indexes */
    two_arrays_bubble_sort(eigenvalues, v.rows, indexes_arr);
    double_array_to_int_array(indexes_arr, int_indexes_arr, v.rows);

    /* create matrix U from eigenvectors corresponding to k smallest eigenvalues */
    for (i = 0; i < k; i++) {
        for (j = 0; j < v.rows; j++) {
            mat_U.data[j][i] = v.data[j][int_indexes_arr[i]];
        }
    }

    /* create matrix T from matrix U */
    for (i = 0; i < mat_U.rows; i++) {
        sq_sum = sum_square_row(mat_U, i);
        if (sq_sum != 0) {
            for (j = 0; j < k; j++) {
                mat_T.data[i][j] = (mat_U.data[i][j] / sq_sum);
            }
        }
    }

    free(int_indexes_arr);
    free(indexes_arr);
    free_matrix_memory(&mat_U);
    return mat_T;
}

/* create a matrix object from a file */
matrix build_input_mat(FILE *f, int rows, int columns) {
    input_mat.rows = rows;
    input_mat.columns = columns;
    input_mat.data =(double **) calloc(rows, sizeof(double*));
    assert(input_mat.data != NULL);
    for(i = 0; i < rows; i++) {
        input_mat.data[i] = (double*) calloc(columns, sizeof(double));
        assert(input_mat.data[i] != NULL);
    }
    delim = ",";
    line = (char *) calloc(BUFSIZ ,sizeof(char));
    assert(line != NULL);
    for(i = 0; i < rows; i++) {
        fgets(line, BUFSIZ, f);
        separated_line = strtok(line, delim);

        j = 0;
        while (separated_line != NULL) {
            input_mat.data[i][j] = strtod(separated_line, NULL);
            j++;
            separated_line = strtok(NULL, delim);
        }
    }
    free(separated_line);
    free(line);
    return input_mat;
}

/* ----------------------------- helper functions ----------------------------- */

int get_file_columns(FILE *f) {
    counter = 1;
    line = (char*) calloc(BUFSIZ, sizeof(char));
    assert(line != NULL);
    fgets(line, BUFSIZ, f);
    for(i = 0; i < (int)strlen(line); ++i) {
        if (line[i] == ',') {
            counter++;
        }
    }
    rewind(f);
    free(line);
    return counter;
}

int get_file_rows(FILE *f) {
    counter = 1;
    line = (char*) calloc( BUFSIZ, sizeof(char));
    assert(line != NULL);
    fgets(line, BUFSIZ, f);
    while (fgets(line, BUFSIZ, f) != NULL) {
        counter++;
    }
    rewind(f);
    free(line);
    return counter;
}

void copy_matrix_values(matrix *source, matrix *destination) {
    assert(source->rows == destination->rows);
    assert(source->columns == destination->columns);
    for (i = 0; i < source->rows; i++) {
        for (j = 0; j < source->columns; j++) {
            destination->data[i][j] = source->data[i][j];
        }
    }
}

void double_array_to_int_array(double *double_arr, int *int_arr, int arr_size) {
    for (i = 0; i < arr_size; i++) {
        int_arr[i] = (int)double_arr[i];
    }
}

/* receive two vectors and returns their euclidean distance */
double calculate_euclidean_distance(double *x_i, double *x_j, int len) {
    distance = 0;
    for (index_1 = 0; index_1 < len; index_1++) {
        distance += (x_i[index_1] - x_j[index_1]) * (x_i[index_1] - x_j[index_1]);
    }
    return sqrt(distance);
}

/* receive two vectors and returns weight in the affinity matrix */
double calculate_weight(double *x_i, double *x_j, int len) {
    euclidean_weight = calculate_euclidean_distance(x_i, x_j, len);
    return exp(-1 * (euclidean_weight * 0.5));
}

/* print a double array to the screen */
void print_double_array(double *arr, int len) {
    for (i = 0; i < len-1; i++) {
        if (arr[i] < 0 && arr[i] > -0.0001) {
            arr[i] = 0.0;
        }
        printf("%.4f ", arr[i]);
    }
    printf("%.4f ", arr[len-1]);
    printf("\n");
}

/* print an int array to the screen */
void print_int_array(int *arr, int len) {
    for (i = 0; i < len-1; i++) {
        if (arr[i] < 0 && arr[i] > -0.0001) {
            arr[i] = 0.0;
        }
        printf("%d ", arr[i]);
    }
    printf("%d ", arr[len-1]);
    printf("\n");
}

/* ----------------------------- matrices functions ----------------------------- */

/* receive a matrix, return 1 if it is symmetric or 0 otherwise */
int is_symmetric(matrix mat) {
    if (mat.rows == 1) {
        return 1;
    }
    for (i = 0; i < mat.rows; i++) {
        for (j = i+1; j < mat.columns; j++) {
            if (mat.data[i][j] != mat.data[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

/* receive a matrix, return 1 if it is diagonal or 0 otherwise */
int is_diagonal(matrix mat) {
    for (i = 0; i < mat.rows; i++) {
        for (j = 0; j < mat.columns; j++) {
            if ((i != j) && (mat.data[i][j] != 0)) {
                return 0;
            }
        }
    }
    return 1;
}

/* receive a matrix and print it (for debugging purposes) */
void print_matrix(matrix mat) {
    for (i = 0; i < mat.rows; i++) {
        for (j = 0; j < mat.columns-1; j++) {
            if (mat.data[i][j] < 0 && mat.data[i][j] == 0.0) {
                printf("0.0000,");
            } else {
                printf("%.4f,", mat.data[i][j]);
            }
        }
        if (mat.data[i][mat.columns-1] < 0 && mat.data[i][mat.columns-1] == 0.0) {
            printf("0.0000\n");
        } else {
            printf("%.4f\n", mat.data[i][mat.columns-1]);
        }
    }
}

/* calculate the subtraction of two matrices */
void subtract_matrices(matrix *mat1, matrix *mat2, matrix *mat3) {
    assert(mat1->rows == mat2->rows);
    assert(mat1->columns == mat2->columns);

    for (i = 0; i < mat1->rows; i++) {
        for (j = 0; j < mat1->columns; j++) {
            mat3->data[i][j] = mat1->data[i][j] - mat2->data[i][j];
        }
    }
}

void create_diagonal_matrix(matrix *mat, int size, double diag_val) {

    /* initialize an array of row double arrays (the rows) */
    mat->data = (double **)calloc(size, sizeof(double*));
    assert(mat->data != NULL);

    /* make each row to point to an array of doubles of size column */
    for (i = 0; i < size; i++) {
        mat->data[i] = (double*)calloc(size, sizeof(double));
        assert(mat->data[i] != NULL);
    }

    /* update matrix attributes and data */
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            mat->data[i][j] = i == j ? diag_val : 0;
        }
    }
    mat->rows = size;
    mat->columns = size;
}

/* receive a matrix object and return its transposed form */
void transpose_matrix(matrix *original, matrix *transposed) {
    for (i = 0; i < original->rows; i++) {
        for (j = 0; j < original->columns; j++) {
            transposed->data[j][i] = original->data[i][j];
        }
    }
}

/* free matrix memory */
void free_matrix_memory(matrix *mat) {
    if (mat->data == NULL) {
        return;
    }
    for (i = 0; i < mat->rows; i++) {
        free(mat->data[i]);
    }
    free(mat->data);
    mat->data = NULL;
}

/* receive two matrices and return its multiplication result */
void multiply_matrices(matrix *mat1, matrix *mat2, matrix *mat3) {
    matrix temp;
    assert(mat1->columns == mat2->rows);
    assert(mat1->rows == mat3->rows);
    assert(mat2->columns == mat3->columns);
    create_diagonal_matrix(&temp, mat1->rows, 0);
    for (n = 0; n < mat1->rows; n++) {
        for (i = 0; i < mat2->columns; i++) {
            for (j = 0; j < mat2->rows; j++) {
                temp.data[n][i] += mat1->data[n][j] * mat2->data[j][i];
            }
        }
    }
    copy_matrix_values(&temp, mat3);
    free_matrix_memory(&temp);
}

/* ----------------------------- kmpeans++ ----------------------------- */

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

/* ----------------------------- main ----------------------------- */
int main(int argc, char* argv[]) {

    if (argc < 3 || argc > 4) {
        printf("Invalid Input!\n");
        return 0;
    } else {

        if (argc == 3) {
            k = 0;
            goal = argv[1];
            file = argv[2];
        } else {
            k = atoi(argv[1]);
            if (k < 0) {
                printf("Invalid Input!\n");
                return 0;
            }
            goal = argv[2];
            file = argv[3];
        }
    }

    f = fopen(file, "r");
    input_rows = get_file_rows(f);

    if (k > input_rows) {
        printf("Invalid Input!\n");
        return 0;
    }

    input_columns = get_file_columns(f);
    input_mat = build_input_mat(f, input_rows, input_columns);

    if (strcmp(goal, "wam") == 0) {
        handle_wam(input_mat);
    } else if (strcmp(goal, "ddg") == 0) {
        handle_ddg(input_mat);
    } else if (strcmp(goal, "lnorm") == 0) {
        handle_lnorm(input_mat);
    } else if (strcmp(goal, "jacobi") == 0) {
        handle_jacobi(input_mat);
    } else {
        printf("Invalid Input!\n");
    }

    /* free memory */
    free_matrix_memory(&input_mat);
    fclose(f);
    return 0;
}
