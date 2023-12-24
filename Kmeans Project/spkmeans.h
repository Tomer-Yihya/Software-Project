#ifndef FINALPROJECT_SPKMEANS_H
#define FINALPROJECT_SPKMEANS_H

#endif

/* data structures */

typedef struct {
    int rows;
    int columns;
    double ** data;
} matrix;

typedef struct {
    double* coordinates; /* (x1,x2,..,xn) */
    int index;           /* the index of closest centroid */
} data_point;

typedef struct {
    double* curr_centroid;
    double* new_centroid;
    int num_of_related_points;
    double epsilon;
} centroid;

/* variables declarations */

FILE *f;
data_point *points;
centroid *centroids;
char *separated_line, *delim, *goal, *line, *file;
int *int_arr, *int_indexes_arr, input_rows, input_columns, counter;
int has_changed, data_point_memory_size, centroid_memory_size, idx;
int i, j, index_1, k, n, r, pivot_i, pivot_j, iteration_number, max_k_diff, k_index;
double **matrix_data, **d1, **d2;
double *eigen_values, *k_diff, *arr, *indexes_arr, *sorted_eigenvalues, *point;
double distance, min_distance, diff, dist, euclidean_weight, pivot_value, sign_theta;
double s, c, t, theta, squares_sum, temp_num, sum, sq_sum, weight, off_a, off_a_tag, distance;
matrix unity, p, p_t, v, jacobi, a, a_tag, mat_T, mat_U, transposed, ddg;
matrix multiplication_result, subtraction_result, temp_mat, w, lnorm, dp_matrix, input_mat;

/* functions signatures */

int get_file_rows(FILE *f);
int is_diagonal(matrix mat);
int is_symmetric(matrix mat);
int get_file_columns(FILE *f);
int get_max_k_diff(double *k_diff, int size);
int calculate_k(double *eigen_values, int size);
int get_centroid_memory_size(int data_point_size);
int get_data_point_memory_size(int data_point_size);
int does_converge(matrix a, matrix a_tag, double epsilon);

void get_pivot(matrix *a);
void print_matrix(matrix mat);
void handle_wam(matrix input_mat);
void handle_ddg(matrix input_mat);
void sqrt_diag_matrix(matrix *mat);
void handle_lnorm(matrix input_mat);
void handle_jacobi(matrix input_mat);
void free_matrix_memory(matrix *mat);
void swap(double *num1, double *num2);
void get_p_matrix(matrix *a, matrix *p);
void print_int_array(int *arr, int len);
void bubble_sort(double *main_arr, int size);
void print_double_array(double *arr, int len);
void get_diagonal_degree_matrix(matrix w, matrix *ddg);
void transpose_matrix(matrix *original, matrix *transposed);
void copy_matrix_values(matrix *source, matrix *destination);
void get_weighted_adjacency_matrix(matrix input_mat, matrix *w);
void subtract_matrices(matrix *mat1, matrix *mat2, matrix *mat3);
void multiply_matrices(matrix *mat1, matrix *mat2, matrix *mat3);
void create_diagonal_matrix(matrix *mat, int size, double diag_val);
void double_array_to_int_array(double *double_arr, int *int_arr, int arr_size);
void two_arrays_bubble_sort(double *main_arr, int size, double *secondary_arr);

double calculate_matrix_off(matrix mat);
double sum_matrix_row(matrix mat, int row);
double sum_square_row(matrix mat, int row);
double calculate_weight(double *x_i, double *x_j, int len);
double calculate_euclidean_distance(double *x_i, double *x_j, int len);
double calculate_distance(double* point1, double* point2, int data_point_size);

double * get_eigen_values(matrix jacobi_mat);

centroid update_centroid(centroid c, int data_point_size);

data_point update_min_distance(data_point point, int data_point_size, int k);

matrix get_jacobi(matrix a);
matrix get_lnorm_matrix(matrix w, matrix ddg);
matrix build_input_mat(FILE *f, int rows, int columns);
matrix get_T_mat(double *eigenvalues, matrix eigenvectors_mat, int k);
