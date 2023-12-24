#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* structures definition */
typedef struct data_point {
    double* coordinates; /* (x1,x2,..,xn) */
    int index;           /* the index of closest centroid */
} data_point;

typedef struct centroid {
    data_point curr_centroid;
    data_point new_centroid;
    int num_of_related_points;
    double epsilon;
} centroid;

/* functions signatures*/
int main(int argc, char* argv[]);
void get_num_of_data_points(char *input_filename);
int get_data_point_size(char *line);
void generate_data_points(); 
void generate_k_centroids(); 
double* get_coordinates(char *s);
double calculate_distance(data_point point1, data_point point2);
data_point update_min_distance(data_point point);
centroid update_centroid(centroid c);
double* get_zeros_array();
int initialize_parameters(int argc, char *argv[]);
void write_output_to_file(centroid *centroids);
int sizeof_data_point();
int sizeof_centroid();


/* variables declartions */
const double EPSILON = 0.001;
FILE *file;
centroid *centroids;
data_point *points; 
data_point point, temp_point;
static int ERROR;
int i, j, k, r, max_iter, data_point_size, num_of_data_points, iteration_number, has_changed, counter;
int data_point_memory_size, centroid_memory_size, line_number, index;
double *coordinates, *temp_coordinates, *zeros_coordinates;
double temp_num, distance, min_distance, diff, dist, temp_coordinate;
char *line, *line_separated, *string_point, *delim, *input_filename, *output_filename;

int main(int argc, char* argv[]) {

    /* get k and max_iter or define max_iter == 200 */
    if (initialize_parameters(argc, argv) == 1) {
        return 1;
    }

    /* updates the number of data points and the size of each data point (use get_data_point_size) */
    get_num_of_data_points(input_filename);
    if (ERROR == 1) {
        printf("An Error Has Occurred\n");
        return ERROR;
    }
    free(line);

    /* calculate memory for each struct */
    data_point_memory_size = sizeof_data_point();
    centroid_memory_size = sizeof_centroid();

    /* memory for k centroids */
    centroids = (centroid*)calloc(k, centroid_memory_size);
    if (centroids == NULL) { 
        ERROR = 1;
        printf("An Error Has Occurred\n");
        return ERROR;
    }

    /* memory for all points */
    points = (data_point*)calloc(num_of_data_points, data_point_memory_size);
    if (points == NULL) {
        ERROR = 1;
        printf("An Error Has Occurred\n");
        return ERROR;
    }

    /* create data points array */
    generate_data_points();
    if (ERROR == 1) {
        printf("An Error Has Occurred\n");
        return ERROR;
    }

    /* create centroids array */
    generate_k_centroids();
    if (ERROR == 1) {
        printf("An Error Has Occurred\n");
        return ERROR;
    }

    /* main loop */
    iteration_number = 0;
    has_changed = 1; /* initiated to True */
    while(iteration_number < max_iter && has_changed == 1) {
        
        /* add data point to closest centroid and update closest centroid index in data point */
        for (i = 0; i < num_of_data_points; i++) {
            points[i] = update_min_distance(points[i]);
        }
        
        /* process the new information */
        for (i = 0; i < k; i++) {
            centroids[i] = update_centroid(centroids[i]);
            if (ERROR == 1) { 
                printf("An Error Has Occurred\n");
                return ERROR; 
            }
        }        

        /* check if all centroids epsilons are smaller than EPSILON = 0.001 */
        has_changed = 0;
        for (j = 0; j < k; j++) {
            if (centroids[j].epsilon >= EPSILON) {
                has_changed = 1;
            }
        }
        iteration_number++;
    }
    
    /* write all centroids to output file */
    write_output_to_file(centroids);

    free(coordinates);
    free(centroids);
    free(points);

    return ERROR;
}

/* get the number of data_points */
void get_num_of_data_points(char *input_filename) {
    counter = 0;
    file = fopen(input_filename, "r");
    if (file == NULL) {
        ERROR = 1;
        return;
    }
    
    line = (char *)(calloc(BUFSIZ, sizeof(char)));
    if (line == NULL) {
        ERROR = 1;
        return;
    }

    while (fgets(line, BUFSIZ, file) != NULL) {
        counter++;
    }
    
    num_of_data_points = counter;
    data_point_size = get_data_point_size(line);
    
    fclose(file);
}

/* get the length of a data_point */ 
int get_data_point_size(char *line) {
    counter = 1;
    for (i = 0; i < (int)strlen(line); i++) {
        if (line[i] == ',') {
            counter++;
        }
    }
    return counter;
}

/* create data_points array */
void generate_data_points() {
    
    file = fopen(input_filename, "r");
    if (file == NULL) { 
        ERROR = 1;
        return;
    }

    line = (char*)(calloc(BUFSIZ, sizeof(char)));
    if (line == NULL) {
        ERROR = 1;
        return;
    }
    
    line_number = 0;
    while (fgets(line, BUFSIZ, file) != NULL) {
        temp_coordinates = get_coordinates(line);
        if (ERROR == 1) {
            break;
        }
        points[line_number].coordinates = temp_coordinates;
        points[line_number].index = 0;
        line_number++;
    }
    free(delim);
    free(line);
    free(line_separated);
    fclose(file);
}

/* create centroids array */
void generate_k_centroids() {
    for (i = 0; i < k ; i++) {
        centroids[i].curr_centroid = points[i];
        centroids[i].epsilon = 2*EPSILON;
        centroids[i].num_of_related_points = 0;
        centroids[i].new_centroid.coordinates = get_zeros_array();
        if (ERROR == 1) { 
            return; 
        }
    }
}

/* create coodrinates for data point*/
double* get_coordinates(char* s) {
    
    coordinates = (double*)calloc(data_point_size, sizeof(double));
    if (coordinates == NULL) {
        ERROR = 1;
        return coordinates;
    }

    delim = (char*)(calloc(1, sizeof(char)));
    if (delim == NULL) { 
        ERROR = 1;
        return coordinates;
    }
    delim[0] = ',';
    
    line_separated = strtok(s, delim);
    j = 0;
    while (line_separated != NULL) {
        coordinates[j] = strtod(line_separated, NULL);
        j++;
        line_separated = strtok(NULL, delim);
    }
    
    return coordinates;
}

/* set all centroids.new centroid[i].coordinates as (0,0...0) */
double* get_zeros_array() {
    coordinates = (double*)calloc(data_point_size, sizeof(double));
    if (coordinates == NULL) { 
        ERROR = 1; 
    }
    return coordinates;
}

/* calculate the distance between points */
double calculate_distance(data_point point1, data_point point2) {
    dist = 0;
    for (r = 0; r < data_point_size; r++) {
		diff = (point1.coordinates[r] - point2.coordinates[r]);
        dist += diff * diff;
	}
	return dist;
}

/* 
   Calculates the distance between a data point and each of the centroids
   updates the centroid (num of releted and new center) and the data point (closest centroid)
*/
data_point update_min_distance(data_point point) {
	
    index = 0;
    min_distance = calculate_distance(point, centroids[0].curr_centroid);
    distance = min_distance;
    
    /* find the closest centroid index */ 
    for (j = 1; j < k; j++) {
        distance = calculate_distance(point, centroids[j].curr_centroid);
        if (distance < min_distance) {
			min_distance = distance;
            index = j;
		}
	}
    point.index = index;
    
    /* add point coordinates to new centroid coordinates */
    for (j = 0; j < data_point_size; j++) {
        centroids[index].new_centroid.coordinates[j] += point.coordinates[j];
    }
    
    /* increase the num of related points */
    centroids[index].num_of_related_points += 1;
    return point;
}

/*
    divide each centroid coordinates by the num of releted points,
    calculate epsilon for every centroid (difference between new center to old center),
    copy the values from new_centroid to curr_centroid,
    set new_centroid values to 0,
    set the centroid num_of_related_points to 0.
*/
centroid update_centroid(centroid c) {
    
    /* update each centroid / num of releted points */
    for (j = 0; j < data_point_size; j++) {
        c.new_centroid.coordinates[j] = c.new_centroid.coordinates[j] / c.num_of_related_points;
    }

    /* update epsilon for every centroid */ 
    c.epsilon = calculate_distance(c.new_centroid, c.curr_centroid);
    
    /* set new centroid values (curr_centroid = new_centroid) */         
    c.curr_centroid = c.new_centroid;
    c.new_centroid.coordinates = get_zeros_array();
    c.num_of_related_points = 0;
    
    return c;
}

/* get k, input filename, output filename and max_iter or define max_iter == 200 */
int initialize_parameters(int argc, char *argv[]) {
    
    if (argc == 4 || argc == 5) {
        k = atoi(argv[1]);
        if (argc == 4) {
            max_iter = 200;
            input_filename = argv[2];
            output_filename = argv[3];
        }
        /* argc == 5 */
        else {
            max_iter = atoi(argv[2]);
            input_filename = argv[3];
            output_filename = argv[4];
        }

        if (k < 0 || max_iter < 0) {
            printf("An Error Has Occurred\n");
            return 1;
        }
    }
    else {
        printf("Invalid Input!\n");
        return 1;
    }
    return 0;
}

/* write centroids to a new file */
void write_output_to_file(centroid* centroids) {

    file = fopen(output_filename, "ab+");
    if (file == NULL) {
        ERROR = 1;
        return;
    }
    
    for (i = 0; i < k; i++) {
        centroid temp_centroid = centroids[i];
        for (j = 0; j < data_point_size-1; j++) {
            temp_num = temp_centroid.curr_centroid.coordinates[j];
            fprintf(file, "%.4f,", temp_num);
        }
        fprintf(file, "%.4f\n", temp_centroid.curr_centroid.coordinates[data_point_size-1]);
    }
    fclose(file);
}

int sizeof_data_point() {
    int size = data_point_size * sizeof(double);    /* one point size */
    size += sizeof(int);                            /* index size */
    return size;
}

int sizeof_centroid() {
    int size = data_point_size * sizeof(double);    /* one point size */
    size = size * 2;                                /* we have two points */
    size += sizeof(double);                         /* epsilon size */
    size += sizeof(int);                            /* num_of_related_points size */
    return size;
}
