#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <math.h>

/* Author: Johan Jacob */

static int dim;
static  int * distances;
static double * pheromone;
static double * eta_values;
static double * nodes;
static double t;
static int * tour;
static int tour_cost;


/*The structure of an ant
* Holds, the starting node, the index to next empty slot in array,
* the array trail which contains the nodes/ cities visited  in order by this ant
* and the cost of the tour made by this ant
 */
typedef struct {
    int start_node;
    int index;
    int * trail;
    double cost;
} ANT;

// function declarations
int main(int argc, char const *argv[]);
FILE * read_words(FILE * f);
FILE * read_coords(FILE * f);
void calculate_distance();
int calculate_solution_cost(int * sol, int size);
void swap_places(int * sol, int i, int j);
void nearest_neighbor(int * sol,int  * new_sol);
int check_valid(int*sol);
void acs(int swarm_count, int * shortest_path);
void global_pheromone_update(double alpha);
int * two_opt(int * sol);
int * two_opt_first(int * sol, int  * distances );
void  two_opt_swap(int *s, int i, int k);
int * hill_climbing(int * sol, int  * distances );
int compute_gain(int * distances, int * sol ,int i, int j);
void exchange(int * sol, int best_i, int best_j);


FILE * read_words(FILE * f){
    char string[100];
    fscanf(f, "%s",string);
    while (strcmp(string,"NODE_COORD_SECTION")!=0) {
        if (strcmp(string,"DIMENSION:")==0 || strcmp(string,"DIMENSION")==0 ) {
            fscanf(f, "%s",string);
            if (strcmp(string,":")==0 ) {
                fscanf(f, "%s",string);
            }
            dim = atoi(string);
        }
        fscanf(f, "%s", string);
    }
    return f;
}
FILE * read_coords(FILE * f){
    int index[1];
    double x[1];
    double y[1];
    for (size_t i = 0; i < dim; i++) {
        fscanf(f, "%d", index);
        fscanf(f, "%lf", x);
        fscanf(f, "%lf", y);
        *(nodes+i) = x[0];
        *(nodes+dim+i) = y[0];
    }
    return f;
}

void calculate_distance(){
    double x_sq, y_sq, value;

    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            if (j == i) {
                *(distances + (i*dim + j)) = INT_MAX;
            }else{
                value = *(nodes+i) - *(nodes+j);
                x_sq = value * value;
                value = *(nodes+dim+i) - *(nodes+dim+j);
                y_sq = value * value;
                value = round(sqrt(x_sq + y_sq));
                *(distances + (i*dim + j)) = value;
                *(distances+(j*dim + i)) = value;
            }
        }
    }
}
int calculate_solution_cost(int * sol, int size){
    int sum = 0;
    for (size_t i = 0; i < dim - 1; i++) {
        sum += *(distances+(sol[i]*dim + sol[i+1]));
    }
    sum += *(distances+(sol[dim-1]*dim + sol[0]));
    return sum;
}

void swap_places(int * sol, int i, int j){
    int temp = sol[j];
    sol[j] = sol[i];
    sol[i] = temp;
}


void nearest_neighbor(int * sol,int  * new_sol){
    int start_node = rand()%dim; // index of randomly selected node
    int visited[dim];
    for (int i = 0; i < dim; i++) {
        visited[i] = 0;
    }
    visited[sol[start_node]] = 1;
    new_sol[0] = sol[start_node];
    int nearest;
    int value;
    int index = sol[start_node];
    for (int j = 1; j < dim; j++) {
        value = *distances;
        for (size_t k = 0; k < dim; k++) {
            if (*(distances+(index*dim + k)) <= value && visited[k] == 0) {
                nearest = k;
                value = *(distances+(index*dim + nearest));
            }
        }
        index = nearest;
        visited[index] = 1;
        new_sol[j] = index;
    }
}

int check_valid(int*sol){
    int visited[dim];
    for (size_t i = 0; i < dim; i++) {
        visited[i] = 0;
    }
    for (size_t j = 0; j < dim; j++) {
        int node = sol[j];
        visited[node] = visited[node] += 1;
    }
    for (int k = 0; k < dim; k++) {
        if (visited[k] != 1) {
            return 1;
        }
    }
    return 0;
}
void global_pheromone_update(double alpha){
    int node1,node2;
    double phi, delta;
    delta = 1.0/tour_cost;
    for (size_t i = 0; i < dim-1; i++) {
        node1 = tour[i];
        node2 = tour[i+1]%dim;
        phi =  *(pheromone+ (node1*dim + node2));
        *(pheromone+ (node1*dim + node2)) = ((1 - alpha)*phi)+(alpha*delta);
        *(pheromone+ (node2*dim + node1))= *(pheromone+ (node1*dim + node2));
        *(eta_values +(node1*dim + node2)) = *(pheromone+ (node1*dim + node2)) * (1.0/ *(distances + (node1 *dim + node2)))* 1.0/ *(distances + (node1*dim + node2));
        *(eta_values +(node2*dim + node1)) = *(eta_values +(node1*dim + node2));
    }
}

int * two_opt(int * sol){
    int best_gain = 1;
    int gain = 0;
    int i, j;
    int best_i, best_j = -1;
    while (best_gain !=0) {
        best_gain = 0;
        for (i = 0; i < dim -1 ; i++) {
            for (j = i+2; j < dim; j++) {
                gain = compute_gain(distances, sol ,i,j);
                if (gain < best_gain) {
                    best_gain = gain;
                    best_i = i;
                    best_j = j;
                }
            }
        }
        if (best_gain < 0) {
            two_opt_swap(sol, best_i, best_j);
        }
    }
    return sol;
}
void two_opt_swap(int *s, int i, int k){
    int temp;
    int a = (i + 1)%dim;
    int b = k;
    while (a < b) {
        temp = s[a];
        s[a] = s[b];
        s[b] = temp;
        a = (a + 1)%dim;
        b = (b - 1)%dim;
    }
}
int compute_gain(int * distances, int * sol ,int i, int j){
    if (i == j || i == (j+1)%dim || (i+1)%dim == j || (i+1)%dim == (j+1)%dim) {
        return INT_MAX;
    }
    int gain = *(distances + (sol[i]*dim + sol[j])) + *(distances + (sol[(i+1)%dim]*dim + sol[(j+1)%dim])) ;
    int old_gain = *(distances + (sol[j]*dim + sol[(j+1)%dim]))  + *(distances + (sol[i]*dim + sol[(i+1)%dim]));
    return (gain - old_gain);
}

void acs(int swarm_count, int * shortest_path){
    time_t start, end;
    double time_duration = 0;
    start = time(NULL);
    int added[dim];
    ANT* nest[swarm_count];
    double S[dim];
    ANT * a1;
    int i,j,k, node,pos,q,valid;
    double value,beta,sigma,max,min,alpha,d, rnd;

    /* We adjust the evaporation alpha depending on the size of the tour*/
    if (dim > 101) {
        alpha = 0.01;
    }else{
        alpha = 0.1;
    }
    beta  = 2.0;
    min = tour_cost;
    for (i = 0; i < swarm_count; i++) {
        // initialise ant
        nest[i] = (ANT*)malloc(sizeof(ANT));
        a1= nest[i];
        a1->trail = (int*)malloc(dim*sizeof(int));
    }
    // loop for 3 minutes
    while (time_duration < 180) {
        end = time(NULL);
        time_duration = difftime(end, start);
        for (i = 0; i < swarm_count; i++) {
            // initialise ant
            a1= nest[i];
            node = rand()%dim;
            a1->start_node = node;
            a1->index = 0;
            a1->cost = INT_MAX;
        }
        for ( i = 0; i < swarm_count; i++) {
            a1 = nest[i];
            pos = a1->start_node;
            sigma = 0;
            for ( j = 0; j < dim; j++) {
                added[j] = 0;
                a1->trail[j] = 0;
            }
            a1->trail[a1->index] = a1->start_node;
            added[a1->start_node] = 1;
            a1->index++;
            // Exploitation
            for (j = 0; j < dim; j++) {
                q = rand()%100;
                max = -1;
                if (q < 95) { // 95% of the time
                    // find next node
                    for ( k = 0; k < dim; k++) {
                        value = *(eta_values +(pos*dim + k));

                        if (value > max && added[k] ==0 && pos != k) {
                            max = value;
                            node = k;
                        }
                    }
                }
                // Exploration
                else{
                    sigma = 0;
                    // calculate probability
                    for (k = 0; k < dim; k++) {
                        if (added[k] == 0) {
                            value = *(eta_values +(pos*dim + k));
                            S[k] = value;
                            sigma += value;
                        }else{
                            S[k] = 0;
                        }
                    }
                    for (k = 0; k < dim; k++) {
                        S[k] = S[k]/sigma;
                    }
                    rnd = (double)rand()/((double)RAND_MAX);
                    d = 0;
                    for (k = 0; k < dim; k++){
                       d += S[k];
                       if (rnd < d && S[k] != 0){
                           node = k;
                           k = dim;
                       }
                   }
                }

                if (added[node] == 0) {
                    a1->trail[a1->index] = node;
                    // local pheromone update (symmetric)
                    *(pheromone+ (pos*dim + node))= ((1 - (alpha))*( *(pheromone+ (pos*dim + node))) + (alpha * t));
                    *(pheromone+ (node*dim + pos)) = *(pheromone+ (pos*dim + node));
                    *(eta_values +(pos*dim + node)) = *(pheromone+ (pos*dim + node)) * (1.0/ *(distances + (pos*dim + node)))* 1.0/ *(distances + (pos*dim + node));
                    *(eta_values +(node*dim + pos)) = *(eta_values +(pos*dim + node));
                    a1->index++;
                    added[node] = 1;
                    pos = node;
                }
            }
            // Check validity of solution
            valid = check_valid(a1->trail);
            if (valid == 0) {
                a1->cost = calculate_solution_cost(a1->trail, dim);
            }else{
                a1->cost = INT_MAX;
            }
        }

        for (i = 0; i < swarm_count; i++) {
            if (nest[i]->cost < min) {
                min = nest[i]->cost;
                tour_cost = min;
                for (size_t j = 0; j < dim; j++) {
                    tour[j] = nest[i]->trail[j];
                }
            }
        }
        tour = two_opt(tour);
        tour_cost = calculate_solution_cost(tour, dim);
        global_pheromone_update(alpha);
    }
    for ( i = 0; i < swarm_count; i++) {
        free(nest[i]->trail);
        free(nest[i]);
    }
}

int main(int argc, char const *argv[]) {
    if (argc < 2) {
        printf("%s\n", "please provide a filename");
        exit(EXIT_FAILURE);
    }else{
        FILE * fp = fopen(argv[1], "r");
        if (!fp) {
            printf("%s\n", "invalid filename");
            exit(EXIT_FAILURE);
        }
        read_words(fp);
        nodes = (double *) malloc (2*dim * sizeof(double *));
        read_coords(fp);
        distances = (int *) malloc (dim*dim * sizeof(int *));
        calculate_distance();
        pheromone = (double *) malloc (dim*dim * sizeof(double *));
        eta_values = (double *) malloc (dim*dim * sizeof(double *));
        tour = (int *) malloc (dim * sizeof(int *));
        int sol[dim];
        for (size_t i = 0; i < dim; i++) {
            sol[i] = i;
        }
        if (argc >= 3) {
            int seed  = atoi(argv[2]);
            srand(seed);
        }else{
            srand(time(NULL));
        }
        for (size_t j = dim -1; j > 0 ; j--) {
            int r = rand()%dim;
            swap_places(sol, j, r);
        }
        nearest_neighbor(sol, tour);
        tour_cost = calculate_solution_cost(tour, dim);
        t = 1.0/(tour_cost*dim);
        for (size_t m = 0; m < dim; m++) {
            for (size_t k = 0; k < dim; k++) {
                *(pheromone +(m*dim + k)) =  t;
                *(pheromone +(k*dim + m))= t;
                *(eta_values +(m*dim + k)) = t * (1.0/ *(distances + (m*dim + k)))* 1.0/ *(distances + (m*dim + k));
                *(eta_values +(k*dim + m)) = *(eta_values +(m*dim + k));
            }
        }
        acs(10, tour);
        for (int i = 0; i < dim; i++) {
            printf("%d%s", tour[i]+1, " " );
        }
        putchar('\n');
        // write the results to the output file
        FILE * fw;
        if (argc >= 4) {
            fw = fopen(argv[3], "w");
        }else{
            fw = fopen("results.txt", "w");
        }
        if (fw) {
            for (size_t k = 0; k < dim; k++) {
                fprintf(fw, "%d\n", tour[k]+1);
            }
            fclose(fw);
        }
        printf("\nTour length:%d\n\n",tour_cost);
        fclose(fp);
        free(nodes);
        free(distances);
        free(pheromone);
        free(eta_values);
        free(tour);
    }

    return 0;
}
