#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>


static int dim;
static int best;

typedef struct {
    int start_node;
    int index;
    int * trail;
    int * visited;
    // int trail[130];
    double cost;
} ANT;

// function declarations
int main(int argc, char const *argv[]);
FILE * read_words(FILE * f);
FILE * read_coords(FILE * f,double nodes[2][dim]);
void calculate_distance(double nodes[2][dim], int distances[dim][dim]);
int calculate_solution_cost(int * sol, int size,int distances[dim][dim]);
void swap_places(int * sol, int i, int j);
void nearest_neighbor(int * sol,int distances[dim][dim],int  * new_sol);
int check_valid(int*sol);
void acs( ANT** nest, int swarm_count, int distances[dim][dim], double pheromone [dim][dim], int * shortest_path);
void global_pheromone_update(ANT* best,double pheromone [dim][dim],double alpha);
void clean_ants(ANT ** nest, int swarm_count);
int exploit(ANT * a, int pos, int * visited, double pheromone [dim][dim],int distances[dim][dim]);
int explore(ANT * a, int pos, double pheromone [dim][dim],int distances[dim][dim]);

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
        if (strcmp(string,"BEST_KNOWN:")==0 || strcmp(string,"BEST_KNOWN")==0 ) {
            fscanf(f, "%s",string);
            if (strcmp(string,":")==0 ) {
                fscanf(f, "%s",string);
            }
            best = atoi(string);
        }
        fscanf(f, "%s", string);
    }
    return f;
}
FILE * read_coords(FILE * f,double nodes[2][dim]){
    int index[1];
    double x[1];
    double y[1];
    for (size_t i = 0; i < dim; i++) {
        fscanf(f, "%d", index);
        fscanf(f, "%lf", x);
        fscanf(f, "%lf", y);
        nodes[0][index[0]-1] = x[0];
        nodes[1][index[0]-1] = y[0];
    }
    return f;
}

void calculate_distance(double nodes[2][dim], int distances[dim][dim]){
    for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < dim; j++) {
            if (j == i) {
                distances[i][j] = 9999;
                distances[j][i] = 9999;
            }else{
                double x_sq = pow(nodes[0][i] - nodes[0][j], 2);
                double y_sq = pow(nodes[1][i] - nodes[1][j], 2);
                distances[i][j] = sqrt(x_sq + y_sq);
                distances[j][i] = sqrt(x_sq + y_sq);
            }
        }
    }
}
int calculate_solution_cost(int * sol, int size, int distances[dim][dim]){
    int sum = 0;
    for (size_t i = 0; i < dim - 1; i++) {
        sum+= distances[sol[i]][sol[i+1]];
    }
    sum += distances[sol[dim-1]][sol[0]];
    return sum;
}
void swap_places(int * sol, int i, int j){
    int temp = sol[j];
    sol[j] = sol[i];
    sol[i] = temp;
}
int find_index(int * sol, int val){
    for (size_t i = 0; i < dim; i++) {
        if (sol[i]==val) {
            return i;
        }
    }
    return -1;
}

void nearest_neighbor(int * sol,int distances[dim][dim],int  * new_sol){
    int start_node = rand()%dim; // index of randomly selected node
    // int start_node = find_index(sol,0);
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
        value = distances[0][0];
        for (size_t k = 0; k < dim; k++) {
            if (distances[index][k] <= value && visited[k] == 0) {
                nearest = k;
                value = distances[index][nearest];
            }
        }
        // printf("%d\n", index);
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
            // printf("%d\t%d\t%d\n",sol[0],k,visited[k] );
            // printf("%s\n", "Invalid Tour");
            return 1;
        }
    }
    return 0;
}
void global_pheromone_update(ANT* best,double pheromone [dim][dim],double alpha){
    int node1,node2;
    double phi, delta;
    delta = 1.0/best->cost;
    for (size_t i = 0; i < dim-1; i++) {
        node1 = best->trail[i];
        node2 = best->trail[i+1];
        phi = pheromone[node1][node2];
        pheromone[node1][node2] = ((1 - alpha)*phi)+(alpha*delta);
        pheromone[node2][node1] = pheromone[node1][node2];
    }
    node1 = best->trail[0];
    node2 = best->trail[dim-1];
    phi = pheromone[node1][node2];
    pheromone[node1][node2] = ((1 - alpha)*phi)+(alpha*delta);
    pheromone[node2][node1] = pheromone[node1][node2];


}

// void acs(int swarm_count, int distances[dim][dim], double pheromone[dim][dim], int * shortest_path){
//     ANT* nest[swarm_count];
//     double S[dim];
//     ANT * a1;
//     ANT best;
//     int i,j,k, node,pos,q,valid;
//     double tau,eta,value,beta,sigma,max,min,alpha,t,d, sum, rnd;
//     alpha = 0.1;
//     t = 0.1;
//     int test = 1000;
//     min = 99999999;
//     for (i = 0; i < swarm_count; i++) {
//         // initialise ant
//         nest[i] = (ANT*)malloc(sizeof(ANT));
//         a1= nest[i];
//         a1->trail = (int*)malloc(dim*sizeof(int));
//         a1->visited = (int*)malloc(dim*sizeof(int));
//     }
//     // loop until condition is met
//     while (test != 0) {
//         test--;
//         for (i = 0; i < swarm_count; i++) {
//             // initialise ant
//             a1= nest[i];
//             node = rand()%dim;
//             a1->start_node = node;
//             // a1->trail[0] = a1->start_node;
//             a1->index = 0;
//             a1->cost = 9999;
//         }
//         for ( i = 0; i < swarm_count; i++) {
//             a1 = nest[i];
//             pos = a1->start_node;
//             sigma = 0;
//             for ( j = 0; j < dim; j++) {
//                 a1->visited[j] = 0;
//                 a1->trail[j] = 0;
//             }
//             a1->trail[a1->index] = a1->start_node;
//             a1->visited[a1->start_node] = 1;
//             a1->index++;
//             // Exploitation
//             for (j = 0; j < dim; j++) {
//                 q = rand()%100;
//                 max = -1;
//                 if (q <= 94) { // 95% of the time
//                     // find next node
//                     for ( k = 0; k < dim; k++) {
//                         tau = pheromone[pos][k];
//                         beta = 2.0;
//                         eta = pow((1.0/(double)distances[pos][k]), beta);
//                         value = tau * eta;
//                         if (value > max && a1->visited[k] ==0 && pos != k) {
//                             max = value;
//                             node = k;
//                         }
//                     }
//                 }
//                 // Exploration
//                 else{
//                     sigma = 0;
//                     sum = 0; // should be 1
//                     // calculate probability
//                     for (k = 0; k < dim; k++) {
//                         if (a1->visited[k] == 0) {
//                             tau = pheromone[pos][k];
//                             // beta = 1.0/95.0;
//                             beta = 2.0;
//                             eta = pow((1.0/(double)distances[pos][k]), beta);
//                             // printf("%lf\n", eta);
//                             value = tau * eta;
//                             S[k] = value;
//                             sigma += value;
//                         }else{
//                             S[k] = 0;
//                         }
//                     }
//                     for (k = 0; k < dim; k++) {
//                         S[k] = S[k]/sigma;
//                         sum += S[k];
//                     }
//
//                     // Pick a random edge
//                     rnd = (double)rand()/((double)RAND_MAX);
//                     d = 0;
//                     for (k = 0; k < dim; k++)
//                        {
//                                d += S[k];
//                                if (rnd < d && S[k] != 0)
//                                {
//                                    node = k;
//                                    k = dim;
//                                }
//                        }
//                 }
//                 if (a1->visited[node] == 0) {
//                     a1->trail[a1->index] = node;
//                     tau = pheromone[pos][node];
//                     // local pheromone update (symmetric)
//                     pheromone[pos][node] = (1 - alpha)*tau + (alpha * t);
//                     pheromone[node][pos] = (1 - alpha)*tau + (alpha * t);
//                     a1->index++;
//                     a1->visited[node] = 1;
//                     pos = node;
//                 }
//             }
//             valid = check_valid(a1->trail);
//             if (valid == 0) {
//                 a1->cost = calculate_solution_cost(a1->trail, dim, distances);
//             }
//         }
//
//         for (i = 0; i < swarm_count; i++) {
//             if (nest[i]->cost < min) {
//                 min = nest[i]->cost;
//                 best = *nest[i];
//                 for (size_t j = 0; j < dim; j++) {
//                     shortest_path[j] = nest[i]->trail[j];
//                 }
//             }
//         }
//     }
//     global_pheromone_update(&best,pheromone,alpha);
//     for ( i = 0; i < swarm_count; i++) {
//         free(nest[i]->trail);
//         free(nest[i]->visited);
//         free(nest[i]);
//     }
// }

void acs( ANT** nest, int swarm_count, int distances[dim][dim], double pheromone [dim][dim], int * shortest_path){
    // declare variables
    ANT best;
    int i, j; // loop indices
    int node,q,valid, pos;
    double tau,min,alpha,t;
    alpha = 0.1;
    t = 0.1; // evaporation
    min = 99999999;
    ANT * a ;
    int flag = 0; // down
    /*
    * Initialise ants.
    * Put m ants on randomly selected cities
    */
    clean_ants(nest, swarm_count);
    // make tour for each ant
    for( i = 0; i < swarm_count; i++){ // for the current ant
        flag = 0;
        a = nest[i];
        pos = a->start_node;
        printf("%d\n", pos);
        for (j = 0; j < dim; j++) { // pick cities
            q = rand()%100;
            if (q <= 94){
                node = exploit(a,pos,a->visited,pheromone,distances);
                // printf("%d\t", node);

            }else{
                node = explore(a,pos,pheromone,distances);

            }
            // add the node
            if (a->visited[node] == 0 && node != -1 && a->index < dim) {
                a->trail[a->index] = node;
                tau = pheromone[pos][node];
                // local pheromone update (symmetric)
                pheromone[pos][node] = (1 - alpha)*tau + (alpha * t);
                pheromone[node][pos] = (1 - alpha)*tau + (alpha * t);
                a->index++;
                a->visited[node] = 1;
                pos = node;
                node = -1;
            }else{
                flag = 1;
            }
        }
        // update ant details
        valid = check_valid(a->trail);
        if (valid == 0 && flag != 1) {
            a->cost = calculate_solution_cost(a->trail, dim, distances);
        }else{
            a->cost = 99999999;
        }
    }
    for (i = 0; i < swarm_count; i++) {
        if (nest[i]->cost < min) {
            min = nest[i]->cost;
            // best = *nest[i];
            for (size_t j = 0; j < dim; j++) {
                best.trail[j] = nest[i]->trail[j];
                shortest_path[j] = nest[i]->trail[j];
            }
        }
    }
    global_pheromone_update(&best,pheromone,alpha);
}

int exploit(ANT * a, int pos, int * visited, double pheromone [dim][dim],int distances[dim][dim]){
    double tau, eta, value;
    // printf("%d\n", a->index);
    // int pos = a->trail[a->index-1];
    double beta = 2.0;
    double max = -1;
    int node = -1;
    for (int k = 0; k < dim; k++) {
        tau = pheromone[pos][k];
        eta = pow((1.0/(double)distances[pos][k]), beta);
        value = tau * eta;
        if (value > max && visited[k] ==0 && pos != k) {
            max = value;
            node = k;
        }
    }
    return node;
}
int explore(ANT * a, int pos, double pheromone [dim][dim],int distances[dim][dim]){
    double tau, eta, value, rnd,d;
    double sigma = 0;
    double S[dim];
    double beta = 2.0;
    int node = -1;
    // int pos = a->trail[a->index-1];
    // calculate probability
    for (int i = 0; i < dim; i++) {
        if (a->visited[i] == 0) {
            tau = pheromone[pos][i];
            eta = pow((1.0/(double)distances[pos][i]), beta);
            value = tau * eta;
            S[i] = value;
            sigma += value;
        }else{
            S[i] = 0;
        }
    }
    for (int i = 0; i < dim; i++) {
        S[i] = S[i]/sigma;
    }

    // Pick a random edge
    rnd = (double)rand()/((double)RAND_MAX);
    // printf("%lf\n", rnd);
    d = 0;
    for (int i = 0; i < dim; i++){
       d += S[i];
       if (rnd < d && S[i] != 0){
        //    printf("%lf\t %lf\t %d\n", rnd,d,i);
           node = i;
           return node;
       }
   }
   return node;
}
void clean_ants(ANT ** nest, int swarm_count){
    ANT * a1;
    for (int i = 0; i < swarm_count; i++) {
        a1= nest[i];
        a1->start_node =  rand()%dim;
        a1->index = 0;
        a1->cost = 9999;
        for (int  j = 0; j < dim; j++) {
            a1->visited[j] = 0;
            a1->trail[j] = 0;
        }
        a1->trail[a1->index] = a1->start_node;
        a1->visited[a1->start_node] = 1;
        a1->index++;
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
        double nodes[2][dim];
        read_coords(fp,nodes);
        int distances[dim][dim];
        calculate_distance(nodes,distances);
        double pheromone[dim][dim];
        for (size_t m = 0; m < dim; m++) {
            for (size_t k = 0; k < dim; k++) {
                pheromone[m][k] = 1;
                pheromone[k][m] = 1;
            }
        }
        int sol[dim];
        for (size_t i = 0; i < dim; i++) {
            sol[i] = i;
        }
        srand(time(NULL));
        for (size_t j = dim -1; j > 0 ; j--) {
            int r = rand()%dim;
            swap_places(sol, j, r);
        }
        int new_sol[dim];
        printf("%s\t%d\n", "Best: ",best);

        /*
        * Make ants
        */
        int swarm_count = 3;
        ANT * nest[swarm_count];
        ANT * a1;
        int test = 1;
        for (int i = 0; i < swarm_count; i++) {
            // initialise ant
            nest[i] = (ANT*)malloc(sizeof(ANT));
            a1= nest[i];
            a1->trail = (int*)malloc(dim*sizeof(int));
            a1->visited = (int*)malloc(dim*sizeof(int));
        }
        while (test !=0 ) {
            acs(nest,swarm_count, distances, pheromone, new_sol);
        }


        int sol_cost = calculate_solution_cost(new_sol, sizeof(new_sol)/sizeof(int), distances);
        printf("%s\t%d\n","ACS: " ,sol_cost);

        for (int i = 0; i < swarm_count; i++) {
            free(nest[i]->trail);
            free(nest[i]->visited);
            free(nest[i]);
        }
        // nearest_neighbor(sol, distances,new_sol);
        // sol_cost = calculate_solution_cost(new_sol, sizeof(new_sol)/sizeof(int), distances);
        // printf("%s\t%d\n","NN: " ,sol_cost);
        printf("%s\t%f\n","Err: " ,((double)(sol_cost - best)/(double)best)*100);
        int valid;
        valid = check_valid(new_sol);
        if (valid == 0) {
            FILE * fw = fopen("near.txt", "w");
            if (fw) {
                for (size_t k = 0; k < dim; k++) {
                    // printf("%d\t", new_sol[k]);
                    fprintf(fw, "%d\n", new_sol[k]+1);
                }
                fclose(fw);
            }
        }
        fclose(fp);
    }
    return 0;
}
