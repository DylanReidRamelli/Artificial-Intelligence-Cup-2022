void read_words(FILE *f) {
  char string[100];
  fscanf(f, "%s", string);
  while (strcmp(string, "NODE_COORD_SECTION") != 0) {
    if (strcmp(string, "DIMENSION:") == 0 || strcmp(string, "DIMENSION") == 0) {
      fscanf(f, "%s", string);
      if (strcmp(string, ":") == 0) {
        fscanf(f, "%s", string);
      }
      dim = atoi(string);
    }
    if (strcmp(string, "BEST_KNOWN:") == 0 ||
        strcmp(string, "BEST_KNOWN") == 0) {
      fscanf(f, "%s", string);
      if (strcmp(string, ":") == 0) {
        fscanf(f, "%s", string);
      }
      best = atoi(string);
    }
    fscanf(f, "%s", string);
  }
}
void read_coords(FILE *f, double nodes[2][dim]) {
  int index[1];
  double x[1];
  double y[1];
  for (size_t i = 0; i < dim; i++) {
    fscanf(f, "%d", index);
    fscanf(f, "%lf", x);
    fscanf(f, "%lf", y);
    nodes[0][index[0] - 1] = x[0];
    nodes[1][index[0] - 1] = y[0];
  }
}
