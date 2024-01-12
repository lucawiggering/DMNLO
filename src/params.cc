#include "DMNLO.h"



using namespace std;

map<string, string> readInputFile(string filename) {
  map<string, string> config;
  errno = 0;
  FILE *input_file;

  if (filename == "-") {
    input_file = stdin;
  } else {
    input_file = fopen(filename.c_str(), "r");
  }
  if (input_file == NULL) {
    fprintf(stderr, "error: could not open file '%s'.\n", filename.c_str());
    char buff[FILENAME_MAX];
    getcwd(buff, FILENAME_MAX);
    printf("Current working dir: %s\n", buff);
    printf("Error %d \n", errno);
    exit(1);
  }

  char *variable = (char *)malloc(128 + 1);
  char *value = (char *)malloc(128 + 1);
  ;
  for (;;) {
    char *input;
    size_t n;
    int chars_read;
    int elements_read;

    input = NULL;
    n = 0;
    chars_read = getline(&input, &n, input_file);
    if (chars_read == -1) {
      break;
    }
    if (chars_read < 2) {
      continue;
    }
    elements_read = sscanf(input, "%128[a-zA-Z0-9_-] = %128s\n", variable, value);
    if (elements_read == 0) {
      // Blank line.
      continue;
    } else if (elements_read == 1 && variable[0] == '#') {
      // Comment.
      continue;
    } else if (elements_read == 1) {
      fprintf(stderr,
              "error: while reading line '%s' of the input file at line %li. (param_card.dat as "
              "input.in?)\n",
              variable, n);
      exit(1);
    } else if (variable[0] == '#') {
      // Comment.
      continue;
    }
    config[variable] = value;
  }
  free(variable);
  free(value);
  fclose(input_file);
  return config;
}