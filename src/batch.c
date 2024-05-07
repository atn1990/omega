#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

static unsigned int k_bound = 8;
static unsigned int t_bound = 1000000;
static unsigned int n_start = 0;

void unix_error(char* msg) {
  fprintf(stderr, "%s: %s\n", msg, strerror(errno));
  exit(EXIT_FAILURE);
}

int run(char* n, char* k) {
  int status;
  pid_t pid;

  pid = fork();
  if (pid < 0) {
    unix_error("fork() error");
  }

  if (pid == 0) {
    ualarm(t_bound, 0);
    if (execl(
          "./driver", "./driver", "--dfa", n, "--k", k, "-v", "0", NULL) < 0) {
      unix_error("exec() error");
    }
  }

  pid = waitpid(pid, &status, 0);
  if (pid < 0) {
    unix_error("waitpid() error");
  }

  if (WIFEXITED(status)) {
    return 0;
  } else {
    return 1;
  }
}

int main(int argc, char *argv[]) {
  char c, n[4], k[3];

  while ((c = getopt(argc, argv, "k:t:n:")) != -1) {
    switch (c) {
      case 'k':
        k_bound = strtol(optarg, NULL, 10);
        break;

      case 't':
        t_bound = strtol(optarg, NULL, 10);
        break;

      case 'n':
        n_start = strtol(optarg, NULL, 10);
        break;
    }
  }

  for (size_t i = n_start; i < 256; i++) {
    sprintf(n, "%zd", i);
    printf("%3zd  ", i);
    fflush(stdout);
    for (size_t j = 1; j < k_bound; j++) {
      sprintf(k, "%zd", j);
      if (run(n, k) > 0) {
        break;
      }
    }
    printf("\n");
  }

  return EXIT_SUCCESS;
}
