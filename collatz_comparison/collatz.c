#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>

#define int_t long long unsigned int

/**
 * Information about a Collatz sequence:
 *
 *   start_value: initial value of the sequence
 *   current_iterate: current iterate
 *   iterations: number of iterations required to reach
 *       the current iterate from the initial value
 */
typedef struct sequence {
    int_t start_value;
    int_t current_iterate;
    int_t iterations;
} Sequence;

/**
 * Information about a search for a long Collatz sequence
 *
 *    start_value: initial value of the longest sequence found so far
 *    iterations: length of the longest sequence found so far
 */
typedef struct longest {
    int_t start_value;
    int_t iterations;
} LongestSequence;

/**
 * Advance a sequence by one iteration and increment the iteration count
 */
void iterate(Sequence *sequence) {

    if (sequence->current_iterate % 2 == 0) {
        sequence->current_iterate = sequence->current_iterate / 2;
    } else {
        sequence->current_iterate = 3*sequence->current_iterate + 1;
    }

    sequence->iterations += 1;

}

/**
 * Test if a sequence has reached one. If it has, update the search results
 * (stored in longest_sequence) if this sequence is the longest found so 
 * far.
 */
bool is_at_one(Sequence *sequence, LongestSequence *longest_sequence) {

    if (sequence->current_iterate == 1) {
        if (sequence->iterations > longest_sequence->iterations) {
            longest_sequence->start_value = sequence->start_value;
            longest_sequence->iterations = sequence->iterations;
        }
        return true;
    } else {
        return false;
    }
}

/**
 * Search for the longest sequence with an initial value at or below
 * max_start_value
 */
int main(int argc, char *argv[]) {

    int_t max_start_value = atoi(argv[1]);
    LongestSequence longest = {1, 0};

    for (int_t start_value = 2; 
         start_value <= max_start_value;
         start_value ++ ) {
        Sequence sequence = { start_value, start_value, 0 };
        while (!is_at_one(&sequence, &longest)) {
            iterate(&sequence);
        }
    }

    printf("Longest sequence: %llu iterations starting from %llu\n",
            longest.iterations, longest.start_value);
}
