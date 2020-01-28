public class Collatz {

    /**
     * Information about a Collatz sequence.
     *
     *   start_value: initial value of the sequence
     *   current_iterate: current iterate
     *   iterations: number of iterations required to reach
     *       the current iterate from the initial value
     */
    public static class Sequence {
        long start_value;
        long current_iterate;
        long iterations;
        public Sequence(long start_value) {
            this.start_value = start_value;
            this.current_iterate = start_value;
            this.iterations = iterations;
        }
    }

    /**
     * Information about a search for a long Collatz sequence.
     *
     *    start_value: initial value of the longest sequence found so far
     *    iterations: length of the longest sequence found so far
     */
    public static class LongestSequence {
        long start_value;
        long iterations;
        public LongestSequence() {
            this.start_value = 1;
            this.iterations = 0;
        }
    }

    /**
     * Advance a sequence by one iteration 
     * and increment the iteration count
     */
    public static void iterate(Sequence sequence) {
    
        if (sequence.current_iterate % 2 == 0) {
            sequence.current_iterate = sequence.current_iterate / 2;
        } else {
            sequence.current_iterate = 3*sequence.current_iterate + 1;
        }
        sequence.iterations += 1;
    }
    
    /**
     * Test if a sequence has reached one. 
     * If it has, update the search results
     * (stored in longest_sequence) 
     * if this sequence is the longest found so far.
     */
    public static boolean is_at_one(
            Sequence sequence, LongestSequence longest_sequence) {
    
        if (sequence.current_iterate == 1) {
            if (sequence.iterations > longest_sequence.iterations) {
                longest_sequence.start_value = sequence.start_value;
                longest_sequence.iterations = sequence.iterations;
            }
            return true;
        } else {
            return false;
        }
    }
    
    /**
     * Search for the longest sequence with an initial value at or below
     * max_start_value (given as first argument)
     */
    public static void main(String[] args) {
    
        long max_start_value = Long.parseLong(args[0]);
        LongestSequence longest = new LongestSequence();
    
        for (long start_value = 2; 
             start_value <= max_start_value;
             start_value++ ) {
            Sequence sequence = new Sequence(start_value);
            while (!is_at_one(sequence, longest)) {
                iterate(sequence);
            }
        }
    
        System.out.println("Longest sequence: " +
                longest.iterations + " iterations starting from " +
                longest.start_value);
    }
}
