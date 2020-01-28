"""
Information about a Collatz sequence

    start_value: initial value of the sequence
    current_iterate: current iterate
    iterations: number of iterations required to reach
        the current iterate from the initial value
"""
class Sequence:
    def __init__(self, start_value):
        self.start_value = start_value
        self.current_iterate = start_value
        self.iterations = 0

"""
Information about a search for a long Collatz sequence

    start_value: initial value of the longest sequence found so far
    iterations: length of the longest sequence found so far
"""
class LongestSequence:
    def __init__(self):
        self.start_value = 1
        self.iterations = 0

"""
Advance a sequence by one iteration and increment the iteration count
"""
def iterate(sequence):

    if sequence.current_iterate % 2 == 0:
        sequence.current_iterate = sequence.current_iterate/2
    else:
        sequence.current_iterate = 3*sequence.current_iterate + 1

    sequence.iterations += 1

"""
Test if a sequence has reached one. If it has, update the search results
(stored in longest_sequence) if this sequence is the longest found so 
far.
"""
def is_at_one(sequence, longest_sequence):

    if sequence.current_iterate == 1:
        if sequence.iterations > longest_sequence.iterations:
            longest_sequence.start_value = sequence.start_value
            longest_sequence.iterations = sequence.iterations
        return True
    else:
        return False

"""
Search for the longest sequence with an initial value at or below
max_start_value
"""
def find_longest_sequence(max_start_value):

    longest = LongestSequence()

    for start_value in range(2, max_start_value+1):
        sequence = Sequence(start_value)
        while not is_at_one(sequence, longest):
            iterate(sequence)
    
    print("Longest sequence: %d iterations starting from %d\n" % ( 
            longest.iterations, longest.start_value))
