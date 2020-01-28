using Printf

"""
Information about a Collatz sequence

    start_value: initial value of the sequence
    current_iterate: current iterate
    iterations: number of iterations required to reach
        the current iterate from the initial value
"""
mutable struct Sequence
    start_value
    current_iterate
    iterations
end

"""
Information about a search for a long Collatz sequence

    start_value: initial value of the longest sequence found so far
    iterations: length of the longest sequence found so far
"""
mutable struct LongestSequence
    start_value
    iterations
end

"""
Advance a sequence by one iteration and increment the iteration count
"""
function iterate!(sequence)

    if sequence.current_iterate % 2 == 0
        sequence.current_iterate = sequence.current_iterate/2
    else
        sequence.current_iterate = 3*sequence.current_iterate + 1
    end

    sequence.iterations += 1

end

"""
Test if a sequence has reached one. If it has, update the search results
(stored in longest_sequence) if this sequence is the longest found so 
far.
"""
function is_at_one!(sequence, longest_sequence)

    if sequence.current_iterate == 1
        if sequence.iterations > longest_sequence.iterations
            longest_sequence.start_value = sequence.start_value
            longest_sequence.iterations = sequence.iterations
        end
        return true
    else
        return false
    end

end

"""
Search for the longest sequence with an initial value at or below
max_start_value
"""
function find_longest_sequence(max_start_value)

    longest = LongestSequence(1,0)

    for start_value = 2:max_start_value
        sequence = Sequence(start_value, start_value, 0)
        while !is_at_one!(sequence, longest)
            iterate!(sequence)
        end
    end
    
    @printf("Longest sequence: %d iterations starting from %d\n", 
            longest.iterations, longest.start_value)

end
