#!/usr/bin/env python3
# Roger Volden

from libc.stdlib cimport calloc, free

def conk(str splint, str read, int penalty):
    '''
    Inputs:
        splint (python string)
        read (python string)
        penalty (int)

    Outputs:
        sums (python list)

    Aligns the splint to the read. The diagonal
    scores are then summed up and returned.
    '''

    # define all of my c variables
    cdef int n = len(splint) + 1
    cdef int m = len(read) + 1
    cdef int i, j, aScore, cMax
    cdef int left, above

    # allocate memory for the alignment matrix
    cdef int* alnMat = <int*>calloc(n * m, sizeof(int))
    if not alnMat:
        raise MemoryError()

    for i in range(1, n):
        for j in range(1, m):
            aScore = baseScore(ord(splint[i-1]), ord(read[j-1]))
            cMax = alnMat[(m*(i-1)) + (j-1)] + aScore
            left = alnMat[(m*(i-1)) + j] - penalty 
            if left > cMax:
                cMax = left
            above = alnMat[(m*i) + (j-1)] - penalty
            if above > cMax:
                cMax = above
            if 0 > cMax:
                cMax = 0
            alnMat[(m*i) + j] = cMax

    # allocate the memory for the summed score array
    cdef int* sums = <int*>calloc(m + 1, sizeof(int))
    if not sums:
        raise MemoryError()

    # sum up the alignment scores
    for i in range(0, n):
        for j in range(i, m):
            sums[j-i] += alnMat[(m*i) + j]

    free(alnMat)
    sumsList = [x for x in sums[:m]]

    return sumsList

cdef baseScore(char a, char b):
    '''
    Inputs:
        a (c char)
        b (c char)

    Outputs:
        int

    Compares the two characters given and returns the
    corresponding score. Ns are counted as matches.
    '''
    if a == b or a == 'N' or b == 'N':
        return 5
    else:
        return -4
