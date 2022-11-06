import csv
import sys
import numpy as np
import os


# Biological Sequence Project #1

# Axel Negron Vega  CIIC4025-016

#-----Global------#
gap_penalty = -2
#-----Global------#


def sequence_matrix(in_file):
    """Goes through csv file and obtains labels from first row and sequences from second row

        :param in_file: csv file comma delimited

    """

    sequences = list()
    with open(in_file) as csv_file:
        reader = csv.reader(csv_file, delimiter="\t")
        for i, line in enumerate(reader):
            if i > 0:
                try:
                    # Remove empty cells in csv file using split
                    dummy = line[0].split(',')
                    x = 0
                    while x < len(dummy):
                        if dummy[x] == '':
                            dummy.pop(x)
                            x -= 1
                        x += 1
                    sequences += dummy
                except:
                    pass

        # Removing spaces in sequence
        sequences = [item.replace(' ', '') for item in sequences]

        ##Creating sequence matrix##

        try:
            sequence_arrs = list()
            for _ in range(int(len(sequences)/2)):
                sequence_arr = np.zeros(
                    (len(sequences[0])+2, len(sequences[1])+2), dtype='U25')  # dtype='|S1')
                temp1 = np.array(list('##'+sequences[1]))
                temp2 = np.array(list(sequences[0]))

                # Set first row in sequence matrix
                sequence_arr[0] = temp1

                # Set first row in sequence matrix
                sequence_arr[1][0] = '#'
                for x in range(2, sequence_arr.shape[0]):
                    sequence_arr[x][0] = temp2[x-2]
                sequence_arrs.append(sequence_arr)
                sequences.pop(0)
                sequences.pop(0)

        except:

            exit()

        ###Matrix form, uncomment to view###
        # y = sequence_arr.view(dtype='U25', type=np.matrix)
        # print(y)
        # return y

        return sequence_arrs


def set_score(matrix, gap_penalty: int):
    """Function construct table values.
        :param matrix: 2-D numpy array or type U25
        :param gap_penalty: int

    """
    ##Setting top and left side##

    for num in range(1, matrix.shape[1]):
        matrix[1][num] = str(gap_penalty*(num-1))

    for num in range(2, matrix.shape[0]):
        matrix[num][1] = str(gap_penalty*(num-1))

    # Sets value of cell based on the values around it
    for row in range(2, matrix.shape[0]):
        for col in range(2, matrix.shape[1]):
            # Checks if protein type is equal
            A, B = (matrix[row][0], matrix[0][col])
            if A == B:
                scoring_matrix = 1
            else:
                scoring_matrix = -1

            # Defines values for diagonal, side and top cell
            diag_val = int(matrix[row-1][col-1])+scoring_matrix
            side_val = int(matrix[row][col-1])+gap_penalty
            top_val = int(matrix[row-1][col])+gap_penalty

            # Grabs max value and sets cell value to that max value
            value = max(side_val, top_val, diag_val)
            matrix[row][col] = str(value)

    return


def backtrack(matrix, total: int, row: int, col: int, sq: str, sq2: str):
    """Recursive method that backtracks through the matrix to grab total points and dna protein sequence
        :param matrix: 2-D numpy array or type U25
        :param total: Likeliness score for given strands
        :param row: Row of matrix
        param col: Column of matrix
        param sq: Sequence 1
        param sq2: Sequence 2

    """
    if (col == 1 and row == 1):
        return (sq, total, sq2)

    A, B = (matrix[row][0], matrix[0][col])

    #If row==1 we can only move left
    if (row == 1):
        return backtrack(matrix, total-2, row, col-1, sq+'-', sq2+B)
    #If col==1 we can only move up#
    elif (col == 1):
        return backtrack(matrix, total-2, row-1, col, sq+A, sq2+'-')
    equal = False
    #Getting left,top and diagonal value
    left_val = int(matrix[row][col-1])+gap_penalty
    top_val = int(matrix[row-1][col])+gap_penalty
    diag_val = int(matrix[row-1][col-1])

    #Updating diagonal value based on if A=B
    if A == B:
        diag_val = int(matrix[row-1][col-1])+1
        equal = True
    else:
        diag_val = int(matrix[row-1][col-1])-1

    #Recursive calls, first check if we move diagonally,then check if we move left, else we move top.
    if (diag_val > left_val and diag_val > top_val):
        if (equal):
            return backtrack(matrix, total+1, row-1, col-1, sq+A, sq2+B)
        else:
            return backtrack(matrix, total-1, row-1, col-1, sq+A, sq2+B)
    else:
        #Left traversal has priority
        if (left_val >= diag_val and left_val >= top_val):
            return backtrack(matrix, total-2, row, col-1, sq+'-', sq2+B)
        else:
            return backtrack(matrix, total-2, row-1, col, sq+A, sq2+'-')


def get_lastpos(matrixes):
    """Function to find last positions in all matrixes
        :param matrixes: list of matrixes
    """
    pos_list = list()
    for x in range(len(matrixes)):
        pos_list.append(matrixes[x].shape)

    return pos_list


def disp_results(results, matrixes):
    """Function to parse results after backtracking, positions the sequences correctly, gets score and displays output
        :param matrixes: list of matrixes
        :param results: list of sequences and scores
    """
    sequences = list()
    scores = list()
    for x in range(len(matrixes)):
        sequences.append((results[x][0][::-1] + ' '+results[x][2][::-1]))
        scores.append(str(results[x][1]))

    for x in range(len(matrixes)):
        print(sequences[x] + ' '+scores[x])


def main(debug: bool = False):
    """Main function takes debug conditional. If true a default path is used for csv. Change path to csv name to test
        :param debug:boolean
    """
    try:
        ##If in debug comment line below##
        in_file = sys.argv[1]
        pass
        if debug:

            # Change input path to debug
            in_file = os.getcwd()+"\\Sequence_file_input.csv"

        ###Creates sequence matrix###
        matrixes = sequence_matrix(in_file)
        for m in matrixes:

            #Sets score inside matrix
            set_score(m, gap_penalty)

        #Gets last positions in matrix
        pos_markers = get_lastpos(matrixes)

        #Get results and display
        results = list()
        for x in range(len(matrixes)):
            sequence, total_score, sequence2 = backtrack(
                matrixes[x], 0, pos_markers[x][0]-1, pos_markers[x][1]-1, '', '')
            results.append((sequence, total_score, sequence2))

        ###Reverse Sequence/Display Results###
        disp_results(results, matrixes)
    except:
        ##If no file was provided##
        pass


main(debug=False)
