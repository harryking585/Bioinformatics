import argparse
#!/usr/bin/env python3


# read_fna(file)
# purpose: reads the input file and separates the headers from the sequences
# @return list
def read_fna(input_file):
    line_list = input_file.readlines()
    header_list = []
    seq_list = []

    seq_str = ''
    for x in line_list:
        if x[0] == '>':
            header_list.append(x[:-1])
            seq_str = ''
        else:
            seq_str += x[:-1]
            try:
                if(line_list[(line_list.index(x)+1)][0] == '>'):
                    seq_list.append(seq_str)
                    seq_str = ''
            except:
                continue

    seq_list.append(seq_str)
    input_file.close()
    return [header_list, seq_list]


# read_subtiution_mtx(file)
# purpose: reads the score matrix file in order to return a dictionary that associates each base pairing with their appropriate score
# @return dictionary
def read_substitution_mtx(mtx_file):
    mtx = []
    line_list = mtx_file.readlines()
    index = 0
    temp_list = []

    for x in line_list:
        if index == 0:
            split_list = x.split(' ')
            temp_list.append(' ')

            for x in split_list:
                if len(x) != 0:
                    temp_list.append(x)

            split_list = temp_list.copy()
            mtx.append(split_list)
            index += 1
            temp_list = []
        else:
            split_list = x.split(' ')

            for x in split_list:
                if len(x) != 0:
                    temp_list.append(x)

            split_list = temp_list.copy()
            mtx.append(split_list)
            temp_list = []

    sub_dict = {}

    for i in range(1, len(mtx[0])):
        bases = mtx[0][i]
        sub_dict[bases] = {}

        for j in range(1, len(mtx)):
            base2 = mtx[j][0]
            score = mtx[j][i]
            sub_dict[bases][base2] = score

    return sub_dict


# write_fna(file, list, list, int)
# purpose: writes the aligned sequences and headers into the specified output file then closes the file
# @return nil
def write_fna(output_file, headers, aligned_sequences, score):

    for i in range(0, len(headers)):
        output_file.write(headers[i] + '; score='+str(score)+'\n')
        output_file.write(aligned_sequences[i] + '\n')

    output_file.close()


# needleman_wunsch_align(list, dictionary)
# purpose: initialize a score matrix that will be filled in accordance with the needleman-wunsch algorithm
#          after the score matrix initialization this function traverses through the score matrix backwards 
#          in order to align the sequences.
# @return list
def needleman_wunsch_align(seq_list, sub_dict):
    score_mtx = []
    seq1 = [*seq_list[0]]
    seq2 = [*seq_list[1]]
    
    # note:
    # num_row = len(seq2) + 1 so all of seq2 can be accounted for as well as the first row being seq1
    # num_cols = len(seq1) + 1 similar to num_row, so that all scores for seq1 can be accounted plus the first col being seq2
    num_rows = len(seq2) + 1 
    num_cols = len(seq1) + 1

    # Used as a 'template' for the size of the matrix before populating with sequences and scores
    score_mtx = [[0] * num_cols for _ in range(num_rows)]
    
    # score_mtx[0][0] = ' ' for myself for easier indexing when solving
    # i.e. I know that all values that populate score_mtx[0][i] are seq1 bases, and all values that populate score_mtx[i][0] are seq2 bases.
    score_mtx[0][0] = ' '
    score_mtx[0][1:] = seq1
    
    # populates the first index of each row with the seq2 bases
    for i in range(1, num_rows):
        score_mtx[i][0] = seq2[i-1]

    # initializes the rows of the score matrix that are multiples of the gap score
    gap_pen = int(sub_dict['A']['-'])
    gap_row = [' ']
    gap_score = gap_pen

    for i in range(1, len(score_mtx[0])):
        gap_row.append(gap_score)
        gap_score += gap_pen

    score_mtx.insert(1, gap_row)
    score_mtx[0].insert(1, ' ')
    score_mtx[1].insert(1, 0)

    gap_score = gap_pen

    for j in range(2, len(score_mtx)):
        score_mtx[j].insert(1, gap_score)
        gap_score += gap_pen
    
    # populates the score matrix with their proper scores in accoradance with needleman-wunsch
    for i in range(2, len(score_mtx)):

        for j in range(2, len(score_mtx[0])):
            score = int(sub_dict[score_mtx[0][j]][score_mtx[i][0]]) # the score of the two bases for each comparison
            neighbour_list = [] # len(3):  up = [0], left = [1], diag = [2]

            if i == 2 and j == 2:
                score_mtx[i][j] = score
            else:
                diag = score_mtx[i-1][j-1] + score
                up = score_mtx[i-1][j] + gap_pen
                left = score_mtx[i][j-1] + gap_pen

                neighbour_list.append(up)
                neighbour_list.append(left)
                neighbour_list.append(diag)
                
                score_mtx[i][j] = max(neighbour_list)


    total_score = score_mtx[len(score_mtx)-1][len(score_mtx[0])-1]

    # trace back begins 
    seq1_len = len(seq1)-1        # length of the first sequence
    seq2_len = len(seq2)-1        # length of the second sequence        
    i_sc = len(score_mtx) - 1       # index of row for score mtx
    j_sc = len(score_mtx[0]) - 1    # index of col for score mtx
    seq1_aligned = ''           # the aligned sequence 1
    seq2_aligned = ''           # the aligned sequence 2
    gap = '-'                   # gap character
    seq1_status = False         # bool to determine whether all of sequence 1 has been aligned
    seq2_status = False         # bool to determine wheteher all of sequence 2 has been aligned

    while(not seq1_status or not seq2_status):
        reference_value = score_mtx[i_sc][j_sc]                                  # reference for the current base value the matrix is at, changes according to traversal
        current_score = int(sub_dict[score_mtx[i_sc][0]][score_mtx[0][j_sc]])    # expected score of the two bases based upon the reference value
        diag_cell = score_mtx[i_sc-1][j_sc-1] + current_score                   
        up_cell = score_mtx[i_sc-1][j_sc]     + gap_pen
        left_cell = score_mtx[i_sc][j_sc-1]   + gap_pen

        if reference_value == diag_cell:        # if the score came from the diag cell
            if not seq1_status:
                seq1_aligned += seq1[seq1_len]
                seq1_len -= 1

            if not seq2_status:
                seq2_aligned += seq2[seq2_len]
                seq2_len -= 1

            i_sc -= 1
            j_sc -= 1
        elif reference_value == up_cell:        # if the score came from the up cell
            if not seq1_status:
                seq1_aligned += gap
                seq1_aligned += seq1[seq1_len]
                seq1_len -= 1

            if not seq2_status:
                seq2_aligned += seq2[seq2_len]
                seq2_len -= 1

            i_sc -= 1
        elif reference_value == left_cell:      # if the score came from the left cell
            if not seq1_status:
                seq1_aligned += seq1[seq1_len]
                seq1_len -= 1

            if not seq2_status:
                seq2_aligned += gap
                seq2_aligned += seq2[seq2_len]
                seq2_len -= 1

            j_sc -= 1

        # traceback terminators
        if seq1_len < 0:        # has sequence 1 finished aligning?
            seq1_status = True

        if seq2_len < 0:        # has sequence 2 finished aligning?
            seq2_status = True

    # reverse both of the aligned sequences as they were built backwards
    seq1_aligned = seq1_aligned[::-1]
    seq2_aligned = seq2_aligned[::-1]

    return [[seq1_aligned, seq2_aligned], total_score]


print('Assignment 3 :: R11616503')

# parser = argparse.ArgumentParser()
# parser.add_argument("-i", "--input", required=True, type=argparse.FileType('r'), dest='input', help="This argument is used to specify the input file")
# parser.add_argument("-o", "--output", required=True, type=argparse.FileType('w'), dest='output', help="This argument is used to specify the output file")
# parser.add_argument("-s", "--scores", required=True, type=argparse.FileType('r'), dest='scores',help="This argument is used to specify the scores file")
# args = parser.parse_args()

# inp_file = args.input
# out_file = args.output
# mtx_file = args.scores

# Used for testing against some of the input files when running into problems
inp_file = open('global_example.fna', 'r')
out_file = open('test_output.fna' , 'w')
mtx_file = open('nucleotide.mtx', 'r')

input = read_fna(inp_file)
sub_dict = read_substitution_mtx(mtx_file)

headers = input[0]
sequences = input[1]


solution_list = needleman_wunsch_align(sequences, sub_dict)

aligned_sequences = solution_list[0]
score = solution_list[1]

write_fna(out_file,headers,aligned_sequences,score)
