# Smith-Waterman Algorithm Implementation for Local Alignment

import sys

def read_fasta(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def smith_waterman_affine(seq1, seq2, gap_open=-5, gap_ext=-1):
    n, m = len(seq1), len(seq2)
    H = [[0]*(m+1) for _ in range(n+1)]
    E = [[float('-inf')] * (m+1) for _ in range(n+1)]
    F = [[float('-inf')] * (m+1) for _ in range(n+1)]
    
    pointer_H = [[None]*(m+1) for _ in range(n+1)]
    pointer_E = [[None]*(m+1) for _ in range(n+1)]
    pointer_F = [[None]*(m+1) for _ in range(n+1)]
    
    max_score = 0
    max_pos = (0, 0, 'H')
    
    for i in range(1, n+1):
        for j in range(1, m+1):

            E[i][j] = max(H[i][j-1] + gap_open + gap_ext, E[i][j-1] + gap_ext)
            pointer_E[i][j] = ('H', i, j-1) if H[i][j-1] + gap_open + gap_ext >= E[i][j-1] + gap_ext else ('E', i, j-1)
            
            F[i][j] = max(H[i-1][j] + gap_open + gap_ext, F[i-1][j] + gap_ext)
            pointer_F[i][j] = ('H', i-1, j) if H[i-1][j] + gap_open + gap_ext >= F[i-1][j] + gap_ext else ('F', i-1, j)
            
            diag = H[i-1][j-1] + (1 if seq1[i-1] == seq2[j-1] else -2)
            H[i][j] = max(0, diag, E[i][j], F[i][j])
            if H[i][j] == 0:
                pointer_H[i][j] = None
            elif H[i][j] == diag:
                pointer_H[i][j] = ('H', i-1, j-1)
            elif H[i][j] == E[i][j]:
                pointer_H[i][j] = ('E', i, j)
            else:
                pointer_H[i][j] = ('F', i, j)
            
            if H[i][j] > max_score:
                max_score = H[i][j]
                max_pos = (i, j, 'H')
            if E[i][j] > max_score:
                max_score = E[i][j]
                max_pos = (i, j, 'E')
            if F[i][j] > max_score:
                max_score = F[i][j]
                max_pos = (i, j, 'F')
                
    return H, E, F, pointer_H, pointer_E, pointer_F, max_score, max_pos

def traceback(seq1, seq2, H, E, F, pointer_H, pointer_E, pointer_F, start_pos, start_matrix):
    aligned_seq1, aligned_seq2, midline = [], [], []
    i, j = start_pos
    matrix_tag = start_matrix

    while True:
        if matrix_tag == 'H':
            if H[i][j] == 0 or pointer_H[i][j] is None:
                break
            prev = pointer_H[i][j]
            if prev[0] == 'H':
                aligned_seq1.append(seq1[i-1])
                aligned_seq2.append(seq2[j-1])
                midline.append('|' if seq1[i-1] == seq2[j-1] else ' ')
                i, j = prev[1], prev[2]
                matrix_tag = 'H'
            else:
                matrix_tag = prev[0]
        elif matrix_tag == 'E':
            if E[i][j] == 0 or pointer_E[i][j] is None:
                break
            prev = pointer_E[i][j]
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j-1])
            midline.append(' ')
            i, j = prev[1], prev[2]
            matrix_tag = prev[0]
        elif matrix_tag == 'F':
            if F[i][j] == 0 or pointer_F[i][j] is None:
                break
            prev = pointer_F[i][j]
            aligned_seq1.append(seq1[i-1])
            aligned_seq2.append('-')
            midline.append(' ')
            i, j = prev[1], prev[2]
            matrix_tag = prev[0]
        else:
            break

    aligned_seq1.reverse()
    aligned_seq2.reverse()
    midline.reverse()

    return ''.join(aligned_seq1), ''.join(midline), ''.join(aligned_seq2)

def local_alignment(seq_file1, seq_file2, gap_open=-5, gap_ext=-1):
    seq1 = read_fasta(seq_file1)
    seq2 = read_fasta(seq_file2)
    H, E, F, pointer_H, pointer_E, pointer_F, max_score, max_info = smith_waterman_affine(seq1, seq2, gap_open, gap_ext)
    i, j, tag = max_info
    aligned_seq1, midline, aligned_seq2 = traceback(seq1, seq2, H, E, F, pointer_H, pointer_E, pointer_F, (i, j), tag)
    print("Score:", max_score)
    print(aligned_seq1)
    print(midline)
    print(aligned_seq2)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python AlignLocal.py <sequence_file1> <sequence_file2>")
        sys.exit(1)
    local_alignment(sys.argv[1], sys.argv[2])
