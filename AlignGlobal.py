# Needleman-Wuncsh Algorithm Implementation for Global Alignment

import sys

def read_fasta(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith('>'))
    return list(sequence)

def initialize_matrices(seq1, seq2, gap):
    n, m = len(seq1), len(seq2)
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    traceback = [["" for _ in range(m + 1)] for _ in range(n + 1)]
    for i in range(1, n + 1):
        score_matrix[i][0] = i * gap
        traceback[i][0] = "U"
    for j in range(1, m + 1):
        score_matrix[0][j] = j * gap
        traceback[0][j] = "L"
    traceback[0][0] = "0"
    return score_matrix, traceback

def match_score(a, b):
    return 1 if a == b else -2

def fill_matrices(seq1, seq2, score_matrix, traceback, gap):
    n, m = len(seq1), len(seq2)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = score_matrix[i - 1][j - 1] + match_score(seq1[i - 1], seq2[j - 1])
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap
            best = max(diagonal, up, left)
            score_matrix[i][j] = best
            if best == diagonal:
                traceback[i][j] = "D"
            elif best == up:
                traceback[i][j] = "U"
            else:
                traceback[i][j] = "L"
    return score_matrix, traceback

def traceback_alignment(seq1, seq2, traceback):
    alignment1, alignment2, midline = [], [], []
    i, j = len(seq1), len(seq2)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i][j] == "D":
            alignment1.append(seq1[i - 1])
            alignment2.append(seq2[j - 1])
            midline.append("|" if seq1[i - 1] == seq2[j - 1] else " ")
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or traceback[i][j] == "U"):
            alignment1.append(seq1[i - 1])
            alignment2.append("-")
            midline.append(" ")
            i -= 1
        elif j > 0 and (i == 0 or traceback[i][j] == "L"):
            alignment1.append("-")
            alignment2.append(seq2[j - 1])
            midline.append(" ")
            j -= 1
    alignment1.reverse()
    alignment2.reverse()
    midline.reverse()
    return alignment1, midline, alignment2

def needleman_wunsch(seq_file1, seq_file2, gap=-3):
    seq1 = read_fasta(seq_file1)
    seq2 = read_fasta(seq_file2)
    score_matrix, traceback_matrix = initialize_matrices(seq1, seq2, gap)
    score_matrix, traceback_matrix = fill_matrices(seq1, seq2, score_matrix, traceback_matrix, gap)
    aligned_seq1, midline, aligned_seq2 = traceback_alignment(seq1, seq2, traceback_matrix)
    print("Alignment Score:", score_matrix[len(seq1)][len(seq2)])
    print("Alignment:")
    print(" ".join(aligned_seq1))
    print(" ".join(midline))
    print(" ".join(aligned_seq2))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python AlignGlobal.py <sequence_file1> <sequence_file2>")
        sys.exit(1)
    needleman_wunsch(sys.argv[1], sys.argv[2])
